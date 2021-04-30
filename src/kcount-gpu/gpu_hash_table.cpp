/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved."

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 (1) Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 (2) Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 (3) Neither the name of the University of California, Lawrence Berkeley National
 Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to
 endorse or promote products derived from this software without specific prior
 written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades
 to the features, functionality or performance of the source code ("Enhancements") to
 anyone; however, if you choose to make your Enhancements available either publicly,
 or directly to Lawrence Berkeley National Laboratory, without imposing a separate
 written license agreement for such Enhancements, then you hereby grant the following
 license: a  non-exclusive, royalty-free perpetual license to install, use, modify,
 prepare derivative works, incorporate into other computer software, distribute, and
 sublicense such enhancements or derivative works thereof, in binary and source code
 form.
*/

#include <iostream>
#include <sstream>
#include <chrono>
#include <tuple>
#include <assert.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

#include "upcxx_utils/colors.h"
#include "gpu_common.hpp"
#include "gpu_hash_table.hpp"

#include "gpu_hash_funcs.cpp"

using namespace std;
using namespace gpu_utils;
using namespace kcount_gpu;

template <int MAX_K>
struct HashTableGPUDriver<MAX_K>::HashTableDriverState {
  cudaEvent_t event;
  QuickTimer ht_timer, kernel_timer;
};

template <int MAX_K>
KmerArray<MAX_K>::KmerArray(const uint64_t *kmer) {
  memcpy(longs, kmer, N_LONGS * sizeof(cu_uint64_t));
}

template <int MAX_K>
__device__ bool kmers_equal(const KmerArray<MAX_K> &kmer1, const KmerArray<MAX_K> &kmer2) {
  int n_longs = kmer1.N_LONGS;  // get_N_LONGS();
  for (int i = 0; i < n_longs; i++) {
    if (kmer1.longs[i] != kmer2.longs[i]) return false;
  }
  return true;
}

template <int MAX_K>
__device__ size_t kmer_hash(const KmerArray<MAX_K> &kmer) {
  return gpu_murmurhash3_64(reinterpret_cast<const void *>(kmer.longs), kmer.N_LONGS * sizeof(cu_uint64_t));
}

template <int MAX_K>
HashTableGPUDriver<MAX_K>::HashTableGPUDriver() {}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::init(int upcxx_rank_me, int upcxx_rank_n, int kmer_len, int max_elems, size_t gpu_avail_mem,
                                     double &init_time, size_t &gpu_bytes_reqd) {
  QuickTimer init_timer;
  init_timer.start();
  this->upcxx_rank_me = upcxx_rank_me;
  this->upcxx_rank_n = upcxx_rank_n;
  this->kmer_len = kmer_len;
  this->gpu_thread = nullptr;
  int device_count = 0;
  cudaErrchk(cudaGetDeviceCount(&device_count));
  int my_gpu_id = upcxx_rank_me % device_count;
  cudaErrchk(cudaSetDevice(my_gpu_id));

  // now check that we have sufficient memory for the required capacity
  size_t elem_size = sizeof(KeyValue<MAX_K>) + sizeof(int);
  size_t elem_buff_size = KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * sizeof(KmerAndExts<MAX_K>);

  // set capacity to max avail from gpu memory - to reduce hash table load
  prime.set((gpu_avail_mem - elem_buff_size) / elem_size, false);
  ht_capacity = prime.get();
  gpu_bytes_reqd = max_elems * elem_size + elem_buff_size;
  if (!upcxx_rank_me)
    cout << KLMAGENTA << "Selecting GPU hash table capacity per rank of " << ht_capacity << " for " << max_elems << " elements\n";

  cudaErrchk(cudaMalloc(&elems_dev, ht_capacity * sizeof(KeyValue<MAX_K>)));
  cudaErrchk(cudaMemset(elems_dev, 0, ht_capacity * sizeof(KeyValue<MAX_K>)));

  // for transferring elements from host to gpu
  elem_buff_host = new KmerAndExts<MAX_K>[KCOUNT_GPU_HASHTABLE_BLOCK_SIZE];

  dstate = new HashTableDriverState();
  init_timer.stop();
  init_time = init_timer.get_elapsed();
}

template <int MAX_K>
HashTableGPUDriver<MAX_K>::~HashTableGPUDriver() {
  if (dstate) delete dstate;
}

template <int MAX_K>
__global__ void gpu_insert_kmer_block(KeyValue<MAX_K> *elems, const KmerAndExts<MAX_K> *elem_buff, uint32_t num_buff_entries,
                                      cu_uint64_t ht_capacity) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;

  int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  int num_inserts = 0, num_dropped = 0;
  if (threadid < num_buff_entries) {
    KmerArray<MAX_K> kmer = elem_buff[threadid].kmer;
    count_t kmer_count = elem_buff[threadid].count;
    char left_ext = elem_buff[threadid].left;
    char right_ext = elem_buff[threadid].right;
    cu_uint64_t slot = kmer_hash(kmer) % ht_capacity;
    auto start_slot = slot;
    num_inserts++;
    int j;
    // restricting the number of probes reduces computation (and imbalance) at the cost of some loss of completeness
    // for a good hash function, this should only kick in when the hash table load is getting really high
    // then we'll start to see the loss of inserts
    const int MAX_PROBE = (ht_capacity < 100 ? ht_capacity : 100);
    for (j = 0; j < MAX_PROBE; j++) {
      cu_uint64_t old_key = atomicCAS(&(elems[slot].key.longs[0]), 0, kmer.longs[0]);
      if (old_key == 0 || old_key == kmer.longs[0]) {
        bool found = true;
        for (int long_i = 1; long_i < N_LONGS; long_i++) {
          cu_uint64_t old_key = atomicCAS(&(elems[slot].key.longs[long_i]), 0, kmer.longs[long_i]);
          if (old_key != 0 && old_key != kmer.longs[long_i]) {
            found = false;
            break;
          }
        }
        if (found) {
          atomicAdd(&(elems[slot].val[0]), kmer_count);
          switch (left_ext) {
            case 'A': atomicAdd(&(elems[slot].val[1]), kmer_count); break;
            case 'C': atomicAdd(&(elems[slot].val[2]), kmer_count); break;
            case 'G': atomicAdd(&(elems[slot].val[3]), kmer_count); break;
            case 'T': atomicAdd(&(elems[slot].val[4]), kmer_count); break;
          }
          switch (right_ext) {
            case 'A': atomicAdd(&(elems[slot].val[5]), kmer_count); break;
            case 'C': atomicAdd(&(elems[slot].val[6]), kmer_count); break;
            case 'G': atomicAdd(&(elems[slot].val[7]), kmer_count); break;
            case 'T': atomicAdd(&(elems[slot].val[8]), kmer_count); break;
          }
          break;
        }
      }
      // linear probing
      // slot = (start_slot + j) % ht_capacity;
      // quadratic probing - worse cache but reduced clustering
      slot = (start_slot + (j + 1) * (j + 1)) % ht_capacity;
    }
    // this entry didn't get inserted because we ran out of probing time (and probably space)
    if (j == MAX_PROBE) num_dropped++;
  }
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::insert_kmer_block(int64_t &num_inserts, int64_t &num_dropped) {
  KmerAndExts<MAX_K> *elem_buff_dev;
  // copy across outside of thread so that we can reuse the elem_buff_host to carry on with inserts while the gpu is running
  cudaErrchk(cudaMalloc(&elem_buff_dev, num_buff_entries * sizeof(KmerAndExts<MAX_K>)));
  cudaErrchk(cudaMemcpy(elem_buff_dev, elem_buff_host, num_buff_entries * sizeof(KmerAndExts<MAX_K>), cudaMemcpyHostToDevice));
  int mingridsize = 0;
  int threadblocksize = 0;
  cudaErrchk(cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_insert_kmer_block<MAX_K>, 0, 0));
  int gridsize = ((uint32_t)num_buff_entries + threadblocksize - 1) / threadblocksize;

  cudaEvent_t start_evt, stop_evt;
  cudaErrchk(cudaEventCreate(&start_evt));
  cudaErrchk(cudaEventCreate(&stop_evt));
  cudaErrchk(cudaEventRecord(start_evt, 0));

  gpu_insert_kmer_block<<<gridsize, threadblocksize>>>(elems_dev, elem_buff_dev, num_buff_entries, ht_capacity);

  cudaErrchk(cudaEventRecord(stop_evt, 0));
  cudaErrchk(cudaEventSynchronize(stop_evt));

  float elapsed_time = 0;
  cudaErrchk(cudaEventElapsedTime(&elapsed_time, start_evt, stop_evt));

  cudaErrchk(cudaEventDestroy(start_evt));
  cudaErrchk(cudaEventDestroy(stop_evt));

  dstate->kernel_timer.inc(elapsed_time / 1000.0);

  cudaFree(elem_buff_dev);
  num_gpu_calls++;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::insert_kmer(const uint64_t *kmer, count_t kmer_count, char left, char right) {
  dstate->ht_timer.start();
  elem_buff_host[num_buff_entries].kmer = kmer;
  elem_buff_host[num_buff_entries].count = kmer_count;
  elem_buff_host[num_buff_entries].left = left;
  elem_buff_host[num_buff_entries].right = right;
  num_buff_entries++;
  if (num_buff_entries == KCOUNT_GPU_HASHTABLE_BLOCK_SIZE) {
    // cp to dev and run kernel
    insert_kmer_block(num_attempted_inserts, num_dropped_entries);
    num_buff_entries = 0;
  }
  dstate->ht_timer.stop();
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::done_inserts() {
  dstate->ht_timer.start();
  if (num_buff_entries) {
    insert_kmer_block(num_attempted_inserts, num_dropped_entries);
    num_buff_entries = 0;
  }
  // delete to make space before returning the hash table entries
  if (elem_buff_host) delete[] elem_buff_host;
  // now copy the gpu hash table values across to the host
  // We only do this once, which requires enough memory on the host to store the full GPU hash table, but since the GPU memory
  // is generally a lot less than the host memory, it should be fine.
  output_elems.resize(ht_capacity);
  output_index = 0;

  cudaErrchk(cudaMemcpy(output_elems.data(), elems_dev, ht_capacity * sizeof(KeyValue<MAX_K>), cudaMemcpyDeviceToHost));
  cudaFree(elems_dev);
  dstate->ht_timer.stop();
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::get_elapsed_time(double &gpu_time, double &gpu_kernel_time) {
  gpu_time = dstate->ht_timer.get_elapsed();
  gpu_kernel_time = dstate->kernel_timer.get_elapsed();
}

template <int MAX_K>
KeyValue<MAX_K> *HashTableGPUDriver<MAX_K>::get_next_entry() {
  if (output_elems.empty() || output_index == output_elems.size()) return nullptr;
  output_index++;
  return &(output_elems[output_index]);
}

template <int MAX_K>
int64_t HashTableGPUDriver<MAX_K>::get_capacity() {
  return ht_capacity;
}

template <int MAX_K>
int64_t HashTableGPUDriver<MAX_K>::get_num_dropped() {
  return num_dropped_entries;
}

template <int MAX_K>
int64_t HashTableGPUDriver<MAX_K>::get_num_inserts() {
  return num_attempted_inserts;
}

template <int MAX_K>
int HashTableGPUDriver<MAX_K>::get_num_gpu_calls() {
  return num_gpu_calls;
}

template class kcount_gpu::HashTableGPUDriver<32>;
#if MAX_BUILD_KMER >= 64
template class kcount_gpu::HashTableGPUDriver<64>;
#endif
#if MAX_BUILD_KMER >= 96
template class kcount_gpu::HashTableGPUDriver<96>;
#endif
#if MAX_BUILD_KMER >= 128
template class kcount_gpu::HashTableGPUDriver<128>;
#endif
#if MAX_BUILD_KMER >= 160
template class kcount_gpu::HashTableGPUDriver<160>;
#endif
