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
#include <fstream>
#include <sstream>
#include <chrono>
#include <tuple>

#include "upcxx_utils/colors.h"
#include "gpu-utils/gpu_compatibility.hpp"
#include "gpu-utils/gpu_common.hpp"
#include "gpu-utils/gpu_utils.hpp"
#include "parse_and_pack.hpp"

using namespace std;
using namespace gpu_common;

__constant__ uint64_t GPU_0_MASK[32] = {
    0x0000000000000000, 0xC000000000000000, 0xF000000000000000, 0xFC00000000000000, 0xFF00000000000000, 0xFFC0000000000000,
    0xFFF0000000000000, 0xFFFC000000000000, 0xFFFF000000000000, 0xFFFFC00000000000, 0xFFFFF00000000000, 0xFFFFFC0000000000,
    0xFFFFFF0000000000, 0xFFFFFFC000000000, 0xFFFFFFF000000000, 0xFFFFFFFC00000000, 0xFFFFFFFF00000000, 0xFFFFFFFFC0000000,
    0xFFFFFFFFF0000000, 0xFFFFFFFFFC000000, 0xFFFFFFFFFF000000, 0xFFFFFFFFFFC00000, 0xFFFFFFFFFFF00000, 0xFFFFFFFFFFFC0000,
    0xFFFFFFFFFFFF0000, 0xFFFFFFFFFFFFC000, 0xFFFFFFFFFFFFF000, 0xFFFFFFFFFFFFFC00, 0xFFFFFFFFFFFFFF00, 0xFFFFFFFFFFFFFFC0,
    0xFFFFFFFFFFFFFFF0, 0xFFFFFFFFFFFFFFFC};

struct kcount_gpu::ParseAndPackDriverState {
  Event_t event;
  int rank_me;
};

inline __device__ uint64_t quick_hash(uint64_t v) {
  v = v * 3935559000370003845 + 2691343689449507681;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}

__device__ uint64_t gpu_minimizer_hash_fast(int m, int kmer_len, int num_longs, uint64_t *longs, uint64_t *rc_longs) {
  const int chunk_step = 32 - ((m + 3) / 4) * 4;  // chunk_step is a multiple of 4

  int base;
  int num_candidates = kmer_len - m + 1;
  const int max_candidates = MAX_BUILD_KMER;
  uint64_t rc_candidates[max_candidates];

  // calculate and temporarily store all revcomp minimizer candidates on the stack
  for (base = 0; base <= kmer_len - m; base += chunk_step) {
    int shift = base % 32;
    int l = base / 32;
    uint64_t tmp = rc_longs[l];
    if (shift) {
      tmp = (tmp << (shift * 2));
      if (l < num_longs - 1) tmp |= rc_longs[l + 1] >> (64 - shift * 2);
    }
    for (int j = 0; j < chunk_step; j++) {
      if (base + j + m > kmer_len) break;
      rc_candidates[base + j] = ((tmp << (j * 2)) & GPU_0_MASK[m]);
    }
  }

  uint64_t minimizer = 0;
  // calculate and compare minimizers from revcomp
  for (base = 0; base <= kmer_len - m; base += chunk_step) {
    int shift = base % 32;
    int l = base / 32;
    uint64_t tmp = longs[l];
    if (shift) {
      tmp = (tmp << (shift * 2));
      if (l < num_longs - 1) tmp |= longs[l + 1] >> (64 - shift * 2);
    }
    for (int j = 0; j < chunk_step; j++) {
      if (base + j + m > kmer_len) break;
      uint64_t fwd_candidate = ((tmp << (j * 2)) & GPU_0_MASK[m]);
      auto &rc_candidate = rc_candidates[num_candidates - base - j - 1];
      uint64_t &least_candidate = (fwd_candidate < rc_candidate) ? fwd_candidate : rc_candidate;
      if (least_candidate > minimizer) minimizer = least_candidate;
    }
  }
  return quick_hash(minimizer);
}

__global__ void parse_and_pack(char *seqs, int minimizer_len, int kmer_len, int num_longs, int seqs_len, int *kmer_targets,
                               int num_ranks) {
  int num_kmers = seqs_len - kmer_len + 1;
  const int MAX_LONGS = (MAX_BUILD_KMER + 31) / 32;
  uint64_t kmer[MAX_LONGS];
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadid < num_kmers) {
    if (pack_seq_to_kmer(&(seqs[threadid]), kmer_len, num_longs, kmer)) {
      uint64_t kmer_rc[MAX_LONGS];
      revcomp(kmer, kmer_rc, kmer_len, num_longs);
      kmer_targets[threadid] = gpu_minimizer_hash_fast(minimizer_len, kmer_len, num_longs, kmer, kmer_rc) % num_ranks;
    } else {
      // indicate invalid with -1
      kmer_targets[threadid] = -1;
    }
  }
}

inline __device__ bool is_valid_base(char base) { return (base != '_' && base != 'N'); }

__global__ void build_supermers(char *seqs, int *kmer_targets, int num_kmers, int kmer_len, int seqs_len,
                                kcount_gpu::SupermerInfo *supermers, unsigned int *num_supermers, unsigned int *num_valid_kmers,
                                int rank_me) {
  // builds a single supermer starting at a given kmer, but only if the kmer is a valid start to a supermer
  int my_valid_kmers = 0;
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadid == 0 && kmer_targets[threadid] != -1) my_valid_kmers++;
  if (threadid > 0 && threadid < num_kmers) {
    int target = kmer_targets[threadid];
    if (target != -1) {
      my_valid_kmers++;
      bool prev_target_ok = false;
      if (threadid == 1) {
        prev_target_ok = true;
      } else {
        if (kmer_targets[threadid - 1] != target) {
          // prev kmer was a different or invalid target
          prev_target_ok = true;
        } else {
          // prev kmer was the same target, but was not a valid start to a supermer
          if (!is_valid_base(seqs[threadid - 2]) || !is_valid_base(seqs[threadid - 1 + kmer_len])) prev_target_ok = true;
        }
      }
      // make sure this is the first kmer for this target
      if (prev_target_ok && is_valid_base(seqs[threadid - 1]) && is_valid_base(seqs[threadid + kmer_len])) {
        int supermer_start_i = threadid - 1;
        int supermer_len = kmer_len + 2;
        // build the supermer
        for (int i = threadid + 1; i < num_kmers - 1; i++) {
          auto next_target = kmer_targets[i];
          int end_pos = supermer_start_i + supermer_len;  // i + kmer_len;
          if (next_target == target && end_pos < seqs_len && is_valid_base(seqs[end_pos]))
            supermer_len++;
          else
            break;
        }
        // get a slot for the supermer
        int slot = atomicAdd(num_supermers, 1);
        supermers[slot].target = target;
        supermers[slot].offset = supermer_start_i;
        supermers[slot].len = supermer_len;
      }
    }
  }
  reduce(my_valid_kmers, num_kmers, num_valid_kmers);
}

static const uint8_t INVALID = 11;

inline __device__ uint8_t get_packed_val(char base) {
  switch (base) {
    case 'a': return 1;
    case 'c': return 2;
    case 'g': return 3;
    case 't': return 4;
    case 'A': return 5;
    case 'C': return 6;
    case 'G': return 7;
    case 'T': return 8;
    case 'N':
    case 'n': return 9;
    case '_':
    case 0: return 0;
    default: return INVALID;  // printf("Invalid value encountered when packing: %d\n", (int)base);
  };
  return 0;
}

__global__ void pack_seqs(char *dev_seqs, char *dev_packed_seqs, int seqs_len) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  int packed_seqs_len = (seqs_len + 1) / 2;
  if (threadid < packed_seqs_len) {
    int seqs_i = threadid * 2;
    char packed = get_packed_val(dev_seqs[seqs_i]);
    if ((int)packed == INVALID) {
      printf("INVALID dev_seqs[%d]=%d, after shifting:%d, packed:%d, seqs_len:%d, packed_seq_len:%d\n", seqs_i, dev_seqs[seqs_i],
             packed << 4, packed, seqs_len, packed_seqs_len);
    }
    packed = packed << 4;
    if (seqs_i + 1 < seqs_len) {
      // do not overflow as each thread handles 1-2 characters in the sequence
      char packed_ = get_packed_val(dev_seqs[seqs_i + 1]);
      if ((int)packed_ == INVALID) {
        printf("INVALID dev_seqs[%d]=%d, dev_seqs[%d]=%d, after shifting:%d, packed:%d, seqs_len:%d, packed_seq_len:%d\n", seqs_i,
               dev_seqs[seqs_i], seqs_i + 1, dev_seqs[seqs_i + 1], packed_ << 4, packed_, seqs_len, packed_seqs_len);
      }
      packed |= packed_;
    }
    dev_packed_seqs[threadid] = packed;
  }
}

kcount_gpu::ParseAndPackGPUDriver::ParseAndPackGPUDriver(int upcxx_rank_me, int upcxx_rank_n, int qual_offset, int kmer_len,
                                                         int num_kmer_longs, int minimizer_len, double &init_time)
    : upcxx_rank_me(upcxx_rank_me)
    , upcxx_rank_n(upcxx_rank_n)
    , kmer_len(kmer_len)
    , qual_offset(qual_offset)
    , num_kmer_longs(num_kmer_longs)
    , minimizer_len(minimizer_len)
    , t_func(0)
    , t_kernel(0) {
  QuickTimer init_timer;
  init_timer.start();
  gpu_utils::set_gpu_device(upcxx_rank_me);
  max_kmers = KCOUNT_SEQ_BLOCK_SIZE - kmer_len + 1;

  ERROR_CHECK(Malloc((void **)&dev_seqs, KCOUNT_SEQ_BLOCK_SIZE));
  ERROR_CHECK(Malloc((void **)&dev_kmer_targets, max_kmers * sizeof(int)));

  ERROR_CHECK(Malloc((void **)&dev_supermers, max_kmers * sizeof(SupermerInfo)));
  ERROR_CHECK(Malloc((void **)&dev_packed_seqs, (KCOUNT_SEQ_BLOCK_SIZE + 1) / 2));
  ERROR_CHECK(Malloc((void **)&dev_num_supermers, sizeof(int)));
  ERROR_CHECK(Malloc((void **)&dev_num_valid_kmers, sizeof(int)));

  // total storage required is approx KCOUNT_SEQ_BLOCK_SIZE * (1 + num_kmers_longs * sizeof(uint64_t) + sizeof(int) + 1)
  dstate = new ParseAndPackDriverState();
  dstate->rank_me = upcxx_rank_me;
  init_timer.stop();
  init_time = init_timer.get_elapsed();
}

kcount_gpu::ParseAndPackGPUDriver::~ParseAndPackGPUDriver() {
  ERROR_CHECK(Free(dev_seqs));
  ERROR_CHECK(Free(dev_kmer_targets));

  ERROR_CHECK(Free(dev_supermers));
  ERROR_CHECK(Free(dev_packed_seqs));
  ERROR_CHECK(Free(dev_num_supermers));
  ERROR_CHECK(Free(dev_num_valid_kmers));

  delete dstate;
}

bool kcount_gpu::ParseAndPackGPUDriver::process_seq_block(const string &seqs, unsigned int &num_valid_kmers) {
  QuickTimer func_timer, kernel_timer;

  if (seqs.length() >= KCOUNT_SEQ_BLOCK_SIZE) return false;
  if (seqs.length() == 0) return false;
  if (seqs.length() < (unsigned int)kmer_len) return false;

  func_timer.start();
  gpu_utils::set_gpu_device(dstate->rank_me);
  ERROR_CHECK(EventCreateWithFlags(&dstate->event, EventDisableTiming | EventBlockingSync));

  int num_kmers = seqs.length() - kmer_len + 1;
  ERROR_CHECK(Memcpy(dev_seqs, &seqs[0], seqs.length(), MemcpyHostToDevice));

  int gridsize, threadblocksize;
  get_kernel_config(seqs.length(), parse_and_pack, gridsize, threadblocksize);
  kernel_timer.start();
  LaunchKernel(parse_and_pack, gridsize, threadblocksize, dev_seqs, minimizer_len, kmer_len, num_kmer_longs, seqs.length(),
               dev_kmer_targets, upcxx_rank_n);

  ERROR_CHECK(Memset(dev_num_supermers, 0, sizeof(int)));
  ERROR_CHECK(Memset(dev_num_valid_kmers, 0, sizeof(int)));
  get_kernel_config(num_kmers, build_supermers, gridsize, threadblocksize);
  LaunchKernel(build_supermers, gridsize, threadblocksize, dev_seqs, dev_kmer_targets, num_kmers, kmer_len, seqs.length(),
               dev_supermers, dev_num_supermers, dev_num_valid_kmers, upcxx_rank_me);
  ERROR_CHECK(Memcpy(&num_valid_kmers, dev_num_valid_kmers, sizeof(unsigned int), MemcpyDeviceToHost));
  unsigned int num_supermers;
  ERROR_CHECK(Memcpy(&num_supermers, dev_num_supermers, sizeof(unsigned int), MemcpyDeviceToHost));
  supermers.resize(num_supermers);
  ERROR_CHECK(Memcpy(&(supermers[0]), dev_supermers, num_supermers * sizeof(SupermerInfo), MemcpyDeviceToHost));
  ERROR_CHECK(EventSynchronize(dstate->event));
  ERROR_CHECK(EventDestroy(dstate->event));
  kernel_timer.stop();
  t_kernel += kernel_timer.get_elapsed();
  func_timer.stop();
  t_func += func_timer.get_elapsed();
  return true;
}

void kcount_gpu::ParseAndPackGPUDriver::pack_seq_block(const string &seqs) {
  gpu_utils::set_gpu_device(dstate->rank_me);
  int packed_seqs_len = (seqs.length() + 1) / 2;
  ERROR_CHECK(Memcpy(dev_seqs, &seqs[0], seqs.length(), MemcpyHostToDevice));
  ERROR_CHECK(Memset(dev_packed_seqs, 0, packed_seqs_len));
  int gridsize, threadblocksize;
  get_kernel_config(packed_seqs_len, pack_seqs, gridsize, threadblocksize);
  GPUTimer t;
  t.start();
  LaunchKernel(pack_seqs, gridsize, threadblocksize, dev_seqs, dev_packed_seqs, seqs.length());
  // this GPUTimer forces a wait for the GPU kernel to complete
  t.stop();
  t_kernel += t.get_elapsed();
  packed_seqs.resize(packed_seqs_len);
  ERROR_CHECK(Memcpy(&(packed_seqs[0]), dev_packed_seqs, packed_seqs_len, MemcpyDeviceToHost));
}

tuple<double, double> kcount_gpu::ParseAndPackGPUDriver::get_elapsed_times() { return {t_func, t_kernel}; }
