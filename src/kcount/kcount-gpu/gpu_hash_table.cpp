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
#include <fstream>
#include <chrono>
#include <tuple>
#include <iomanip>
#include <assert.h>

#include "upcxx_utils/colors.h"
#include "gpu-utils/gpu_compatibility.hpp"
#include "gpu-utils/gpu_common.hpp"
#include "gpu-utils/gpu_utils.hpp"
#include "gpu_hash_table.hpp"
#include "prime.hpp"
#include "gpu_hash_funcs.hpp"
#ifdef USE_TCF
#include "tcf_wrapper.hpp"
#else
// quotient filter calls stubbed out
namespace two_choice_filter {
#define TCF_RESULT uint8_t
struct TCF {
  static TCF *generate_on_device(bool *, int) { return nullptr; }
  static void free_on_device(TCF *) {}
  __device__ bool get_my_tile() { return false; }
  __device__ bool insert_if_not_exists(bool b, uint64_t, bool, TCF_RESULT, bool) { return false; }
  int get_fill() { return 0; }
  int get_num_slots() { return 0; }
};
__device__ uint8_t pack_extensions(char left, char right) { return 0; }
__device__ bool unpack_extensions(uint8_t storage, char &left, char &right) { return false; }
static uint64_t estimate_memory(uint64_t max_num_kmers) { return 0; }
static bool get_tcf_sizing_from_mem(uint64_t available_bytes) { return false; }
}  // namespace two_choice_filter
#endif

using namespace std;
using namespace gpu_common;
using namespace kcount_gpu;

// convenience functions
#define SDBG(fmt, ...) \
  if (!upcxx_rank_me) printf(KLMAGENTA "GPU kcount: " fmt KNORM "\n", ##__VA_ARGS__)

#define SWARN(fmt, ...) \
  if (!upcxx_rank_me) printf(KLRED "WARN GPU kcount: " fmt KNORM "\n", ##__VA_ARGS__)

#define WARN(fmt, ...) printf(KLRED "WARN GPU kcount: " fmt KNORM "\n", ##__VA_ARGS__)

const uint64_t KEY_EMPTY = 0xffffffffffffffff;
const uint64_t KEY_TRANSITION = 0xfffffffffffffffe;
const uint8_t KEY_EMPTY_BYTE = 0xff;

template <int MAX_K>
__device__ void kmer_set(KmerArray<MAX_K> &kmer1, const KmerArray<MAX_K> &kmer2) {
  int N_LONGS = kmer1.N_LONGS;
  uint64_t old_key;
  for (int i = 0; i < N_LONGS - 1; i++) {
    old_key = atomicExch((unsigned long long *)&(kmer1.longs[i]), kmer2.longs[i]);
    if (old_key != KEY_EMPTY) WARN("old key should be KEY_EMPTY");
  }
  old_key = atomicExch((unsigned long long *)&(kmer1.longs[N_LONGS - 1]), kmer2.longs[N_LONGS - 1]);
  if (old_key != KEY_TRANSITION) WARN("old key should be KEY_TRANSITION");
}

template <int MAX_K>
__device__ bool kmers_equal(const KmerArray<MAX_K> &kmer1, const KmerArray<MAX_K> &kmer2) {
  int n_longs = kmer1.N_LONGS;
  for (int i = 0; i < n_longs; i++) {
    uint64_t old_key = atomicAdd((unsigned long long *)&(kmer1.longs[i]), 0ULL);
    if (old_key != kmer2.longs[i]) return false;
  }
  return true;
}

template <int MAX_K>
__device__ size_t kmer_hash(const KmerArray<MAX_K> &kmer) {
  return gpu_murmurhash3_64(reinterpret_cast<const void *>(kmer.longs), kmer.N_LONGS * sizeof(uint64_t));
}

__device__ int8_t get_ext(CountsArray &counts, int pos, int8_t *ext_map) {
  count_t top_count = 0, runner_up_count = 0;
  int top_ext_pos = 0;
  count_t kmer_count = counts.kmer_count;
  for (int i = pos; i < pos + 4; i++) {
    if (counts.ext_counts[i] >= top_count) {
      runner_up_count = top_count;
      top_count = counts.ext_counts[i];
      top_ext_pos = i;
    } else if (counts.ext_counts[i] > runner_up_count) {
      runner_up_count = counts.ext_counts[i];
    }
  }
  int dmin_dyn = (1.0 - DYN_MIN_DEPTH) * kmer_count;
  if (dmin_dyn < 2.0) dmin_dyn = 2.0;
  if (top_count < dmin_dyn) return 'X';
  if (runner_up_count >= dmin_dyn) return 'F';
  return ext_map[top_ext_pos - pos];
}

__device__ bool ext_conflict(ext_count_t *ext_counts, int start_idx) {
  int idx = -1;
  for (int i = start_idx; i < start_idx + 4; i++) {
    if (ext_counts[i]) {
      // conflict
      if (idx != -1) return true;
      idx = i;
    }
  }
  return false;
}

template <int MAX_K>
__global__ void gpu_merge_ctg_kmers(KmerCountsMap<MAX_K> read_kmers, const KmerCountsMap<MAX_K> ctg_kmers,
                                    uint64_t *insert_counts) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  int8_t ext_map[4] = {'A', 'C', 'G', 'T'};
  int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  uint64_t attempted_inserts = 0;
  uint64_t dropped_inserts = 0;
  uint64_t new_inserts = 0;
  if (threadid < ctg_kmers.capacity) {
    count_t kmer_count = ctg_kmers.vals[threadid].kmer_count;
    ext_count_t *ext_counts = ctg_kmers.vals[threadid].ext_counts;
    if (kmer_count && !ext_conflict(ext_counts, 0) && !ext_conflict(ext_counts, 4)) {
      KmerArray<MAX_K> kmer = ctg_kmers.keys[threadid];
      uint64_t slot = kmer_hash(kmer) % read_kmers.capacity;
      auto start_slot = slot;
      attempted_inserts++;
      const int MAX_PROBE = (read_kmers.capacity < KCOUNT_HT_MAX_PROBE ? read_kmers.capacity : KCOUNT_HT_MAX_PROBE);
      for (int j = 0; j < MAX_PROBE; j++) {
        uint64_t old_key = atomicCAS((unsigned long long *)&(read_kmers.keys[slot].longs[N_LONGS - 1]), KEY_EMPTY, KEY_TRANSITION);
        if (old_key == KEY_EMPTY) {
          new_inserts++;
          memcpy(&read_kmers.vals[slot], &ctg_kmers.vals[threadid], sizeof(CountsArray));
          kmer_set(read_kmers.keys[slot], kmer);
          break;
        } else if (old_key == kmer.longs[N_LONGS - 1]) {
          if (kmers_equal(read_kmers.keys[slot], kmer)) {
            // existing kmer from reads - only replace if the kmer is non-UU
            // there is no need for atomics here because all ctg kmers are unique; hence only one thread will ever match this kmer
            int8_t left_ext = get_ext(read_kmers.vals[slot], 0, ext_map);
            int8_t right_ext = get_ext(read_kmers.vals[slot], 4, ext_map);
            if (left_ext == 'X' || left_ext == 'F' || right_ext == 'X' || right_ext == 'F')
              memcpy(&read_kmers.vals[slot], &ctg_kmers.vals[threadid], sizeof(CountsArray));
            break;
          }
        }
        // quadratic probing - worse cache but reduced clustering
        slot = (start_slot + (j + 1) * (j + 1)) % read_kmers.capacity;
        if (j == MAX_PROBE - 1) dropped_inserts++;
      }
    }
  }
  reduce(attempted_inserts, ctg_kmers.capacity, &(insert_counts[0]));
  reduce(dropped_inserts, ctg_kmers.capacity, &(insert_counts[1]));
  reduce(new_inserts, ctg_kmers.capacity, &(insert_counts[2]));
}

template <int MAX_K>
__global__ void gpu_compact_ht(KmerCountsMap<MAX_K> elems, KmerExtsMap<MAX_K> compact_elems, uint64_t *elem_counts) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  uint64_t dropped_inserts = 0;
  uint64_t unique_inserts = 0;
  int8_t ext_map[4] = {'A', 'C', 'G', 'T'};
  if (threadid < elems.capacity) {
    if (elems.vals[threadid].kmer_count) {
      KmerArray<MAX_K> kmer = elems.keys[threadid];
      uint64_t slot = kmer_hash(kmer) % compact_elems.capacity;
      auto start_slot = slot;
      // we set a constraint on the max probe to track whether we are getting excessive collisions and need a bigger default
      // compact table
      const int MAX_PROBE = (compact_elems.capacity < KCOUNT_HT_MAX_PROBE ? compact_elems.capacity : KCOUNT_HT_MAX_PROBE);
      // look for empty slot in compact hash table
      for (int j = 0; j < MAX_PROBE; j++) {
        uint64_t old_key =
            atomicCAS((unsigned long long *)&(compact_elems.keys[slot].longs[N_LONGS - 1]), KEY_EMPTY, kmer.longs[N_LONGS - 1]);
        if (old_key == KEY_EMPTY) {
          // found empty slot - there will be no duplicate keys since we're copying across from another hash table
          unique_inserts++;
          memcpy((void *)compact_elems.keys[slot].longs, kmer.longs, sizeof(uint64_t) * (N_LONGS - 1));
          // compute exts
          int8_t left_ext = get_ext(elems.vals[threadid], 0, ext_map);
          int8_t right_ext = get_ext(elems.vals[threadid], 4, ext_map);
          if (elems.vals[threadid].kmer_count < 2) WARN("elem should have been purged, count %d", elems.vals[threadid].kmer_count);
          compact_elems.vals[slot].count = elems.vals[threadid].kmer_count;
          compact_elems.vals[slot].left = left_ext;
          compact_elems.vals[slot].right = right_ext;
          break;
        }
        // quadratic probing - worse cache but reduced clustering
        slot = (start_slot + (j + 1) * (j + 1)) % compact_elems.capacity;
        if (j == MAX_PROBE - 1) dropped_inserts++;
      }
    }
  }
  reduce(dropped_inserts, compact_elems.capacity, &(elem_counts[0]));
  reduce(unique_inserts, compact_elems.capacity, &(elem_counts[1]));
}

template <int MAX_K>
__global__ void gpu_purge_invalid(KmerCountsMap<MAX_K> elems, uint64_t *elem_counts) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  uint64_t num_purged = 0;
  uint64_t num_elems = 0;
  if (threadid < elems.capacity) {
    if (elems.vals[threadid].kmer_count) {
      uint64_t ext_sum = 0;
      for (int j = 0; j < 8; j++) ext_sum += elems.vals[threadid].ext_counts[j];
      if (elems.vals[threadid].kmer_count < 2 || !ext_sum) {
        memset(&elems.vals[threadid], 0, sizeof(CountsArray));
        memset((void *)elems.keys[threadid].longs, KEY_EMPTY_BYTE, N_LONGS * sizeof(uint64_t));
        num_purged++;
      } else {
        num_elems++;
      }
    }
  }
  reduce(num_purged, elems.capacity, &(elem_counts[0]));
  reduce(num_elems, elems.capacity, &(elem_counts[1]));
}

static __constant__ char to_base[] = {'0', 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T', 'N'};

inline __device__ char to_base_func(int index, int pp) {
  if (index > 9) {
    WARN("index out of range for to_base: %d, packed seq pos %d", index, pp);
    return 0;
  }
  if (index == 0) return '_';
  return to_base[index];
}

__global__ void gpu_unpack_supermer_block(SupermerBuff unpacked_supermer_buff, SupermerBuff packed_supermer_buff, int buff_len) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadid >= buff_len) return;
  uint8_t packed = packed_supermer_buff.seqs[threadid];
  if (packed == '_') return;
  uint8_t left_side = (packed & 240) >> 4;
  unpacked_supermer_buff.seqs[threadid * 2] = to_base_func(left_side, packed);
  if (packed_supermer_buff.counts) unpacked_supermer_buff.counts[threadid * 2] = packed_supermer_buff.counts[threadid];
  uint8_t right_side = packed & 15;
  unpacked_supermer_buff.seqs[threadid * 2 + 1] = to_base_func(right_side, packed);
  if (packed_supermer_buff.counts) unpacked_supermer_buff.counts[threadid * 2 + 1] = packed_supermer_buff.counts[threadid];
}

inline __device__ bool is_valid_base(char base) {
  return (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == '0' || base == 'N');
}

inline __device__ bool bad_qual(char base) { return (base == 'a' || base == 'c' || base == 'g' || base == 't'); }

inline __device__ void inc_ext(char ext, ext_count_t kmer_count, ext_count_t *ext_counts) {
  switch (ext) {
    case 'A': atomicAddUint16_thres(&(ext_counts[0]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
    case 'C': atomicAddUint16_thres(&(ext_counts[1]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
    case 'G': atomicAddUint16_thres(&(ext_counts[2]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
    case 'T': atomicAddUint16_thres(&(ext_counts[3]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
  }
}

template <int MAX_K>
__device__ bool get_kmer_from_supermer(SupermerBuff supermer_buff, uint32_t buff_len, int kmer_len, uint64_t *kmer, char &left_ext,
                                       char &right_ext, count_t &count) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  int num_kmers = buff_len - kmer_len + 1;
  if (threadid >= num_kmers) return false;
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  if (!pack_seq_to_kmer(&(supermer_buff.seqs[threadid]), kmer_len, N_LONGS, kmer)) return false;
  if (threadid + kmer_len >= buff_len) return false;  // printf("out of bounds %d >= %d\n", threadid + kmer_len, buff_len);
  left_ext = supermer_buff.seqs[threadid - 1];
  right_ext = supermer_buff.seqs[threadid + kmer_len];
  if (left_ext == '_' || right_ext == '_') return false;
  if (!left_ext || !right_ext) return false;
  if (supermer_buff.counts) {
    count = supermer_buff.counts[threadid];
  } else {
    count = 1;
    if (bad_qual(left_ext)) left_ext = '0';
    if (bad_qual(right_ext)) right_ext = '0';
  }
  if (!is_valid_base(left_ext)) {
    WARN("threadid %d, invalid char for left nucleotide %d", threadid, (uint8_t)left_ext);
    return false;
  }
  if (!is_valid_base(right_ext)) {
    WARN("threadid %d, invalid char for right nucleotide %d", threadid, (uint8_t)right_ext);
    return false;
  }
  uint64_t kmer_rc[N_LONGS];
  revcomp(kmer, kmer_rc, kmer_len, N_LONGS);
  for (int l = 0; l < N_LONGS; l++) {
    if (kmer_rc[l] == kmer[l]) continue;
    if (kmer_rc[l] < kmer[l]) {
      // swap
      char tmp = left_ext;
      left_ext = comp_nucleotide(right_ext);
      right_ext = comp_nucleotide(tmp);

      // FIXME: we should be able to have a 0 extension even for revcomp - we do for non-revcomp
      // if (!left_ext || !right_ext) return false;

      memcpy(kmer, kmer_rc, N_LONGS * sizeof(uint64_t));
    }
    break;
  }
  return true;
}

template <int MAX_K>
__device__ bool gpu_insert_kmer(KmerCountsMap<MAX_K> elems, uint64_t hash_val, KmerArray<MAX_K> &kmer, char left_ext,
                                char right_ext, char prev_left_ext, char prev_right_ext, count_t kmer_count, uint64_t &new_inserts,
                                uint64_t &dropped_inserts, bool ctg_kmers, bool use_qf, bool update_only) {
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  uint64_t slot = hash_val % elems.capacity;
  auto start_slot = slot;
  const int MAX_PROBE = (elems.capacity < 200 ? elems.capacity : 200);
  bool found_slot = false;
  bool kmer_found_in_ht = false;
  uint64_t old_key = KEY_TRANSITION;
  for (int j = 0; j < MAX_PROBE; j++) {
    // we have to be careful here not to end up with multiple threads on the same warp accessing the same slot, because
    // that will cause a deadlock. So we loop over all statements in each CAS spin to ensure that all threads get a
    // chance to execute
    do {
      old_key = atomicCAS((unsigned long long *)&(elems.keys[slot].longs[N_LONGS - 1]), KEY_EMPTY, KEY_TRANSITION);
      if (old_key != KEY_TRANSITION) {
        if (old_key == KEY_EMPTY) {
          if (update_only) {
            old_key = atomicExch((unsigned long long *)&(elems.keys[slot].longs[N_LONGS - 1]), KEY_EMPTY);
            if (old_key != KEY_TRANSITION) WARN("old key should be KEY_TRANSITION");
            return false;
          }
          kmer_set(elems.keys[slot], kmer);
          found_slot = true;
        } else if (old_key == kmer.longs[N_LONGS - 1]) {
          if (kmers_equal(elems.keys[slot], kmer)) {
            found_slot = true;
            kmer_found_in_ht = true;
          }
        }
      }
    } while (old_key == KEY_TRANSITION);
    if (found_slot) break;
    // quadratic probing - worse cache but reduced clustering
    slot = (start_slot + j * j) % elems.capacity;
    // this entry didn't get inserted because we ran out of probing time (and probably space)
    if (j == MAX_PROBE - 1) dropped_inserts++;
  }
  if (found_slot) {
    ext_count_t *ext_counts = elems.vals[slot].ext_counts;
    if (ctg_kmers) {
      // the count is the min of all counts. Use CAS to deal with the initial zero value
      int prev_count = atomicCAS(&elems.vals[slot].kmer_count, 0, kmer_count);
      if (prev_count)
        atomicMin(&elems.vals[slot].kmer_count, kmer_count);
      else
        new_inserts++;
    } else {
      assert(kmer_count == 1);
      int prev_count = atomicAdd(&elems.vals[slot].kmer_count, kmer_count);
      if (!prev_count) new_inserts++;
    }
    ext_count_t kmer_count_uint16 = min(kmer_count, UINT16_MAX);
    inc_ext(left_ext, kmer_count_uint16, ext_counts);
    inc_ext(right_ext, kmer_count_uint16, ext_counts + 4);
    if (use_qf && !update_only && !kmer_found_in_ht && !ctg_kmers) {
      // kmer was not in hash table, so it must have been found in the qf
      // add the extensions from the previous entry stored in the qf
      inc_ext(prev_left_ext, 1, ext_counts);
      inc_ext(prev_right_ext, 1, ext_counts + 4);
      // inc the overall kmer count
      atomicAdd(&elems.vals[slot].kmer_count, 1);
    }
  }
  return true;
}

template <int MAX_K>
__global__ void gpu_insert_supermer_block(KmerCountsMap<MAX_K> elems, SupermerBuff supermer_buff, uint32_t buff_len, int kmer_len,
                                          bool ctg_kmers, InsertStats *insert_stats, two_choice_filter::TCF *tcf) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  uint64_t attempted_inserts = 0, dropped_inserts = 0, new_inserts = 0, num_unique_qf = 0, dropped_inserts_qf = 0;
  if (threadid > 0 && threadid < buff_len) {
    attempted_inserts++;
    KmerArray<MAX_K> kmer;
    char left_ext, right_ext;
    count_t kmer_count;
    if (get_kmer_from_supermer<MAX_K>(supermer_buff, buff_len, kmer_len, kmer.longs, left_ext, right_ext, kmer_count)) {
      if (kmer.longs[N_LONGS - 1] == KEY_EMPTY) WARN("block equal to KEY_EMPTY");
      if (kmer.longs[N_LONGS - 1] == KEY_TRANSITION) WARN("block equal to KEY_TRANSITION");
      auto hash_val = kmer_hash(kmer);
      char prev_left_ext = '0', prev_right_ext = '0';
      bool use_qf = (tcf != nullptr);
      bool update_only = (use_qf && !ctg_kmers);
      bool updated = gpu_insert_kmer(elems, hash_val, kmer, left_ext, right_ext, prev_left_ext, prev_right_ext, kmer_count,
                                     new_inserts, dropped_inserts, ctg_kmers, use_qf, update_only);

      if (update_only && !updated) {
        // not found in the hash table - look in the qf
        bool found = false;

        auto packed = two_choice_filter::pack_extensions(left_ext, right_ext);

        TCF_RESULT result = 0;
        bool success = tcf->insert_if_not_exists(tcf->get_my_tile(), hash_val, packed, result, found);

        if (success) {
          if (!found) {
            // inserted successfully
            num_unique_qf++;

            // does this need to be asserted?
            assert(prev_left_ext == '0' && prev_right_ext == '0');

          } else {
            // found successfully
            two_choice_filter::unpack_extensions(result, prev_left_ext, prev_right_ext);
            gpu_insert_kmer(elems, hash_val, kmer, left_ext, right_ext, prev_left_ext, prev_right_ext, kmer_count, new_inserts,
                            dropped_inserts, ctg_kmers, use_qf, false);
          }

        } else {
          // dropped
          dropped_inserts_qf++;
          // now insert it into the main hash table - this will be purged later if it's a singleton
          gpu_insert_kmer(elems, hash_val, kmer, left_ext, right_ext, prev_left_ext, prev_right_ext, kmer_count, new_inserts,
                          dropped_inserts, ctg_kmers, false, false);
        }
      }
    }
  }
  reduce(attempted_inserts, buff_len, &insert_stats->attempted);
  reduce(dropped_inserts, buff_len, &insert_stats->dropped);
  reduce(dropped_inserts_qf, buff_len, &insert_stats->dropped_qf);
  reduce(new_inserts, buff_len, &insert_stats->new_inserts);
  reduce(num_unique_qf, buff_len, &insert_stats->num_unique_qf);
}

template <int MAX_K>
struct HashTableGPUDriver<MAX_K>::HashTableDriverState {
  Event_t event;
  QuickTimer insert_timer, kernel_timer;
  two_choice_filter::TCF *tcf = nullptr;
};

template <int MAX_K>
void KmerArray<MAX_K>::set(const uint64_t *kmer) {
  memcpy(longs, kmer, N_LONGS * sizeof(uint64_t));
}

template <int MAX_K>
void KmerCountsMap<MAX_K>::init(int64_t ht_capacity) {
  capacity = ht_capacity;
  ERROR_CHECK(Malloc(&keys, capacity * sizeof(KmerArray<MAX_K>)));
  ERROR_CHECK(Memset((void *)keys, KEY_EMPTY_BYTE, capacity * sizeof(KmerArray<MAX_K>)));
  ERROR_CHECK(Malloc(&vals, capacity * sizeof(CountsArray)));
  ERROR_CHECK(Memset(vals, 0, capacity * sizeof(CountsArray)));
}

template <int MAX_K>
void KmerCountsMap<MAX_K>::clear() {
  ERROR_CHECK(Free((void *)keys));
  ERROR_CHECK(Free(vals));
}

template <int MAX_K>
void KmerExtsMap<MAX_K>::init(int64_t ht_capacity) {
  capacity = ht_capacity;
  ERROR_CHECK(Malloc(&keys, capacity * sizeof(KmerArray<MAX_K>)));
  ERROR_CHECK(Memset((void *)keys, KEY_EMPTY_BYTE, capacity * sizeof(KmerArray<MAX_K>)));
  ERROR_CHECK(Malloc(&vals, capacity * sizeof(CountExts)));
  ERROR_CHECK(Memset(vals, 0, capacity * sizeof(CountExts)));
}

template <int MAX_K>
void KmerExtsMap<MAX_K>::clear() {
  ERROR_CHECK(Free((void *)keys));
  ERROR_CHECK(Free(vals));
}

template <int MAX_K>
HashTableGPUDriver<MAX_K>::HashTableGPUDriver() {}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::init(int upcxx_rank_me, int upcxx_rank_n, int kmer_len, size_t max_elems, size_t max_ctg_elems,
                                     size_t num_errors, size_t gpu_avail_mem, string &msgs, string &warnings, bool use_qf) {
  this->upcxx_rank_me = upcxx_rank_me;
  this->upcxx_rank_n = upcxx_rank_n;
  this->kmer_len = kmer_len;
  pass_type = READ_KMERS_PASS;
  gpu_utils::set_gpu_device(upcxx_rank_me);
  dstate = new HashTableDriverState();

  // reserve space for the fixed size buffer for passing data to the GPU
  size_t elem_buff_size = KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * (3 + sizeof(count_t));
  gpu_avail_mem -= elem_buff_size;
  ostringstream log_msgs, log_warnings;
  log_msgs << "Elem buff size " << elem_buff_size << " (avail mem now " << gpu_avail_mem << ")\n";
  size_t elem_size = sizeof(KmerArray<MAX_K>) + sizeof(CountsArray);
  // expected size of compact hash table
  size_t compact_elem_size = sizeof(KmerArray<MAX_K>) + sizeof(CountExts);

  double elem_size_ratio = (double)compact_elem_size / (double)elem_size;
  log_msgs << "Element size for main HT " << elem_size << " and for compact HT " << compact_elem_size << " (ratio " << fixed
           << setprecision(3) << elem_size_ratio << ")\n";

  double target_load_factor = 0.66;
  double load_multiplier = 1.0 / target_load_factor;
  // There are several different structures that all have to fit in the GPU memory. We first compute the
  // memory required by all of them at the target load factor, and then reduce uniformly if there is insufficient
  // 1. The read kmers hash table. With the QF, this is the size of the number of unique kmers. Without the QF,
  //    it is that size plus the size of the errors. In addition, this hash table needs to be big enough to have
  //    all the ctg kmers added too,
  size_t max_read_kmers = load_multiplier * (max_elems + max_ctg_elems + (use_qf ? 0 : num_errors));
  size_t read_kmers_size = max_read_kmers * elem_size;
  // 2. The QF, if used. This is the size of all the unique read kmers plus the errors, plus some wiggle room. The
  //    QF uses so little memory that we can afford to oversize some
  size_t max_qf_kmers = load_multiplier * (use_qf ? max_elems + num_errors : 0) * 1.3;
  size_t qf_size = use_qf ? two_choice_filter::estimate_memory(max(1.0, log2(max_qf_kmers))) : 0;
  // 3. The ctg kmers hash table (only present if this is not the first contigging round)
  size_t max_ctg_kmers = load_multiplier * max_ctg_elems;
  size_t ctg_kmers_size = max_ctg_kmers * elem_size;
  // 4. The final compact hash table, which is the size needed to store all the unique kmers from both the reads and contigs.
  size_t max_compact_kmers = load_multiplier * (max_elems + max_ctg_elems);
  size_t compact_kmers_size = max_compact_kmers * compact_elem_size;

  log_msgs << "Element counts: read kmers " << max_read_kmers << ", qf " << max_qf_kmers << ", ctg kmer " << max_ctg_kmers
           << ", compact ht " << max_compact_kmers << "\n";
  //  for the total size, the read kmer hash table must exist with just the QF, then just the ctg kmers, then just the compact kmers
  //  so we choose the largest of these options
  size_t tot_size = read_kmers_size + max(qf_size, max(ctg_kmers_size, compact_kmers_size));
  log_msgs << "Hash table sizes: read kmers " << read_kmers_size << ", qf " << qf_size << ", ctg kmers " << ctg_kmers_size
           << ", compact ht " << compact_kmers_size << ", total " << tot_size << "\n";

  // keep some in reserve as a buffer
  double mem_ratio = (double)(0.8 * gpu_avail_mem) / tot_size;
  if (mem_ratio < 0.9)
    log_warnings << "Insufficent memory for " << fixed << setprecision(3) << target_load_factor
                 << " load factor across all data structures; reducing by a factor of " << mem_ratio
                 << "; this could result in an OOM";
  max_read_kmers *= mem_ratio;
  max_qf_kmers *= mem_ratio;
  max_ctg_kmers *= mem_ratio;
  max_compact_kmers *= mem_ratio;
  log_msgs << "Adjusted element counts by " << fixed << setprecision(3) << mem_ratio << ": read kmers " << max_read_kmers << ", qf "
           << max_qf_kmers << ", ctg kmers " << max_ctg_kmers << ", compact ht " << max_compact_kmers << "\n";

  size_t qf_bytes_used = 0;
  if (use_qf) {
    qf_bytes_used = two_choice_filter::estimate_memory(max_qf_kmers);
    if (qf_bytes_used == 0) {
      use_qf = false;
    } else {
      auto sizing_controller = two_choice_filter::get_tcf_sizing_from_mem(qf_bytes_used);
      dstate->tcf = two_choice_filter::TCF::generate_on_device(&sizing_controller, 42);
    }
  }

  // find the first prime number lower than the available slots, and no more than 3x the max number of elements
  primes::Prime prime;
  prime.set(max_read_kmers, false);
  auto ht_capacity = prime.get();
  auto ht_bytes_used = ht_capacity * elem_size;

  log_msgs << "GPU read kmers hash table has capacity per rank of " << ht_capacity << " and uses " << ht_bytes_used << " (QF uses "
           << qf_bytes_used << ")\n";

  // uncomment to debug OOMs
  // cout << "ht bytes used " << (ht_bytes_used / 1024 / 1024) << "MB\n";
  read_kmers_dev.init(ht_capacity);
  // for transferring packed elements from host to gpu
  elem_buff_host.seqs = new char[KCOUNT_GPU_HASHTABLE_BLOCK_SIZE];
  // these are not used for kmers from reads
  elem_buff_host.counts = nullptr;
  // buffer on the device
  ERROR_CHECK(Malloc(&packed_elem_buff_dev.seqs, KCOUNT_GPU_HASHTABLE_BLOCK_SIZE));
  ERROR_CHECK(Malloc(&unpacked_elem_buff_dev.seqs, KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * 2));
  packed_elem_buff_dev.counts = nullptr;
  unpacked_elem_buff_dev.counts = nullptr;

  ERROR_CHECK(Malloc(&gpu_insert_stats, sizeof(InsertStats)));
  ERROR_CHECK(Memset(gpu_insert_stats, 0, sizeof(InsertStats)));

  msgs = log_msgs.str();
  warnings = log_warnings.str();
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::init_ctg_kmers(uint64_t max_elems, size_t gpu_avail_mem) {
  pass_type = CTG_KMERS_PASS;
  // free up space
  if (dstate->tcf) two_choice_filter::TCF::free_on_device(dstate->tcf);
  dstate->tcf = nullptr;

  size_t elem_buff_size = KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * (1 + sizeof(count_t)) * 3;
  size_t elem_size = sizeof(KmerArray<MAX_K>) + sizeof(CountsArray);
  size_t max_slots = 0.97 * (gpu_avail_mem - elem_buff_size) / elem_size;
  primes::Prime prime;
  prime.set(min(max_slots, (size_t)(max_elems * 3)), false);
  auto ht_capacity = prime.get();
  ctg_kmers_dev.init(ht_capacity);
  elem_buff_host.counts = new count_t[KCOUNT_GPU_HASHTABLE_BLOCK_SIZE];
  ERROR_CHECK(Malloc(&packed_elem_buff_dev.counts, KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * sizeof(count_t)));
  ERROR_CHECK(Malloc(&unpacked_elem_buff_dev.counts, 2 * KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * sizeof(count_t)));
  ERROR_CHECK(Memset(gpu_insert_stats, 0, sizeof(InsertStats)));
}

template <int MAX_K>
HashTableGPUDriver<MAX_K>::~HashTableGPUDriver() {
  if (dstate) {
    // this happens when there is no ctg kmers pass
    if (dstate->tcf) two_choice_filter::TCF::free_on_device(dstate->tcf);
    delete dstate;
  }
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::insert_supermer_block() {
  dstate->insert_timer.start();
  bool is_ctg_kmers = (pass_type == CTG_KMERS_PASS);
  ERROR_CHECK(Memcpy(packed_elem_buff_dev.seqs, elem_buff_host.seqs, buff_len, MemcpyHostToDevice));
  ERROR_CHECK(Memset(unpacked_elem_buff_dev.seqs, 0, buff_len * 2));
  if (is_ctg_kmers)
    ERROR_CHECK(Memcpy(packed_elem_buff_dev.counts, elem_buff_host.counts, buff_len * sizeof(count_t), MemcpyHostToDevice));

  int gridsize, threadblocksize;
  dstate->kernel_timer.start();
  get_kernel_config(buff_len, gpu_unpack_supermer_block, gridsize, threadblocksize);
  LaunchKernel(gpu_unpack_supermer_block, gridsize, threadblocksize, unpacked_elem_buff_dev, packed_elem_buff_dev, buff_len);
  get_kernel_config(buff_len * 2, gpu_insert_supermer_block<MAX_K>, gridsize, threadblocksize);
  // gridsize = gridsize * threadblocksize;
  // threadblocksize = 1;
  LaunchKernel(gpu_insert_supermer_block, gridsize, threadblocksize, is_ctg_kmers ? ctg_kmers_dev : read_kmers_dev,
               unpacked_elem_buff_dev, buff_len * 2, kmer_len, is_ctg_kmers, gpu_insert_stats, dstate->tcf);
  // the kernel time is not going to be accurate, because we are not waiting for the kernel to complete
  // need to uncomment the line below, which will decrease performance by preventing the overlap of GPU and CPU execution
  ERROR_CHECK(DeviceSynchronize());
  dstate->kernel_timer.stop();
  num_gpu_calls++;
  dstate->insert_timer.stop();
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::insert_supermer(const string &supermer_seq, count_t supermer_count) {
  if (buff_len + supermer_seq.length() + 1 >= KCOUNT_GPU_HASHTABLE_BLOCK_SIZE) {
    insert_supermer_block();
    buff_len = 0;
  }
  memcpy(&(elem_buff_host.seqs[buff_len]), supermer_seq.c_str(), supermer_seq.length());
  if (pass_type == CTG_KMERS_PASS) {
    for (int i = 0; i < (int)supermer_seq.length(); i++) elem_buff_host.counts[buff_len + i] = supermer_count;
  }
  buff_len += supermer_seq.length();
  elem_buff_host.seqs[buff_len] = '_';
  if (pass_type == CTG_KMERS_PASS) elem_buff_host.counts[buff_len] = 0;
  buff_len++;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::purge_invalid(uint64_t &num_purged, uint64_t &num_entries) {
  num_purged = num_entries = 0;
  uint64_t *counts_gpu;
  int NUM_COUNTS = 2;
  ERROR_CHECK(Malloc(&counts_gpu, NUM_COUNTS * sizeof(uint64_t)));
  ERROR_CHECK(Memset(counts_gpu, 0, NUM_COUNTS * sizeof(uint64_t)));
  GPUTimer t;
  int gridsize, threadblocksize;
  get_kernel_config(read_kmers_dev.capacity, gpu_purge_invalid<MAX_K>, gridsize, threadblocksize);
  t.start();
  // now purge all invalid kmers (do it on the gpu)
  LaunchKernel(gpu_purge_invalid, gridsize, threadblocksize, read_kmers_dev, counts_gpu);
  t.stop();
  dstate->kernel_timer.inc(t.get_elapsed());

  uint64_t counts_host[NUM_COUNTS];
  ERROR_CHECK(Memcpy(&counts_host, counts_gpu, NUM_COUNTS * sizeof(uint64_t), MemcpyDeviceToHost));
  num_purged = counts_host[0];
  num_entries = counts_host[1];
#ifdef DEBUG
  auto expected_num_entries = read_kmers_stats.new_inserts - num_purged;
  if (num_entries != (int)expected_num_entries)
    WARN("mismatch %lu != %lu diff %lu new inserts %lu num purged %lu", num_entries, expected_num_entries,
         (num_entries - (int)expected_num_entries), read_kmers_stats.new_inserts, num_purged);
#endif
  read_kmers_dev.num = num_entries;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::flush_inserts() {
  if (buff_len) {
    insert_supermer_block();
    buff_len = 0;
  }
  ERROR_CHECK(Memcpy(pass_type == READ_KMERS_PASS ? &read_kmers_stats : &ctg_kmers_stats, gpu_insert_stats, sizeof(InsertStats),
                     MemcpyDeviceToHost));
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::done_all_inserts(uint64_t &num_dropped, uint64_t &num_unique, uint64_t &num_purged) {
  uint64_t num_entries = 0;
  purge_invalid(num_purged, num_entries);
  read_kmers_dev.num = num_entries;
  if (elem_buff_host.seqs) delete[] elem_buff_host.seqs;
  if (elem_buff_host.counts) delete[] elem_buff_host.counts;
  ERROR_CHECK(Free(packed_elem_buff_dev.seqs));
  ERROR_CHECK(Free(unpacked_elem_buff_dev.seqs));
  if (packed_elem_buff_dev.counts) ERROR_CHECK(Free(packed_elem_buff_dev.counts));
  if (unpacked_elem_buff_dev.counts) ERROR_CHECK(Free(unpacked_elem_buff_dev.counts));
  ERROR_CHECK(Free(gpu_insert_stats));
  // overallocate to reduce collisions
  num_entries *= 1.3;
  // now compact the hash table entries
  uint64_t *counts_gpu;
  int NUM_COUNTS = 2;
  ERROR_CHECK(Malloc(&counts_gpu, NUM_COUNTS * sizeof(uint64_t)));
  ERROR_CHECK(Memset(counts_gpu, 0, NUM_COUNTS * sizeof(uint64_t)));
  KmerExtsMap<MAX_K> compact_read_kmers_dev;
  compact_read_kmers_dev.init(num_entries);
  GPUTimer t;
  int gridsize, threadblocksize;
  get_kernel_config(read_kmers_dev.capacity, gpu_compact_ht<MAX_K>, gridsize, threadblocksize);
  t.start();
  LaunchKernel(gpu_compact_ht, gridsize, threadblocksize, read_kmers_dev, compact_read_kmers_dev, counts_gpu);
  t.stop();
  dstate->kernel_timer.inc(t.get_elapsed());
  read_kmers_dev.clear();
  uint64_t counts_host[NUM_COUNTS];
  ERROR_CHECK(Memcpy(&counts_host, counts_gpu, NUM_COUNTS * sizeof(uint64_t), MemcpyDeviceToHost));
  ERROR_CHECK(Free(counts_gpu));
  num_dropped = counts_host[0];
  num_unique = counts_host[1];
#ifdef DEBUG
  if (num_unique != read_kmers_dev.num) WARN("mismatch in expected entries %lu != %lu", num_unique, read_kmers_dev.num);
#endif
  // now copy the gpu hash table values across to the host
  // We only do this once, which requires enough memory on the host to store the full GPU hash table, but since the GPU memory
  // is generally a lot less than the host memory, it should be fine.
  output_keys.resize(num_entries);
  output_vals.resize(num_entries);
  output_index = 0;
  ERROR_CHECK(Memcpy(output_keys.data(), (void *)compact_read_kmers_dev.keys,
                     compact_read_kmers_dev.capacity * sizeof(KmerArray<MAX_K>), MemcpyDeviceToHost));
  ERROR_CHECK(Memcpy(output_vals.data(), compact_read_kmers_dev.vals, compact_read_kmers_dev.capacity * sizeof(CountExts),
                     MemcpyDeviceToHost));
  compact_read_kmers_dev.clear();
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::done_ctg_kmer_inserts(uint64_t &attempted_inserts, uint64_t &dropped_inserts, uint64_t &new_inserts) {
  uint64_t *counts_gpu;
  int NUM_COUNTS = 3;
  ERROR_CHECK(Malloc(&counts_gpu, NUM_COUNTS * sizeof(uint64_t)));
  ERROR_CHECK(Memset(counts_gpu, 0, NUM_COUNTS * sizeof(uint64_t)));
  GPUTimer t;
  int gridsize, threadblocksize;
  get_kernel_config(ctg_kmers_dev.capacity, gpu_merge_ctg_kmers<MAX_K>, gridsize, threadblocksize);
  t.start();
  LaunchKernel(gpu_merge_ctg_kmers, gridsize, threadblocksize, read_kmers_dev, ctg_kmers_dev, counts_gpu);
  t.stop();
  dstate->kernel_timer.inc(t.get_elapsed());
  ctg_kmers_dev.clear();
  uint64_t counts_host[NUM_COUNTS];
  ERROR_CHECK(Memcpy(&counts_host, counts_gpu, NUM_COUNTS * sizeof(uint64_t), MemcpyDeviceToHost));
  ERROR_CHECK(Free(counts_gpu));
  attempted_inserts = counts_host[0];
  dropped_inserts = counts_host[1];
  new_inserts = counts_host[2];
  read_kmers_dev.num += new_inserts;
  read_kmers_stats.new_inserts += new_inserts;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::get_elapsed_time(double &insert_time, double &kernel_time) {
  insert_time = dstate->insert_timer.get_elapsed();
  kernel_time = dstate->kernel_timer.get_elapsed();
}

template <int MAX_K>
pair<KmerArray<MAX_K> *, CountExts *> HashTableGPUDriver<MAX_K>::get_next_entry() {
  if (output_keys.empty() || output_index == output_keys.size()) return {nullptr, nullptr};
  output_index++;
  return {&(output_keys[output_index - 1]), &(output_vals[output_index - 1])};
}

template <int MAX_K>
int64_t HashTableGPUDriver<MAX_K>::get_capacity() {
  if (pass_type == READ_KMERS_PASS)
    return read_kmers_dev.capacity;
  else
    return ctg_kmers_dev.capacity;
}

template <int MAX_K>
int64_t HashTableGPUDriver<MAX_K>::get_final_capacity() {
  return read_kmers_dev.capacity;
}

template <int MAX_K>
InsertStats &HashTableGPUDriver<MAX_K>::get_stats() {
  if (pass_type == READ_KMERS_PASS)
    return read_kmers_stats;
  else
    return ctg_kmers_stats;
}

template <int MAX_K>
int HashTableGPUDriver<MAX_K>::get_num_gpu_calls() {
  return num_gpu_calls;
}

template <int MAX_K>
double HashTableGPUDriver<MAX_K>::get_qf_load_factor() {
  if (dstate->tcf) return (double)dstate->tcf->get_fill() / dstate->tcf->get_num_slots();
  return 0;
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
