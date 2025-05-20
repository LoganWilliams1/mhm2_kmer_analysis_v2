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

#include <upcxx/upcxx.hpp>


#include "upcxx_utils/colors.h"
// #include "gpu-utils/gpu_common.hpp"

//#include <Kokkos_Core.hpp>
#include "kokkos_gpu_ht.hpp"
#include "prime.hpp"
#include "kokkos_gpu_hash_funcs.hpp"
#ifdef USE_TCF
#include "tcf_wrapper.hpp"
#else
// two choice filter calls stubbed out
namespace two_choice_filter {
#define TCF_RESULT uint8_t
struct TCF {
//   static TCF *generate_on_device(bool *, int) { return nullptr; }
//   static void free_on_device(TCF *) {}
//   __device__ bool get_my_tile() { return false; }
//   __device__ bool insert_with_delete(bool, uint64_t, uint8_t) { return false; }
//   __device__ bool query(bool, uint64_t, TCF_RESULT &) { return false; }
//   __device__ bool remove(bool, uint64_t) { return false; }
//   int get_fill() { return 0; }
//   int get_num_slots() { return 0; }
};
// __device__ uint8_t pack_extensions(char left, char right) { return 0; }
// __device__ bool unpack_extensions(uint8_t storage, char &left, char &right) { return false; }
// static uint64_t estimate_memory(uint64_t max_num_kmers) { return 0; }
// static bool get_tcf_sizing_from_mem(uint64_t available_bytes) { return false; }

}  // namespace two_choice_filter
#endif

using namespace std;
// using namespace gpu_common;
using namespace kcount_gpu;

// convenience functions
#define SDBG(fmt, ...) \
  if (!upcxx_rank_me) printf(KLMAGENTA "GPU kcount: " fmt KNORM "\n", ##__VA_ARGS__)

#define SWARN(fmt, ...) \
  if (!upcxx_rank_me) printf(KLRED "WARN GPU kcount %d: " fmt KNORM "\n", __LINE__, ##__VA_ARGS__)

#define WARN(fmt, ...) printf(KLRED "WARN GPU kcount %d:" fmt KNORM "\n", __LINE__, ##__VA_ARGS__)

template <int MAX_K>
KOKKOS_FUNCTION void kmer_set(KmerArray<MAX_K> &kmer1, const KmerArray<MAX_K> &kmer2) {
  const uint64_t KEY_EMPTY = 0xffffffffffffffff;
  const uint64_t KEY_TRANSITION = 0xfffffffffffffffe;  
  int N_LONGS = kmer1.N_LONGS;
  uint64_t old_key;
  for (int i = 0; i < N_LONGS - 1; i++) {
    old_key = Kokkos::atomic_exchange((unsigned long long *)&(kmer1.longs[i]), kmer2.longs[i]);
    //if (old_key != KEY_EMPTY) WARN("old key should be KEY_EMPTY");
    if (old_key != KEY_EMPTY) Kokkos::printf("old key should be KEY_EMPTY");
  }
  old_key = Kokkos::atomic_exchange((unsigned long long *)&(kmer1.longs[N_LONGS - 1]), kmer2.longs[N_LONGS - 1]);
  //if (old_key != KEY_TRANSITION) WARN("old key should be KEY_TRANSITION");
  if (old_key != KEY_TRANSITION) Kokkos::printf("old key should be KEY_TRANSITION");
}

template <int MAX_K>
KOKKOS_FUNCTION bool kmers_equal(const KmerArray<MAX_K> &kmer1, const KmerArray<MAX_K> &kmer2) {
  int n_longs = kmer1.N_LONGS;
  for (int i = 0; i < n_longs; i++) {
    uint64_t old_key = Kokkos::atomic_fetch_add((unsigned long long *)&(kmer1.longs[i]), 0ULL);
    if (old_key != kmer2.longs[i]) return false;
  }
  return true;
}

template <int MAX_K>
KOKKOS_FUNCTION size_t kmer_hash(const KmerArray<MAX_K> &kmer) {
  return gpu_murmurhash3_64(reinterpret_cast<const void *>(kmer.longs), kmer.N_LONGS * sizeof(uint64_t));
}

KOKKOS_FUNCTION int8_t get_ext(CountsArray &counts, int pos, int8_t *ext_map) {
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

KOKKOS_FUNCTION bool ext_conflict(ext_count_t *ext_counts, int start_idx) {
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
void gpu_merge_ctg_kmers(KmerCountsMap<MAX_K> read_kmers, const KmerCountsMap<MAX_K> ctg_kmers,
                                    uint64_t *insert_counts) {
  // unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
  // int8_t ext_map[4] = {'A', 'C', 'G', 'T'};
  // int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  // uint64_t attempted_inserts = 0;
  // uint64_t dropped_inserts = 0;
  // uint64_t new_inserts = 0;
  // if (threadid < ctg_kmers.capacity) {
  //   count_t kmer_count = ctg_kmers.vals[threadid].kmer_count;
  //   ext_count_t *ext_counts = ctg_kmers.vals[threadid].ext_counts;
  //   if (kmer_count && !ext_conflict(ext_counts, 0) && !ext_conflict(ext_counts, 4)) {
  //     KmerArray<MAX_K> kmer = ctg_kmers.keys[threadid];
  //     uint64_t slot = kmer_hash(kmer) % read_kmers.capacity;
  //     auto start_slot = slot;
  //     attempted_inserts++;
  //     const int MAX_PROBE = (read_kmers.capacity < KCOUNT_HT_MAX_PROBE ? read_kmers.capacity : KCOUNT_HT_MAX_PROBE);
  //     for (int j = 0; j < MAX_PROBE; j++) {
  //       uint64_t old_key = atomicCAS((unsigned long long *)&(read_kmers.keys[slot].longs[N_LONGS - 1]), KEY_EMPTY, KEY_TRANSITION);
  //       if (old_key == KEY_EMPTY) {
  //         new_inserts++;
  //         memcpy(&read_kmers.vals[slot], &ctg_kmers.vals[threadid], sizeof(CountsArray));
  //         kmer_set(read_kmers.keys[slot], kmer);
  //         break;
  //       } else if (old_key == kmer.longs[N_LONGS - 1]) {
  //         if (kmers_equal(read_kmers.keys[slot], kmer)) {
  //           // existing kmer from reads - only replace if the kmer is non-UU
  //           // there is no need for atomics here because all ctg kmers are unique; hence only one thread will ever match this kmer
  //           int8_t left_ext = get_ext(read_kmers.vals[slot], 0, ext_map);
  //           int8_t right_ext = get_ext(read_kmers.vals[slot], 4, ext_map);
  //           if (left_ext == 'X' || left_ext == 'F' || right_ext == 'X' || right_ext == 'F')
  //             memcpy(&read_kmers.vals[slot], &ctg_kmers.vals[threadid], sizeof(CountsArray));
  //           break;
  //         }
  //       }
  //       // quadratic probing - worse cache but reduced clustering
  //       slot = (start_slot + (j + 1) * (j + 1)) % read_kmers.capacity;
  //       if (j == MAX_PROBE - 1) dropped_inserts++;
  //     }
  //   }
  // }
  // reduce(attempted_inserts, ctg_kmers.capacity, &(insert_counts[0]));
  // reduce(dropped_inserts, ctg_kmers.capacity, &(insert_counts[1]));
  // reduce(new_inserts, ctg_kmers.capacity, &(insert_counts[2]));
}

template <int MAX_K>
void gpu_compact_ht(KmerCountsMap<MAX_K> elems, KmerExtsMap<MAX_K> compact_elems, uint64_t& num_purged, uint64_t& num_entries) {
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  uint64_t dropped_inserts = 0;
  uint64_t unique_inserts = 0;
  const uint64_t KEY_EMPTY = 0xffffffffffffffff;

  Kokkos::View<KmerArray<MAX_K>*> keys_view = elems.keys_v;

  Kokkos::View<CountsArray*> vals_view = elems.vals_v;

  Kokkos::View<KmerArray<MAX_K>*> compact_keys_view = compact_elems.keys_v;

  Kokkos::View<CountExts*> compact_vals_view = compact_elems.vals_v;    

  Kokkos::parallel_reduce("gpu_compact_ht", elems.capacity, KOKKOS_LAMBDA (int i, uint64_t& dropped_inserts, uint64_t& unique_inserts) {
    if (vals_view(i).kmer_count) {
      int8_t ext_map[4] = {'A', 'C', 'G', 'T'};
      KmerArray<MAX_K> kmer = keys_view(i);
      uint64_t slot = kmer_hash(kmer) % compact_elems.capacity;
      auto start_slot = slot;
      // we set a constraint on the max probe to track whether we are getting excessive collisions and need a bigger default
      // compact table      
      const int MAX_PROBE = (compact_elems.capacity < KCOUNT_HT_MAX_PROBE ? compact_elems.capacity : KCOUNT_HT_MAX_PROBE);
      // look for empty slot in compact hash table
      for (int j = 0; j < MAX_PROBE; j++) {
        uint64_t old_key = 
            Kokkos::atomic_compare_exchange((unsigned long long *)&(compact_keys_view(slot).longs[N_LONGS - 1]), KEY_EMPTY, kmer.longs[N_LONGS - 1]);
        if (old_key == KEY_EMPTY) {
          // found empty slot - there will be no duplicate keys since we're copying across from another hash table
          unique_inserts++;
          for (int k = 0; k < N_LONGS - 1; k++) {
            compact_keys_view(slot).longs[k] = kmer.longs[k];
          }
          // compute exts
          int8_t left_ext = get_ext(vals_view(i), 0, ext_map);
          int8_t right_ext = get_ext(vals_view(i), 4, ext_map);
          //if (vals_view(i).kmer_count < 2) WARN("elem should have been purged, count %d", vals_view(i).kmer_count);
          if (vals_view(i).kmer_count < 2) Kokkos::printf("elem should have been purged, count %d", vals_view(i).kmer_count);
          compact_vals_view(slot).count = vals_view(i).kmer_count;
          compact_vals_view(slot).left = left_ext;
          compact_vals_view(slot).right = right_ext;
          break;
        }
        // quadratic probing - worse cache but reduced clustering
        slot = (start_slot + (j + 1) * (j+ 1)) % compact_elems.capacity;
        if (j == MAX_PROBE - 1) dropped_inserts++;
      }
    }
  }, num_purged, num_entries);
}

template <int MAX_K>
void gpu_purge_invalid(KmerCountsMap<MAX_K> elems, uint64_t& num_purged, uint64_t& num_entries) {
  int N_LONGS = KmerArray<MAX_K>::N_LONGS;

  const uint64_t KEY_EMPTY = 0xffffffffffffffff;

  Kokkos::View<KmerArray<MAX_K>*> keys_view = elems.keys_v;

  Kokkos::View<CountsArray*> vals_view = elems.vals_v;

  uint64_t capacity = elems.capacity;

  Kokkos::parallel_reduce("gpu_purge_invalid", capacity, KOKKOS_LAMBDA (int i, uint64_t& local_purged, uint64_t& local_elems) {
    if (vals_view(i).kmer_count) {
      uint64_t ext_sum = 0;
      for (int j = 0; j < 8; j++) ext_sum += vals_view(i).ext_counts[j];
      if (vals_view(i).kmer_count < 2 || !ext_sum) {
        vals_view(i).kmer_count = 0;
        for (int j = 0; j < 8; j++) vals_view(i).ext_counts[j] = 0;
        for (int j = 0; j < N_LONGS; j++) keys_view(i).longs[j] = KEY_EMPTY;
        local_purged++;
      } else {
        local_elems++;
      }
    }
  }, num_purged, num_entries);
}


KOKKOS_INLINE_FUNCTION char to_base_func(int index, int pp) {
  if (index > 9) {
    //WARN("index out of range for to_base: %d, packed seq pos %d", index, pp);
    Kokkos::printf("index out of range for to_base: %d, packed seq pos %d", index, pp);
    return 0;
  }
  if (index == 0) return '_';
  switch (index) {
    case 1: return 'a';
    case 2: return 'c';
    case 3: return 'g';
    case 4: return 't';
    case 5: return 'A';
    case 6: return 'C';
    case 7: return 'G';
    case 8: return 'T';
    case 9: return 'N';
    default: return '0';
  }
  return '0';

}

void gpu_unpack_supermer_block(SupermerBuff unpacked_supermer_buff, SupermerBuff packed_supermer_buff, int buff_len) {

  Kokkos::parallel_for("gpu_unpack_supermer_block", buff_len - 1, KOKKOS_LAMBDA (int i) {
    uint8_t packed = packed_supermer_buff.seqs_v(i);
    if (packed == '_') return;
    uint8_t left_side = (packed & 240) >> 4;
    unpacked_supermer_buff.seqs_v(i * 2) = to_base_func(left_side, packed);
    if (packed_supermer_buff.counts_v.size() > 1) unpacked_supermer_buff.counts_v(i * 2) = packed_supermer_buff.counts_v(i);
    uint8_t right_side = packed & 15;
    unpacked_supermer_buff.seqs_v(i * 2 + 1) = to_base_func(right_side, packed);
    if (packed_supermer_buff.counts_v.size() > 1) unpacked_supermer_buff.counts_v(i * 2 + 1) = packed_supermer_buff.counts_v(i);
  });  

   Kokkos::fence();

}

KOKKOS_INLINE_FUNCTION bool is_valid_base(char base) {
  return (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == '0' || base == 'N');
}

KOKKOS_INLINE_FUNCTION bool bad_qual(char base) { return (base == 'a' || base == 'c' || base == 'g' || base == 't'); }

KOKKOS_INLINE_FUNCTION uint16_t atomicAddUint16(uint16_t *address, uint16_t val) {
  unsigned int *base_address = (unsigned int *)((size_t)address & ~2);
  unsigned int long_val = ((size_t)address & 2) ? ((unsigned int)val << 16) : val;
  unsigned int long_old = Kokkos::atomic_fetch_add(base_address, long_val);
  return ((size_t)address & 2) ? (uint16_t)(long_old >> 16) : (uint16_t)(long_old & 0xffff);
}

KOKKOS_INLINE_FUNCTION void atomicAddUint16_thres(uint16_t *address, uint16_t val, uint16_t thres) {
  if (atomicAddUint16(address, 0) < thres - val) atomicAddUint16(address, val);
}

KOKKOS_INLINE_FUNCTION void inc_ext(char ext, ext_count_t kmer_count, ext_count_t *ext_counts) {
  switch (ext) {
    case 'A': atomicAddUint16_thres(&(ext_counts[0]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
    case 'C': atomicAddUint16_thres(&(ext_counts[1]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
    case 'G': atomicAddUint16_thres(&(ext_counts[2]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
    case 'T': atomicAddUint16_thres(&(ext_counts[3]), kmer_count, KCOUNT_MAX_KMER_COUNT); return;
  }
}

KOKKOS_INLINE_FUNCTION bool pack_seq_to_kmer(char *seqs, int kmer_len, int num_longs, uint64_t *kmer) {
  int l = 0, prev_l = 0;
  uint64_t longs = 0;
  memset(kmer, 0, sizeof(uint64_t) * num_longs);
  // each thread extracts one kmer
  for (int k = 0; k < kmer_len; k++) {
    char s = seqs[k];
    switch (s) {
      case 'a': s = 'A'; break;
      case 'c': s = 'C'; break;
      case 'g': s = 'G'; break;
      case 't': s = 'T'; break;
      case 'A':
      case 'C':
      case 'G':
      case 'T': break;
      default: return false;
    }
    int j = k % 32;
    prev_l = l;
    l = k / 32;
    // we do it this way so we can operate on the variable longs in a register, rather than local memory in the array
    if (l > prev_l) {
      kmer[prev_l] = longs;
      longs = 0;
      prev_l = l;
    }
    uint64_t x = (s & 4) >> 1;
    longs |= ((x + ((x ^ (s & 2)) >> 1)) << (2 * (31 - j)));
  }
  kmer[l] = longs;
  return true;
}

const uint64_t KOKKOS_TWINS[256] = {
    0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F, 0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F, 0xFB, 0xBB, 0x7B, 0x3B,
    0xEB, 0xAB, 0x6B, 0x2B, 0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B, 0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
    0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07, 0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23, 0xD3, 0x93, 0x53, 0x13,
    0xC3, 0x83, 0x43, 0x03, 0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E, 0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
    0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A, 0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A, 0xF6, 0xB6, 0x76, 0x36,
    0xE6, 0xA6, 0x66, 0x26, 0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06, 0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
    0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02, 0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D, 0xDD, 0x9D, 0x5D, 0x1D,
    0xCD, 0x8D, 0x4D, 0x0D, 0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29, 0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
    0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25, 0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05, 0xF1, 0xB1, 0x71, 0x31,
    0xE1, 0xA1, 0x61, 0x21, 0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01, 0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
    0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C, 0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28, 0xD8, 0x98, 0x58, 0x18,
    0xC8, 0x88, 0x48, 0x08, 0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24, 0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
    0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00};

KOKKOS_INLINE_FUNCTION void revcomp(uint64_t *longs, uint64_t *rc_longs, int kmer_len, int num_longs, Kokkos::View<uint64_t[256]> twins_v) {
  int last_long = (kmer_len + 31) / 32;
  for (size_t i = 0; i < last_long; i++) {
    uint64_t v = longs[i];
    rc_longs[last_long - 1 - i] = (twins_v(v & 0xFF) << 56) | (twins_v((v >> 8) & 0xFF) << 48) |
                                  (twins_v((v >> 16) & 0xFF) << 40) | (twins_v((v >> 24) & 0xFF) << 32) |
                                  (twins_v((v >> 32) & 0xFF) << 24) | (twins_v((v >> 40) & 0xFF) << 16) |
                                  (twins_v((v >> 48) & 0xFF) << 8) | (twins_v((v >> 56)));
  }
  uint64_t shift = (kmer_len % 32) ? 2 * (32 - (kmer_len % 32)) : 0;
  uint64_t shiftmask = (kmer_len % 32) ? (((((uint64_t)1) << shift) - 1) << (64 - shift)) : ((uint64_t)0);
  rc_longs[0] = rc_longs[0] << shift;
  for (size_t i = 1; i < last_long; i++) {
    rc_longs[i - 1] |= (rc_longs[i] & shiftmask) >> (64 - shift);
    rc_longs[i] = rc_longs[i] << shift;
  }
}

KOKKOS_INLINE_FUNCTION char comp_nucleotide(char ch) {
  switch (ch) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    case '0': return '0';
    case 'U':
    case 'R':
    case 'Y':
    case 'K':
    case 'M':
    case 'S':
    case 'W':
    case 'B':
    case 'D':
    case 'H':
    case 'V': return 'N';
    default: return 0;
  }
}

template <int MAX_K>
KOKKOS_FUNCTION bool get_kmer_from_supermer(Kokkos::View<char*> seqs_view, Kokkos::View<count_t*> counts_view, uint32_t buff_len, int kmer_len, uint64_t *kmer, char &left_ext,
                                       char &right_ext, count_t &count, unsigned int kokkos_index, Kokkos::View<uint64_t[256]> twins_v) {


  int num_kmers = buff_len - kmer_len + 1;
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  if (!pack_seq_to_kmer(&(seqs_view(kokkos_index)), kmer_len, N_LONGS, kmer)) return false;
  if (kokkos_index + kmer_len >= buff_len) return false;
  left_ext = seqs_view(kokkos_index - 1);
  right_ext = seqs_view(kokkos_index + kmer_len);
  if (left_ext == '_' || right_ext == '_') return false;
  if (!left_ext || !right_ext) return false;
  if (counts_view.size() > 1) {
    count = counts_view(kokkos_index);
  } else {
    count = 1;
    if (bad_qual(left_ext)) left_ext = '0';
    if (bad_qual(right_ext)) right_ext = '0';
  }
  if (!is_valid_base(left_ext)) {
    //WARN("threadid %d, invalid char for left nucleotide %d", kokkos_index, (uint8_t)left_ext);
    Kokkos::printf("threadid %d, invalid char for left nucleotide %d", kokkos_index, (uint8_t)left_ext);
    return false;
  }
  if (!is_valid_base(right_ext)) {
    //WARN("threadid %d, invalid char for right nucleotide %d", kokkos_index, (uint8_t)right_ext);
    Kokkos::printf("threadid %d, invalid char for right nucleotide %d", kokkos_index, (uint8_t)right_ext);
    return false;
  }
  uint64_t kmer_rc[N_LONGS];
  revcomp(kmer, kmer_rc, kmer_len, N_LONGS, twins_v);
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
KOKKOS_FUNCTION bool gpu_insert_kmer(Kokkos::View<KmerArray<MAX_K>*> keys_view, Kokkos::View<CountsArray*> vals_view, uint64_t capacity, uint64_t hash_val, KmerArray<MAX_K> &kmer, char left_ext,
                                char right_ext, char prev_left_ext, char prev_right_ext, count_t kmer_count, uint64_t &new_inserts,
                                uint64_t &dropped_inserts, bool ctg_kmers, bool use_qf, bool update_only) {
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;
  const uint64_t KEY_EMPTY = 0xffffffffffffffff;
  const uint64_t KEY_TRANSITION = 0xfffffffffffffffe;
  uint64_t slot = hash_val % capacity;
  auto start_slot = slot;
  const int MAX_PROBE = (capacity < 200 ? capacity : 200);
  bool found_slot = false;
  bool kmer_found_in_ht = false;
  uint64_t old_key = KEY_TRANSITION;
  for (int j = 0; j < MAX_PROBE; j++) {
    // we have to be careful here not to end up with multiple threads on the same warp accessing the same slot, because
    // that will cause a deadlock. So we loop over all statements in each CAS spin to ensure that all threads get a
    // chance to execute
    do {
      old_key = Kokkos::atomic_compare_exchange((unsigned long long *)&(keys_view(slot).longs[N_LONGS - 1]), KEY_EMPTY, KEY_TRANSITION);
      if (old_key != KEY_TRANSITION) {
        if (old_key == KEY_EMPTY) {
          if (update_only) {
            old_key = Kokkos::atomic_exchange((unsigned long long *)&(keys_view(slot).longs[N_LONGS - 1]), KEY_EMPTY);
            //if (old_key != KEY_TRANSITION) WARN("old key should be KEY_TRANSITION");
            if (old_key != KEY_TRANSITION) Kokkos::printf("old key should be KEY_TRANSITION");
            return false;
          }
          kmer_set(keys_view(slot), kmer);
          found_slot = true;
        } else if (old_key == kmer.longs[N_LONGS - 1]) {
          if (kmers_equal(keys_view(slot), kmer)) {
            found_slot = true;
            kmer_found_in_ht = true;
          }
        }
      }
    } while (old_key == KEY_TRANSITION);
    if (found_slot) break;
    // quadratic probing - worse cache but reduced clustering
    slot = (start_slot + j * j) % capacity;
    // this entry didn't get inserted because we ran out of probing time (and probably space)
    if (j == MAX_PROBE - 1) dropped_inserts++;
  }
  if (found_slot) {
    ext_count_t *ext_counts = vals_view(slot).ext_counts;
    if (ctg_kmers) {
      // the count is the min of all counts. Use CAS to deal with the initial zero value
      int prev_count = Kokkos::atomic_compare_exchange(&vals_view(slot).kmer_count, 0, kmer_count);
      if (prev_count)
        Kokkos::atomic_min(&vals_view(slot).kmer_count, kmer_count);
      else
        new_inserts++;
    } else {
      assert(kmer_count == 1);
      int prev_count = Kokkos::atomic_fetch_add(&vals_view(slot).kmer_count, kmer_count);
      if (!prev_count) new_inserts++;
    }
    //FIXME: error: no matching function for call to 'min' w/ IntelOneAPI on GR CPU
    ext_count_t kmer_count_uint16 = min(kmer_count, static_cast<uint32_t>(UINT16_MAX));
    inc_ext(left_ext, kmer_count_uint16, ext_counts);
    inc_ext(right_ext, kmer_count_uint16, ext_counts + 4);
    if (use_qf && !update_only && !kmer_found_in_ht && !ctg_kmers) {
      // kmer was not in hash table, so it must have been found in the qf
      // add the extensions from the previous entry stored in the qf
      inc_ext(prev_left_ext, 1, ext_counts);
      inc_ext(prev_right_ext, 1, ext_counts + 4);
      // inc the overall kmer count
      Kokkos::atomic_add(&vals_view(slot).kmer_count, 1);
    }
  }
  return true;
}

template <int MAX_K>
void gpu_insert_supermer_block(KmerCountsMap<MAX_K> elems, SupermerBuff supermer_buff, uint32_t buff_len, int kmer_len, bool ctg_kmers, InsertStats& insert_stats, two_choice_filter::TCF *tcf, Kokkos::View<uint64_t[256]> twins_v) {

  const uint64_t KEY_EMPTY = 0xffffffffffffffff;
  const uint64_t KEY_TRANSITION = 0xfffffffffffffffe;
  const int N_LONGS = KmerArray<MAX_K>::N_LONGS;

  uint64_t att, dr, dr_qf, new_ins, u_qf;

  Kokkos::View<char*> seqs_view = supermer_buff.seqs_v;
  Kokkos::View<count_t*> counts_view = supermer_buff.counts_v;

  Kokkos::View<KmerArray<MAX_K>*> keys_view = elems.keys_v;
  Kokkos::View<CountsArray*> vals_view = elems.vals_v;
  uint64_t capacity = elems.capacity;

  Kokkos::parallel_reduce("gpu_insert_supermer_block", buff_len, KOKKOS_LAMBDA(const int& kokkos_index, uint64_t& attempted_inserts, uint64_t& dropped_inserts, uint64_t& dropped_qf, uint64_t& num_new_inserts, uint64_t& unique_qf) {
    if (kokkos_index > 0) {
      attempted_inserts++;
      KmerArray<MAX_K> kmer;
      char left_ext, right_ext;
      uint32_t kmer_count;
      if (get_kmer_from_supermer<MAX_K>(seqs_view, counts_view, buff_len, kmer_len, kmer.longs, left_ext, right_ext, kmer_count, kokkos_index, twins_v)) {
        //if (kmer.longs[N_LONGS - 1] == KEY_EMPTY) WARN("block equal to KEY_EMPTY");
        if (kmer.longs[N_LONGS - 1] == KEY_EMPTY) Kokkos::printf("block equal to KEY_EMPTY");
        //if (kmer.longs[N_LONGS - 1] == KEY_TRANSITION) WARN("block equal to KEY_TRANSITION");
        if (kmer.longs[N_LONGS - 1] == KEY_TRANSITION) Kokkos::printf("block equal to KEY_TRANSITION");
        auto hash_val = kmer_hash(kmer);
        char prev_left_ext = '0', prev_right_ext = '0';
        // bool use_qf = (tcf != nullptr);
        bool use_qf = false;
        bool update_only = (use_qf && !ctg_kmers);
        bool updated = gpu_insert_kmer(keys_view, vals_view, capacity, hash_val, kmer, left_ext, right_ext, prev_left_ext, prev_right_ext, kmer_count,
                                      num_new_inserts, dropped_inserts, ctg_kmers, use_qf, update_only);
    //   //   // if (update_only && !updated) {
    //   //   //   auto packed = two_choice_filter::pack_extensions(left_ext, right_ext);
    //   //   //   TCF_RESULT result = 0;
    //   //   //   if (tcf->query(tcf->get_my_tile(), hash_val, result)) {
    //   //   //     // found successfully
    //   //   //     tcf->remove(tcf->get_my_tile(), hash_val);
    //   //   //     two_choice_filter::unpack_extensions(result, prev_left_ext, prev_right_ext);
    //   //   //     gpu_insert_kmer(elems, hash_val, kmer, left_ext, right_ext, prev_left_ext, prev_right_ext, kmer_count, new_inserts,
    //   //   //                     dropped_inserts, ctg_kmers, use_qf, false);
    //   //   //   } else {
    //   //   //     if (tcf->insert_with_delete(tcf->get_my_tile(), hash_val, packed)) {
    //   //   //       // inserted successfully
    //   //   //       num_unique_qf++;
    //   //   //     } else {
    //   //   //       // dropped
    //   //   //       dropped_inserts_qf++;
    //   //   //       // now insert it into the main hash table - this will be purged later if it's a singleton
    //   //   //       gpu_insert_kmer(elems, hash_val, kmer, left_ext, right_ext, prev_left_ext, prev_right_ext, kmer_count, new_inserts,
    //   //   //                       dropped_inserts, ctg_kmers, false, false);
    //   //   //     }
    //   //   //   }
    //   //   // }

      }
    }


  }, att, dr, dr_qf, new_ins, u_qf);

  insert_stats.attempted = att;
  insert_stats.dropped = dr;
  insert_stats.dropped_qf = dr_qf;
  insert_stats.new_inserts = new_ins;
  insert_stats.num_unique_qf = u_qf;

}

template <int MAX_K>
struct HashTableGPUDriver<MAX_K>::HashTableDriverState {
  two_choice_filter::TCF *tcf = nullptr;
};


// Function to create & initialize Kokkos keys & values Views from KmerArray & CountsArray structs
template <int MAX_K>
void KmerCountsMap<MAX_K>::init(int64_t ht_capacity) {
  capacity = ht_capacity;
  
  printf("\n\n ... DEBUG_0 - entering KmerCountsMap<MAX_K>::init function ... \n\n");
  // ht_capacity TIOGA:  1 248 138 943 
  printf("\n\n DEBUG_1 - ht_capacity: %d\n\n", ht_capacity);
  // MAX_K TIOGA:  32
  printf("\n\n DEBUG_2 - MAX_K:  %d\n\n", MAX_K);
   
   // Create 1D Kokkos View containing keys view with sized as "capacity"
   keys_v = Kokkos::View<KmerArray<MAX_K>*>("CountsMap keys from capacity", capacity);
   
   // create 1D host space view of keys (suggested by Jan)
   Kokkos::View<KmerArray<MAX_K>*, Kokkos::HostSpace> h_keys_v("CountsMap keys - HostSpace", capacity);

  // initialize to hexadecimal integer
  for (int64_t i = 0; i < capacity; i++) {                 // rows
    for (int j = 0; j < KmerArray<MAX_K>::N_LONGS; j++) {  // cols
      // FIXME: square brackets for "j" cols indexing variable needed b/c range specifier, "KmerArray<MAX_K>::N_LONGS" \
      // is not a pointer 
      h_keys_v(i).longs[j] = 0xffffffffffffffff;
    }
  }
  
  // Copy h_keys_v vals to keys_v
  Kokkos::deep_copy(keys_v, h_keys_v);

  vals_v = Kokkos::View<CountsArray*>("CountsMap vals", capacity);
  
  Kokkos::Profiling::popRegion();

}


template <int MAX_K>
void KmerCountsMap<MAX_K>::clear() {
  // unsure if this will be needed
}

template <int MAX_K>
void KmerExtsMap<MAX_K>::init(int64_t ht_capacity) {
  const uint64_t KEY_EMPTY = 0xffffffffffffffff;
  capacity = ht_capacity;

  keys_v = Kokkos::View<KmerArray<MAX_K>*>("ExtsMap keys", capacity);
  typename Kokkos::View<KmerArray<MAX_K>*>::HostMirror h_keys_v = Kokkos::create_mirror_view(keys_v);
  for (int64_t i = 0; i < capacity; i++) {
    for (int j = 0; j < KmerArray<MAX_K>::N_LONGS; j++) {
      h_keys_v(i).longs[j] = 0xffffffffffffffff;
    }
  }
  Kokkos::deep_copy(keys_v, h_keys_v);  

  vals_v = Kokkos::View<CountExts*>("ExtsMap vals", capacity);
}

template <int MAX_K>
void KmerExtsMap<MAX_K>::clear() {
  // unsure is this will be needed
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
  // gpu_utils::set_gpu_device(upcxx_rank_me);
  dstate = new HashTableDriverState();

  // reserve space for the fixed size buffer for passing data to the GPU
  size_t elem_buff_size = KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * (3 + sizeof(count_t));
  if (gpu_avail_mem > 0) {
    gpu_avail_mem -= elem_buff_size;
  }
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

  // size_t qf_size = use_qf ? two_choice_filter::estimate_memory(max(1.0, log2(max_qf_kmers))) : 0;
  size_t qf_size = 0;

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
  double mem_ratio = (gpu_avail_mem > 0) ? (double)(0.8 * gpu_avail_mem) / tot_size : 1;
  if (mem_ratio < 0.9)
    log_warnings << "Insufficent memory for " << fixed << setprecision(3) << target_load_factor
                 << " load factor across all data structures; reducing to " << (mem_ratio * target_load_factor)
                 << "; this could result in an OOM or dropped kmers";
  max_read_kmers *= mem_ratio;
  max_qf_kmers *= mem_ratio;
  max_ctg_kmers *= mem_ratio;
  max_compact_kmers *= mem_ratio;
  log_msgs << "Adjusted element counts by " << fixed << setprecision(3) << mem_ratio << ": read kmers " << max_read_kmers << ", qf "
           << max_qf_kmers << ", ctg kmers " << max_ctg_kmers << ", compact ht " << max_compact_kmers << "\n";

  size_t qf_bytes_used = 0;
  if (use_qf) {
    // qf_bytes_used = two_choice_filter::estimate_memory(max_qf_kmers);
    // if (qf_bytes_used == 0) {
    //   use_qf = false;
    // } else {
    //   auto sizing_controller = two_choice_filter::get_tcf_sizing_from_mem(qf_bytes_used);
    //   dstate->tcf = two_choice_filter::TCF::generate_on_device(&sizing_controller, 42);
    // }
  }

  // find the first prime number lower than the available slots, and no more than 3x the max number of elements
  primes::Prime prime;
  prime.set(max_read_kmers, false);
  // JAN - max_read_kmers is 1118881598
  printf("JAN - max_read_kmers is %zu\n\n", max_read_kmers);
  auto ht_capacity = prime.get();
  printf("JAN - ht_capacity is %lu\n\n", ht_capacity);

  auto ht_bytes_used = ht_capacity * elem_size;

  log_msgs << "GPU read kmers hash table has capacity per rank of " << ht_capacity << " and uses " << ht_bytes_used << " (QF uses "
           << qf_bytes_used << ")\n";

  // uncomment to debug OOMs
  // cout << "ht bytes used " << (ht_bytes_used / 1024 / 1024) << "MB\n";

  read_kmers_dev.init(ht_capacity);

  // buffer on the device

  packed_elem_buff_dev.seqs_v = Kokkos::View<char*>("packed_elem_buff_dev_seqs", KCOUNT_GPU_HASHTABLE_BLOCK_SIZE);
  unpacked_elem_buff_dev.seqs_v = Kokkos::View<char*>("unpacked_elem_buff_dev_seqs", KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * 2);

  packed_elem_buff_dev.counts_v = Kokkos::View<count_t*>("packed_elem_buff_dev_counts");
  unpacked_elem_buff_dev.counts_v = Kokkos::View<count_t*>("unpacked_elem_buff_dev_counts");

  // for transferring packed elements from host to gpu
  elem_buff_host.h_seqs_v = Kokkos::create_mirror_view(packed_elem_buff_dev.seqs_v);
  // these are not used for kmers from reads
  elem_buff_host.h_counts_v = Kokkos::create_mirror_view(packed_elem_buff_dev.counts_v);   

  msgs = log_msgs.str();
  warnings = log_warnings.str();

  twins_v = Kokkos::View<uint64_t[256]>("twins_array");
  Kokkos::View<uint64_t[256]>::HostMirror h_twins_v = Kokkos::create_mirror_view(twins_v);
  for (int i = 0; i < 256; i++) {
    h_twins_v(i) = KOKKOS_TWINS[i];
  }
  Kokkos::deep_copy(twins_v, h_twins_v);  

}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::init_ctg_kmers(uint64_t max_elems, size_t gpu_avail_mem) {
  // pass_type = CTG_KMERS_PASS;
  // // free up space
  // if (dstate->tcf) two_choice_filter::TCF::free_on_device(dstate->tcf);
  // dstate->tcf = nullptr;

  // size_t elem_buff_size = KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * (1 + sizeof(count_t)) * 3;
  // size_t elem_size = sizeof(KmerArray<MAX_K>) + sizeof(CountsArray);
  // size_t max_slots = 0.97 * (gpu_avail_mem - elem_buff_size) / elem_size;
  // primes::Prime prime;
  // prime.set(min(max_slots, (size_t)(max_elems * 3)), false);
  // auto ht_capacity = prime.get();
  // ctg_kmers_dev.init(ht_capacity);
  // elem_buff_host.counts = new count_t[KCOUNT_GPU_HASHTABLE_BLOCK_SIZE];
  // ERROR_CHECK(Malloc(&packed_elem_buff_dev.counts, KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * sizeof(count_t)));
  // ERROR_CHECK(Malloc(&unpacked_elem_buff_dev.counts, 2 * KCOUNT_GPU_HASHTABLE_BLOCK_SIZE * sizeof(count_t)));
  // ERROR_CHECK(Memset(gpu_insert_stats, 0, sizeof(InsertStats)));
}

template <int MAX_K>
HashTableGPUDriver<MAX_K>::~HashTableGPUDriver() {
  if (dstate) {
    // this happens when there is no ctg kmers pass
    // if (dstate->tcf) two_choice_filter::TCF::free_on_device(dstate->tcf);
    delete dstate;
  }
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::insert_supermer_block() {
  bool is_ctg_kmers = (pass_type == CTG_KMERS_PASS);

  Kokkos::deep_copy(packed_elem_buff_dev.seqs_v, elem_buff_host.h_seqs_v);

  Kokkos::deep_copy(unpacked_elem_buff_dev.seqs_v, 0);
  // if (is_ctg_kmers)
  //   ERROR_CHECK(Memcpy(packed_elem_buff_dev.counts, elem_buff_host.counts, buff_len * sizeof(count_t), MemcpyHostToDevice));


  gpu_unpack_supermer_block(unpacked_elem_buff_dev, packed_elem_buff_dev, buff_len);


  const uint32_t unpacked_buff_len = buff_len * 2;

  // gpu_insert_supermer_block(is_ctg_kmers ? ctg_kmers_dev : read_kmers_dev, unpacked_elem_buff_dev, buff_len * 2, kmer_len, is_ctg_kmers, gpu_insert_stats_v, dstate->tcf, twins_v);
  gpu_insert_supermer_block(read_kmers_dev, unpacked_elem_buff_dev, unpacked_buff_len, kmer_len, is_ctg_kmers, read_kmers_stats, dstate->tcf, twins_v);

  num_gpu_calls++;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::insert_supermer(const string &supermer_seq, count_t supermer_count) {
  if (buff_len + supermer_seq.length() + 1 >= KCOUNT_GPU_HASHTABLE_BLOCK_SIZE) {
    insert_supermer_block();
    buff_len = 0;
  }
  for (int i = 0; i < supermer_seq.length(); i++) {
    elem_buff_host.h_seqs_v(i + buff_len) = supermer_seq[i];
  }

  // if (pass_type == CTG_KMERS_PASS) {
  //   for (int i = 0; i < (int)supermer_seq.length(); i++) elem_buff_host.counts[buff_len + i] = supermer_count;
  // }
  buff_len += supermer_seq.length();
  elem_buff_host.h_seqs_v(buff_len) = '_';
  // if (pass_type == CTG_KMERS_PASS) elem_buff_host.counts[buff_len] = 0;
  buff_len++;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::purge_invalid(uint64_t &num_purged, uint64_t &num_entries) {
  num_purged = num_entries = 0;

  // now purge all invalid kmers (do it on the gpu)
  gpu_purge_invalid(read_kmers_dev, num_purged, num_entries);

#ifdef DEBUG
  auto expected_num_entries = read_kmers_stats.new_inserts - num_purged;
  if (num_entries != expected_num_entries)
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
  // ERROR_CHECK(Memcpy(pass_type == READ_KMERS_PASS ? &read_kmers_stats : &ctg_kmers_stats, gpu_insert_stats, sizeof(InsertStats),
  //                    MemcpyDeviceToHost));

}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::done_all_inserts(uint64_t &num_dropped, uint64_t &num_unique, uint64_t &num_purged) {
  uint64_t num_entries = 0;
  purge_invalid(num_purged, num_entries);
  read_kmers_dev.num = num_entries;

  elem_buff_host.h_seqs_v = Kokkos::View<char*>::HostMirror();
  elem_buff_host.h_counts_v = Kokkos::View<count_t*>::HostMirror();
  packed_elem_buff_dev.seqs_v = Kokkos::View<char*>();
  packed_elem_buff_dev.counts_v = Kokkos::View<count_t*>();
  unpacked_elem_buff_dev.seqs_v = Kokkos::View<char*>();
  unpacked_elem_buff_dev.counts_v = Kokkos::View<count_t*>();    
  // overallocate to reduce collisions
  num_entries *= 1.3;
  // now compact the hash table entries
  KmerExtsMap<MAX_K> compact_read_kmers_dev;
  compact_read_kmers_dev.init(num_entries);

  gpu_compact_ht(read_kmers_dev, compact_read_kmers_dev, num_dropped, num_unique);

  read_kmers_dev.clear();
 
#ifdef DEBUG
  if (num_unique != read_kmers_dev.num) WARN("mismatch in expected entries %lu != %lu", num_unique, read_kmers_dev.num);
#endif
  // now copy the gpu hash table values across to the host
  // We only do this once, which requires enough memory on the host to store the full GPU hash table, but since the GPU memory
  // is generally a lot less than the host memory, it should be fine.

  output_keys_v = Kokkos::create_mirror_view(compact_read_kmers_dev.keys_v);
  output_vals_v = Kokkos::create_mirror_view(compact_read_kmers_dev.vals_v);
  begin_iterate();

  Kokkos::deep_copy(output_keys_v, compact_read_kmers_dev.keys_v);
  Kokkos::deep_copy(output_vals_v, compact_read_kmers_dev.vals_v);
  compact_read_kmers_dev.clear();
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::done_ctg_kmer_inserts(uint64_t &attempted_inserts, uint64_t &dropped_inserts,
                                                      uint64_t &new_inserts) {
  // uint64_t *counts_gpu;
  // int NUM_COUNTS = 3;
  // ERROR_CHECK(Malloc(&counts_gpu, NUM_COUNTS * sizeof(uint64_t)));
  // ERROR_CHECK(Memset(counts_gpu, 0, NUM_COUNTS * sizeof(uint64_t)));
  // GPUTimer t;
  // int gridsize, threadblocksize;
  // get_kernel_config(ctg_kmers_dev.capacity, gpu_merge_ctg_kmers<MAX_K>, gridsize, threadblocksize);
  // t.start();
  // LaunchKernel(gpu_merge_ctg_kmers, gridsize, threadblocksize, read_kmers_dev, ctg_kmers_dev, counts_gpu);
  // t.stop();
  // dstate->kernel_timer.inc(t.get_elapsed());
  // ctg_kmers_dev.clear();
  // uint64_t counts_host[NUM_COUNTS];
  // ERROR_CHECK(Memcpy(&counts_host, counts_gpu, NUM_COUNTS * sizeof(uint64_t), MemcpyDeviceToHost));
  // ERROR_CHECK(Free(counts_gpu));
  // attempted_inserts = counts_host[0];
  // dropped_inserts = counts_host[1];
  // new_inserts = counts_host[2];
  // read_kmers_dev.num += new_inserts;
  // read_kmers_stats.new_inserts += new_inserts;
}

template <int MAX_K>
void HashTableGPUDriver<MAX_K>::get_elapsed_time(double &insert_time, double &kernel_time) {

}

template<int MAX_K>
void HashTableGPUDriver<MAX_K>::begin_iterate() {
  output_index = 0;
}

template <int MAX_K>
pair<KmerArray<MAX_K> *, CountExts *> HashTableGPUDriver<MAX_K>::get_next_entry() {

  if (!output_keys_v.size() || output_index == output_keys_v.size()) return {nullptr, nullptr};
  output_index++;
  return {&(output_keys_v(output_index - 1)), &(output_vals_v(output_index - 1))};
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
  // if (dstate->tcf) return (double)dstate->tcf->get_fill() / dstate->tcf->get_num_slots();
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
