#pragma once

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

#include <vector>
#include <array>
#include <unordered_map>
#include <thread>

#include "hash_funcs.h"

#include <Kokkos_Core.hpp>

namespace kcount_gpu {

enum PASS_TYPE { READ_KMERS_PASS = 0, CTG_KMERS_PASS = 1 };

using count_t = uint32_t;
using ext_count_t = uint16_t;

struct CountsArray {
  static const int n_exts = 8;
  uint32_t kmer_count;
  uint16_t ext_counts[n_exts];
};

struct CountExts {
  count_t count;
  int8_t left, right;
};

template <int MAX_K>
struct KmerArray {
  static const int N_LONGS = (MAX_K + 31) / 32;
  uint64_t longs[N_LONGS];

  // void set(const uint64_t *x);
};

struct SupermerBuff {
  // char *seqs;
  // count_t *counts;
  Kokkos::View<char*> seqs_v;
  Kokkos::View<char*>::HostMirror h_seqs_v;

  Kokkos::View<count_t*> counts_v;
  Kokkos::View<count_t*>::HostMirror h_counts_v;
};

// Bytes used per element:
// k = 21: 8+20 = 28
// k = 33, 55: 16+20 = 36
// k = 77: 24+20 = 44
// k = 99: 32+20 = 52
template <int MAX_K>
struct KmerCountsMap {
  // Arrays for keys and values. They are separate because the keys get initialized with max number and the vals with zero
  // KmerArray<MAX_K> *keys = nullptr;
  Kokkos::View<KmerArray<MAX_K>*> keys_v;

  // CountsArray *vals = nullptr;
  Kokkos::View<CountsArray*> vals_v;
  
  uint64_t capacity = 0;
  uint64_t num = 0;

  void init(int64_t ht_capacity);
  void clear();
};

template <int MAX_K>
struct KmerExtsMap {
  // KmerArray<MAX_K> *keys = nullptr;
  Kokkos::View<KmerArray<MAX_K>*> keys_v;

  // CountExts *vals = nullptr;
  Kokkos::View<CountExts*> vals_v;

  int64_t capacity = 0;

  void init(int64_t ht_capacity);
  void clear();
};

struct InsertStats {
  uint64_t dropped = 0;
  uint64_t dropped_qf = 0;
  uint64_t attempted = 0;
  uint64_t new_inserts = 0;
  uint64_t num_unique_qf = 0;
};

template <int MAX_K>
class HashTableGPUDriver {
  static const int N_LONGS = (MAX_K + 31) / 32;
  struct HashTableDriverState;
  // stores CUDA specific variables
  HashTableDriverState *dstate = nullptr;

  int upcxx_rank_me;
  int upcxx_rank_n;
  int kmer_len;
  uint64_t buff_len = 0;
  // std::vector<KmerArray<MAX_K>> output_keys;
  // std::vector<CountExts> output_vals;
  typename Kokkos::View<KmerArray<MAX_K>*>::HostMirror output_keys_v;
  typename Kokkos::View<CountExts*>::HostMirror output_vals_v;
  size_t output_index = 0;

  KmerCountsMap<MAX_K> read_kmers_dev;
  KmerCountsMap<MAX_K> ctg_kmers_dev;

  // for buffering elements in the host memory
  SupermerBuff elem_buff_host;

  // for transferring host memory buffer to device
  SupermerBuff unpacked_elem_buff_dev;
  SupermerBuff packed_elem_buff_dev;

  InsertStats read_kmers_stats;
  InsertStats ctg_kmers_stats;
  // InsertStats *gpu_insert_stats;
  typename Kokkos::View<InsertStats>::HostMirror read_kmers_stats_v;
  typename Kokkos::View<InsertStats>::HostMirror ctg_kmers_stats_v;  
  typename Kokkos::View<InsertStats> gpu_insert_stats_v;
  int num_gpu_calls = 0;

  const uint64_t KEY_EMPTY = 0xffffffffffffffff;
  const uint64_t KEY_TRANSITION = 0xfffffffffffffffe;
  const uint8_t KEY_EMPTY_BYTE = 0xff;  

  void insert_supermer_block();
  void purge_invalid(uint64_t &num_purged, uint64_t &num_entries);

  typename Kokkos::View<uint64_t[256]> twins_v;


 public:
  PASS_TYPE pass_type;

  HashTableGPUDriver();
  ~HashTableGPUDriver();

  void init(int upcxx_rank_me, int upcxx_rank_n, int kmer_len, size_t max_elems, size_t max_ctg_elems, size_t num_errors,
            size_t gpu_avail_mem, std::string &msgs, std::string &warnings, bool use_qf);

  void init_ctg_kmers(uint64_t max_elems, size_t gpu_avail_mem);

  void insert_supermer(const std::string &supermer_seq, count_t supermer_count);

  void flush_inserts();

  void done_ctg_kmer_inserts(uint64_t &attempted_inserts, uint64_t &dropped_inserts, uint64_t &new_inserts);

  void done_all_inserts(uint64_t &num_dropped, uint64_t &num_unique, uint64_t &num_purged);

  void begin_iterate();
  
  std::pair<KmerArray<MAX_K> *, CountExts *> get_next_entry();

  void get_elapsed_time(double &insert_time, double &kernel_time);

  int64_t get_capacity();

  int64_t get_final_capacity();

  InsertStats &get_stats();

  int get_num_gpu_calls();

  double get_qf_load_factor();
};

}  // namespace kcount_gpu
