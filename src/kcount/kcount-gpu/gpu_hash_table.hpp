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

namespace kcount_gpu {

enum PASS_TYPE { READ_KMERS_PASS = 0, CTG_KMERS_PASS = 1 };

using count_t = uint32_t;
using ext_count_t = uint16_t;

struct CountsArray {
  count_t kmer_count;
  ext_count_t ext_counts[8];
};

struct CountExts {
  count_t count;
  int8_t left, right;
};

template <int MAX_K>
struct KmerArray {
  static const int N_LONGS = (MAX_K + 31) / 32;
  uint64_t longs[N_LONGS];

  void set(const uint64_t *x);
};

struct SupermerBuff {
  char *seqs;
  count_t *counts;
};

// Bytes used per element:
// k = 21: 8+20 = 28
// k = 33, 55: 16+20 = 36
// k = 77: 24+20 = 44
// k = 99: 32+20 = 52
template <int MAX_K>
struct KmerCountsMap {
  // Arrays for keys and values. They are separate because the keys get initialized with max number and the vals with zero
  KmerArray<MAX_K> *keys = nullptr;
  CountsArray *vals = nullptr;
  int64_t capacity = 0;
  int num = 0;

  void init(int64_t ht_capacity);
  void clear();
};

template <int MAX_K>
struct KmerExtsMap {
  KmerArray<MAX_K> *keys = nullptr;
  CountExts *vals = nullptr;
  int64_t capacity = 0;

  void init(int64_t ht_capacity);
  void clear();
};

struct InsertStats {
  unsigned int dropped = 0;
  unsigned int attempted = 0;
  unsigned int new_inserts = 0;
  unsigned int num_unique_qf = 0;
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
  int buff_len = 0;
  std::vector<KmerArray<MAX_K>> output_keys;
  std::vector<CountExts> output_vals;
  size_t output_index = 0;

  KmerCountsMap<MAX_K> read_kmers_dev;
  KmerCountsMap<MAX_K> ctg_kmers_dev;

  // for buffering elements in the host memory
  SupermerBuff elem_buff_host = {0};
  // for transferring host memory buffer to device
  SupermerBuff unpacked_elem_buff_dev = {0};
  SupermerBuff packed_elem_buff_dev = {0};

  InsertStats read_kmers_stats;
  InsertStats ctg_kmers_stats;
  InsertStats *gpu_insert_stats;
  int num_gpu_calls = 0;

  void insert_supermer_block();
  void purge_invalid(int &num_purged, int &num_entries);

 public:
  PASS_TYPE pass_type;

  HashTableGPUDriver();
  ~HashTableGPUDriver();

  void init(int upcxx_rank_me, int upcxx_rank_n, int kmer_len, int max_elems, size_t gpu_avail_mem, double &init_time,
            size_t &gpu_bytes_reqd, size_t &ht_bytes_used, size_t &qf_bytes_used, bool use_qf);

  void init_ctg_kmers(int max_elems, size_t gpu_avail_mem);

  void insert_supermer(const std::string &supermer_seq, count_t supermer_count);

  void flush_inserts();

  void done_ctg_kmer_inserts(int &attempted_inserts, int &dropped_inserts, int &new_inserts);

  void done_all_inserts(int &num_dropped, int &num_unique, int &num_purged);

  std::pair<KmerArray<MAX_K> *, CountExts *> get_next_entry();

  void get_elapsed_time(double &insert_time, double &kernel_time);

  int64_t get_capacity();

  int64_t get_final_capacity();

  InsertStats &get_stats();

  int get_num_gpu_calls();

  double get_qf_load_factor();

  uint64_t get_qf_failures();
};

}  // namespace kcount_gpu
