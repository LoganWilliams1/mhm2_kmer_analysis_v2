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

#include <stdio.h>
#include <iostream>

#include "gpu-utils/gpu_compatibility.hpp"
#include "gpu-utils/gpu_common.hpp"

#define FULL 0xFFFFFFFF
#define EMPTY 0xFFFFFFFE
#define FULL_MASK 0xffffffff

#ifdef CUDA_GPU
#define WARP_SIZE 32
#endif
#ifdef HIP_GPU
#define WARP_SIZE 64
#endif

struct cstr_type {
  char* start_ptr;
  int length;
  __device__ cstr_type() {
    start_ptr = nullptr;
    length = 0;
  }
  __device__ cstr_type(char* ptr, int len) {
    start_ptr = ptr;
    length = len;
  }

  __device__ bool operator==(const cstr_type& in2) {
    bool str_eq = (length == in2.length) & (length != EMPTY) & (length != FULL);
    if (str_eq)
      for (int i = 0; i < in2.length; i++) {
        if (start_ptr[i] != in2.start_ptr[i]) {
          str_eq = false;
          break;
        }
      }
    return str_eq;
  }

};

__device__ void cstr_copy(cstr_type& str1, cstr_type& str2);

namespace gpu_loc_assem {

struct ExtCounts {
  uint32_t count_A;
  uint32_t count_C;
  uint32_t count_G;
  uint32_t count_T;

  __device__ void print() {  // for debugging
    printf("count_A:%d, count_C:%d, count_G:%d, count_T:%d\n", count_A, count_C, count_G, count_T);
  }

  // TODO: replace  numeric_limits by either a suitable constant/predefined value or find a device alternate
  __device__ void inc(char ext, int count) {
    switch (ext) {
      case 'A': atomicAdd(&count_A, count); break;
      case 'C': atomicAdd(&count_C, count); break;
      case 'G': atomicAdd(&count_G, count); break;
      case 'T': atomicAdd(&count_T, count); break;
    }
  }
};

struct MerBase {
  char base;
  uint32_t nvotes_hi_q, nvotes, rating;
  __device__ void print() {  // for debuggin
    printf("base:%c, nvotes_hiq_q:%d, nvotes:%d, rating:%d\n", base, nvotes_hi_q, nvotes, rating);
  }

  __device__ uint16_t get_base_rating(int depth) {
    double min_viable = max(LASSM_MIN_VIABLE_DEPTH * depth, 2.0);
    double min_expected_depth = max(LASSM_MIN_EXPECTED_DEPTH * depth, 2.0);
    if (nvotes == 0) return 0;
    if (nvotes == 1) return 1;
    if (nvotes < min_viable) return 2;
    if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q < min_viable) return 3;
    if (min_expected_depth > nvotes && nvotes >= min_viable && nvotes_hi_q >= min_viable) return 4;
    if (nvotes >= min_expected_depth && nvotes_hi_q < min_viable) return 5;
    if (nvotes >= min_expected_depth && min_viable < nvotes_hi_q && nvotes_hi_q < min_expected_depth) return 6;
    return 7;
  }
};

struct MerFreqs {
  // how many times this kmer has occurred: don't need to count beyond 65536
  // count of high quality extensions and low quality extensions - structure comes from kmer_dht.hpp
  ExtCounts hi_q_exts, low_q_exts;
  // the final extensions chosen - A,C,G,T, or F,X
  char ext;
  // the count of the final extension
  int count;
  __device__ bool comp_merbase(MerBase& elem1, MerBase& elem2) {
    if (elem1.rating != elem2.rating) return elem1.rating > elem2.rating;
    if (elem1.nvotes_hi_q != elem2.nvotes_hi_q) return elem1.nvotes_hi_q > elem2.nvotes_hi_q;
    if (elem1.nvotes != elem2.nvotes) return elem1.nvotes > elem2.nvotes;

    return true;
  }

  __device__ void sort_merbase(MerBase (&merbases)[4]) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        if (comp_merbase(merbases[i], merbases[j])) {
          MerBase temp = merbases[i];
          merbases[i] = merbases[j];
          merbases[j] = temp;
        }
      }
    }
  }
  __device__ void set_ext(int seq_depth) {
    // set extension similarly to how it is done with localassm in mhm
    MerBase mer_bases[4] = {{.base = 'A', .nvotes_hi_q = hi_q_exts.count_A, .nvotes = low_q_exts.count_A},
                            {.base = 'C', .nvotes_hi_q = hi_q_exts.count_C, .nvotes = low_q_exts.count_C},
                            {.base = 'G', .nvotes_hi_q = hi_q_exts.count_G, .nvotes = low_q_exts.count_G},
                            {.base = 'T', .nvotes_hi_q = hi_q_exts.count_T, .nvotes = low_q_exts.count_T}};
    for (int i = 0; i < 4; i++) {
      mer_bases[i].rating = mer_bases[i].get_base_rating(seq_depth);
    }

    // sort bases in descending order of quality
    sort_merbase(mer_bases);

    int top_rating = mer_bases[0].rating;
    int runner_up_rating = mer_bases[1].rating;
    // assert(top_rating >= runner_up_rating);// for now finding a way around for assertion
    if (top_rating < runner_up_rating) printf("******* POSSIBLE ERROR IN sort_merbase************");
    int top_rated_base = mer_bases[0].base;
    ext = 'X';
    count = 0;
    // no extension (base = 0) if the runner up is close to the top rating
    // except, if rating is 7 (best quality), then all bases of rating 7 are forks
    if (top_rating > LASSM_RATING_THRES) {  // must have at least minViable bases
      if (top_rating <= 3) {                // must be uncontested
        if (runner_up_rating == 0) ext = top_rated_base;
      } else if (top_rating < 6) {
        if (runner_up_rating < 3) ext = top_rated_base;
      } else if (top_rating == 6) {  // viable and fair hiQ support
        if (runner_up_rating < 4) ext = top_rated_base;
      } else {  // strongest rating trumps
        if (runner_up_rating < 7) {
          ext = top_rated_base;
        } else {
          if (mer_bases[2].rating == 7 || mer_bases[0].nvotes == mer_bases[1].nvotes)
            ext = 'F';
          else if (mer_bases[0].nvotes > mer_bases[1].nvotes)
            ext = mer_bases[0].base;
          else if (mer_bases[1].nvotes > mer_bases[0].nvotes)
            ext = mer_bases[1].base;
        }
      }
    }
    for (int i = 0; i < 4; i++) {
      if (mer_bases[i].base == ext) {
        count = mer_bases[i].nvotes;
        break;
      }
    }
  }
};

}  // end of namespace gpu_loc_assem
struct loc_ht {
  cstr_type key;
  gpu_loc_assem::MerFreqs val;
  __device__ loc_ht() : key{}, val{} {}
  __device__ loc_ht(cstr_type in_key, gpu_loc_assem::MerFreqs in_val) {
    key = in_key;
    val = in_val;
  }
  __device__ static bool is_valid(const loc_ht& x) { return x.key.length != FULL; } // EMPTY is valid
};

struct loc_ht_bool {
  cstr_type key;
  bool val;
  __device__ loc_ht_bool() : key{}, val{} {}
  __device__ loc_ht_bool(cstr_type in_key, bool in_val) {
    key = in_key;
    val = in_val;
  }
  __device__ static bool is_valid(const loc_ht_bool& x) { return x.key.length != FULL; } // EMPTY is valid
};

__device__ void print_mer(cstr_type& mer);
__global__ void ht_kernel(loc_ht* ht, char* contigs, int* offset_sum, int kmer_size);
__device__ bool ht_insert(loc_ht* thread_ht, cstr_type kmer_key, cstr_type ctg_val, uint32_t max_size);
__device__ bool ht_insert(loc_ht_bool* thread_ht, cstr_type kmer_key, bool bool_val, uint32_t max_size);
__device__ void ht_delete(loc_ht* thread_ht, cstr_type kmer_key, uint32_t max_size);
__device__ loc_ht& ht_get(loc_ht* thread_ht, cstr_type kmer_key, uint32_t max_size);
__device__ unsigned hash_func(cstr_type key, uint32_t max_size);
__device__ void count_mers(loc_ht* thrd_loc_ht, char* loc_r_reads, uint32_t max_ht_size, char* loc_r_quals, int32_t* reads_r_offset,
                           int32_t& r_rds_cnt, int32_t* rds_count_r_sum, double& loc_ctg_depth, int& mer_len, uint32_t& qual_offset,
                           int64_t& excess_reads, const long int idx);
__device__ char walk_mers(loc_ht* thrd_loc_ht, loc_ht_bool* thrd_ht_bool, uint32_t max_ht_size, int& mer_len,
                          cstr_type& mer_walk_temp, cstr_type& longest_walk, cstr_type& walk, const int idx, int max_walk_len);
__global__ void iterative_walks_kernel(uint64_t* cid, uint32_t* ctg_offsets, char* contigs, char* reads_r, char* quals_r,
                                       uint32_t* reads_r_offset, uint32_t* rds_count_r_sum, double* ctg_depth, loc_ht* global_ht,
                                       uint32_t* prefix_ht, loc_ht_bool* global_ht_bool, int kmer_len, uint32_t max_mer_len_off,
                                       uint32_t* term_counts, int64_t num_walks, int64_t max_walk_len, int64_t sum_ext,
                                       int32_t max_read_size, int32_t max_read_count, uint32_t qual_offset, char* longest_walks,
                                       char* mer_walk_temp, uint32_t* final_walk_lens, int tot_ctgs);
