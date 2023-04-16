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

#include "localassm_struct.hpp"

namespace localassm_driver {

struct accum_data {
  std::vector<uint64_t> ht_sizes;
  std::vector<uint64_t> l_reads_count;
  std::vector<uint64_t> r_reads_count;
  std::vector<uint64_t> ctg_sizes;
  void clear() {
    ht_sizes.clear();
    l_reads_count.clear();
    r_reads_count.clear();
    ctg_sizes.clear();
  }
};

struct ctg_bucket {
  std::vector<CtgWithReads> ctg_vec;
  accum_data sizes_vec;
  uint32_t l_max, r_max, max_contig_sz, max_read_sz;
  uint64_t tot_ht, tot_ctg, tot_l_reads, tot_r_reads, count;
  ctg_bucket(uint32_t max_read_sz)
      : ctg_vec{}
      , sizes_vec{}
      , l_max{0}
      , r_max{0}
      , max_contig_sz{0}
      , max_read_sz{max_read_sz} 
      , tot_ht{0}
      , tot_ctg{0}
      , tot_l_reads{0}
      , tot_r_reads{0}
      , count{0}
      {}
  ctg_bucket(const ctg_bucket &copy) = default;
  ctg_bucket(ctg_bucket &&move) = default;
  void add(CtgWithReads &&cwr) {
    assert(max_read_sz >= cwr.get_max_read_size() && "CtgWithRead max_read_size exceeds ctg_bucket max_read_sz");
    auto max_reads = cwr.max_reads;
    auto seq_size = cwr.seq.size();
    auto reads_left_size = cwr.reads_left.size();
    auto reads_right_size = cwr.reads_right.size();
    ctg_vec.emplace_back(std::move(cwr));
    uint32_t temp_ht_size = max_reads * max_read_sz;
    sizes_vec.ht_sizes.push_back(temp_ht_size);
    tot_ht += temp_ht_size;
    sizes_vec.ctg_sizes.push_back(seq_size);
    tot_ctg += seq_size;
    sizes_vec.l_reads_count.push_back(reads_left_size);
    tot_l_reads += reads_left_size;
    sizes_vec.r_reads_count.push_back(reads_right_size);
    tot_r_reads += reads_right_size;
    count++;
  
    if (l_max < reads_left_size) l_max = reads_left_size;
    if (r_max < reads_right_size) r_max = reads_right_size;
    if (max_contig_sz < seq_size) max_contig_sz = seq_size;
  }
  void clear() {
    ctg_vec.clear();
    sizes_vec.clear();
    l_max = r_max = max_contig_sz = max_read_sz = 0;
  }
};

void localassm_driver(std::vector<CtgWithReads>& data_in, uint32_t max_ctg_size, uint32_t max_read_size, uint32_t max_r_count,
                      uint32_t max_l_count, int mer_len, int max_kmer_len, accum_data& sizes_outliers, int walk_len_limit,
                      int qual_offset, int my_rank, size_t gpu_mem_avail, int debug_line);

}  // namespace localassm_driver
