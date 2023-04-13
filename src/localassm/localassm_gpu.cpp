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

#include <fstream>
#include <iostream>
#include <regex>
#include <upcxx/upcxx.hpp>
#include <memory>

#include "alignments.hpp"
#include "contigs.hpp"
#include "kmer_dht.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"
#include "localassm_core.hpp"
#include "devices_gpu.hpp"

#include "localassm-gpu/driver.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;
using namespace localassm_core;

static vector<ReadSeq> reads_to_reads(vector<ReadSeq> read_in) {
  vector<ReadSeq> reads_out;
  for (int i = 0; i < min((int)read_in.size(), (int)LASSM_MAX_COUNT_MERS_READS); i++) {
    ReadSeq temp_seq_in = read_in[i];
    ReadSeq temp_seq_out;

    temp_seq_out.read_id = temp_seq_in.read_id;
    temp_seq_out.seq = temp_seq_in.seq;
    temp_seq_out.quals = temp_seq_in.quals;
    reads_out.push_back(temp_seq_out);
  }
  return reads_out;
}

static CtgWithReads ctgs_to_ctgs(CtgWithReads ctg_in) {
  CtgWithReads ctg_out;
  ctg_out.cid = ctg_in.cid;
  ctg_out.seq = ctg_in.seq;
  ctg_out.depth = ctg_in.depth;
  ctg_out.max_reads = ctg_in.max_reads;
  vector<ReadSeq> temp_reads = reads_to_reads(ctg_in.reads_left);
  ctg_out.reads_left = temp_reads;
  temp_reads = reads_to_reads(ctg_in.reads_right);
  ctg_out.reads_right = temp_reads;
  return ctg_out;
}

static void bucket_ctgs(localassm_driver::ctg_bucket &zero_slice, localassm_driver::ctg_bucket &mid_slice,
                        localassm_driver::ctg_bucket &outlier_slice, CtgsWithReadsDHT &ctgs_dht,
                        unsigned max_read_size, IntermittentTimer &ctg_buckets_timer) {
  ctg_buckets_timer.start();
  for (auto ctg = ctgs_dht.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_dht.get_next_local_ctg()) {
    CtgWithReads temp_in = ctgs_to_ctgs(*ctg);
    temp_in.max_reads =
        temp_in.reads_left.size() > temp_in.reads_right.size() ? temp_in.reads_left.size() : temp_in.reads_right.size();
    if (temp_in.max_reads == 0) {
      zero_slice.ctg_vec.push_back(temp_in);
    } else if (temp_in.max_reads > 0 && temp_in.max_reads < 10) {
      mid_slice.ctg_vec.push_back(temp_in);
      uint32_t temp_ht_size = temp_in.max_reads * max_read_size;
      mid_slice.sizes_vec.ht_sizes.push_back(temp_ht_size);
      mid_slice.sizes_vec.ctg_sizes.push_back(temp_in.seq.size());
      mid_slice.sizes_vec.l_reads_count.push_back(temp_in.reads_left.size());
      mid_slice.sizes_vec.r_reads_count.push_back(temp_in.reads_right.size());
      if (mid_slice.l_max < temp_in.reads_left.size()) mid_slice.l_max = temp_in.reads_left.size();
      if (mid_slice.r_max < temp_in.reads_right.size()) mid_slice.r_max = temp_in.reads_right.size();
      if (mid_slice.max_contig_sz < temp_in.seq.size()) mid_slice.max_contig_sz = temp_in.seq.size();
    } else {
      outlier_slice.ctg_vec.push_back(temp_in);
      uint32_t temp_ht_size = temp_in.max_reads * max_read_size;
      outlier_slice.sizes_vec.ht_sizes.push_back(temp_ht_size);
      outlier_slice.sizes_vec.ctg_sizes.push_back(temp_in.seq.size());
      outlier_slice.sizes_vec.l_reads_count.push_back(temp_in.reads_left.size());
      outlier_slice.sizes_vec.r_reads_count.push_back(temp_in.reads_right.size());
      if (outlier_slice.l_max < temp_in.reads_left.size()) outlier_slice.l_max = temp_in.reads_left.size();
      if (outlier_slice.r_max < temp_in.reads_right.size()) outlier_slice.r_max = temp_in.reads_right.size();
      if (outlier_slice.max_contig_sz < temp_in.seq.size()) outlier_slice.max_contig_sz = temp_in.seq.size();
    }
  }
  ctg_buckets_timer.stop();
}

void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                 int qual_offset, unsigned max_read_size) {
  BarrierTimer timer(__FILEFUNC__);
  // walk should never be more than this. Note we use the maximum insert size from all libraries
  int walk_len_limit = insert_avg + 2 * insert_stddev;
  WalkMetrics wm;

  IntermittentTimer count_mers_timer(__FILENAME__ + string(":") + "count_mers"),
      walk_mers_timer(__FILENAME__ + string(":") + "walk_mers"), ctg_buckets_timer(__FILENAME__ + string(":") + "bucket_ctgs"),
      loc_assem_kernel_timer(__FILENAME__ + string(":") + "GPU_locassem");

  ProgressBar progbar(ctgs_dht.get_local_num_ctgs(), "Extending contigs");

  localassm_driver::ctg_bucket zero_slice, mid_slice, outlier_slice;
  bucket_ctgs(zero_slice, mid_slice, outlier_slice, ctgs_dht, max_read_size, ctg_buckets_timer);
  ctg_buckets_timer.done_all();

  auto gpu_avail_mem_per_rank = get_gpu_avail_mem_per_rank();  // implicit local_team barier
  future<> fut_outlier = make_future();
  if (outlier_slice.ctg_vec.size() > 0)
    fut_outlier = upcxx_utils::execute_in_thread_pool(
        [&outlier_slice, max_read_size, walk_len_limit, qual_offset, max_kmer_len, kmer_len, gpu_avail_mem_per_rank, &loc_assem_kernel_timer]() {
          loc_assem_kernel_timer.start();
          localassm_driver::localassm_driver(outlier_slice.ctg_vec, outlier_slice.max_contig_sz, max_read_size, outlier_slice.r_max,
                                             outlier_slice.l_max, kmer_len, max_kmer_len, outlier_slice.sizes_vec, walk_len_limit,
                                             qual_offset, local_team().rank_me(), gpu_avail_mem_per_rank);
          loc_assem_kernel_timer.stop();
        });

  // work steal while either:
  //    my outliers are running on the GPU
  // OR other members of the local_team are still processing outliers (and consuming GPU memory)
  // OR if there are less than 100 mid_slice contigs to localassm
  future<> fut_local_barrier = upcxx::barrier_async(upcxx::local_team());
  upcxx::discharge();
  auto tot_mids{mid_slice.ctg_vec.size()};
  while ((!fut_outlier.ready() && mid_slice.ctg_vec.size() > 0) ||
         (!fut_local_barrier.ready() && mid_slice.ctg_vec.size() > 0) ||
         (mid_slice.ctg_vec.size() <= 100 && mid_slice.ctg_vec.size() > 0)) {
    auto ctg = &mid_slice.ctg_vec.back();
    extend_ctg(ctg, wm, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset, walk_len_limit, count_mers_timer,
               walk_mers_timer);
    ctgs.add_contig({.id = ctg->cid, .seq = ctg->seq, .depth = ctg->depth});
    mid_slice.ctg_vec.pop_back();
    upcxx::progress();
  }
  fut_outlier.wait();
  fut_local_barrier.wait();
  auto cpu_exts{tot_mids - mid_slice.ctg_vec.size()};
  LOG("Number of Local Contig Extensions processed on CPU:", cpu_exts, "\n");

  auto remaining_gpu_avail_mem_per_rank = get_gpu_avail_mem_per_rank();  // implicit local_team barier
  if (mid_slice.ctg_vec.size() > 0) {
    loc_assem_kernel_timer.start();
    localassm_driver::localassm_driver(mid_slice.ctg_vec, mid_slice.max_contig_sz, max_read_size, mid_slice.r_max, mid_slice.l_max,
                                       kmer_len, max_kmer_len, mid_slice.sizes_vec, walk_len_limit, qual_offset,
                                       local_team().rank_me(), remaining_gpu_avail_mem_per_rank);
    loc_assem_kernel_timer.stop();
  }
  
  for (int j = 0; j < zero_slice.ctg_vec.size(); j++) {
    CtgWithReads temp_ctg = zero_slice.ctg_vec[j];
    ctgs.add_contig({.id = temp_ctg.cid, .seq = temp_ctg.seq, .depth = temp_ctg.depth});
  }
  for (int j = 0; j < mid_slice.ctg_vec.size(); j++) {
    CtgWithReads temp_ctg = mid_slice.ctg_vec[j];
    ctgs.add_contig({.id = temp_ctg.cid, .seq = temp_ctg.seq, .depth = temp_ctg.depth});
  }
  for (int j = 0; j < outlier_slice.ctg_vec.size(); j++) {
    CtgWithReads temp_ctg = outlier_slice.ctg_vec[j];
    ctgs.add_contig({.id = temp_ctg.cid, .seq = temp_ctg.seq, .depth = temp_ctg.depth});
  }

  count_mers_timer.done_all();
  walk_mers_timer.done_all();
  loc_assem_kernel_timer.done_all();
  // implicit barrier from BarrierTimer
}
