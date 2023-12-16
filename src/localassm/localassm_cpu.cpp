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
#include <memory>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "contigs.hpp"
#include "kmer_dht.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"
#include "localassm_core.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;
using namespace localassm_core;

void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                 int qual_offset, unsigned max_read_size) {
  BarrierTimer timer(__FILEFUNC__);
  // walk should never be more than this. Note we use the maximum insert size from all libraries
  int walk_len_limit = insert_avg + 2 * insert_stddev;
  WalkMetrics wm;

  IntermittentTimer count_mers_timer(__FILENAME__ + string(":") + "count_mers"),
      walk_mers_timer(__FILENAME__ + string(":") + "walk_mers");
  ProgressBar progbar(ctgs_dht.get_local_num_ctgs(), "Extending contigs");
  ofstream ctgs_outf("ctgs-locassm-" + to_string(rank_me()) + "-k" + to_string(kmer_len));
  for (auto ctg = ctgs_dht.get_first_local_ctg(); ctg != nullptr; ctg = ctgs_dht.get_next_local_ctg()) {
    progbar.update();
    ctgs_outf << ctg->seq << " " << ctg->depth << " " << ctg->reads_left.size() << " " << ctg->reads_right.size() << endl;
    Contig ext_contig;
    extend_ctg(ctg, wm, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset, walk_len_limit, count_mers_timer,
               walk_mers_timer);
    ctgs.add_contig({.id = ctg->cid, .seq = ctg->seq, .depth = ctg->depth});
  }
  progbar.done();
  ctgs_outf.close();
  count_mers_timer.done_all();
  walk_mers_timer.done_all();
  barrier();
  SLOG_VERBOSE("Walk terminations: ", reduce_one(wm.term_counts[0], op_fast_add, 0).wait(), " X, ",
               reduce_one(wm.term_counts[1], op_fast_add, 0).wait(), " F, ", reduce_one(wm.term_counts[2], op_fast_add, 0).wait(),
               " R\n");
  auto tot_num_reads = reduce_one(wm.num_reads, op_fast_add, 0).wait();
  auto tot_num_walks = reduce_one(wm.num_walks, op_fast_add, 0).wait();
  auto tot_sum_ext = reduce_one(wm.sum_ext, op_fast_add, 0).wait();
  auto tot_sum_clen = reduce_one(wm.sum_clen, op_fast_add, 0).wait();
  auto tot_max_walk_len = reduce_one(wm.max_walk_len, op_fast_max, 0).wait();
  auto tot_excess_reads = reduce_one(wm.excess_reads, op_fast_add, 0).wait();
  auto num_ctgs = ctgs_dht.get_num_ctgs();
  auto tot_num_sides = reduce_one(wm.num_sides, op_fast_add, 0).wait();
  SLOG_VERBOSE("Used a total of ", tot_num_reads, " reads, max per ctg ", wm.max_num_reads, " avg per ctg ",
               (num_ctgs > 0 ? (tot_num_reads / num_ctgs) : 0), ", dropped ", perc_str(tot_excess_reads, tot_num_reads),
               " excess reads\n");
  SLOG_VERBOSE("Could walk ", perc_str(tot_num_sides, num_ctgs * 2), " contig sides\n");
  if (tot_sum_clen)
    SLOG_VERBOSE("Found ", tot_num_walks, " walks, total extension length ", tot_sum_ext, " extended ",
                 (double)(tot_sum_ext + tot_sum_clen) / tot_sum_clen, "\n");
  if (tot_num_walks)
    SLOG_VERBOSE("Average walk length ", tot_sum_ext / tot_num_walks, ", max walk length ", tot_max_walk_len, "\n");
}
