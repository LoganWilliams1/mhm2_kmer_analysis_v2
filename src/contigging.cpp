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

#include "contigging.hpp"

#include "gasnet_stats.hpp"
#include "kcount/kcount.hpp"
#include "kmer_dht.hpp"
#include "stage_timers.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"

using namespace upcxx;
using namespace upcxx_utils;

using std::fixed;
using std::setprecision;
using std::shared_ptr;
using std::string;
using std::tie;
using std::vector;

// template <int MAX_K>
// void traverse_debruijn_graph(unsigned kmer_len, dist_object<KmerDHT<MAX_K>> &kmer_dht, Contigs &my_uutigs);
// void localassm(int max_kmer_len, int kmer_len, PackedReadsList &packed_reads_list, int insert_avg, int insert_stddev,
//                int qual_offset, Contigs &ctgs, const Alns &alns);

template <int MAX_K>
void contigging(int kmer_len, int prev_kmer_len, int &rlen_limit, PackedReadsList &packed_reads_list, Contigs &ctgs,
                Histogrammer &histogrammer, Options &options) {
  auto loop_start_t = clock_now();
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Contig generation k = ", kmer_len, KNORM, "\n");
  SLOG("\n");
  LOG_MEM("Starting contigging k=" + to_string(kmer_len));
  bool is_debug = false;
#ifdef DEBUG
  is_debug = true;
#endif

  auto max_kmer_store = options.max_kmer_store_mb * ONE_MB;

  string uutigs_fname("uutigs-" + to_string(kmer_len) + ".fasta");
  if (options.ctgs_fname != uutigs_fname) {
    Kmer<MAX_K>::set_k(kmer_len);
    // duration of kmer_dht
    stage_timers.analyze_kmers->start();
    auto my_num_kmers = reduce_all(PackedReads::estimate_num_kmers(kmer_len, packed_reads_list), op_fast_add).wait() / rank_n();
    auto my_num_ctg_kmers = reduce_all(ctgs.get_num_ctg_kmers(kmer_len), op_fast_add).wait() / rank_n();
    dist_object<KmerDHT<MAX_K>> kmer_dht(world(), my_num_kmers, my_num_ctg_kmers, max_kmer_store, options.max_rpcs_in_flight,
                                         options.use_qf, options.sequencing_depth);
    LOG_MEM("Allocated kmer_dht");
    barrier();
    begin_gasnet_stats("kmer_analysis k = " + to_string(kmer_len));
    analyze_kmers(kmer_len, prev_kmer_len, options.qual_offset, packed_reads_list, options.dmin_thres, ctgs, kmer_dht,
                  options.dump_kmers);
    LOG_MEM("Analyzed kmers");
    end_gasnet_stats();
    stage_timers.analyze_kmers->stop();
    auto avg_kmer_count = kmer_dht->get_avg_kmer_count();
    SLOG_VERBOSE("Changing sequencing depth from ", options.sequencing_depth, " to ", (int)avg_kmer_count, "\n");
    options.sequencing_depth = (int)avg_kmer_count;
    barrier();
    LOG_MEM("Analyzed kmers");
    // stage_timers.dbjg_traversal->start();
    // begin_gasnet_stats("dbjg_traversal k = " + to_string(kmer_len));
    // traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
    // end_gasnet_stats();
    // LOG_MEM("Traversed dbg");
    // stage_timers.dbjg_traversal->stop();
    // if (is_debug) {
    //   stage_timers.dump_ctgs->start();
    //   ctgs.dump_contigs(uutigs_fname, 0, "uutig_");
    //   stage_timers.dump_ctgs->stop();
    // }
  }
  LOG_MEM("Generated contigs k=" + to_string(kmer_len));

  if (kmer_len < options.kmer_lens.back()) {
    if (kmer_len == options.kmer_lens.front()) {
      size_t num_reads = PackedReads::get_total_local_num_reads(packed_reads_list);
      auto avg_num_reads = reduce_one(num_reads, op_fast_add, 0).wait() / rank_n();
      auto max_num_reads = reduce_one(num_reads, op_fast_max, 0).wait();
      SLOG_VERBOSE("Avg reads per rank ", avg_num_reads, " max ", max_num_reads, " (balance ",
                   (double)avg_num_reads / max_num_reads, ")\n");
      // if (options.shuffle_reads) {
      //   stage_timers.shuffle_reads->start();
      //   begin_gasnet_stats("shuffle_reads k = " + to_string(kmer_len));
      //   shuffle_reads(options.qual_offset, packed_reads_list, ctgs);
      //   end_gasnet_stats();
      //   stage_timers.shuffle_reads->stop();
      //   LOG_MEM("Shuffled reads");
      //   num_reads = 0;
      //   for (auto packed_reads : packed_reads_list) {
      //     num_reads += packed_reads->get_local_num_reads();
      //   }
      //   avg_num_reads = reduce_one(num_reads, op_fast_add, 0).wait() / rank_n();
      //   max_num_reads = reduce_one(num_reads, op_fast_max, 0).wait();
      //   SLOG_VERBOSE("After shuffle: avg reads per rank ", avg_num_reads, " max ", max_num_reads, " (load balance ",
      //                (double)avg_num_reads / max_num_reads, ")\n");
      //   rlen_limit = 0;
      //   for (auto packed_reads : packed_reads_list) {
      //     rlen_limit = max(rlen_limit, (int)packed_reads->get_max_read_len());
      //   }
      // }
    }
    barrier();
    // Alns alns;
    // stage_timers.alignments->start();
    // begin_gasnet_stats("alignment k = " + to_string(kmer_len));
    // bool first_ctg_round = (kmer_len == options.kmer_lens[0]);
    // auto [kernel_elapsed, aln_comms_elapsed] =
    //     find_alignments<MAX_K>(kmer_len, packed_reads_list, max_kmer_store, options.max_rpcs_in_flight, ctgs, alns,
    //                            KLIGN_SEED_SPACE, rlen_limit, options.optimize_for == "contiguity", 0, options.klign_rget_buf_size);
    // end_gasnet_stats();
    // stage_timers.kernel_alns->inc_elapsed(kernel_elapsed);
    // stage_timers.aln_comms->inc_elapsed(aln_comms_elapsed);
    // stage_timers.alignments->stop();
    // barrier();
    // LOG_MEM("Aligned reads to contigs");
    // if (is_debug) alns.dump_single_file("ctg-alns-" + to_string(kmer_len) + ".blast", Alns::Format::BLAST);
    // histogrammer.calculate_insert_size(alns);
    // // insert size should never be larger than this; if it is that signals some error in the assembly
    // barrier();
    // stage_timers.localassm->start();
    // begin_gasnet_stats("local_assembly k = " + to_string(kmer_len));
    // localassm(LASSM_MAX_KMER_LEN, kmer_len, packed_reads_list, histogrammer.ins_avg, histogrammer.ins_stddev, options.qual_offset,
    //           ctgs, alns);
    // end_gasnet_stats();
    // stage_timers.localassm->stop();
    // LOG_MEM("Local assembly completed");
  }
  Timings::wait_pending();
  barrier();
  if (is_debug || options.checkpoint) {
    stage_timers.dump_ctgs->start();
    string contigs_fname("contigs-" + to_string(kmer_len) + ".fasta");
    ctgs.dump_contigs(contigs_fname, 0, "contig_");
    stage_timers.dump_ctgs->stop();
  }
  SLOG(KBLUE "_________________________", KNORM, "\n");
  ctgs.print_stats(500);
  std::chrono::duration<double> loop_t_elapsed = clock_now() - loop_start_t;
  SLOG("\n");
  SLOG(KBLUE, "Completed contig round k = ", kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), " s at ",
       get_current_time(), " (", get_size_str(get_free_mem()), " free memory on node 0)", KNORM, "\n");
  LOG_MEM("After contig round k = " + to_string(kmer_len));
  Timings::wait_pending();
  barrier();
}
