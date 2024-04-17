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

#include "post_assembly.hpp"

#include "aln_depths.hpp"
#include "fastq.hpp"
#include "gasnet_stats.hpp"
#include "klign.hpp"
#include "packed_reads.hpp"
#include "stage_timers.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"

using namespace upcxx_utils;

using std::shared_ptr;
using std::string;
using std::vector;

void post_assembly(Contigs &ctgs, Options &options) {
  SLOG(KBLUE, "_________________________", KNORM, "\n");
  SLOG(KBLUE, "Post processing", KNORM, "\n\n");
  LOG_MEM("Starting Post Assembly");

  // set up output files
  SLOG_VERBOSE("Writing SAM headers\n");
  // shared_ptr<dist_ofstream> sh_sam_of = make_shared<dist_ofstream>("final_assembly.sam");
  dist_ofstream sam_of("final_assembly.sam");
  upcxx::future<> fut_sam = make_future();
  fut_sam = Alns::write_sam_header(sam_of, options.reads_fnames, ctgs, options.min_ctg_print_len);
  auto num_read_groups = options.reads_fnames.size();
  SLOG_VERBOSE("Preparing aln depths for post assembly abundance\n");
  AlnDepths aln_depths(ctgs, options.min_ctg_print_len, num_read_groups);
  LOG_MEM("After Post Assembly Ctgs Depths");

  int num_subsets = 5;
  int subset_size = ctgs.size() / num_subsets;
  size_t num_ctgs = ctgs.size();

  for (int subset_i = 0; subset_i < num_subsets; subset_i++) {
    size_t begin_idx = subset_i * subset_size;
    size_t end_idx = (subset_i == num_subsets - 1 ? num_ctgs : (subset_i + 1) * subset_size);
    SLOG("num_ctgs ", num_ctgs, " subset begin ", begin_idx, " subset end ", end_idx, "\n");
    ctgs.set_range(begin_idx, end_idx);

    // build kmer_ctg_dht
    auto max_kmer_store = options.max_kmer_store_mb * ONE_MB;
    int64_t all_num_ctgs = reduce_all(ctgs.size(), op_fast_add).wait();
    const int MAX_K = (POST_ASM_ALN_K + 31) / 32 * 32;
    auto sh_kmer_ctg_dht = build_kmer_ctg_dht<MAX_K>(POST_ASM_ALN_K, max_kmer_store, options.max_rpcs_in_flight, ctgs,
                                                     options.min_ctg_print_len, true);
    auto &kmer_ctg_dht = *sh_kmer_ctg_dht;
    LOG_MEM("After Post Assembly Built Kmer Seeds");

    KlignTimers timers;
    unsigned rlen_limit = 0;
    int read_group_id = 0;
    for (auto &reads_fname : options.reads_fnames) {
      SLOG(KBLUE, "_________________________", KNORM, "\n");
      SLOG(KBLUE, "Processing file ", reads_fname, KNORM, "\n");
      SLOG("\n");
      vector<string> one_file_list;
      one_file_list.push_back(reads_fname);
      FastqReaders::open_all_file_blocking(one_file_list);
      PackedReadsList packed_reads_list;
      packed_reads_list.push_back(new PackedReads(options.qual_offset, reads_fname, true));
      auto packed_reads = packed_reads_list[0];
      auto short_name = get_basename(packed_reads->get_fname());

      stage_timers.cache_reads->start();
      packed_reads->load_reads();
      auto max_read_len = packed_reads->get_max_read_len();
      rlen_limit = max(rlen_limit, max_read_len);
      DBG("max_read_len=", max_read_len, " rlen_limit=", rlen_limit, "\n");
      stage_timers.cache_reads->stop();
      max_read_len = reduce_all(max_read_len, op_fast_max).wait();
      LOG_MEM("Read " + short_name);

      stage_timers.alignments->start();

      bool report_cigar = true;
      int kmer_len = POST_ASM_ALN_K;

      Alns alns;
      vector<ReadRecord> read_records(packed_reads->get_local_num_reads());
      fetch_ctg_maps(kmer_ctg_dht, packed_reads, read_records, KLIGN_SEED_SPACE, timers);
      int max_clen = reduce_all(ctgs.get_max_clen(), op_fast_max).wait();
      compute_alns<MAX_K>(packed_reads, read_records, alns, read_group_id, rlen_limit, max_clen, report_cigar, true, all_num_ctgs,
                          options.klign_rget_buf_size, timers);
      stage_timers.alignments->stop();
      LOG_MEM("Aligned Post Assembly Reads " + short_name);

      delete packed_reads;
      LOG_MEM("Purged Post Assembly Reads" + short_name);

#ifdef PAF_OUTPUT_FORMAT
      string aln_name("final_assembly-" + short_name + ".paf");
      alns.dump_single_file(aln_name);
      SLOG("\n", KBLUE, "PAF alignments can be found at ", options.output_dir, "/", aln_name, KNORM, "\n");
#elif BLAST6_OUTPUT_FORMAT
      string aln_name("final_assembly-" + short_name + ".b6");
      alns.dump_single_file(aln_name);
      SLOG("\n", KBLUE, "Blast alignments can be found at ", options.output_dir, "/", aln_name, KNORM, "\n");
#endif
      LOG_MEM("After Post Assembly Alignments Saved");
      // Dump 1 file at a time with proper read groups
      auto fut_flush = alns.write_sam_alignments(sam_of, options.min_ctg_print_len);
      fut_sam = when_all(fut_sam, fut_flush);

      LOG_MEM("After Post Assembly SAM Saved");

      stage_timers.compute_ctg_depths->start();
      // compute depths 1 column at a time
      aln_depths.compute_for_read_group(alns, read_group_id);
      stage_timers.compute_ctg_depths->stop();
      LOG_MEM("After Post Assembly Depths Saved");

      read_group_id++;
    }
    stage_timers.kernel_alns->inc_elapsed(timers.aln_kernel.get_elapsed());
    Timings::wait_pending();
  }

  fut_sam.wait();
  sam_of.close();
  SLOG("\n", KBLUE, "Aligned unmerged reads to final assembly: SAM file can be found at ", options.output_dir,
       "/final_assembly.sam", KNORM, "\n");

  string fname("final_assembly_depths.txt");
  SLOG_VERBOSE("Writing ", fname, "\n");
  aln_depths.done_computing();
  aln_depths.dump_depths(fname, options.reads_fnames);
  SLOG(KBLUE, "Contig depths (abundances) can be found at ", options.output_dir, "/", fname, KNORM, "\n");

  SLOG(KBLUE, "_________________________", KNORM, "\n");
}
