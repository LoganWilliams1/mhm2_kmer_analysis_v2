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
#include "shuffle_reads.hpp"
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

  auto start_t = clock_now();
  // set up output files
  SLOG_VERBOSE("Writing SAM headers\n");
  dist_ofstream sam_header_ofs("final_assembly.header.sam");
  stage_timers.dump_alns->start();
  Alns::write_sam_header(sam_header_ofs, options.reads_fnames, ctgs, options.min_ctg_print_len).wait();
  sam_header_ofs.close();
  stage_timers.dump_alns->stop();
  auto num_read_groups = options.reads_fnames.size();
  SLOG_VERBOSE("Preparing aln depths for post assembly abundance\n");
  AlnDepths aln_depths(ctgs, options.min_ctg_print_len, num_read_groups);
  LOG_MEM("After Post Assembly Ctgs Depths");
  auto max_kmer_store = options.max_kmer_store_mb * ONE_MB;
  const int MAX_K = (POST_ASM_ALN_K + 31) / 32 * 32;
  const bool REPORT_CIGAR = true;
  const bool USE_BLASTN_SCORES = true;
  size_t tot_num_reads = 0;
  size_t tot_num_bases = 0;
  size_t tot_num_ctgs = ctgs.size();
  size_t tot_num_ctg_bases = ctgs.get_length();
  SLOG(KBLUE, "Processing contigs in ", options.post_assm_subsets, " subsets", KNORM, "\n");
  for (int read_group_id = 0; read_group_id < options.reads_fnames.size(); read_group_id++) {
    string &reads_fname = options.reads_fnames[read_group_id];
    SLOG(KBLUE, "_________________________", KNORM, "\n");
    SLOG(KBLUE, "Processing file ", reads_fname, KNORM, "\n");
    vector<string> one_file_list = {reads_fname};
    FastqReaders::open_all_file_blocking(one_file_list);
    PackedReads packed_reads(options.qual_offset, reads_fname, true);
    auto short_name = get_basename(packed_reads.get_fname());
    stage_timers.cache_reads->start();
    packed_reads.load_reads(options.adapter_fname);
    unsigned rlen_limit = packed_reads.get_max_read_len();
    stage_timers.cache_reads->stop();
    LOG_MEM("Read " + short_name);
    barrier();
    ctgs.clear_slices();
    for (int subset_i = 0; subset_i < options.post_assm_subsets; subset_i++) {
      SLOG(KBLUE, "\nContig subset ", subset_i, KNORM, ":\n");
      ctgs.set_next_slice(options.post_assm_subsets);
      int64_t all_num_ctgs = reduce_all(ctgs.size(), op_fast_add).wait();
      stage_timers.build_aln_seed_index->start();
      auto sh_kmer_ctg_dht = build_kmer_ctg_dht<MAX_K>(POST_ASM_ALN_K, max_kmer_store, options.max_rpcs_in_flight, ctgs,
                                                       options.min_ctg_print_len, true);
      stage_timers.build_aln_seed_index->stop();
      auto &kmer_ctg_dht = *sh_kmer_ctg_dht;
      LOG_MEM("After Post Assembly Built Kmer Seeds");
      Alns alns;
      stage_timers.alignments->start();
      KlignTimers aln_timers;
      vector<ReadRecord> read_records(packed_reads.get_local_num_reads());
      fetch_ctg_maps(kmer_ctg_dht, &packed_reads, read_records, KLIGN_SEED_SPACE, aln_timers);
      compute_alns<MAX_K>(&packed_reads, read_records, alns, read_group_id, rlen_limit, REPORT_CIGAR, USE_BLASTN_SCORES,
                          all_num_ctgs, options.klign_rget_buf_size, aln_timers);
      stage_timers.kernel_alns->inc_elapsed(aln_timers.aln_kernel.get_elapsed());
      stage_timers.aln_comms->inc_elapsed(aln_timers.fetch_ctg_maps.get_elapsed() + aln_timers.rget_ctg_seqs.get_elapsed());
      stage_timers.alignments->stop();
      LOG_MEM("Aligned Post Assembly Reads " + short_name);

#ifdef PAF_OUTPUT_FORMAT
      string aln_name("final_assembly-" + short_name + ".paf");
      alns.dump_single_file(aln_name, Alns::Format::PAF);
      SLOG("\n", KBLUE, "PAF alignments can be found at ", options.output_dir, "/", aln_name, KNORM, "\n");
#elif BLAST6_OUTPUT_FORMAT
      string aln_name("final_assembly-" + short_name + ".b6");
      alns.dump_single_file(aln_name, Alns::Format::BLAST);
      SLOG("\n", KBLUE, "Blast alignments can be found at ", options.output_dir, "/", aln_name, KNORM, "\n");
#endif
      LOG_MEM("After Post Assembly Alignments Saved");
      // Dump 1 file at a time with proper read groups
      stage_timers.dump_alns->start();
      dist_ofstream sam_ofs(short_name + ".sam");
      alns.write_sam_alignments(sam_ofs, options.min_ctg_print_len).wait();
      sam_ofs.close();
      stage_timers.dump_alns->stop();

      LOG_MEM("After Post Assembly SAM Saved");

      stage_timers.compute_ctg_depths->start();
      // compute depths 1 column at a time
      aln_depths.compute_for_read_group(alns, read_group_id);
      stage_timers.compute_ctg_depths->stop();
      LOG_MEM("After Post Assembly Depths Saved");
    }
    tot_num_reads += packed_reads.get_local_num_reads();
    tot_num_bases += packed_reads.get_local_bases();
    LOG_MEM("Purged Post Assembly Reads" + short_name);
  }
  Timings::wait_pending();

  auto all_tot_num_reads = reduce_one(tot_num_reads, op_fast_add, 0).wait();
  auto all_tot_num_bases = reduce_one(tot_num_bases, op_fast_add, 0).wait();
  auto all_num_ctgs = reduce_one(tot_num_ctgs, op_fast_add, 0).wait();
  auto all_ctgs_len = reduce_one(tot_num_ctg_bases, op_fast_add, 0).wait();

  SLOG(KBLUE "_________________________", KNORM, "\n");
  SLOG("Alignment statistics\n");
  SLOG("  Reads: ", all_tot_num_reads, "\n");
  SLOG("  Bases: ", all_tot_num_bases, "\n");
  SLOG("  Mapped reads*: ", 0, "\n");
  SLOG("  Mapped bases*: ", 0, "\n");
  SLOG("  Ref scaffolds: ", all_num_ctgs, "\n");
  SLOG("  Ref bases: ", all_ctgs_len, "\n");
  /*
  Reads: 21470321354
  Mapped reads: 20362450458
  Mapped bases: 3030227153794
  Ref scaffolds: 55342847
  Ref bases: 74970251022

  Percent mapped: 94.840
  Percent proper pairs: 70.095
  Average coverage: 40.419
  Average coverage with deletions: 40.412
  Standard deviation: 148.923
  Percent scaffolds with any coverage: 100.00
  Percent of reference bases covered: 99.76
  */

  stage_timers.alignments->inc_elapsed(stage_timers.build_aln_seed_index->get_elapsed());

  SLOG(KBLUE "_________________________", KNORM, "\n");
  SLOG("Stage timing:\n");
  SLOG("    ", stage_timers.cache_reads->get_final(), "\n");
  SLOG("    ", stage_timers.alignments->get_final(), "\n");
  SLOG("      -> ", stage_timers.build_aln_seed_index->get_final(), "\n");
  SLOG("      -> ", stage_timers.kernel_alns->get_final(), "\n");
  SLOG("      -> ", stage_timers.aln_comms->get_final(), "\n");
  SLOG("    ", stage_timers.compute_ctg_depths->get_final(), "\n");
  SLOG("    ", stage_timers.dump_alns->get_final(), "\n");
  // if (options.shuffle_reads) SLOG("    ", stage_timers.shuffle_reads->get_final(), "\n");
  SLOG("    FASTQ total read time: ", FastqReader::get_io_time(), "\n");
  SLOG(KBLUE "_________________________", KNORM, "\n");
  std::chrono::duration<double> t_elapsed = clock_now() - start_t;
  SLOG("Finished in ", std::setprecision(2), std::fixed, t_elapsed.count(), " s at ", get_current_time(), " for ", MHM2_VERSION,
       "\n");

  SLOG("\n", KBLUE, "Aligned unmerged reads to final assembly. Files can be found in directory \"", options.output_dir,
       "\":\n  \"final_assembly.header.sam\" contains header information",
       "\n  \"*.sam\" files contain alignments per input/read file", "\n  \"final_assembly_depths.text\" contains scaffold depths",
       KNORM, "\n");
  string fname("final_assembly_depths.txt");
  SLOG_VERBOSE("Writing ", fname, "\n");
  aln_depths.done_computing();
  aln_depths.dump_depths(fname, options.reads_fnames);

  SLOG(KBLUE, "_________________________", KNORM, "\n");
}
