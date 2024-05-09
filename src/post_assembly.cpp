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

#include <bitset>

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
#include "upcxx_utils/three_tier_aggr_store.hpp"

using namespace upcxx_utils;
using namespace std;

struct CtgBaseRange {
  cid_t cid;
  int clen;
  int cstart;
  int cstop;
};

#define BITSET_SIZE 256

static size_t get_target_rank(cid_t cid) { return std::hash<cid_t>{}(cid) % upcxx::rank_n(); }

struct CtgBasesCovered {
  int clen;
  vector<bitset<BITSET_SIZE>> bases_covered;

  int num_bits_set() {
    int num_bits = 0;
    for (auto &b : bases_covered) num_bits += b.count();
    return num_bits;
  }
};

class CtgsCovered {
  using local_ctgs_covered_map_t = HASH_TABLE<cid_t, CtgBasesCovered>;
  using ctgs_covered_map_t = upcxx::dist_object<local_ctgs_covered_map_t>;
  ctgs_covered_map_t ctgs_covered;
  ThreeTierAggrStore<CtgBaseRange> ctg_bases_covered_store;

  void update_ctg_bases(const CtgBaseRange &ctg_base_range) {
    cid_t cid = ctg_base_range.cid;
    int clen = ctg_base_range.clen;
    auto elem = ctgs_covered->find(cid);
    if (elem == ctgs_covered->end()) {
      CtgBasesCovered new_elem = {.clen = clen, .bases_covered = {}};
      new_elem.bases_covered.resize(clen / BITSET_SIZE + (clen % BITSET_SIZE != 0), 0);
      elem = ctgs_covered->insert({cid, new_elem}).first;
    }
    if (elem->second.clen != ctg_base_range.clen) DIE("clens don't match ", elem->second.clen, " != ", ctg_base_range.clen);
    // update the ranges of covered bases
    for (int i = ctg_base_range.cstart; i < ctg_base_range.cstop; i++) {
      elem->second.bases_covered[i / BITSET_SIZE][i % BITSET_SIZE] = 1;
    }
  }

 public:
  CtgsCovered(size_t num_ctgs)
      : ctgs_covered(local_ctgs_covered_map_t{})
      , ctg_bases_covered_store() {
    ctgs_covered->reserve(num_ctgs * 2 + 2000);  // entries for self + distributed calcs
    int64_t mem_to_use = 0.1 * get_free_mem(true) / local_team().rank_n();
    size_t max_store_bytes = std::max(mem_to_use, (int64_t)sizeof(CtgBasesCovered) * 100);
    ctg_bases_covered_store.set_size("Ctg Bases Covered", max_store_bytes);
    ctg_bases_covered_store.set_update_func(
        [&self = *this](CtgBaseRange ctg_base_range) { self.update_ctg_bases(ctg_base_range); });
  }

  void add_ctg_range(cid_t cid, int clen, int cstart, int cstop) {
    CtgBaseRange ctg_base_range = {.cid = cid, .clen = clen, .cstart = cstart, .cstop = cstop};
    auto tgt = get_target_rank(cid);
    if (tgt != rank_me()) {
      ctg_bases_covered_store.update(tgt, ctg_base_range);
    } else {
      update_ctg_bases(ctg_base_range);
    }
  }

  void done_adding() { ctg_bases_covered_store.flush_updates(); }

  size_t get_covered_bases() {
    size_t num_covered_bases = 0;
    for (auto &[key, val] : *ctgs_covered) num_covered_bases += val.num_bits_set();
    return num_covered_bases;
  }
};

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
  CtgsCovered ctgs_covered(ctgs.size());
  auto max_kmer_store = options.max_kmer_store_mb * ONE_MB;
  const int MAX_K = (POST_ASM_ALN_K + 31) / 32 * 32;
  const bool REPORT_CIGAR = true;
  const bool USE_BLASTN_SCORES = true;
  size_t tot_num_reads = 0;
  size_t tot_num_bases = 0;
  size_t tot_reads_aligned = 0;
  size_t tot_bases_aligned = 0;
  size_t tot_proper_pairs = 0;
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
    dist_ofstream sam_ofs(short_name + ".sam");
    KlignTimers aln_timers;
    Alns alns;
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
      stage_timers.alignments->start();
      vector<ReadRecord> read_records(packed_reads.get_local_num_reads());
      fetch_ctg_maps(kmer_ctg_dht, &packed_reads, read_records, KLIGN_SEED_SPACE, aln_timers);
      compute_alns<MAX_K>(&packed_reads, read_records, alns, read_group_id, rlen_limit, REPORT_CIGAR, USE_BLASTN_SCORES,
                          all_num_ctgs, options.klign_rget_buf_size, aln_timers);
      stage_timers.kernel_alns->inc_elapsed(aln_timers.aln_kernel.get_elapsed());
      stage_timers.aln_comms->inc_elapsed(aln_timers.fetch_ctg_maps.get_elapsed() + aln_timers.rget_ctg_seqs.get_elapsed());
      stage_timers.alignments->stop();
      LOG_MEM("Aligned Post Assembly Reads " + short_name);
    }
#ifdef PAF_OUTPUT_FORMAT
    string aln_name("final_assembly-" + short_name + ".paf");
    alns.dump_single_file(aln_name, Alns::Format::PAF);
    SLOG("\n", KBLUE, "PAF alignments can be found at ", options.output_dir, "/", aln_name, KNORM, "\n");
    LOG_MEM("After Post Assembly Alignments Saved");
#elif BLAST6_OUTPUT_FORMAT
    string aln_name("final_assembly-" + short_name + ".b6");
    alns.dump_single_file(aln_name, Alns::Format::BLAST);
    SLOG("\n", KBLUE, "Blast alignments can be found at ", options.output_dir, "/", aln_name, KNORM, "\n");
    LOG_MEM("After Post Assembly Alignments Saved");
#endif
    // the alignments have to be accumulated per read so they can be sorted to keep alignments to each read together
    sort_alns<MAX_K>(alns, aln_timers, packed_reads.get_fname()).wait();
    size_t num_reads_aligned, num_bases_aligned, num_proper_pairs;
    alns.compute_stats(num_reads_aligned, num_bases_aligned, num_proper_pairs);
    tot_reads_aligned += num_reads_aligned;
    tot_bases_aligned += num_bases_aligned;
    tot_proper_pairs += num_proper_pairs;
    for (auto &aln : alns) ctgs_covered.add_ctg_range(aln.cid, aln.clen, aln.cstart, aln.cstop);
    //  Dump 1 file at a time with proper read groups
    stage_timers.dump_alns->start();
    alns.write_sam_alignments(sam_ofs, options.min_ctg_print_len).wait();
    stage_timers.dump_alns->stop();
    LOG_MEM("After Post Assembly SAM Saved");
    stage_timers.compute_ctg_depths->start();
    // compute depths 1 column at a time
    aln_depths.compute_for_read_group(alns, read_group_id);
    stage_timers.compute_ctg_depths->stop();
    LOG_MEM("After Post Assembly Depths Saved");
    sam_ofs.close();
    tot_num_reads += packed_reads.get_local_num_reads();
    tot_num_bases += packed_reads.get_local_bases();
    LOG_MEM("Purged Post Assembly Reads" + short_name);
  }
  Timings::wait_pending();
  ctgs.clear_slices();
  aln_depths.done_computing();
  string fname("final_assembly_depths.txt");
  SLOG_VERBOSE("Writing ", fname, "\n");
  aln_depths.dump_depths(fname, options.reads_fnames);
  ctgs_covered.done_adding();
  auto covered_bases = ctgs_covered.get_covered_bases();

  auto all_covered_bases = reduce_one(covered_bases, op_fast_add, 0).wait();
  auto all_tot_num_reads = reduce_one(tot_num_reads, op_fast_add, 0).wait();
  auto all_tot_num_bases = reduce_one(tot_num_bases, op_fast_add, 0).wait();
  auto all_num_ctgs = reduce_one(ctgs.size(), op_fast_add, 0).wait();
  auto all_ctgs_len = reduce_all(ctgs.get_length(), op_fast_add).wait();
  auto all_reads_aligned = reduce_one(tot_reads_aligned, op_fast_add, 0).wait();
  auto all_bases_aligned = reduce_all(tot_bases_aligned, op_fast_add).wait();
  auto all_proper_pairs = reduce_one(tot_proper_pairs, op_fast_add, 0).wait();
  size_t all_unmapped_bases = 0;

  // every rank needs the avg coverage to calculate the std deviation
  double avg_coverage = (double)all_bases_aligned / all_ctgs_len;
  double sum_diffs = 0;
  for (auto &ctg : ctgs) sum_diffs += pow((double)ctg.depth - avg_coverage, 2.0) * (double)ctg.seq.length();
  auto all_sum_diffs = reduce_one(sum_diffs, op_fast_add, 0).wait();
  auto std_dev_coverage = sqrt(all_sum_diffs / all_ctgs_len);

  SLOG(KBLUE "_________________________", KNORM, "\n");
  SLOG("Alignment statistics\n");
  SLOG("  Reads: ", all_tot_num_reads, "\n");
  SLOG("  Bases: ", all_tot_num_bases, "\n");
  SLOG("  Mapped reads: ", all_reads_aligned, "\n");
  SLOG("  Mapped bases: ", all_bases_aligned, "\n");
  SLOG("  Ref scaffolds: ", all_num_ctgs, "\n");
  SLOG("  Ref bases: ", all_ctgs_len, "\n");
  SLOG("  Percent mapped: ", 100.0 * (double)all_reads_aligned / all_tot_num_reads, "\n");
  SLOG("  Percent bases mapped: ", 100.0 * (double)all_bases_aligned / all_tot_num_bases, "\n");
  //  a proper pair is where both sides of the pair map to the same contig in the correct orientation, less than 32kbp apart
  SLOG("  Percent proper pairs: ", 100.0 * (double)all_proper_pairs * 2.0 / all_tot_num_reads, "\n");
  // average depth per base
  SLOG("  Average coverage: ", avg_coverage, "\n");
  SLOG("  Average coverage with deletions: N/A\n");
  // standard deviation of depth per base
  SLOG("  Standard deviation: ", std_dev_coverage, "\n");
  // this will always be 100% because the contigs are created from the reads in the first place
  SLOG("  Percent scaffolds with any coverage: 100.0\n");
  SLOG("  Percent of reference bases covered: ", 100.0 * (double)all_covered_bases / all_ctgs_len, "\n");

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
  chrono::duration<double> t_elapsed = clock_now() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(), " for ", MHM2_VERSION, "\n");

  SLOG("\n", KBLUE, "Aligned unmerged reads to final assembly. Files can be found in directory \"", options.output_dir,
       "\":\n  \"final_assembly.header.sam\" contains header information",
       "\n  \"*.sam\" files contain alignments per input/read file", "\n  \"final_assembly_depths.text\" contains scaffold depths",
       KNORM, "\n");
  SLOG(KBLUE, "_________________________", KNORM, "\n");
}
