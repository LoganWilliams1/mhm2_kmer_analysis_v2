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

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <string_view>
#include <unordered_set>
#include <deque>

#include <algorithm>
#include <forward_list>
#include <iterator>
#include <iostream>
#include <random>
#include <thread>
#include <upcxx/upcxx.hpp>

#include "klign.hpp"
#include "kmer.hpp"
#include "contigs.hpp"
#include "ssw.hpp"
#include "utils.hpp"
#include "zstr.hpp"
#include "aligner_cpu.hpp"

#include "upcxx_utils/limit_outstanding.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/thread_pool.hpp"
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "upcxx_utils/timers.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

void init_aligner(int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty, int ambiguity_penalty,
                  int rlen_limit, bool compute_cigar);
void cleanup_aligner();
void kernel_align_block(CPUAligner &cpu_aligner, vector<Aln> &kernel_alns, vector<string> &ctg_seqs, vector<string> &read_seqs,
                        Alns *alns, future<> &active_kernel_fut, int read_group_id, int max_clen, int max_rlen,
                        KlignTimers &klign_timers);

template <int MAX_K>
struct KmersReadsBuffer {
  vector<Kmer<MAX_K>> kmers;
  vector<ReadRecordPtr> read_records;

  void add(const Kmer<MAX_K> &kmer, ReadRecord *read_record, int read_offset, bool read_is_rc) {
    kmers.push_back(kmer);
    read_records.push_back({read_record, read_offset, read_is_rc});
  }

  size_t size() const { return kmers.size(); }

  static size_t get_size_per() { return sizeof(Kmer<MAX_K>) + sizeof(ReadRecordPtr); }

  void clear() {
    kmers.clear();
    read_records.clear();
  }

  bool empty() const { return kmers.empty(); }
  void reserve(size_t sz) {
    kmers.reserve(sz);
    read_records.reserve(sz);
  }
};  // struct KmersReadsBuffer

struct RgetRequest {
  CtgLoc ctg_loc;
  string rname;
  string read_seq;
  int rstart;
  int cstart;
  char orient;
  int overlap_len;
  int read_group_id;
};

template <int MAX_K>
class KmerCtgDHT {
  using local_kmer_map_t = HASH_TABLE<Kmer<MAX_K>, vector<CtgLoc>>;
  using kmer_map_t = dist_object<local_kmer_map_t>;
  kmer_map_t kmer_map;
  vector<global_ptr<char>> global_ctg_seqs;
  ThreeTierAggrStore<KmerAndCtgLoc<MAX_K>> kmer_store;
  dist_object<int64_t> num_dropped_seed_to_ctgs;

 public:
  size_t kmer_seed_lookups = 0;
  size_t unique_kmer_seed_lookups = 0;
  unsigned kmer_len;
  bool allow_multi_kmers = false;

  KmerCtgDHT(int max_store_size, int max_rpcs_in_flight, bool allow_multi_kmers, uint64_t num_contig_kmers = 0)
      : kmer_map(local_kmer_map_t{})
      , global_ctg_seqs({})
      , kmer_store()
      , allow_multi_kmers(allow_multi_kmers)
      , num_dropped_seed_to_ctgs(0) {
    if (allow_multi_kmers) SLOG_VERBOSE("Finding multiple kmer to contig mappings\n");
    kmer_len = Kmer<MAX_K>::get_k();
    kmer_store.set_size("insert ctg seeds", max_store_size, max_rpcs_in_flight, num_contig_kmers);
    kmer_store.set_update_func([&kmer_map = this->kmer_map, &num_dropped_seed_to_ctgs = this->num_dropped_seed_to_ctgs,
                                allow_multi_kmers](KmerAndCtgLoc<MAX_K> kmer_and_ctg_loc) {
      CtgLoc ctg_loc = kmer_and_ctg_loc.ctg_loc;
      auto it = kmer_map->find(kmer_and_ctg_loc.kmer);
      if (it == kmer_map->end()) {
        // always insert if not found
        it = kmer_map->insert({kmer_and_ctg_loc.kmer, {ctg_loc}}).first;
      } else {
        if (allow_multi_kmers) {
          // only add if we haven't hit the threshold to prevent excess memory usage and load imbalance here
          if (it->second.size() < KLIGN_MAX_CTGS_PER_KMER) it->second.push_back(ctg_loc);
        } else {
          // there are conflicts so don't allow any kmer mappings. This improves the assembly when scaffolding k is smaller than
          // the final contigging k, e.g. sk=33
          it->second.clear();
          (*num_dropped_seed_to_ctgs)++;
        }
      }
    });
  }

  void clear() {
    LOG("aggregated kmer seed lookups ", perc_str(kmer_seed_lookups - unique_kmer_seed_lookups, kmer_seed_lookups), ", total ",
        kmer_seed_lookups, "\n");
    for (auto &gptr : global_ctg_seqs) upcxx::deallocate(gptr);
    local_kmer_map_t().swap(*kmer_map);  // release all memory
    kmer_store.clear();
  }

  ~KmerCtgDHT() { clear(); }

  global_ptr<char> add_ctg_seq(string seq) {
    auto seq_gptr = upcxx::allocate<char>(seq.length() + 1);
    global_ctg_seqs.push_back(seq_gptr);  // remember to dealloc!
    strcpy(seq_gptr.local(), seq.c_str());
    return seq_gptr;
  }

  intrank_t get_target_rank(const Kmer<MAX_K> &kmer) const { return std::hash<Kmer<MAX_K>>{}(kmer) % rank_n(); }

  future<size_t> fut_get_num_kmers(bool all = false) {
    auto &pr = Timings::get_promise_reduce();
    if (!all) return pr.reduce_one(kmer_map->size(), op_fast_add, 0);
    return pr.reduce_all(kmer_map->size(), op_fast_add);
  }

  future<int64_t> fut_get_num_dropped_seed_to_ctgs(bool all = false) {
    auto &pr = Timings::get_promise_reduce();
    if (!all) return pr.reduce_one(*num_dropped_seed_to_ctgs, op_fast_add, 0);
    return pr.reduce_all(*num_dropped_seed_to_ctgs, op_fast_add);
  }

  void add_kmer(const Kmer<MAX_K> &kmer_fw, CtgLoc &ctg_loc) {
    Kmer<MAX_K> kmer_rc = kmer_fw.revcomp();
    ctg_loc.is_rc = false;
    const Kmer<MAX_K> *kmer_lc = &kmer_fw;
    if (kmer_rc < kmer_fw) {
      kmer_lc = &kmer_rc;
      ctg_loc.is_rc = true;
    }
    KmerAndCtgLoc<MAX_K> kmer_and_ctg_loc = {*kmer_lc, ctg_loc};
    kmer_store.update(get_target_rank(*kmer_lc), kmer_and_ctg_loc);
  }

  void flush_add_kmers() {
    BarrierTimer timer(__FILEFUNC__, false);  // barrier on exit, not entrance
    kmer_store.flush_updates();
    kmer_store.clear();
    size_t max_ctgs = 0;
    // determine max number of ctgs mapped to by a single kmer
    for (auto &elem : *kmer_map) {
      auto num_ctgs_for_kmer = elem.second.size();
      max_ctgs = ::max(max_ctgs, num_ctgs_for_kmer);
    }
    auto &pr = Timings::get_promise_reduce();
    auto fut_reduce = when_all(pr.reduce_one(max_ctgs, op_fast_max, 0));
    auto fut_report = fut_reduce.then([](auto all_max_ctgs) {
      if (all_max_ctgs > 1) SLOG_VERBOSE("Max contigs mapped by a single kmer: ", all_max_ctgs, "\n");
    });
    Timings::set_pending(fut_report);
  }

  future<vector<CtgLocAndKmerIdx>> get_ctgs_with_kmers(int target_rank, const vector<Kmer<MAX_K>> &kmers) {
    DBG_VERBOSE("Sending request for ", kmers.size(), " to ", target_rank, "\n");
    auto fut_rpc = rpc(
        target_rank,
        [](const vector<Kmer<MAX_K>> &kmers, kmer_map_t &kmer_map) {
          BaseTimer t("with get_ctgs_with_kmers rpc");
          t.start();
          vector<CtgLocAndKmerIdx> ctg_locs;
          ctg_locs.reserve(kmers.size());
          for (int i = 0; i < kmers.size(); i++) {
            auto &kmer = kmers[i];
            assert(kmer.is_least());
            assert(kmer.is_valid());
            const auto it = kmer_map->find(kmer);
            if (it == kmer_map->end()) continue;
            for (auto &ctg_loc : it->second) ctg_locs.push_back({ctg_loc, i});
          }
          t.stop();
          // For Issue137 tracking
          static int num_calls = 0;
          if (ctg_locs.size() * sizeof(CtgLocAndKmerIdx) > KLIGN_MAX_RPC_MESSAGE_SIZE || rand() % 1000 == 0)
            LOG("Received ", kmers.size(), " kmers ", get_size_str(kmers.size() * sizeof(Kmer<MAX_K>)), " responding with ",
                ctg_locs.size(), " ctg_locs ", get_size_str(ctg_locs.size() * sizeof(CtgLocAndKmerIdx)), " in ", t.get_elapsed(),
                " s num_calls=", num_calls, "\n");
          num_calls++;
          return ctg_locs;
        },
        kmers, kmer_map);
    return fut_rpc;
  }

  void dump_ctg_kmers() {
    BarrierTimer timer(__FILEFUNC__, false);  // barrier on exit not entrance
    string dump_fname = "ctg_kmers-" + to_string(kmer_len) + ".txt.gz";
    get_rank_path(dump_fname, rank_me());
    zstr::ofstream dump_file(dump_fname);
    ostringstream out_buf;
    ProgressBar progbar(kmer_map->size(), "Dumping kmers to " + dump_fname);
    int64_t i = 0;
    for (auto &elem : *kmer_map) {
      for (auto &ctg_loc : elem.second) {
        out_buf << elem.first << " " << ctg_loc.cid << " " << ctg_loc.clen << " " << ctg_loc.depth << " " << ctg_loc.pos << " "
                << ctg_loc.is_rc << "\n";
        i++;
        if (!(i % 1000)) {
          dump_file << out_buf.str();
          out_buf = ostringstream();
        }
      }
      progbar.update();
    }
    if (!out_buf.str().empty()) dump_file << out_buf.str();
    dump_file.close();
    auto fut_rep = when_all(progbar.set_done(), this->fut_get_num_kmers()).then([](auto tot_num_kmers) {
      SLOG_VERBOSE("Dumped ", tot_num_kmers, " kmers\n");
    });
    Timings::set_pending(fut_rep);
  }

  void build(Contigs &ctgs, unsigned min_ctg_len) {
    BarrierTimer timer(__FILEFUNC__);
    int64_t num_kmers = 0;
    ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
    // estimate and reserve room in the local map to avoid excessive reallocations
    int64_t est_num_kmers = 0;
    size_t ctg_seq_lengths = 0, min_len_ctgs = 0;
    for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
      auto ctg = it;
      auto len = ctg->seq.length();
      ctg_seq_lengths += len;
      min_len_ctgs++;
      if (len < min_ctg_len) continue;
      est_num_kmers += len - kmer_len + 1;
    }
    est_num_kmers = upcxx::reduce_all(est_num_kmers, upcxx::op_fast_add).wait();
    auto tot_ctg_lengths = reduce_one(ctg_seq_lengths, op_fast_add, 0).wait();
    auto tot_num_ctgs = reduce_one(min_len_ctgs, op_fast_add, 0).wait();
    auto my_reserve = 1.2 * est_num_kmers / rank_n() + 2000;  // 120% to keep the map fast
    SLOG_VERBOSE("Estimated ", est_num_kmers, " contig ", kmer_len, "-kmers for ", tot_num_ctgs, " contigs with total len ",
                 tot_ctg_lengths, "\n");
    auto my_required_mem = my_reserve * (sizeof(typename local_kmer_map_t::value_type) +
                                         sizeof(CtgLoc));  // 1 map key/value entry + 1 CtgLoc within the vector
    LOG("Reserving ", my_reserve, " for my entries at approx ", get_size_str(my_required_mem * upcxx::local_team().rank_n()),
        " per node for the local hashtable\n");
    LOG_MEM("Before reserving entries for ctg kmers");
    kmer_map->reserve(my_reserve);
    global_ctg_seqs.reserve(ctgs.size());
    LOG_MEM("After reserving entries for ctg kmers");
    vector<Kmer<MAX_K>> kmers;
    for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
      auto ctg = it;
      progbar.update();
      if (ctg->seq.length() < min_ctg_len) continue;
      global_ptr<char> seq_gptr = add_ctg_seq(ctg->seq);
      CtgLoc ctg_loc = {.cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length(), .depth = (float)ctg->depth};
      Kmer<MAX_K>::get_kmers(kmer_len, string_view(ctg->seq.data(), ctg->seq.size()), kmers, true);
      num_kmers += kmers.size();
      for (unsigned i = 0; i < kmers.size(); i++) {
        ctg_loc.pos = i;
        if (ctg_loc.pos >= ctg_loc.clen) DIE("ctg_loc.pos ", ctg_loc.pos, " ctg_loc.clen ", ctg_loc.clen);
        if (!kmers[i].is_valid()) continue;
        add_kmer(kmers[i], ctg_loc);
      }
      progress();
    }
    flush_add_kmers();
    auto fut = progbar.set_done();
    auto &pr = Timings::get_promise_reduce();

    auto fut_reduce =
        when_all(fut, pr.reduce_one(num_kmers, op_fast_add, 0), fut_get_num_kmers(), fut_get_num_dropped_seed_to_ctgs());

    auto fut_report = fut_reduce.then([=](auto tot_num_kmers, auto num_kmers_in_ht, auto num_dropped_seed_to_ctgs) {
      LOG("Estimated room for ", my_reserve, " my final count ", kmer_map->size(), "\n");
      SLOG_VERBOSE("Total contigs >= ", min_ctg_len, ": ", tot_num_ctgs, " seq_length: ", tot_ctg_lengths, "\n");
      SLOG_VERBOSE("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");

      if (num_dropped_seed_to_ctgs)
        SLOG_VERBOSE("For k = ", kmer_len, " dropped ", num_dropped_seed_to_ctgs, " non-unique seed-to-contig mappings (",
                     setprecision(2), fixed, (100.0 * num_dropped_seed_to_ctgs / tot_num_kmers), "%)\n");
    });
    Timings::set_pending(fut_report);
  }

  void purge_high_count_seeds() {
    BarrierTimer timer(__FILEFUNC__);
    // FIXME: Hack for Issue137 to be fixed robustly later
    size_t num_purged = 0;
    size_t num_ctg_locs_purged = 0;
    size_t max_ctg_locs = 0;
    double total_depth = 0.0;
    double total_purged_depth = 0.0;
    double max_depth = 0.0;
    vector<Kmer<MAX_K>> to_be_purged;
    for (auto &elem : *kmer_map) {
      auto &[kmer, ctg_locs] = elem;
      double depth = 0.0;
      for (auto &ctg_loc : ctg_locs) {
        depth += ctg_loc.depth;
      }
      max_ctg_locs = max(max_ctg_locs, ctg_locs.size());
      max_depth = max(max_depth, depth);
      total_depth += depth;
      // we should never actually exceed this threshold because we limit the number of insertions
      if (KLIGN_MAX_CTGS_PER_KMER > 0 && ctg_locs.size() >= KLIGN_MAX_CTGS_PER_KMER) {
        DBG("Removing kmer with ", ctg_locs.size(), " ctg_loc hits aggregated depth=", depth, ": ", kmer.to_string(), "\n");
        to_be_purged.push_back(kmer);
        num_purged++;
        num_ctg_locs_purged += ctg_locs.size();
        total_purged_depth += depth;
      }
    }

    for (auto &kmer : to_be_purged) {
      auto it = kmer_map->find(kmer);
      if (it == kmer_map->end()) DIE("Could not find kmer with too many ctg_locs!");
      kmer_map->erase(it);
    }

    auto &pr = Timings::get_promise_reduce();
    auto fut_reduce = when_all(pr.reduce_one(num_purged, op_fast_add, 0), pr.reduce_one(num_ctg_locs_purged, op_fast_add, 0),
                               pr.reduce_one(max_ctg_locs, op_fast_max, 0), pr.reduce_one(total_depth, op_fast_add, 0),
                               pr.reduce_one(total_purged_depth, op_fast_add, 0), pr.msm_reduce_one(kmer_map->size(), 0),
                               pr.msm_reduce_one(num_purged, 0), pr.msm_reduce_one(max_ctg_locs, 0),
                               pr.msm_reduce_one(total_purged_depth), pr.msm_reduce_one(max_depth, 0));
    auto fut_report = fut_reduce.then([=](auto all_num_purged, auto all_ctg_locs_purged, auto global_max_ctg_locs,
                                          auto all_total_depth, auto all_total_purged_depth, auto msm_kmers_remaining,
                                          auto msm_kmers_purged, auto msm_max_ctg_locs, auto msm_purged_depth, auto msm_max_depth) {
      LOG("Purge of high count ctg kmer seeds: num_purged=", num_purged, " num_ctg_locs_purged=", num_ctg_locs_purged,
          " max_ctg_locs=", max_ctg_locs, " total_depth=", total_depth, " total_purged_depth=", total_purged_depth, "\n");
      SLOG_VERBOSE("All purged ", all_num_purged, " high count ctg kmer seeds.  total_ctg_locs_purged=", all_ctg_locs_purged,
                   " max_ctg_locs=", global_max_ctg_locs, " all_total_depth=", all_total_depth,
                   " all_purged_depth=", all_total_purged_depth, " seeds_remaining=", msm_kmers_remaining.to_string(),
                   " seeds_purged=", msm_kmers_purged.to_string(), " max_ctg_locs=", msm_max_ctg_locs.to_string(),
                   " purged_depth=", msm_purged_depth.to_string(), " max_depth=", msm_max_depth.to_string(), "\n");
    });
    Timings::set_pending(fut_report);
  }
};  // class KmerCtgDHT

static tuple<int, int, int> get_start_positions(int kmer_len, int clen, int pos_in_ctg, int pos_in_read, int rlen) {
  // calculate available bases before and after the seeded kmer
  int ctg_bases_left_of_kmer = pos_in_ctg;
  int ctg_bases_right_of_kmer = clen - ctg_bases_left_of_kmer - kmer_len;
  int read_bases_left_of_kmer = pos_in_read;
  int read_bases_right_of_kmer = rlen - kmer_len - pos_in_read;
  int left_of_kmer = min(read_bases_left_of_kmer, ctg_bases_left_of_kmer);
  int right_of_kmer = min(read_bases_right_of_kmer, ctg_bases_right_of_kmer);

  int cstart = pos_in_ctg - left_of_kmer;
  int rstart = pos_in_read - left_of_kmer;
  int overlap_len = left_of_kmer + kmer_len + right_of_kmer;
  return {cstart, rstart, overlap_len};
}

class Aligner {
  int64_t num_alns;
  int64_t num_perfect_alns;
  int kmer_len;

  vector<Aln> kernel_alns;
  vector<string> ctg_seqs;
  vector<string> read_seqs;

  future<> active_kernel_fut;

  int64_t max_clen = 0;
  int64_t max_rlen = 0;
  CPUAligner cpu_aligner;

  int64_t ctg_bytes_fetched = 0;
  int64_t rget_calls = 0;
  int64_t large_rget_calls = 0;
  int64_t large_rget_bytes = 0;
  int64_t local_ctg_fetches = 0;
  int64_t remote_ctg_fetches = 0;

  Alns *alns;

 private:
  vector<vector<RgetRequest>> rget_requests;

  void align_read(const string &rname, int64_t cid, const string_view &rseq, const string_view &ctg_seq, int rstart, int cstart,
                  char orient, int overlap_len, int read_group_id, KlignTimers &timers) {
    num_alns++;
    size_t clen = ctg_seq.length();
    size_t rlen = rseq.length();
    assert(clen > cstart);
    assert(rlen > rstart);
    assert(overlap_len <= clen && overlap_len <= rlen);
    string_view cseq = string_view(ctg_seq.data() + cstart, overlap_len);
    if (cseq.compare(0, overlap_len, rseq, rstart, overlap_len) == 0) {
      num_perfect_alns++;
      int rstop = rstart + overlap_len;
      int cstop = cstart + overlap_len;
      if (orient == '-') switch_orient(rstart, rstop, rlen);
      int score1 = overlap_len * cpu_aligner.ssw_aligner.get_match_score();
      Aln aln(rname, cid, rstart, rstop, rlen, cstart, cstop, clen, orient, score1, 0, 0, read_group_id);
      assert(aln.is_valid());
      if (cpu_aligner.ssw_filter.report_cigar) aln.set_sam_string(to_string(overlap_len) + "=");
      alns->add_aln(aln);
    } else {
      max_clen = max((int64_t)cseq.size(), max_clen);
      max_rlen = max((int64_t)rseq.size(), max_rlen);
      int64_t num_kernel_alns = kernel_alns.size() + 1;
      int64_t tot_mem_est = num_kernel_alns * (max_clen + max_rlen + 2 * sizeof(int) + 5 * sizeof(short));
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      // this is also the way it's done in meraligner
      kernel_alns.emplace_back(rname, cid, rlen, cstart, clen, orient);
      ctg_seqs.emplace_back(cseq);
      read_seqs.emplace_back(rseq);
      if (num_kernel_alns >= KLIGN_GPU_BLOCK_SIZE) {
        kernel_align_block(cpu_aligner, kernel_alns, ctg_seqs, read_seqs, alns, active_kernel_fut, read_group_id, max_clen,
                           max_rlen, timers);
        clear_aln_bufs();
      }
    }
  }

  void clear() {
    if (kernel_alns.size() || !active_kernel_fut.is_ready())
      DIE("clear called with alignments in the buffer or active kernel - was flush_remaining called before destrutor?\n");
    clear_aln_bufs();
  }

  void do_rget_irregular(int target, KlignTimers &timers) {
    assert(!upcxx::in_progress());
    assert(upcxx::master_persona().active_with_caller());
    deque<pair<global_ptr<char>, size_t>> src;
    deque<pair<char *, size_t>> dest;
    HASH_TABLE<cid_t, string> ctgs_fetched;
    auto sz = rget_requests[target].size();
    ctgs_fetched.reserve(sz);
    upcxx::future<> fut_chain = make_future();
    for (auto &req : rget_requests[target]) {
      auto clen = req.ctg_loc.clen;
      auto it = ctgs_fetched.find(req.ctg_loc.cid);
      if (it == ctgs_fetched.end()) {
        it = ctgs_fetched.insert({req.ctg_loc.cid, string(clen, ' ')}).first;
        if (KLIGN_MAX_IRREGULAR_RGET > 0 && clen >= KLIGN_MAX_IRREGULAR_RGET) {
          // issue a normal rget
          auto fut = rget(req.ctg_loc.seq_gptr, const_cast<char *>(it->second.data()), clen);
          fut_chain = when_all(fut, fut_chain);
          large_rget_calls++;
          large_rget_bytes += clen;
        } else {
          src.push_back({req.ctg_loc.seq_gptr, clen});
          dest.push_back({const_cast<char *>(it->second.data()), clen});
          ctg_bytes_fetched += clen;
        }
        progress();
      }
    }
    assert(src.size() == dest.size());
    // SLOG_VERBOSE("Using rget_irregular to fetch ", perc_str(ctgs_fetched.size(), rget_requests[target].size()), " contigs\n");
    timers.rget_ctg_seqs.start();
    if (src.size()) rget_irregular(src.begin(), src.end(), dest.begin(), dest.end()).wait();
    fut_chain.wait();
    timers.rget_ctg_seqs.stop();
    rget_calls++;
    for (auto &req : rget_requests[target]) {
      auto cid = req.ctg_loc.cid;
      auto it = ctgs_fetched.find(cid);
      if (it == ctgs_fetched.end()) DIE("Could not find the sequence for the contig ", cid);
      string &ctg_seq = it->second;
      align_read(req.rname, cid, req.read_seq, ctg_seq, req.rstart, req.cstart, req.orient, req.overlap_len, req.read_group_id,
                 timers);
      progress();
    }
    rget_requests[target].clear();
  }

 public:
  Aligner(int kmer_len, Alns &alns, int rlen_limit, bool report_cigar, bool use_blastn_scores)
      : num_alns(0)
      , num_perfect_alns(0)
      , kmer_len(kmer_len)
      , kernel_alns({})
      , ctg_seqs({})
      , read_seqs({})
      , active_kernel_fut(make_future())
      , cpu_aligner(report_cigar, use_blastn_scores)
      , alns(&alns)
      , rget_requests(rank_n())
      , rget_calls(0)
      , large_rget_calls(0)
      , large_rget_bytes(0)
      , local_ctg_fetches(0)
      , remote_ctg_fetches(0) {
    init_aligner((int)cpu_aligner.ssw_aligner.get_match_score(), (int)cpu_aligner.ssw_aligner.get_mismatch_penalty(),
                 (int)cpu_aligner.ssw_aligner.get_gap_opening_penalty(), (int)cpu_aligner.ssw_aligner.get_gap_extending_penalty(),
                 (int)cpu_aligner.ssw_aligner.get_ambiguity_penalty(), rlen_limit, report_cigar);
  }

  ~Aligner() {
    clear();
    cleanup_aligner();
  }

  future<int64_t> fut_get_num_perfect_alns(bool all = false) {
    auto &pr = Timings::get_promise_reduce();
    if (!all) return pr.reduce_one(num_perfect_alns, op_fast_add, 0);
    return pr.reduce_all(num_perfect_alns, op_fast_add);
  }

  future<int64_t> fut_get_num_alns(bool all = false) {
    auto &pr = Timings::get_promise_reduce();
    if (!all) return pr.reduce_one(num_alns, op_fast_add, 0);
    return pr.reduce_all(num_alns, op_fast_add);
  }

  void clear_aln_bufs() {
    kernel_alns.clear();
    ctg_seqs.clear();
    read_seqs.clear();
    max_clen = 0;
    max_rlen = 0;
  }

  void flush_remaining(int read_group_id, KlignTimers &timers) {
    assert(!upcxx::in_progress());
    assert(upcxx::master_persona().active_with_caller());
    DBG(__FILEFUNC__);
    BaseTimer t(__FILEFUNC__);
    t.start();
    for (auto target : upcxx_utils::foreach_rank_by_node()) {
      if (rget_requests[target].size()) do_rget_irregular(target, timers);
    }
    auto num = kernel_alns.size();
    if (num) {
      kernel_align_block(cpu_aligner, kernel_alns, ctg_seqs, read_seqs, alns, active_kernel_fut, read_group_id, max_clen, max_rlen,
                         timers);
      clear_aln_bufs();
    }
    bool is_ready = active_kernel_fut.is_ready();
    active_kernel_fut.wait();
    t.stop();
  }

  void compute_alns_for_read(CtgAndReadLocsMap &aligned_ctgs_map, const string &rname, const string &rseq_fw, int read_group_id,
                             int rget_buf_size, KlignTimers &timers) {
    assert(!upcxx::in_progress());
    assert(upcxx::master_persona().active_with_caller());
    int rlen = rseq_fw.length();
    string rseq_rc;
    string tmp_ctg;
    string rseq;
    for (auto &ctg_and_read_locs : aligned_ctgs_map) {
      vector<CtgAndReadLoc> &locs = ctg_and_read_locs.second.locs;
      auto seq_gptr = ctg_and_read_locs.second.seq_gptr;
      int clen = ctg_and_read_locs.second.clen;
      auto cid = ctg_and_read_locs.second.cid;
      auto depth = ctg_and_read_locs.second.depth;
      if (locs.empty()) continue;
      // sort list of positions in ctg so overlaps can be filtered when processing alignments
      sort(locs.begin(), locs.end(), [](CtgAndReadLoc c1, CtgAndReadLoc c2) {
        if (c1.pos_in_ctg == c2.pos_in_ctg) return c1.pos_in_read < c2.pos_in_read;
        return c1.pos_in_ctg < c2.pos_in_ctg;
      });
      int prev_cstart = -1;
      for (auto &ctg_and_read_loc : locs) {
        progress();
        int pos_in_ctg = ctg_and_read_loc.pos_in_ctg;
        bool ctg_is_rc = ctg_and_read_loc.ctg_is_rc;
        int pos_in_read = ctg_and_read_loc.pos_in_read;
        bool read_is_rc = ctg_and_read_loc.read_is_rc;
        if (pos_in_ctg >= clen || pos_in_read >= rlen)
          DIE("pos_in_ctg ", pos_in_ctg, " ", clen, " pos_in_read ", pos_in_read, " rlen ", rlen);

        char orient = '+';
        if (ctg_is_rc != read_is_rc) {
          // it's revcomp in either contig or read, but not in both or neither
          orient = '-';
          pos_in_read = rlen - (kmer_len + pos_in_read);
        }
        auto [cstart, rstart, overlap_len] = get_start_positions(kmer_len, clen, pos_in_ctg, pos_in_read, rlen);
        assert(cstart >= 0 && cstart + overlap_len <= clen);
        assert(overlap_len <= 2 * rlen);
        // has this alignment already been calculated?
        if (prev_cstart == cstart) continue;
        prev_cstart = cstart;
        if (ctg_is_rc != read_is_rc) {
          if (rseq_rc.empty()) rseq_rc = revcomp(rseq_fw);
          rseq = rseq_rc;
        } else {
          rseq = rseq_fw;
        }
        bool on_node = seq_gptr.is_local();
#ifdef DEBUG
        //     test both on node and off node on a single node
        if (on_node && local_team().rank_n() == rank_n()) on_node = (seq_gptr.where() % 2) == (rank_me() % 2);
#endif
        if (on_node) {
          // on same node already
          string_view ctg_seq = string_view(seq_gptr.local(), clen);
          align_read(rname, cid, rseq, ctg_seq, rstart, cstart, orient, overlap_len, read_group_id, timers);
          local_ctg_fetches++;
        } else {
          auto target = seq_gptr.where();
          if (pos_in_ctg >= clen) DIE(" pos_in_ctg ", pos_in_ctg, " clen ", clen);
          CtgLoc ctg_loc({cid, seq_gptr, clen, depth, pos_in_ctg, ctg_is_rc});
          rget_requests[target].push_back({ctg_loc, rname, rseq, rstart, cstart, orient, overlap_len, read_group_id});
          if (rget_requests[target].size() == rget_buf_size) do_rget_irregular(target, timers);
          remote_ctg_fetches++;
        }
      }
    }
  }

  void finish_alns() {
    assert(!upcxx::in_progress());
    assert(upcxx::master_persona().active_with_caller());
    DBG(__FILEFUNC__);
    if (!kernel_alns.empty()) DIE("sort_alns called while alignments are still pending to be processed - ", kernel_alns.size());
    if (!active_kernel_fut.is_ready()) SWARN("Waiting for active_kernel - has flush_remaining() been called?\n");
    wait_wrapper(active_kernel_fut);
  }

  void log_ctg_bytes_fetched() {
    assert(!upcxx::in_progress());
    assert(upcxx::master_persona().active_with_caller());
    DBG(__FILEFUNC__);
    auto &pr = Timings::get_promise_reduce();
    auto fut_reduce = when_all(pr.reduce_one(rget_calls, op_fast_add, 0), pr.reduce_one(ctg_bytes_fetched, op_fast_add, 0),
                               pr.reduce_one(ctg_bytes_fetched, op_fast_max, 0), pr.reduce_one(local_ctg_fetches, op_fast_add, 0),
                               pr.reduce_one(remote_ctg_fetches, op_fast_add, 0), pr.reduce_one(large_rget_calls, op_fast_add, 0),
                               pr.reduce_one(large_rget_bytes, op_fast_add, 0));

    auto fut_report =
        fut_reduce.then([](auto all_rget_calls, auto all_ctg_bytes_fetched, auto max_ctg_bytes_fetched, auto all_local_ctg_fetches,
                           auto all_remote_ctg_fetches, auto all_large_rget_calls, auto all_large_rget_bytes) {
          if (all_ctg_bytes_fetched > 0)
            SLOG_VERBOSE("Contig bytes fetched ", get_size_str(all_ctg_bytes_fetched), " balance ",
                         (double)all_ctg_bytes_fetched / (rank_n() * max_ctg_bytes_fetched), " average rget size ",
                         get_size_str(all_ctg_bytes_fetched / all_rget_calls), ", large_rgets ", all_large_rget_calls, " avg ",
                         get_size_str(all_large_rget_calls > 0 ? all_large_rget_bytes / all_large_rget_calls : 0), "\n");

          if (all_local_ctg_fetches > 0)
            SLOG_VERBOSE("Local contig fetches ", perc_str(all_local_ctg_fetches, all_local_ctg_fetches + all_remote_ctg_fetches),
                         "\n");
        });
    Timings::set_pending(fut_report);
    ctg_bytes_fetched = 0;
  }
};  // class Aligner

template <int MAX_K>
static upcxx::future<> fetch_ctg_maps_for_target(int target_rank, KmerCtgDHT<MAX_K> &kmer_ctg_dht,
                                                 KmersReadsBuffer<MAX_K> &kmers_reads_buffer, int64_t &num_alns,
                                                 int64_t &num_excess_alns_reads, int64_t &bytes_sent, int64_t &bytes_received,
                                                 int64_t &max_bytes_sent, int64_t &max_bytes_received, int64_t &num_rpcs) {
  assert(!upcxx::in_progress());
  assert(upcxx::master_persona().active_with_caller());
  assert(kmers_reads_buffer.size());
  int64_t sent_msg_size = (sizeof(Kmer<MAX_K>) * kmers_reads_buffer.kmers.size());
  bytes_sent += sent_msg_size;
  max_bytes_sent = std::max(max_bytes_sent, sent_msg_size);
  num_rpcs++;
  // move and consume kmers_read_buffer.  Keep scope until the future completes.
  auto sh_krb = make_shared<KmersReadsBuffer<MAX_K>>();
  std::swap(*sh_krb, kmers_reads_buffer);
  auto fut_get_ctgs = kmer_ctg_dht.get_ctgs_with_kmers(target_rank, sh_krb->kmers);
  auto fut_rpc_returned =
      fut_get_ctgs
          .then([target_rank, &kmers_reads_buffer = *sh_krb, &num_alns, &num_excess_alns_reads, &bytes_received,
                 &max_bytes_received](const vector<CtgLocAndKmerIdx> &ctg_locs_and_kmers_idx) {
            int64_t received_msg_size = (sizeof(CtgLocAndKmerIdx) * ctg_locs_and_kmers_idx.size());
            bytes_received += received_msg_size;
            max_bytes_received = std::max(max_bytes_received, received_msg_size);
            int kmer_len = Kmer<MAX_K>::get_k();
            // iterate through the vector of ctg locations, each one is associated with an index in the read list
            for (int i = 0; i < ctg_locs_and_kmers_idx.size(); i++) {
              auto &ctg_loc = ctg_locs_and_kmers_idx[i].ctg_loc;
              auto kmer_idx = ctg_locs_and_kmers_idx[i].kmer_i;
              auto read_record = kmers_reads_buffer.read_records[kmer_idx].read_record;
              auto rlen = read_record->rlen;
              auto pos_in_read = kmers_reads_buffer.read_records[kmer_idx].read_offset;
              auto read_is_rc = kmers_reads_buffer.read_records[kmer_idx].is_rc;
              assert(read_record->is_valid());
              if (KLIGN_MAX_ALNS_PER_READ && read_record->aligned_ctgs_map.size() >= KLIGN_MAX_ALNS_PER_READ) {
                // too many mappings for this read, stop adding to it
                num_excess_alns_reads++;
                continue;
              }
              auto it = read_record->aligned_ctgs_map.find(ctg_loc.cid);
              auto new_pos_in_read = (ctg_loc.is_rc == read_is_rc ? pos_in_read : rlen - (kmer_len + pos_in_read));
              if (ctg_loc.pos >= ctg_loc.clen) DIE("ctg_loc.pos ", ctg_loc.pos, " ctg_loc.clen ", ctg_loc.clen);
              auto [cstart, rstart, overlap_len] = get_start_positions(kmer_len, ctg_loc.clen, ctg_loc.pos, new_pos_in_read, rlen);
              if (it == read_record->aligned_ctgs_map.end())
                it = read_record->aligned_ctgs_map
                         .insert({ctg_loc.cid, {ctg_loc.cid, ctg_loc.seq_gptr, ctg_loc.clen, ctg_loc.depth, {}}})
                         .first;
              it->second.locs.push_back({ctg_loc.pos, ctg_loc.is_rc, cstart, pos_in_read, read_is_rc});
              num_alns++;
            }
          })
          .then([sh_krb]() {});
  progress();
  return fut_rpc_returned;
};  // fetch_ctg_maps_for_target

template <int MAX_K>
void fetch_ctg_maps(KmerCtgDHT<MAX_K> &kmer_ctg_dht, PackedReads *packed_reads, vector<ReadRecord> &read_records, int seed_space,
                    KlignTimers &timers) {
  assert(!upcxx::in_progress());
  assert(upcxx::master_persona().active_with_caller());
  DBG("fetch_ctg_maps on ", packed_reads->get_fname(), "\n");
  timers.fetch_ctg_maps.start();
  int64_t bytes_sent = 0;
  int64_t bytes_received = 0;
  int64_t max_bytes_sent = 0;
  int64_t max_bytes_received = 0;
  int64_t num_reads = 0;
  int64_t num_alns = 0;
  int64_t num_excess_alns_reads = 0;
  int64_t num_kmers = 0;
  int64_t num_rpcs = 0;

  auto num_local_reads = packed_reads->get_local_num_reads();
  ProgressBar progbar(num_local_reads,
                      string("Fetching ctg maps for alignments - ") + upcxx_utils::get_basename(packed_reads->get_fname()));

  upcxx::future<> fetch_fut_chain = make_future();
  vector<KmersReadsBuffer<MAX_K>> kmers_reads_buffers(num_local_reads > 0 ? rank_n() : 0);
  vector<Kmer<MAX_K>> kmers;
  int kmer_len = Kmer<MAX_K>::get_k();

  // Do not exceed maximum size for RPC messages in either direction, assuming 1 hit per kmer
  size_t max_rdvz_buffer_size = KLIGN_MAX_RPC_MESSAGE_SIZE / ::max(sizeof(Kmer<MAX_K>), sizeof(CtgLocAndKmerIdx));

  // further limit private memory buffering (per rank) to 10% of memory
  size_t max_mem = KLIGN_MAX_MEM_PCT * upcxx_utils::get_max_free_mem_per_rank() / 100;
  size_t max_mem_buffer_size =
      (max_mem - rank_n() * (sizeof(KmersReadsBuffer<MAX_K>))) / (rank_n() * KmersReadsBuffer<MAX_K>::get_size_per());
  max_rdvz_buffer_size = ::min(max_rdvz_buffer_size, max_mem_buffer_size) + 1;

  LOG("max_mem for buffers=", get_size_str(max_mem), " KLIGN_MAX_RPC_MESSAGE_SIZE=", get_size_str(KLIGN_MAX_RPC_MESSAGE_SIZE),
      " max_mem_buffer_size=", max_mem_buffer_size, " max_rdvz_buffer_size=", max_rdvz_buffer_size, "\n");

  // stagger first batches with a reduced count before flush
  size_t max_kmer_buffer_size = max_rdvz_buffer_size / 3 + rand() % (2 * max_rdvz_buffer_size / 3) + 1;
  int stagger_count = 0;

  for (int ri = 0; ri < num_local_reads; ri++) {
    progress();
    progbar.update();
    string read_seq;
    packed_reads->get_read_seq(ri, read_seq);
    // this happens when a placeholder read with just a single N character is added after merging reads
    if (kmer_len > read_seq.length()) continue;
    num_reads++;
    read_records[ri] = ReadRecord({ri, (int)read_seq.length()});
    Kmer<MAX_K>::get_kmers(kmer_len, string_view(read_seq.data(), read_seq.size()), kmers, true);
    if (!kmers.size()) continue;

    for (int i = 0; i < kmers.size(); i += seed_space) {
      const Kmer<MAX_K> &kmer_fw = kmers[i];
      if (!kmer_fw.is_valid()) continue;
      num_kmers++;
      const Kmer<MAX_K> kmer_rc = kmer_fw.revcomp();
      assert(kmer_fw.is_valid() && kmer_rc.is_valid());
      const Kmer<MAX_K> *kmer_lc = &kmer_fw;
      bool is_rc = false;
      if (kmer_rc < kmer_fw) {
        kmer_lc = &kmer_rc;
        is_rc = true;
      }
      assert(kmer_lc->is_least());
      auto target = kmer_ctg_dht.get_target_rank(*kmer_lc);
      auto &tgt_kr_buf = kmers_reads_buffers[target];
      if (tgt_kr_buf.empty()) tgt_kr_buf.reserve(max_rdvz_buffer_size);  // lazy but avoid excessive resizes
      tgt_kr_buf.add(*kmer_lc, &read_records[ri], i, is_rc);
      if (tgt_kr_buf.size() >= max_kmer_buffer_size) {
        auto fetch_fut = fetch_ctg_maps_for_target(target, kmer_ctg_dht, tgt_kr_buf, num_alns, num_excess_alns_reads, bytes_sent,
                                                   bytes_received, max_bytes_sent, max_bytes_received, num_rpcs);
        upcxx_utils::limit_outstanding_futures(fetch_fut).wait();
        fetch_fut_chain = when_all(fetch_fut_chain, fetch_fut);
        if (++stagger_count > 2 * rank_n() / 3)
          max_kmer_buffer_size = max_rdvz_buffer_size;  // stagger achieved on first 2/3 of ranks, reset to max size
        assert(tgt_kr_buf.empty());
      }
    }
    kmers.clear();
  }

  // flush any buffers
  if (num_local_reads > 0) {
    for (auto target : upcxx_utils::foreach_rank_by_node()) {  // stagger by rank_me, round robin by node
      auto &tgt_kr_buf = kmers_reads_buffers[target];
      if (tgt_kr_buf.size()) {
        auto fetch_fut = fetch_ctg_maps_for_target(target, kmer_ctg_dht, tgt_kr_buf, num_alns, num_excess_alns_reads, bytes_sent,
                                                   bytes_received, max_bytes_sent, max_bytes_received, num_rpcs);
        upcxx_utils::limit_outstanding_futures(fetch_fut).wait();
        fetch_fut_chain = when_all(fetch_fut_chain, fetch_fut);
      }
    }
  }
  upcxx_utils::flush_outstanding_futures();
  fetch_fut_chain.wait();
  auto fut_prog_done = progbar.set_done();
  timers.fetch_ctg_maps.stop();
  upcxx_utils::Timings::set_pending(fut_prog_done);

  // defer reporting
  auto &pr = Timings::get_promise_reduce();
  auto fut_reduce = when_all(pr.reduce_one(num_reads, op_fast_add, 0), pr.reduce_one(num_kmers, op_fast_add, 0),
                             pr.reduce_one(bytes_sent, op_fast_add, 0), pr.reduce_one(bytes_received, op_fast_add, 0),
                             pr.reduce_one(max_bytes_sent, op_fast_max, 0), pr.reduce_one(max_bytes_received, op_fast_max, 0),
                             pr.reduce_one(num_rpcs, op_fast_add, 0), pr.reduce_one(num_excess_alns_reads, op_fast_add, 0));

  auto fut_report = fut_reduce.then([](auto all_num_reads, auto all_num_kmers, auto all_bytes_sent, auto all_bytes_received,
                                       auto all_max_bytes_sent, auto all_max_bytes_received, auto all_num_rpcs,
                                       auto all_excess_alns_reads) {
    if (rank_me() == 0) {
      SLOG_VERBOSE("Parsed ", all_num_reads, " reads and extracted ", all_num_kmers, " kmers\n");
      SLOG_VERBOSE("Sent ", get_size_str(all_bytes_sent), " (", get_size_str(all_num_rpcs > 0 ? all_bytes_sent / all_num_rpcs : 0),
                   " avg msg, ", get_size_str(all_max_bytes_sent), " max msg) of kmers and received ",
                   get_size_str(all_bytes_received), " (", get_size_str(all_num_rpcs > 0 ? all_bytes_received / all_num_rpcs : 0),
                   " avg msg, ", get_size_str(all_max_bytes_received), " max msg )\n");
      if (all_excess_alns_reads)
        SLOG_VERBOSE("Dropped ", all_excess_alns_reads, " alignments in excess of ", KLIGN_MAX_ALNS_PER_READ, " per read\n");
    }
  });
  Timings::set_pending(fut_report);
};  // fetch_ctg_maps

template <int MAX_K>
void compute_alns(PackedReads *packed_reads, vector<ReadRecord> &read_records, Alns &alns, int read_group_id, int rlen_limit,
                  bool report_cigar, bool use_blastn_scores, int64_t all_num_ctgs, int rget_buf_size, KlignTimers &timers) {
  assert(!upcxx::in_progress());
  assert(upcxx::master_persona().active_with_caller());
  auto short_name = get_basename(packed_reads->get_fname());
  DBG("compute_alns on ", short_name, "\n");
  timers.compute_alns.start();
  int kmer_len = Kmer<MAX_K>::get_k();
  int64_t num_reads_aligned = 0;
  int64_t num_reads = 0;
  Aligner aligner(Kmer<MAX_K>::get_k(), alns, rlen_limit, report_cigar, use_blastn_scores);
  string read_seq, read_id, read_quals;
  ProgressBar progbar(packed_reads->get_local_num_reads(), string("Computing alignments on ") + short_name);
  for (auto &read_record : read_records) {
    progress();
    progbar.update();
    if (kmer_len > read_record.rlen) continue;
    num_reads++;
    // compute alignments
    if (read_record.aligned_ctgs_map.size()) {
      num_reads_aligned++;
      packed_reads->get_read(read_record.index, read_id, read_seq, read_quals);
      aligner.compute_alns_for_read(read_record.aligned_ctgs_map, read_id, read_seq, read_group_id, rget_buf_size, timers);
    }
  }
  aligner.flush_remaining(read_group_id, timers);
  aligner.finish_alns();
  Timings::set_pending(progbar.set_done());
  timers.compute_alns.stop();
  read_records.clear();

  // defer reporting
  aligner.log_ctg_bytes_fetched();
  auto &pr = Timings::get_promise_reduce();
  auto fut_reduce =
      when_all(pr.reduce_one(num_reads, op_fast_add, 0), pr.reduce_one(num_reads_aligned, op_fast_add, 0),
               aligner.fut_get_num_alns(), aligner.fut_get_num_perfect_alns(), pr.reduce_one(alns.get_num_dups(), op_fast_add, 0),
               pr.reduce_one(alns.get_num_bad(), op_fast_add, 0), pr.reduce_one(alns.size(), op_fast_add, 0));

  auto fut_report = fut_reduce.then([short_name](auto all_num_reads, auto all_num_reads_aligned, auto all_num_alns,
                                                 auto all_num_perfect, auto all_num_dups, auto all_num_bad, auto all_num_good) {
    SLOG_VERBOSE("Found ", all_num_alns, " alignments in ", short_name, "\n");
    SLOG_VERBOSE("  perfect ", perc_str(all_num_perfect, all_num_alns), "\n");
    SLOG_VERBOSE("  good ", perc_str(all_num_good, all_num_alns), "\n");
    SLOG_VERBOSE("  bad ", perc_str(all_num_bad, all_num_alns), "\n");
    SLOG_VERBOSE("  duplicates ", perc_str(all_num_dups, all_num_alns), "\n");
    SLOG_VERBOSE("Mapped ", perc_str(all_num_reads_aligned, all_num_reads), " reads to contigs, average mappings per read ",
                 (double)all_num_alns / all_num_reads_aligned, "\n");
  });
  Timings::set_pending(fut_report);
};  // compute_alns

template <int MAX_K>
future<> sort_alns(Alns &alns, KlignTimers &timers, const string &reads_fname) {
  // barrier before possibly losing attention for sort & append
  LOG("Entering barrier after all reads have been aligned and before sorting alignments\n");
  BarrierTimer sort_bt("Sorting Alignments for " + get_basename(reads_fname));

  LOG("Sorting ", alns.size(), " alignments in a new thread after discharge\n");
  // wrap expensive sort in thread and keep some attention while waiting
  auto sort_lambda = [&sort_t = timers.sort_t, &alns]() {
    sort_t.start();
    alns.sort_alns();
    sort_t.stop();
    LOG("Finished sorting ", alns.size(), " alignments\n");
  };
  future<> fut_sort = upcxx_utils::ThreadPool::get_single_pool().enqueue_serially(sort_lambda);
  return fut_sort;
}

template <int MAX_K>
shared_ptr<KmerCtgDHT<MAX_K>> build_kmer_ctg_dht(unsigned kmer_len, int max_store_size, int max_rpcs_in_flight, Contigs &ctgs,
                                                 int min_ctg_len, bool allow_multi_kmers) {
  BarrierTimer timer(__FILEFUNC__);
  Kmer<MAX_K>::set_k(kmer_len);
  uint64_t num_ctg_kmers = 0;
  size_t num_ctgs = 0;
  for (auto &ctg : ctgs) {
    if (ctg.seq.length() >= Kmer<MAX_K>::get_k()) {
      num_ctg_kmers += ctg.seq.length() - Kmer<MAX_K>::get_k() + 1;
      num_ctgs++;
    }
  }
  auto sh_kmer_ctg_dht = make_shared<KmerCtgDHT<MAX_K>>(max_store_size, max_rpcs_in_flight, allow_multi_kmers, num_ctg_kmers);
  auto &kmer_ctg_dht = *sh_kmer_ctg_dht;
  kmer_ctg_dht.build(ctgs, min_ctg_len);
  kmer_ctg_dht.purge_high_count_seeds();
  return sh_kmer_ctg_dht;
}

template <int MAX_K>
pair<double, double> find_alignments(unsigned kmer_len, PackedReadsList &packed_reads_list, int max_store_size,
                                     int max_rpcs_in_flight, Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit,
                                     bool use_blastn_scores, int min_ctg_len, int rget_buf_size) {
  BarrierTimer timer(__FILEFUNC__);
  SLOG_VERBOSE("Aligning with seed size of ", kmer_len, " and seed space ", seed_space, "\n");

  const bool ALLOW_MULTI_KMERS = false;
  const bool REPORT_CIGAR = false;
  auto sh_kmer_ctg_dht =
      build_kmer_ctg_dht<MAX_K>(kmer_len, max_store_size, max_rpcs_in_flight, ctgs, min_ctg_len, ALLOW_MULTI_KMERS);
  auto &kmer_ctg_dht = *sh_kmer_ctg_dht;
  KlignTimers timers;
  int64_t all_num_ctgs = reduce_all(ctgs.size(), op_fast_add).wait();

#ifdef DEBUG
// kmer_ctg_dht.dump_ctg_kmers();
#endif
  int read_group_id = 0;
  upcxx::promise prom_gpu(1);
  for (auto packed_reads : packed_reads_list) {
    vector<ReadRecord> read_records(packed_reads->get_local_num_reads());
    fetch_ctg_maps(kmer_ctg_dht, packed_reads, read_records, seed_space, timers);
    auto sh_alns = make_shared<Alns>();
    Alns &alns_for_sample = *sh_alns;
    compute_alns<MAX_K>(packed_reads, read_records, alns_for_sample, read_group_id, rlen_limit, REPORT_CIGAR, use_blastn_scores,
                        all_num_ctgs, rget_buf_size, timers);
    if (!sh_alns->empty()) {
      auto fut_sort = sort_alns<MAX_K>(alns_for_sample, timers, packed_reads->get_fname()).then([sh_alns, &alns]() {
        assert(upcxx::master_persona().active_with_caller());
        alns.append(*sh_alns);
      });
      Timings::set_pending(fut_sort);
    }
    read_group_id++;
  }
  Timings::wait_pending();
  timers.done_all();
  double aln_kernel_elapsed = timers.aln_kernel.get_elapsed();
  double aln_comms_elapsed = timers.fetch_ctg_maps.get_elapsed() + timers.rget_ctg_seqs.get_elapsed();
  timers.clear();
  return {aln_kernel_elapsed, aln_comms_elapsed};
};  // find_alignments
