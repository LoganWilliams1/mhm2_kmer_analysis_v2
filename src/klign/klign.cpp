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

#include <algorithm>
#include <iostream>
#include <thread>
#include <upcxx/upcxx.hpp>

#include "klign.hpp"
#include "kmer.hpp"
#include "ssw.hpp"
#include "upcxx_utils/limit_outstanding.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "zstr.hpp"
#include "aligner_cpu.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

using cid_t = int64_t;

struct KlignTimers {
  upcxx_utils::IntermittentTimer fetch_ctg_maps, compute_alns, rget_ctg_seqs, aln_kernel;

  KlignTimers()
      : fetch_ctg_maps("klign: fetch ctg maps")
      , compute_alns("klign: compute alns")
      , rget_ctg_seqs("klign: rget ctg seqs")
      , aln_kernel("klign: aln kernel") {}

  void done_all() {
    fetch_ctg_maps.done_all();
    compute_alns.done_all();
    rget_ctg_seqs.done_all();
    aln_kernel.done_all();
  }

  void clear() {
    fetch_ctg_maps.clear();
    compute_alns.clear();
    rget_ctg_seqs.clear();
    aln_kernel.clear();
  }
};

void init_aligner(int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty, int ambiguity_penalty,
                  int rlen_limit, bool compute_cigar);
void cleanup_aligner();
void kernel_align_block(CPUAligner &cpu_aligner, vector<Aln> &kernel_alns, vector<string> &ctg_seqs, vector<string> &read_seqs,
                        Alns *alns, future<> &active_kernel_fut, int read_group_id, int max_clen, int max_rlen,
                        IntermittentTimer &aln_kernel_timer);

struct CtgLoc {
  cid_t cid;
  global_ptr<char> seq_gptr;
  int clen;
  float depth;
  int pos;
  bool is_rc;
};

struct CtgAndReadLoc {
  CtgLoc ctg_loc;
  int32_t cstart;
  int32_t pos_in_read : 31;
  int32_t read_is_rc : 1;
};

template <int MAX_K>
struct KmerAndCtgLoc {
  Kmer<MAX_K> kmer;
  CtgLoc ctg_loc;
  UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
};

struct CtgLocAndKmerIdx {
  CtgLoc ctg_loc;
  int kmer_i;
};

using CtgAndReadLocsMap = HASH_TABLE<cid_t, vector<CtgAndReadLoc>>;

struct ReadRecord {
  int64_t index : 40;
  int64_t rlen : 24;

  CtgAndReadLocsMap aligned_ctgs_map;

  ReadRecord()
      : index(-1)
      , rlen(0) {}

  ReadRecord(int index, int rlen)
      : index(index)
      , rlen(rlen)
      , aligned_ctgs_map{} {}

  bool is_valid() const { return index >= 0 && rlen > 0; }
};

struct ReadRecordPtr {
  ReadRecord *read_record;
  int read_offset;
  bool is_rc;
};

template <int MAX_K>
struct KmersReadsBuffer {
  vector<Kmer<MAX_K>> kmers;
  vector<ReadRecordPtr> read_records;

  void add(const Kmer<MAX_K> &kmer, ReadRecord *read_record, int read_offset, bool read_is_rc) {
    kmers.push_back(kmer);
    read_records.push_back({read_record, read_offset, read_is_rc});
  }

  size_t size() { return kmers.size(); }

  void clear() {
    kmers.clear();
    read_records.clear();
  }
};

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
          it->second.push_back(ctg_loc);
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

  int64_t get_num_kmers(bool all = false) {
    if (!all) return reduce_one(kmer_map->size(), op_fast_add, 0).wait();
    return reduce_all(kmer_map->size(), op_fast_add).wait();
  }

  int64_t get_num_dropped_seed_to_ctgs(bool all = false) {
    if (!all) return reduce_one(*num_dropped_seed_to_ctgs, op_fast_add, 0).wait();
    return reduce_all(*num_dropped_seed_to_ctgs, op_fast_add).wait();
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
    size_t max_ctgs = 0;
    // determine max number of ctgs mapped to by a single kmer
    for (auto &elem : *kmer_map) {
      max_ctgs = ::max(max_ctgs, elem.second.size());
    }
    barrier();
    auto all_max_ctgs = reduce_one(max_ctgs, op_fast_max, 0).wait();
    if (all_max_ctgs > 1) SLOG_VERBOSE("Max contigs mapped by a single kmer: ", all_max_ctgs, "\n");
    BarrierTimer timer(__FILEFUNC__, false);  // barrier on exit, not entrance
    kmer_store.flush_updates();
    kmer_store.clear();
  }

  future<vector<CtgLocAndKmerIdx>> get_ctgs_with_kmers(int target_rank, vector<Kmer<MAX_K>> &kmers) {
    DBG_VERBOSE("Sending request for ", kmers.size(), " to ", target_rank, "\n");
    return rpc(
        target_rank,
        [allow_multi_kmers = this->allow_multi_kmers](vector<Kmer<MAX_K>> kmers, kmer_map_t &kmer_map) {
          vector<CtgLocAndKmerIdx> ctg_locs;
          for (int i = 0; i < kmers.size(); i++) {
            auto &kmer = kmers[i];
            assert(kmer.is_least());
            assert(kmer.is_valid());
            const auto it = kmer_map->find(kmer);
            if (it == kmer_map->end()) continue;
            for (auto &ctg_loc : it->second) {
              ctg_locs.push_back({ctg_loc, i});
            }
          }
          return ctg_locs;
        },
        kmers, kmer_map);
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
    progbar.done();
    SLOG_VERBOSE("Dumped ", this->get_num_kmers(), " kmers\n");
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
    SLOG_VERBOSE("Estimated ", est_num_kmers, " contig ", kmer_len, "-kmers for ", tot_num_ctgs, " contigs with total len ", tot_ctg_lengths, "\n");
    auto my_required_mem = my_reserve * (sizeof(typename local_kmer_map_t::value_type) + sizeof(CtgLoc)); // 1 map key/value entry + 1 CtgLoc within the vector
    LOG("Reserving ", my_reserve, " for my entries at approx ", get_size_str(my_required_mem * upcxx::local_team().rank_n()), " per node for the local hashtable\n");
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
        // if (kmers[i].to_string() == "ACATCTACCGCTAGAGGATTA")
        //   WARN("kmer ", kmers[i].to_string(), " found in contig ", ctg->id, " in position ", i);
        ctg_loc.pos = i;
        if (!kmers[i].is_valid()) continue;
        add_kmer(kmers[i], ctg_loc);
      }
      progress();
    }
    auto fut = progbar.set_done();
    flush_add_kmers();
    auto tot_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
    fut.wait();
    auto num_kmers_in_ht = get_num_kmers();
    LOG("Estimated room for ", my_reserve, " my final count ", kmer_map->size(), "\n");
    SLOG_VERBOSE("Total contigs >= ", min_ctg_len, ": ", tot_num_ctgs, " seq_length: ", tot_ctg_lengths, "\n");
    SLOG_VERBOSE("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");
    auto num_dropped_seed_to_ctgs = get_num_dropped_seed_to_ctgs();
    if (num_dropped_seed_to_ctgs)
      SLOG_VERBOSE("For k = ", kmer_len, " dropped ", num_dropped_seed_to_ctgs, " non-unique seed-to-contig mappings (",
                   setprecision(2), fixed, (100.0 * num_dropped_seed_to_ctgs / tot_num_kmers), "%)\n");
  }
};

static tuple<int, int, int> get_start_positions(int kmer_len, const CtgLoc &ctg_loc, int pos_in_read, int rlen) {
  // calculate available bases before and after the seeded kmer
  int ctg_bases_left_of_kmer = ctg_loc.pos;
  int ctg_bases_right_of_kmer = ctg_loc.clen - ctg_bases_left_of_kmer - kmer_len;
  int read_bases_left_of_kmer = pos_in_read;
  int read_bases_right_of_kmer = rlen - kmer_len - pos_in_read;
  int left_of_kmer = min(read_bases_left_of_kmer, ctg_bases_left_of_kmer);
  int right_of_kmer = min(read_bases_right_of_kmer, ctg_bases_right_of_kmer);

  int cstart = ctg_loc.pos - left_of_kmer;
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
  int rget_calls = 0;
  int local_ctg_fetches = 0;
  int remote_ctg_fetches = 0;

  Alns *alns;

 private:
  vector<vector<RgetRequest>> rget_requests;

  void align_read(const string &rname, int64_t cid, const string_view &rseq, const string_view &ctg_seq, int rstart, int cstart,
                  char orient, int overlap_len, int read_group_id, KlignTimers &timers) {
    num_alns++;
    size_t clen = ctg_seq.length();
    size_t rlen = rseq.length();
    string_view cseq = string_view(ctg_seq.data() + cstart, overlap_len);
    if (cseq.compare(0, overlap_len, rseq, rstart, overlap_len) == 0) {
      num_perfect_alns++;
      int rstop = rstart + overlap_len;
      int cstop = cstart + overlap_len;
      if (orient == '-') switch_orient(rstart, rstop, rlen);
      int score1 = overlap_len * cpu_aligner.ssw_aligner.get_match_score();
      Aln aln(rname, cid, rstart, rstop, rlen, cstart, cstop, clen, orient, score1, 0, 0, read_group_id);
      assert(aln.is_valid());
      if (cpu_aligner.ssw_filter.report_cigar) aln.set_sam_string(rseq, to_string(overlap_len) + "=");
      alns->add_aln(aln);
    } else {
      max_clen = max((int64_t)cseq.size(), max_clen);
      max_rlen = max((int64_t)rseq.size(), max_rlen);
      int64_t num_kernel_alns = kernel_alns.size() + 1;
      int64_t tot_mem_est = num_kernel_alns * (max_clen + max_rlen + 2 * sizeof(int) + 5 * sizeof(short));
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      // this is also the way it's done in meraligner
      kernel_alns.emplace_back(rname, cid, 0, 0, rlen, cstart, 0, clen, orient);
      ctg_seqs.emplace_back(cseq);
      read_seqs.emplace_back(rseq);
      if (num_kernel_alns >= KLIGN_GPU_BLOCK_SIZE) {
        kernel_align_block(cpu_aligner, kernel_alns, ctg_seqs, read_seqs, alns, active_kernel_fut, read_group_id, max_clen,
                           max_rlen, timers.aln_kernel);
        clear_aln_bufs();
      }
    }
  }

  void clear() {
    if (kernel_alns.size() || !active_kernel_fut.ready())
      DIE("clear called with alignments in the buffer or active kernel - was flush_remaining called before destrutor?\n");
    clear_aln_bufs();
  }

  void do_rget_irregular(int target, KlignTimers &timers) {
    vector<pair<global_ptr<char>, size_t>> src;
    vector<pair<char *, size_t>> dest;
    HASH_TABLE<cid_t, string> ctgs_fetched;
    for (auto &req : rget_requests[target]) {
      auto clen = req.ctg_loc.clen;
      auto it = ctgs_fetched.find(req.ctg_loc.cid);
      if (it == ctgs_fetched.end()) {
        it = ctgs_fetched.insert({req.ctg_loc.cid, string(clen, ' ')}).first;
        src.push_back({req.ctg_loc.seq_gptr, clen});
        dest.push_back({const_cast<char *>(it->second.data()), clen});
        ctg_bytes_fetched += clen;
      }
    }
    // SLOG_VERBOSE("Using rget_irregular to fetch ", perc_str(ctgs_fetched.size(), rget_requests[target].size()), " contigs\n");
    timers.rget_ctg_seqs.start();
    rget_irregular(src.begin(), src.end(), dest.begin(), dest.end()).wait();
    timers.rget_ctg_seqs.stop();
    rget_calls++;
    for (auto &req : rget_requests[target]) {
      auto cid = req.ctg_loc.cid;
      auto it = ctgs_fetched.find(cid);
      if (it == ctgs_fetched.end()) DIE("Could not find the sequence for the contig ", cid);
      string &ctg_seq = it->second;
      align_read(req.rname, cid, req.read_seq, ctg_seq, req.rstart, req.cstart, req.orient, req.overlap_len, req.read_group_id,
                 timers);
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

  int64_t get_num_perfect_alns(bool all = false) {
    if (!all) return reduce_one(num_perfect_alns, op_fast_add, 0).wait();
    return reduce_all(num_perfect_alns, op_fast_add).wait();
  }

  int64_t get_num_alns(bool all = false) {
    if (!all) return reduce_one(num_alns, op_fast_add, 0).wait();
    return reduce_all(num_alns, op_fast_add).wait();
  }

  void clear_aln_bufs() {
    kernel_alns.clear();
    ctg_seqs.clear();
    read_seqs.clear();
    max_clen = 0;
    max_rlen = 0;
  }

  void flush_remaining(int read_group_id, KlignTimers &timers) {
    BaseTimer t(__FILEFUNC__);
    t.start();
    for (auto target : upcxx_utils::foreach_rank_by_node()) {
      if (rget_requests[target].size()) do_rget_irregular(target, timers);
    }
    auto num = kernel_alns.size();
    if (num) {
      kernel_align_block(cpu_aligner, kernel_alns, ctg_seqs, read_seqs, alns, active_kernel_fut, read_group_id, max_clen, max_rlen,
                         timers.aln_kernel);
      clear_aln_bufs();
    }
    bool is_ready = active_kernel_fut.ready();
    active_kernel_fut.wait();
    t.stop();
  }

  void compute_alns_for_read(CtgAndReadLocsMap &aligned_ctgs_map, const string &rname, const string &rseq_fw, int read_group_id,
                             KlignTimers &timers) {
    int rlen = rseq_fw.length();
    string rseq_rc;
    string tmp_ctg;
    for (auto &ctg_and_read_locs : aligned_ctgs_map) {
      progress();
      for (auto &ctg_and_read_loc : ctg_and_read_locs.second) {
        int pos_in_read = ctg_and_read_loc.pos_in_read;
        bool read_is_rc = ctg_and_read_loc.read_is_rc;
        auto &ctg_loc = ctg_and_read_loc.ctg_loc;
        char orient = '+';
        string rseq;
        if (ctg_loc.is_rc != read_is_rc) {
          // it's revcomp in either contig or read, but not in both or neither
          orient = '-';
          pos_in_read = rlen - (kmer_len + pos_in_read);
          if (rseq_rc.empty()) rseq_rc = revcomp(rseq_fw);
          rseq = rseq_rc;
        } else {
          rseq = rseq_fw;
        }
        auto [cstart, rstart, overlap_len] = get_start_positions(kmer_len, ctg_loc, pos_in_read, rlen);
        assert(cstart >= 0 && cstart + overlap_len <= ctg_loc.clen);
        assert(overlap_len <= 2 * rlen);
        bool on_node = ctg_loc.seq_gptr.is_local();
#ifdef DEBUG
        //    test both on node and off node on a single node
        if (on_node && local_team().rank_n() == rank_n()) on_node = (ctg_loc.seq_gptr.where() % 2) == (rank_me() % 2);
#endif
        if (on_node) {
          // on same node already
          string_view ctg_seq = string_view(ctg_loc.seq_gptr.local(), ctg_loc.clen);
          align_read(rname, ctg_loc.cid, rseq, ctg_seq, rstart, cstart, orient, overlap_len, read_group_id, timers);
          local_ctg_fetches++;
        } else {
          auto target = ctg_loc.seq_gptr.where();
          rget_requests[target].push_back({ctg_loc, rname, rseq, rstart, cstart, orient, overlap_len, read_group_id});
          if (rget_requests[target].size() == KLIGN_RGET_BUF_SIZE) do_rget_irregular(target, timers);
          remote_ctg_fetches++;
        }
      }
    }
  }

  void sort_alns() {
    if (!kernel_alns.empty()) DIE("sort_alns called while alignments are still pending to be processed - ", kernel_alns.size());
    if (!active_kernel_fut.ready()) SWARN("Waiting for active_kernel - has flush_remaining() been called?\n");
    active_kernel_fut.wait();
    alns->sort_alns().wait();
  }

  void log_ctg_bytes_fetched() {
    auto all_rget_calls = reduce_one(rget_calls, op_fast_add, 0).wait();
    auto all_ctg_bytes_fetched = reduce_one(ctg_bytes_fetched, op_fast_add, 0).wait();
    auto max_ctg_bytes_fetched = reduce_one(ctg_bytes_fetched, op_fast_max, 0).wait();
    if (all_ctg_bytes_fetched > 0)
      SLOG_VERBOSE("Contig bytes fetched ", get_size_str(all_ctg_bytes_fetched), " balance ",
                   (double)all_ctg_bytes_fetched / (rank_n() * max_ctg_bytes_fetched), " average rget size ",
                   get_size_str(all_ctg_bytes_fetched / all_rget_calls), "\n");
    ctg_bytes_fetched = 0;
    auto all_local_ctg_fetches = reduce_one(local_ctg_fetches, op_fast_add, 0).wait();
    auto all_remote_ctg_fetches = reduce_one(remote_ctg_fetches, op_fast_add, 0).wait();
    if (all_local_ctg_fetches > 0)
      SLOG_VERBOSE("Local contig fetches ", perc_str(all_local_ctg_fetches, all_local_ctg_fetches + all_remote_ctg_fetches), "\n");
  }
};

template <int MAX_K>
static upcxx::future<> fetch_ctg_maps_for_target(int target_rank, KmerCtgDHT<MAX_K> &kmer_ctg_dht,
                                                 KmersReadsBuffer<MAX_K> &kmers_reads_buffer, int64_t &num_alns,
                                                 int64_t &num_excess_alns_reads, int64_t &bytes_sent, int64_t &bytes_received,
                                                 int64_t &max_bytes_sent, int64_t &max_bytes_received, int64_t &num_rpcs) {
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
                 &max_bytes_received](const vector<CtgLocAndKmerIdx> ctg_locs_and_kmers_idx) {
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
              auto [cstart, rstart, overlap_len] = get_start_positions(kmer_len, ctg_loc, new_pos_in_read, rlen);
              bool overlaps = false;
              if (it == read_record->aligned_ctgs_map.end()) {
                it = read_record->aligned_ctgs_map.insert({ctg_loc.cid, {}}).first;
              } else {
                // check to see if alignment will overlap one from an earlier seed
                for (auto &prev_ctg_loc : it->second) {
                  if (cstart + rlen >= prev_ctg_loc.cstart && cstart < prev_ctg_loc.cstart + rlen) {
                    overlaps = true;
                    break;
                  }
                }
              }
              if (!overlaps) {
                it->second.push_back({ctg_loc, cstart, pos_in_read, read_is_rc});
                num_alns++;
              }
            }
          })
          .then([sh_krb]() {});
  return fut_rpc_returned;
}

template <int MAX_K>
void fetch_ctg_maps(KmerCtgDHT<MAX_K> &kmer_ctg_dht, PackedReads *packed_reads, vector<ReadRecord> &read_records, int seed_space,
                    KlignTimers &timers) {
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

  upcxx::future<> fetch_fut_chain = make_future();
  vector<KmersReadsBuffer<MAX_K>> kmers_reads_buffers(rank_n());
  vector<Kmer<MAX_K>> kmers;
  int kmer_len = Kmer<MAX_K>::get_k();

  // Do not exceed 256KB RPC messages in either direction
  size_t max_rdzv_message_size = 256 * 1024;
  size_t max_rdvz_buffer_size = max_rdzv_message_size / ::max(sizeof(Kmer<MAX_K>), sizeof(CtgLocAndKmerIdx));
  size_t max_kmer_buffer_size = ::min((size_t)KLIGN_KMERS_BUF_SIZE / rank_n(), max_rdvz_buffer_size);

  ProgressBar progbar(packed_reads->get_local_num_reads(), "Fetching ctg maps for alignments");
  for (int ri = 0; ri < packed_reads->get_local_num_reads(); ri++) {
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
      kmers_reads_buffers[target].add(*kmer_lc, &read_records[ri], i, is_rc);
      if (kmers_reads_buffers[target].size() >= max_kmer_buffer_size) {
        auto fetch_fut =
            fetch_ctg_maps_for_target(target, kmer_ctg_dht, kmers_reads_buffers[target], num_alns, num_excess_alns_reads,
                                      bytes_sent, bytes_received, max_bytes_sent, max_bytes_received, num_rpcs);
        fetch_fut_chain = when_all(fetch_fut_chain, fetch_fut);
        upcxx_utils::limit_outstanding_futures(fetch_fut_chain).wait();
      }
    }
    kmers.clear();
  }
  for (auto target : upcxx_utils::foreach_rank_by_node()) {  // stagger by rank_me, round robin by node
    if (kmers_reads_buffers[target].size()) {
      auto fetch_fut = fetch_ctg_maps_for_target(target, kmer_ctg_dht, kmers_reads_buffers[target], num_alns, num_excess_alns_reads,
                                                 bytes_sent, bytes_received, max_bytes_sent, max_bytes_received, num_rpcs);
      fetch_fut_chain = when_all(fetch_fut_chain, fetch_fut);
      upcxx_utils::limit_outstanding_futures(fetch_fut_chain).wait();
    }
  }
  fetch_fut_chain.wait();
  upcxx_utils::flush_outstanding_futures();
  upcxx_utils::Timings::set_pending(progbar.set_done());
  Timings::wait_pending();
  barrier();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto all_bytes_sent = reduce_one(bytes_sent, op_fast_add, 0).wait();
  auto all_bytes_received = reduce_one(bytes_received, op_fast_add, 0).wait();
  auto all_max_bytes_sent = reduce_one(max_bytes_sent, op_fast_max, 0).wait();
  auto all_max_bytes_received = reduce_one(max_bytes_received, op_fast_max, 0).wait();
  auto all_num_rpcs = reduce_one(num_rpcs, op_fast_add, 0).wait();
  auto all_excess_alns_reads = reduce_one(num_excess_alns_reads, op_fast_add, 0).wait();
  if (rank_me() == 0) {
    SLOG_VERBOSE("Parsed ", all_num_reads, " reads and extracted ", all_num_kmers, " kmers\n");
    SLOG_VERBOSE("Sent ", get_size_str(all_bytes_sent), " (", get_size_str(all_num_rpcs > 0 ? all_bytes_sent / all_num_rpcs : 0),
                 " avg msg, ", get_size_str(all_max_bytes_sent), " max msg) of kmers and received ",
                 get_size_str(all_bytes_received), " (", get_size_str(all_num_rpcs > 0 ? all_bytes_received / all_num_rpcs : 0),
                 " avg msg, ", get_size_str(all_max_bytes_received), " max msg )\n");
    if (all_excess_alns_reads)
      SLOG_VERBOSE("Dropped ", all_excess_alns_reads, " alignments in excess of ", KLIGN_MAX_ALNS_PER_READ, " per read\n");
  }
  timers.fetch_ctg_maps.stop();
  barrier();
}

template <int MAX_K>
void compute_alns(PackedReads *packed_reads, vector<ReadRecord> &read_records, Alns &alns, int read_group_id, int rlen_limit,
                  bool report_cigar, bool use_blastn_scores, int64_t all_num_ctgs, KlignTimers &timers) {
  timers.compute_alns.start();
  int kmer_len = Kmer<MAX_K>::get_k();
  int64_t num_reads_aligned = 0;
  int64_t num_reads = 0;
  Alns alns_for_sample;
  Aligner aligner(Kmer<MAX_K>::get_k(), alns_for_sample, rlen_limit, report_cigar, use_blastn_scores);
  string read_seq, read_id, read_quals;
  ProgressBar progbar(packed_reads->get_local_num_reads(), "Computing alignments");
  for (auto &read_record : read_records) {
    progress();
    progbar.update();
    if (kmer_len > read_record.rlen) continue;
    num_reads++;
    // compute alignments
    if (read_record.aligned_ctgs_map.size()) {
      num_reads_aligned++;
      packed_reads->get_read(read_record.index, read_id, read_seq, read_quals);
      aligner.compute_alns_for_read(read_record.aligned_ctgs_map, read_id, read_seq, read_group_id, timers);
    }
  }
  aligner.flush_remaining(read_group_id, timers);
  Timings::set_pending(progbar.set_done());
  Timings::wait_pending();
  barrier();
  read_records.clear();
  aligner.sort_alns();
  aligner.log_ctg_bytes_fetched();

  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_reads_aligned = reduce_one(num_reads_aligned, op_fast_add, 0).wait();
  auto all_num_alns = aligner.get_num_alns();
  auto all_num_perfect = aligner.get_num_perfect_alns();
  auto all_num_dups = reduce_one(alns_for_sample.get_num_dups(), op_fast_add, 0).wait();
  auto all_num_bad = reduce_one(alns_for_sample.get_num_bad(), op_fast_add, 0).wait();
  auto all_num_good = reduce_one(alns_for_sample.size(), op_fast_add, 0).wait();
  SLOG("Found ", all_num_alns, " alignments:\n");
  SLOG("  perfect ", perc_str(all_num_perfect, all_num_alns), "\n");
  SLOG("  good ", perc_str(all_num_good, all_num_alns), "\n");
  SLOG("  bad ", perc_str(all_num_bad, all_num_alns), "\n");
  SLOG("  duplicates ", perc_str(all_num_dups, all_num_alns), "\n");
  SLOG_VERBOSE("Mapped ", perc_str(all_num_reads_aligned, all_num_reads), " reads to contigs, average mappings per read ",
               (double)all_num_alns / all_num_reads_aligned, "\n");
  alns.append(alns_for_sample);
  timers.compute_alns.stop();
  barrier();
}

template <int MAX_K>
double find_alignments(unsigned kmer_len, PackedReadsList &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool report_cigar, bool use_blastn_scores,
                       int min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);
  Kmer<MAX_K>::set_k(kmer_len);
  KlignTimers timers;
  SLOG_VERBOSE("Aligning with seed size of ", kmer_len, " and seed space ", seed_space, "\n");
  int64_t all_num_ctgs = reduce_all(ctgs.size(), op_fast_add).wait();
  uint64_t num_ctg_kmers = 0;
  for (auto &ctg : ctgs)
    if (ctg.seq.length() >= Kmer<MAX_K>::get_k()) num_ctg_kmers += ctg.seq.length() - Kmer<MAX_K>::get_k() + 1;
  bool allow_multi_kmers = report_cigar;
  KmerCtgDHT<MAX_K> kmer_ctg_dht(max_store_size, max_rpcs_in_flight, allow_multi_kmers, num_ctg_kmers);
  barrier();
  kmer_ctg_dht.build(ctgs, min_ctg_len);
#ifdef DEBUG
// kmer_ctg_dht.dump_ctg_kmers();
#endif
  int read_group_id = 0;
  for (auto packed_reads : packed_reads_list) {
    vector<ReadRecord> read_records(packed_reads->get_local_num_reads());
    fetch_ctg_maps(kmer_ctg_dht, packed_reads, read_records, seed_space, timers);
    compute_alns<MAX_K>(packed_reads, read_records, alns, read_group_id, rlen_limit, report_cigar, use_blastn_scores, all_num_ctgs,
                        timers);
    read_group_id++;
  }
  barrier();
  timers.done_all();
  double aln_kernel_elapsed = timers.aln_kernel.get_elapsed();
  timers.clear();
  return aln_kernel_elapsed;
}
