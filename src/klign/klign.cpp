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
#include "upcxx_utils/fixed_size_cache.hpp"
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
  upcxx_utils::IntermittentTimer fetch_ctg_maps, compute_alns, get_ctgs, rget_ctg_seqs, aln_kernel;

  KlignTimers()
      : fetch_ctg_maps("klign: fetch ctg maps")
      , compute_alns("klign: compute alns")
      , get_ctgs("klign: get ctgs with kmers")
      , rget_ctg_seqs("klign: rget ctg seqs")
      , aln_kernel("klign: aln kernel") {}

  void done_all() {
    fetch_ctg_maps.done_all();
    compute_alns.done_all();
    get_ctgs.done_all();
    rget_ctg_seqs.done_all();
    aln_kernel.done_all();
  }

  void clear() {
    fetch_ctg_maps.clear();
    compute_alns.clear();
    get_ctgs.clear();
    rget_ctg_seqs.clear();
    aln_kernel.clear();
  }
};

static KlignTimers timers;

void init_aligner(int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty, int ambiguity_penalty,
                  int rlen_limit);
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
  int pos_in_read;
  bool read_is_rc;
  int cstart;
};

template <int MAX_K>
struct KmerAndCtgLoc {
  Kmer<MAX_K> kmer;
  CtgLoc ctg_loc;
  UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
};

using CtgAndReadLocsMap = HASH_TABLE<cid_t, vector<CtgAndReadLoc>>;

struct ReadRecord {
  int index;
  int rlen;

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

template <int MAX_K>
struct KmerWithTarget {
  Kmer<MAX_K> kmer;
  int target;
  int seq_offset;
  int read_i;
};

template <int MAX_K>
struct KmersReadsBuffer {
  vector<Kmer<MAX_K>> kmers;
  vector<ReadRecord *> read_records;

  void add(const Kmer<MAX_K> &kmer, ReadRecord *read_record) {
    kmers.push_back(kmer);
    read_records.push_back(read_record);
  }

  size_t size() { return kmers.size(); }

  void clear() {
    kmers.clear();
    read_records.clear();
  }
};

struct RgetRequest {
  CtgLoc ctg_loc;
  int get_start;
  int get_len;
  string rname;
  string read_subseq;
  string ctg_seq;
  int rstart;
  int rlen;
  int cstart;
  char orient;
  int overlap_len;
  int read_group_id;
};

static int seed_space = 1;

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
      : kmer_map({})
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
      if (allow_multi_kmers) {
        if (it == kmer_map->end()) it = kmer_map->insert({kmer_and_ctg_loc.kmer, {}}).first;
        it->second.push_back(ctg_loc);
      } else {
        if (it == kmer_map->end()) {
          it = kmer_map->insert({kmer_and_ctg_loc.kmer, {}}).first;
          it->second.push_back(ctg_loc);
        } else {
          // there are conflicts so don't allow any kmer mappings
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

  void reserve_ctg_seqs(size_t sz) { global_ctg_seqs.reserve(sz); }

  global_ptr<char> add_ctg_seq(string seq) {
    auto seq_gptr = upcxx::allocate<char>(seq.length() + 1);
    global_ctg_seqs.push_back(seq_gptr);  // remember to dealloc!
    strcpy(seq_gptr.local(), seq.c_str());
    return seq_gptr;
  }

  void reserve(int64_t mysize) { kmer_map->reserve(mysize); }

  int64_t size() const { return kmer_map->size(); }

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
    BarrierTimer timer(__FILEFUNC__, false);  // barrier on exit, not entrance
    kmer_store.flush_updates();
    kmer_store.clear();
  }

  future<vector<KmerAndCtgLoc<MAX_K>>> get_ctgs_with_kmers(int target_rank, vector<Kmer<MAX_K>> &kmers) {
    DBG_VERBOSE("Sending request for ", kmers.size(), " to ", target_rank, "\n");
    return rpc(
        target_rank,
        [](vector<Kmer<MAX_K>> kmers, kmer_map_t &kmer_map) {
          vector<KmerAndCtgLoc<MAX_K>> kmer_ctg_locs;
          kmer_ctg_locs.reserve(kmers.size());
          for (auto &kmer : kmers) {
            assert(kmer.is_least());
            assert(kmer.is_valid());
            const auto it = kmer_map->find(kmer);
            if (it == kmer_map->end()) continue;
            for (auto &ctg_loc : it->second) {
              kmer_ctg_locs.push_back({kmer, ctg_loc});
            }
          }
          DBG_VERBOSE("processed get_ctgs_with_kmers ", kmers.size(), " ", get_size_str(kmers.size() * sizeof(Kmer<MAX_K>)),
                      ", returning ", kmer_ctg_locs.size(), " ", get_size_str(kmer_ctg_locs.size() * sizeof(KmerAndCtgLoc<MAX_K>)),
                      "\n");
          return kmer_ctg_locs;
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

class CtgCache {
 private:
  std::unordered_set<cid_t> fetched_ctgs;
  int64_t refetches = 0;

 public:
  FixedSizeCache<cid_t, string> cache;
  int64_t hits = 0;
  int64_t lookups = 0;
  int64_t extended_fetches = 0;
  int64_t local_hits = 0;

  void init(size_t size) {
    cache.set_invalid_key(std::numeric_limits<cid_t>::max());
    cache.reserve(size);
  }

  void clear() {
    cache.clear();
    fetched_ctgs.clear();
  }

  void print_stats() {
    auto msm_ctg_local_hits = min_sum_max_reduce_one(local_hits, 0).wait();
    auto msm_ctg_cache_hits = min_sum_max_reduce_one(hits, 0).wait();
    auto msm_ctg_lookups = min_sum_max_reduce_one(lookups + local_hits, 0).wait();
    auto msm_ctg_refetches = min_sum_max_reduce_one(refetches, 0).wait();
    auto msm_ctg_extended_fetches = min_sum_max_reduce_one(extended_fetches, 0).wait();
    auto msm_ctg_clobberings = min_sum_max_reduce_one(cache.get_clobberings(), 0).wait();
    SLOG_VERBOSE("Hits on ctg cache: ", perc_str(msm_ctg_cache_hits.sum, msm_ctg_lookups.sum), " (my) cache size/load ",
                 perc_str(cache.capacity(), cache.size()), " clobberings ", cache.get_clobberings(), "\n");
    SLOG_VERBOSE("Local contig hits bypassing cache: ", perc_str(msm_ctg_local_hits.sum, msm_ctg_lookups.sum), "\n");
    SLOG_VERBOSE("ctg_local_hits: ", msm_ctg_local_hits.to_string(), "\n");
    SLOG_VERBOSE("ctg_cache_hits: ", msm_ctg_cache_hits.to_string(), "\n");
    SLOG_VERBOSE("ctg_lookups: ", msm_ctg_lookups.to_string(), "\n");
    SLOG_VERBOSE("ctg_refetches: ", msm_ctg_refetches.to_string(), "\n");
    SLOG_VERBOSE("ctg_extended_fetches: ", msm_ctg_extended_fetches.to_string(), "\n");
    SLOG_VERBOSE("ctg_cache_clobberings: ", msm_ctg_clobberings.to_string(), "\n");
  }

  void insert(RgetRequest &req) {
    auto itpair = cache.insert({req.ctg_loc.cid, req.ctg_seq});
    // unsuccessful insert - possible because we are fetching many locations at once
    // if (!itpair.second) DIE("Duplicate cid in ctg cache insert: ", req.ctg_loc.cid, "; should not be possible\n");
    if (fetched_ctgs.find(req.ctg_loc.cid) == fetched_ctgs.end())
      fetched_ctgs.insert(req.ctg_loc.cid);
    else
      refetches++;
  }
};

class Aligner {
  int64_t num_alns;
  int64_t num_perfect_alns;
  int64_t num_overlaps;
  int kmer_len;

  vector<Aln> kernel_alns;
  vector<string> ctg_seqs;
  vector<string> read_seqs;

  future<> active_kernel_fut;

  int64_t max_clen = 0;
  int64_t max_rlen = 0;
  CPUAligner cpu_aligner;

  CtgCache ctg_cache;
  int64_t ctg_bytes_fetched = 0;
  int rget_calls = 0;

  Alns *alns;

 private:
  vector<vector<RgetRequest>> rget_requests;

  void align_read(const string &rname, int64_t cid, const string_view &rseq, const string_view &ctg_seq, int rstart, int rlen,
                  int cstart, int clen, char orient, int overlap_len, int read_group_id) {
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
      alns->add_aln(aln, cpu_aligner.allow_multi);
    } else {
      max_clen = max((int64_t)cseq.size(), max_clen);
      max_rlen = max((int64_t)rseq.size(), max_rlen);
      int64_t num_alns = kernel_alns.size() + 1;
      int64_t tot_mem_est = num_alns * (max_clen + max_rlen + 2 * sizeof(int) + 5 * sizeof(short));
      // contig is the ref, read is the query - done this way so that we can potentially do multiple alns to each read
      // this is also the way it's done in meraligner
      kernel_alns.emplace_back(rname, cid, 0, 0, rlen, cstart, 0, clen, orient);
      ctg_seqs.emplace_back(cseq);
      read_seqs.emplace_back(rseq);
      if (num_alns >= KLIGN_GPU_BLOCK_SIZE) {
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
    ctg_cache.clear();
  }

  void do_rget_irregular(int target, bool insert_into_cache) {
    vector<pair<global_ptr<char>, size_t>> src;
    vector<pair<char *, size_t>> dest;
    for (auto &req : rget_requests[target]) {
      src.push_back({req.ctg_loc.seq_gptr + req.get_start, req.get_len});
      dest.push_back({const_cast<char *>(req.ctg_seq.data()) + req.get_start, req.get_len});
      ctg_bytes_fetched += req.get_len;
    }
    // SLOG("flush to ", target, ": rget_irregular with ", src.size(), " elements\n");
    timers.rget_ctg_seqs.start();
    rget_irregular(src.begin(), src.end(), dest.begin(), dest.end()).wait();
    timers.rget_ctg_seqs.stop();
    rget_calls++;
    for (auto &req : rget_requests[target]) {
      align_read(req.rname, req.ctg_loc.cid, req.read_subseq, req.ctg_seq, req.rstart, req.rlen, req.cstart, req.ctg_loc.clen,
                 req.orient, req.overlap_len, req.read_group_id);
      num_alns++;
      if (insert_into_cache) ctg_cache.insert(req);
    }
    rget_requests[target].clear();
  }

 public:
  Aligner(int kmer_len, Alns &alns, int rlen_limit, bool allow_multi, bool compute_cigar, bool use_blastn_scores,
          int64_t all_num_ctgs)
      : num_alns(0)
      , num_perfect_alns(0)
      , num_overlaps(0)
      , kmer_len(kmer_len)
      , kernel_alns({})
      , ctg_seqs({})
      , read_seqs({})
      , active_kernel_fut(make_future())
      , cpu_aligner(allow_multi, compute_cigar, use_blastn_scores)
      , alns(&alns) {
    ctg_cache.init(3 * all_num_ctgs / rank_n() + 1024);
    init_aligner((int)cpu_aligner.ssw_aligner.get_match_score(), (int)cpu_aligner.ssw_aligner.get_mismatch_penalty(),
                 (int)cpu_aligner.ssw_aligner.get_gap_opening_penalty(), (int)cpu_aligner.ssw_aligner.get_gap_extending_penalty(),
                 (int)cpu_aligner.ssw_aligner.get_ambiguity_penalty(), rlen_limit);
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

  int64_t get_num_overlaps(bool all = false) {
    if (!all) return reduce_one(num_overlaps, op_fast_add, 0).wait();
    return reduce_all(num_overlaps, op_fast_add).wait();
  }

  void clear_aln_bufs() {
    kernel_alns.clear();
    ctg_seqs.clear();
    read_seqs.clear();
    max_clen = 0;
    max_rlen = 0;
  }

  void flush_remaining(int read_group_id) {
    BaseTimer t(__FILEFUNC__);
    t.start();
    for (auto target : upcxx_utils::foreach_rank_by_node()) {
      if (rget_requests[target].size()) do_rget_irregular(target, false);
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

  void compute_alns_for_read(CtgAndReadLocsMap *aligned_ctgs_map, const string &rname, const string &rseq_fw, int read_group_id) {
    int rlen = rseq_fw.length();
    string rseq_rc;
    string tmp_ctg;
    for (auto &ctg_and_read_locs : *aligned_ctgs_map) {
      progress();
      for (auto &ctg_and_read_loc : ctg_and_read_locs.second) {
        int pos_in_read = ctg_and_read_loc.pos_in_read;
        bool read_kmer_is_rc = ctg_and_read_loc.read_is_rc;
        auto &ctg_loc = ctg_and_read_loc.ctg_loc;
        char orient = '+';
        string rseq_ptr;
        if (ctg_loc.is_rc != read_kmer_is_rc) {
          // it's revcomp in either contig or read, but not in both or neither
          orient = '-';
          pos_in_read = rlen - (kmer_len + pos_in_read);
          if (rseq_rc.empty()) rseq_rc = revcomp(rseq_fw);
          rseq_ptr = string(rseq_rc);
        } else {
          rseq_ptr = string(rseq_fw);
        }
        auto [cstart, rstart, overlap_len] = get_start_positions(kmer_len, ctg_loc, pos_in_read, rlen);
        // use the whole read, to account for possible indels
        string read_subseq = rseq_ptr.substr(0, rlen);

        assert(cstart >= 0 && cstart + overlap_len <= ctg_loc.clen);
        assert(overlap_len <= 2 * rlen);

        bool on_node = ctg_loc.seq_gptr.is_local();
#ifdef DEBUG
        // test both on node and off node ctg cache
        if (on_node) on_node = (ctg_loc.seq_gptr.where() % 2) == (rank_me() % 2);
#endif
        if (on_node) {
          // on same node already
          ctg_cache.local_hits++;
          string_view ctg_seq = string_view(ctg_loc.seq_gptr.local(), ctg_loc.clen);
          align_read(rname, ctg_loc.cid, read_subseq, ctg_seq, rstart, rlen, cstart, ctg_loc.clen, orient, overlap_len,
                     read_group_id);
          num_alns++;
        } else {
          bool found = false;
          ctg_cache.lookups++;
          auto it = ctg_cache.cache.find(ctg_loc.cid);
          auto get_start = cstart, get_len = overlap_len;
          if (it != ctg_cache.cache.end()) {
            string &ctg_str = it->second;  // the actual underlying string, not view
            string_view ctg_seq = string_view(ctg_str.data(), ctg_str.size());
            found = true;
            // find the first and last blank within the overlap region on cached contig (if any)
            for (int i = 0; i < overlap_len; i++) {
              if (ctg_seq[get_start] == ' ') {
                found = false;
                break;
              }
              get_start++;
              get_len--;
            }
            if (!found) {
              ctg_cache.extended_fetches++;
              while (get_len > 0) {
                if (ctg_seq[get_start + get_len - 1] == ' ') break;
                get_len--;
              }
            } else {
              ctg_cache.hits++;
              align_read(rname, ctg_loc.cid, read_subseq, ctg_seq, rstart, rlen, cstart, ctg_loc.clen, orient, overlap_len,
                         read_group_id);
            }
          }
          if (!found) {
            string ctg_str(ctg_loc.clen, ' ');
            // also get extra bordering blank bases on either side of the contig for negligable extra overhead and likely fewer
            // rgets
            const int extra_bases = 384;
            for (int i = 0; i < extra_bases; i++) {
              if (get_start == 0) break;
              if (ctg_str[get_start - 1] == ' ') {
                get_start--;
                get_len++;
              } else {
                break;
              }
            }
            for (int i = 0; i < extra_bases; i++) {
              if (get_start + get_len >= ctg_str.size()) break;
              if (ctg_str[get_start + get_len] == ' ')
                get_len++;
              else
                break;
            }
            // finally extend this fetch to the end of the contig if a small percentage of the contig remains unfetched
            int edge_bases = ctg_str.size() * 0.1;
            if (get_start < edge_bases) {
              get_len += get_start;
              get_start = 0;
            }
            if (ctg_str.size() - get_start - get_len < edge_bases) {
              get_len = ctg_str.size() - get_start;
            }
            assert(get_start >= 0);
            assert(get_start + get_len <= ctg_str.size());

            auto target = ctg_loc.seq_gptr.where();
            RgetRequest req = {ctg_loc, get_start, get_len, rname,  read_subseq, ctg_str,
                               rstart,  rlen,      cstart,  orient, overlap_len, read_group_id};
            rget_requests[target].push_back(req);
            if (rget_requests[target].size() == KLIGN_RGET_BUF_SIZE) do_rget_irregular(target, true);
          }
        }
      }
    }
  }

  void sort_alns() {
    if (!kernel_alns.empty()) {
      DIE("sort_alns called while alignments are still pending to be processed - ", kernel_alns.size());
    }
    if (!active_kernel_fut.ready()) {
      SWARN("Waiting for active_kernel - has flush_remaining() been called?\n");
    }
    active_kernel_fut.wait();
    alns->sort_alns().wait();
  }

  void log_ctg_bytes_fetched() {
    auto all_rget_calls = reduce_one(rget_calls, op_fast_add, 0).wait();
    auto all_ctg_bytes_fetched = reduce_one(ctg_bytes_fetched, op_fast_add, 0).wait();
    auto max_ctg_bytes_fetched = reduce_one(ctg_bytes_fetched, op_fast_max, 0).wait();
    SLOG_VERBOSE("Contig bytes fetched ", get_size_str(all_ctg_bytes_fetched), " balance ",
                 (double)all_ctg_bytes_fetched / (rank_n() * max_ctg_bytes_fetched), "\n");
    ctg_bytes_fetched = 0;
  }
};

template <int MAX_K>
static void build_alignment_index(KmerCtgDHT<MAX_K> &kmer_ctg_dht, Contigs &ctgs, unsigned min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  ProgressBar progbar(ctgs.size(), "Extracting seeds from contigs");
  // estimate and reserve room in the local map to avoid excessive reallocations
  int64_t est_num_kmers = 0;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    auto len = ctg->seq.length();
    if (len < min_ctg_len) continue;
    est_num_kmers += len - kmer_ctg_dht.kmer_len + 1;
  }
  est_num_kmers = upcxx::reduce_all(est_num_kmers, upcxx::op_fast_add).wait();
  auto my_reserve = 1.2 * est_num_kmers / rank_n() + 2000;  // 120% to keep the map fast
  kmer_ctg_dht.reserve(my_reserve);
  kmer_ctg_dht.reserve_ctg_seqs(ctgs.size());
  size_t ctg_seq_lengths = 0, min_len_ctgs = 0;
  vector<Kmer<MAX_K>> kmers;
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() < min_ctg_len) continue;
    ctg_seq_lengths += ctg->seq.length();
    min_len_ctgs++;
    global_ptr<char> seq_gptr = kmer_ctg_dht.add_ctg_seq(ctg->seq);
    CtgLoc ctg_loc = {.cid = ctg->id, .seq_gptr = seq_gptr, .clen = (int)ctg->seq.length(), .depth = (float)ctg->depth};
    Kmer<MAX_K>::get_kmers(kmer_ctg_dht.kmer_len, string_view(ctg->seq.data(), ctg->seq.size()), kmers, true);
    num_kmers += kmers.size();
    for (unsigned i = 0; i < kmers.size(); i++) {
      // if (kmers[i].to_string() == "ACATCTACCGCTAGAGGATTA")
      //   WARN("kmer ", kmers[i].to_string(), " found in contig ", ctg->id, " in position ", i);
      ctg_loc.pos = i;
      if (!kmers[i].is_valid()) continue;
      kmer_ctg_dht.add_kmer(kmers[i], ctg_loc);
    }
    progress();
  }
  auto fut = progbar.set_done();
  kmer_ctg_dht.flush_add_kmers();
  auto tot_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto tot_num_ctgs = reduce_one(min_len_ctgs, op_fast_add, 0).wait();
  auto tot_ctg_lengths = reduce_one(ctg_seq_lengths, op_fast_add, 0).wait();
  fut.wait();
  auto num_kmers_in_ht = kmer_ctg_dht.get_num_kmers();
  LOG("Estimated room for ", my_reserve, " my final count ", kmer_ctg_dht.size(), "\n");
  SLOG_VERBOSE("Total contigs >= ", min_ctg_len, ": ", tot_num_ctgs, " seq_length: ", tot_ctg_lengths, "\n");
  SLOG_VERBOSE("Processed ", tot_num_kmers, " seeds from contigs, added ", num_kmers_in_ht, "\n");
  auto num_dropped_seed_to_ctgs = kmer_ctg_dht.get_num_dropped_seed_to_ctgs();
  if (num_dropped_seed_to_ctgs)
    SLOG_VERBOSE("For k = ", kmer_ctg_dht.kmer_len, " dropped ", num_dropped_seed_to_ctgs, " non-unique seed-to-contig mappings (",
                 setprecision(2), fixed, (100.0 * num_dropped_seed_to_ctgs / tot_num_kmers), "%)\n");
}

template <int MAX_K>
static upcxx::future<> fetch_ctg_maps_for_target(int target_rank, KmerCtgDHT<MAX_K> &kmer_ctg_dht,
                                                 KmersReadsBuffer<MAX_K> &kmers_reads_buffer, int64_t &num_alns,
                                                 int64_t &num_excess_alns_reads, int64_t &bytes_sent, int64_t &bytes_received,
                                                 int64_t &num_rpcs) {
  return make_future<>();
}

template <int MAX_K>
void fetch_ctg_maps(KmerCtgDHT<MAX_K> &kmer_ctg_dht, PackedReads *packed_reads, vector<ReadRecord> &read_records) {
  timers.fetch_ctg_maps.start();
  int64_t bytes_sent = 0;
  int64_t bytes_received = 0;
  int64_t num_reads = 0;
  int64_t num_alns = 0;
  int64_t num_excess_alns_reads = 0;
  int64_t num_kmers = 0;
  int64_t num_rpcs = 0;

  upcxx::future<> fetch_fut = make_future();
  vector<KmersReadsBuffer<MAX_K>> kmers_reads_buffers(rank_n());
  vector<Kmer<MAX_K>> kmers;
  int kmer_len = Kmer<MAX_K>::get_k();

  ProgressBar progbar(packed_reads->get_local_num_reads(), "Fetching ctg maps for alignments");
  for (int ri = 0; ri < packed_reads->get_local_num_reads(); ri++) {
    progress();
    progbar.update();
    string read_seq;
    packed_reads->get_read_seq(ri, read_seq);
    // this happens when a placeholder read with just a single N character is added after merging reads
    if (kmer_len > read_seq.length()) continue;
    num_reads++;
    num_kmers += read_seq.length() - kmer_len;
    read_records[ri] = ReadRecord({ri, (int)read_seq.length()});
    Kmer<MAX_K>::get_kmers(kmer_len, string_view(read_seq.data(), read_seq.size()), kmers, true);
    if (!kmers.size()) continue;
    for (int i = 0; i < kmers.size(); i += seed_space) {
      const Kmer<MAX_K> &kmer_fw = kmers[i];
      if (!kmer_fw.is_valid()) continue;
      num_kmers++;
      const Kmer<MAX_K> kmer_rc = kmer_fw.revcomp();
      const Kmer<MAX_K> *kmer_lc = &kmer_fw;
      assert(kmer_fw.is_valid() && kmer_rc.is_valid());
      bool is_rc = false;
      if (kmer_rc < kmer_fw) {
        kmer_lc = &kmer_rc;
        is_rc = true;
      }
      assert(kmer_lc->is_least());
      auto target = kmer_ctg_dht.get_target_rank(*kmer_lc);
      kmers_reads_buffers[target].add(*kmer_lc, &read_records[ri]);
      if (kmers_reads_buffers[target].size() >= KLIGN_KMERS_BUF_SIZE / rank_n())
        fetch_fut = when_all(fetch_ctg_maps_for_target(target, kmer_ctg_dht, kmers_reads_buffers[target], num_alns,
                                                       num_excess_alns_reads, bytes_sent, bytes_received, num_rpcs),
                             fetch_fut);
    }
    kmers.clear();
  }
  for (auto target : upcxx_utils::foreach_rank_by_node()) {  // stagger by rank_me, round robin by node
    if (kmers_reads_buffers[target].size())
      fetch_fut = when_all(fetch_ctg_maps_for_target(target, kmer_ctg_dht, kmers_reads_buffers[target], num_alns,
                                                     num_excess_alns_reads, bytes_sent, bytes_received, num_rpcs),
                           fetch_fut);
  }
  when_all(fetch_fut, progbar.set_done()).wait();
  upcxx_utils::flush_outstanding_futures();
  barrier();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_num_kmers = reduce_one(num_kmers, op_fast_add, 0).wait();
  auto kmers_bytes_sent = all_num_kmers * Kmer<MAX_K>::get_N_LONGS() * sizeof(longs_t);
  auto all_bytes_sent = reduce_one(bytes_sent, op_fast_add, 0).wait();
  auto all_bytes_received = reduce_one(bytes_received, op_fast_add, 0).wait();
  auto all_num_rpcs = reduce_one(num_rpcs, op_fast_add, 0).wait();
  auto all_excess_alns_reads = reduce_one(num_excess_alns_reads, op_fast_add, 0).wait();
  if (rank_me() == 0) {
    SLOG_VERBOSE("Parsed ", all_num_reads, " reads and extracted ", all_num_kmers, " kmers\n");
    SLOG_VERBOSE("Sent ", get_size_str(all_bytes_sent), " (", get_size_str(all_bytes_sent / all_num_rpcs),
                 " avg msg) of kmers and received ", get_size_str(all_bytes_received), " (",
                 get_size_str(all_bytes_received / all_num_rpcs), " avg msg)\n");
    if (all_excess_alns_reads)
      SLOG_VERBOSE("Dropped ", all_excess_alns_reads, " alignments in excess of ", KLIGN_MAX_ALNS_PER_READ, " per read\n");
  }
  timers.fetch_ctg_maps.stop();
  barrier();
}

template <int MAX_K>
double find_alignments(unsigned kmer_len, PackedReadsList &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool use_kmer_cache, bool compute_cigar,
                       bool use_blastn_scores, int min_ctg_len) {
  BarrierTimer timer(__FILEFUNC__);
  Kmer<MAX_K>::set_k(kmer_len);
  SLOG_VERBOSE("Aligning with seed size of ", kmer_len, "\n");
  int64_t all_num_ctgs = reduce_all(ctgs.size(), op_fast_add).wait();
  bool allow_multi_kmers = compute_cigar;
  uint64_t num_ctg_kmers = 0;
  for (auto &ctg : ctgs)
    if (ctg.seq.length() >= Kmer<MAX_K>::get_k()) num_ctg_kmers += ctg.seq.length() - Kmer<MAX_K>::get_k() + 1;
  KmerCtgDHT<MAX_K> kmer_ctg_dht(max_store_size, max_rpcs_in_flight, allow_multi_kmers, num_ctg_kmers);
  barrier();
  build_alignment_index(kmer_ctg_dht, ctgs, min_ctg_len);
#ifdef DEBUG
// kmer_ctg_dht.dump_ctg_kmers();
#endif
  int read_group_id = 0;
  for (auto packed_reads : packed_reads_list) {
    vector<ReadRecord> read_records(packed_reads->get_local_num_reads());
    fetch_ctg_maps(kmer_ctg_dht, packed_reads, read_records);
    //    compute_alns<MAX_K>(packed_reads, read_records, alns, read_group_id, rlen_limit, post_asm, compute_cigar,
    //    use_blastn_scores,
    //                        all_num_ctgs);
    read_group_id++;
  }
  barrier();
  timers.done_all();
  double aln_kernel_elapsed = timers.aln_kernel.get_elapsed();
  timers.clear();
  barrier();
  auto num_alns = alns.size();
  auto num_dups = alns.get_num_dups();
  if (num_dups) SLOG_VERBOSE("Number of duplicate alignments ", perc_str(num_dups, num_alns), "\n");
  barrier();
  return aln_kernel_elapsed;
}
