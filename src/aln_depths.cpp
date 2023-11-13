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

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "contigs.hpp"
#include "upcxx_utils/flat_aggr_store.hpp"
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "upcxx_utils/limit_outstanding.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx_utils;

template <typename T>
struct AvgVar {
  T avg, var;
  AvgVar()
      : avg(0)
      , var(0) {}
};  // template struct AvgVar

struct CtgLen {
  cid_t cid;
  int clen;
  int owning_rank;
};

struct CtgStats {
  cid_t cid;
  vector<AvgVar<float>> rg_stats;
  UPCXX_SERIALIZED_FIELDS(cid, rg_stats);
};

struct CtgBaseDepths {
  using base_count_t = uint16_t;
  using read_group_base_count_t = vector<base_count_t *>;  // lazy allocate when read group applies alignments data
  cid_t cid;
  int clen;
  int owning_rank;
  read_group_base_count_t read_group_base_counts;
  vector<AvgVar<float>> rg_stats;  // num_read_groups

  // default empty for serialization
  CtgBaseDepths()
      : cid{}
      , clen{}
      , owning_rank{}
      , read_group_base_counts{}
      , rg_stats{} {}
  // for owning rank's placeholder entry in dist_object
  CtgBaseDepths(cid_t cid)
      : cid{cid}
      , clen{}
      , owning_rank(rank_me())
      , read_group_base_counts{}
      , rg_stats{} {}
  // for target ranks's calculations
  CtgBaseDepths(cid_t cid, int num_read_groups, int clen, int owning_rank)
      : cid(cid)
      , clen(clen)
      , owning_rank(owning_rank)
      , read_group_base_counts(num_read_groups, nullptr)
      , rg_stats(num_read_groups) {}
  // for update function
  CtgBaseDepths(const CtgStats &ctg_stats)
      : cid(ctg_stats.cid)
      , clen{}
      , owning_rank{}
      , read_group_base_counts{}
      , rg_stats(ctg_stats.rg_stats) {}
  ~CtgBaseDepths() { clear_rg_bases(); }
  void clear_rg_bases() {
    for (auto rgbc_ptr : read_group_base_counts) assert(rgbc_ptr == nullptr && "All allocations are freed");
    if (!read_group_base_counts.empty()) read_group_base_count_t().swap(read_group_base_counts);
  }
  base_count_t *get_base_counts(int read_group_id) {
    assert(read_group_id >= 0 && read_group_id < read_group_base_counts.size());
    auto &rgbc_ptr = read_group_base_counts[read_group_id];
    if (!rgbc_ptr) {
      rgbc_ptr = new base_count_t[clen]();  //  allocate per-base counts for this read group
      for (auto i = 0; i < clen; i++) assert(rgbc_ptr[i] == 0);
      DBG("Allocated new read group base counts read_group_id=", read_group_id, " prt=", rgbc_ptr, " clen=", clen, "\n");
    }
    return rgbc_ptr;
  }
  void free_base_counts(int read_group_id) {
    assert(read_group_id >= 0 && read_group_id < read_group_base_counts.size());
    auto &rgbc_ptr = read_group_base_counts[read_group_id];
    if (rgbc_ptr) {
      DBG("Deleted read group base counts read_group_id=", read_group_id, " prt=", rgbc_ptr, "\n");
      delete[] rgbc_ptr;  // free memory
    }
    rgbc_ptr = nullptr;
  }

  void add_alignment_depth(int read_group_id, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
    auto ctg_base_counts = get_base_counts(read_group_id);
    DBG("cid=", cid, " len=", aln_stop - aln_start, " counting aln_start=", aln_start, " aln_stop=", aln_stop,
        " read_group_id=", read_group_id, " ctg_base_counts=", ctg_base_counts, " merge_start=", aln_merge_start,
        " merge_stop=", aln_merge_stop, "\n");
    // DBG_VERBOSE("cid=", cid, " counting aln_start=", aln_start, " aln_stop=", aln_stop, " read_group_id=", read_group_id,
    //    " contig.size()=", ctg_base_counts.size(), "\n");
    assert(aln_start >= 0 && "Align start >= 0");
    assert(aln_stop <= clen && "Align stop <= contig.size()");
    for (int i = aln_start; i < aln_stop; i++) {
      ctg_base_counts[i]++;
    }
    // disabled by default -- counts are "better" if this does not happen
    // instead, count *insert coverage* not *read coverage*
    // insert coverage should be closer to a Poisson distribution

    // this is the merged region in a merged read - will be counted double
    if (aln_merge_start != -1 && aln_merge_stop != -1) {
      assert(aln_merge_start >= 0 && "merge start >= 0");
      assert(aln_merge_stop <= clen && "merge_stop <= size");
      for (int i = aln_merge_start; i < aln_merge_stop; i++) {
        ctg_base_counts[i]++;
      }
    }
  }  // add_alignment_depth

  void calc_stats(int read_group_id, int edge_base_len) {
    DBG("Calc stats on cid=", cid, " rg=", read_group_id, " edge_base_len=", edge_base_len, "\n");
    assert(read_group_id >= 0);
    assert(!read_group_base_counts.empty());
    assert(read_group_id < read_group_base_counts.size());
    auto rgbc_ptr = read_group_base_counts[read_group_id];

    AvgVar<double> stats{};  // use double while calculating.
    size_t adjusted_clen = clen - 2 * edge_base_len;
    if (adjusted_clen <= 0) adjusted_clen = 1;
    if (rgbc_ptr != nullptr) {
      auto &ctg_base_depths = rgbc_ptr;
      for (int i = edge_base_len; i < (int)clen - edge_base_len; i++) {
        stats.avg += ctg_base_depths[i];
      }
      stats.avg /= adjusted_clen;

      for (int i = edge_base_len; i < (int)clen - edge_base_len; i++) {
        stats.var += pow((double)ctg_base_depths[i] - stats.avg, 2.0);
      }
      stats.var /= adjusted_clen;
    }
    rg_stats[read_group_id].avg = stats.avg;
    rg_stats[read_group_id].var = stats.var;
    DBG("Calc on cid=", cid, " read_group_id=", read_group_id, " rgbc_ptr=", rgbc_ptr, " avg=", stats.avg, " var=", stats.var,
        " adjusted_len=", adjusted_clen, " clen=", clen, "\n");
    free_base_counts(read_group_id);
  }  // calc_stats

};   // struct CtgBaseDepths

class CtgsDepths {
 private:
  using local_ctgs_depths_map_t = HASH_TABLE<cid_t, CtgBaseDepths>;
  using ctgs_depths_map_t = upcxx::dist_object<local_ctgs_depths_map_t>;
  ctgs_depths_map_t ctgs_depths;
  Contigs &ctgs;
  int edge_base_len, min_ctg_len, num_read_groups, finished_read_groups;
  size_t max_store_bytes;
  vector<upcxx_utils::PromiseBarrier> all_done_rg_prom_barriers;
  upcxx::future<> fut_done;
  HASH_TABLE<cid_t, CtgBaseDepths>::iterator ctgs_depths_iter;
  struct CtgAlnDepth {
    cid_t cid;
    int read_group_id, aln_start, aln_stop, aln_merge_start, aln_merge_stop;
  };
  ThreeTierAggrStore<CtgAlnDepth> ctg_aln_depth_store;
  ThreeTierAggrStore<CtgLen> add_ctg_store;
  ThreeTierAggrStore<CtgStats> return_stats_store;

  static size_t get_target_rank(cid_t cid) { return std::hash<cid_t>{}(cid) % upcxx::rank_n(); }

  void _add_ctg(const CtgLen &ctg_len) {
    CtgBaseDepths newctg(ctg_len.cid, num_read_groups, ctg_len.clen, ctg_len.owning_rank);
    ctgs_depths->insert({ctg_len.cid, std::move(newctg)});
    DBG("Added cid=", ctg_len.cid, " for calcs\n");
  }

  void _update_ctg_aln_depth(cid_t cid, int read_group_id, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
    const auto it = ctgs_depths->find(cid);
    if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
    CtgBaseDepths &cbd = it->second;
    cbd.add_alignment_depth(read_group_id, aln_start, aln_stop, aln_merge_start, aln_merge_stop);
  }

  void _return_stats(const CtgStats &ctg_stats) {
    auto cid = ctg_stats.cid;
    const auto it = ctgs_depths->find(cid);
    if (it == ctgs_depths->end()) DIE("Expected that returned contig depths already is already in place ", cid, "\n");
    auto &cbd = it->second;
    assert(cbd.rg_stats.empty());
    assert(cbd.owning_rank == rank_me());
    assert(cbd.read_group_base_counts.empty());
    cbd.rg_stats = ctg_stats.rg_stats;
  }

 public:
  CtgsDepths(Contigs &ctgs, int edge_base_len, int min_ctg_len, int num_read_groups)
      : ctgs_depths(local_ctgs_depths_map_t{})
      , ctgs(ctgs)
      , edge_base_len(edge_base_len)
      , min_ctg_len(min_ctg_len)
      , num_read_groups(num_read_groups == 0 ? 1 : num_read_groups)
      , finished_read_groups(0)
      , ctg_aln_depth_store()
      , add_ctg_store()
      , all_done_rg_prom_barriers(num_read_groups)
      , fut_done{} {
    ctgs_depths->reserve(ctgs.size() * 2 + 2000);  // entries for self + distributed calcs
    int64_t mem_to_use = 0.1 * get_free_mem(true) / local_team().rank_n();
    max_store_bytes = std::max(mem_to_use, (int64_t)sizeof(CtgAlnDepth) * 100);

    add_ctg_store.set_size("Add Ctgs", max_store_bytes);
    add_ctg_store.set_update_func([&self = *this](const CtgLen &ctg_len) { self._add_ctg(ctg_len); });

    ctg_aln_depth_store.set_update_func([&self = *this](const CtgAlnDepth &cad) {
      self._update_ctg_aln_depth(cad.cid, cad.read_group_id, cad.aln_start, cad.aln_stop, cad.aln_merge_start, cad.aln_merge_stop);
    });

    return_stats_store.set_update_func([&self = *this](const CtgStats &ctg_stats) { self._return_stats(ctg_stats); });

    // PromiseBarriers for when all ranks have completed a read group
    fut_done = make_future();
    int read_group_id = 0;
    for (auto &rg_prom : all_done_rg_prom_barriers) {
      auto fut_calc = rg_prom.get_future().then([&self = *this, read_group_id]() {
        DBG("Calculating stats on read_group=", read_group_id, "\n");
        for (auto &ctgs_depth : *self.ctgs_depths) {
          auto &[cid, ctg_base_depths] = ctgs_depth;
          if (ctg_base_depths.owning_rank == rank_me() && get_target_rank(cid) != rank_me() &&
              ctg_base_depths.read_group_base_counts.empty())
            continue;  // ignore placeholders here
          ctg_base_depths.calc_stats(read_group_id, self.edge_base_len);
        }
        DBG("Done calculating stats on read_group=", read_group_id, "\n");
      });
      fut_done = when_all(fut_done, fut_calc);
      read_group_id++;
    }
  }
  ~CtgsDepths() { assert(finished_read_groups == num_read_groups); }

  void finish_read_group(int read_group_id) {
    assert(read_group_id >= 0 && read_group_id < num_read_groups && finished_read_groups < num_read_groups);
    DBG("Finishing rg=", read_group_id, "\n");
    all_done_rg_prom_barriers[read_group_id].fulfill();
    all_done_rg_prom_barriers[read_group_id].get_future().wait();
    DBG("Finished rg=", read_group_id, "\n");
    finished_read_groups++;
  }

  void finish_all() {
    BarrierTimer timer(__FILEFUNC__);
    DBG("Finishing all\n");
    ctg_aln_depth_store.clear();
    return_stats_store.set_size("Return stats",
                                max_store_bytes * sizeof(CtgStats) / (sizeof(CtgStats) + num_read_groups * 2 * sizeof(float)));
    fut_done.wait();
    if (num_read_groups != finished_read_groups)
      DIE("finish_all cannot be called before all read groups have called finish_read_group");
    Timings::set_pending(fut_done);

    SLOG_VERBOSE("Sending read group stats back to original rank\n");
    for (auto &key_val : *ctgs_depths) {
      auto &[cid, cbds] = key_val;
      auto tgt = get_target_rank(cid);
      if (tgt != rank_me()) continue;  // placeholder for one of my loaded contigs to be updated
      tgt = cbds.owning_rank;
      if (tgt == rank_me()) {
        // bypass.  All good
      } else {
        CtgStats cs{.cid = cid, .rg_stats = std::move(cbds.rg_stats)};
        return_stats_store.update(tgt, cs);
      }
      cbds.clear_rg_bases();  // free some memory
    }
    return_stats_store.flush_updates();
    return_stats_store.clear();
  }

  void flush_contigs() {
    DBG("Flushing contigs\n");
    add_ctg_store.flush_updates();
    add_ctg_store.clear();
    ctg_aln_depth_store.set_size("Ctg Aln Depths", max_store_bytes);
  }

  void flush_ctg_alns() {
    DBG("Flush ctg alns\n");
    ctg_aln_depth_store.flush_updates();
  }

  upcxx::future<int64_t> fut_get_num_ctgs() { return reduce_one((int64_t)ctgs_depths->size(), upcxx::op_fast_add, 0); }

  void add_new_ctg(cid_t cid, int clen) {
    CtgLen ctg_len{.cid = cid, .clen = clen, .owning_rank = rank_me()};
    auto tgt = get_target_rank(cid);
    if (tgt != rank_me()) {
      add_ctg_store.update(tgt, ctg_len);
      // also store locally for return
      CtgBaseDepths cbd(cid);
      ctgs_depths->insert({cid, cbd});
      DBG("Added cid=", cid, " placeholder\n");
    } else {
      _add_ctg(ctg_len);
      DBG("Added cid=", cid, " bypass for calcs\n");
    }
  }

  void update_ctg_aln_depth(cid_t cid, int read_group_id, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
    if (aln_start >= aln_stop) return;
    CtgAlnDepth cad{cid, read_group_id, aln_start, aln_stop, aln_merge_start, aln_merge_stop};
    ctg_aln_depth_store.update(get_target_rank(cid), cad);
  }

  CtgBaseDepths *get_first_local_ctg() {
    ctgs_depths_iter = ctgs_depths->begin();
    if (ctgs_depths_iter == ctgs_depths->end()) return nullptr;
    auto ctg = &ctgs_depths_iter->second;
    ctgs_depths_iter++;
    return ctg;
  }

  CtgBaseDepths *get_next_local_ctg() {
    if (ctgs_depths_iter == ctgs_depths->end()) return nullptr;
    auto ctg = &ctgs_depths_iter->second;
    ctgs_depths_iter++;
    return ctg;
  }

  // return a vector of pairs of avg,var for total and each read_group
  future<vector<AvgVar<float>>> fut_get_depth(cid_t cid) {
    auto it = ctgs_depths->find(cid);
    if (it != ctgs_depths->end()) {
      return make_future(it->second.rg_stats);
    }
    LOG("Falling back to rpc for cid=", cid, "\n");
    auto target_rank = get_target_rank(cid);
    // DBG_VERBOSE("Sending rpc to ", target_rank, " for cid=", cid, "\n");
    return upcxx::rpc(
        target_rank,
        [](ctgs_depths_map_t &ctgs_depths, cid_t cid, int edge_base_len) -> vector<AvgVar<float>> {
          const auto it = ctgs_depths->find(cid);
          if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
          const auto &ctg_base_depths = it->second;
          auto &read_group_base_counts = ctg_base_depths.read_group_base_counts;
          for (auto &rg_base_counts_ptr : read_group_base_counts) {
            DBG("Testing ", rg_base_counts_ptr, "\n");
            assert(rg_base_counts_ptr == nullptr);
          }
          return ctg_base_depths.rg_stats;
        },
        ctgs_depths, cid, edge_base_len);
  }  // fut_get_depth

  void add_contigs(const Contigs &ctgs) {
    BarrierTimer timer(__FILEFUNC__);
    size_t bases = 0;
    for (auto &ctg : ctgs) {
      int clen = ctg.seq.length();
      if (clen < min_ctg_len) continue;
      bases += clen;
    }

    LOG("Locally processing ", ctgs.size(), " contigs with ", bases,
        " total bases. minimum_mem_required=", get_size_str(bases * 2 * sizeof(CtgBaseDepths::base_count_t)), "\n");

    LOG_MEM("Before allocating per_base ctgs_depths");
    SLOG_VERBOSE("Processing contigs, using an edge base length of ", edge_base_len, " and a min ctg len of ", min_ctg_len, "\n");
    for (const auto &ctg : ctgs) {
      int clen = ctg.seq.length();
      if (clen < min_ctg_len) continue;
      add_new_ctg(ctg.id, clen);
      upcxx::progress();
    }
    flush_contigs();
    LOG_MEM("After allocating per_base ctgs_depths");
  }

  static shared_ptr<CtgsDepths> build_ctgs_depths(Contigs &ctgs, int min_ctg_len, int num_read_groups) {
    BarrierTimer timer(__FILEFUNC__);
    int edge_base_len = (min_ctg_len >= 75 ? 75 : 0);

    size_t bases = 0;
    for (auto &ctg : ctgs) {
      int clen = ctg.seq.length();
      if (clen < min_ctg_len) continue;
      bases += clen;
    }

    auto sh_ctgs_depths = make_shared<CtgsDepths>(ctgs, edge_base_len, min_ctg_len, num_read_groups);
    auto &ctgs_depths = *sh_ctgs_depths;
    ctgs_depths.add_contigs(ctgs);
    return sh_ctgs_depths;
  }

  void compute_aln_depths_by_read_group(const Alns &alns, bool double_count_merged_region, int read_group_id = -1) {
    assert(!upcxx::in_progress());
    DBG("alns=", alns.size(), "\n");
    auto unmerged_rlen = alns.calculate_unmerged_rlen();
    int64_t num_bad_alns = 0;
    SLOG_VERBOSE("Computing aln depths for ctgs\n");
    ProgressBar progbar(alns.size(), "Processing alignments");
    for (auto &aln : alns) {
      progbar.update();
      assert(read_group_id == -1 || aln.read_group_id == read_group_id);
      // aln.check_quality();
      //  this gives abundances more in line with what we see in MetaBAT, which uses a 97% identity cut-off
      if (min_ctg_len && aln.calc_identity() < 97) {
        num_bad_alns++;
        continue;
      }

      // convert to coords for use here
      assert(aln.is_valid());
      // set to -1 if this read is not merged
      int aln_cstart_merge = -1, aln_cstop_merge = -1;
      // FIXME: need to somehow communicate to the update func the range of double counting for a merged read.
      // This is the area > length of read pair that is in the middle of the read
      if (double_count_merged_region && aln.rlen > unmerged_rlen) {
        // merged read
        int merge_offset = (aln.rlen - unmerged_rlen) / 2;
        aln_cstart_merge = (merge_offset > aln.rstart ? merge_offset - aln.rstart : 0) + aln.cstart;
        int stop_merge = aln.rlen - merge_offset;
        aln_cstop_merge = aln.cstop - (stop_merge < aln.rstop ? aln.rstop - stop_merge : 0);
        // the aln may not include the merged region
        if (aln_cstart_merge >= aln_cstop_merge) aln_cstart_merge = -1;
      }
      // as per MetaBAT analysis, ignore the 75 bases at either end because they are likely to be in error
      auto adjusted_start = ::max(aln.cstart, edge_base_len);
      auto adjusted_stop = ::min(aln.cstop, aln.clen - 1 - edge_base_len);
      // DBG_VERBOSE("Sending update for ", aln.to_string(), " st=", adjusted_start, " end=", adjusted_stop, " edge_base_len=",
      // edge_base_len,
      //    "\n");
      update_ctg_aln_depth(aln.cid, read_group_id < 0 ? 0 : aln.read_group_id, adjusted_start, adjusted_stop, aln_cstart_merge,
                           aln_cstop_merge);
      upcxx::progress();
    }
    auto fut_progbar = progbar.set_done();
    flush_ctg_alns();
    finish_read_group(read_group_id < 0 ? 0 : read_group_id);
    auto &pr = Timings::get_promise_reduce();
    auto fut_reduce = when_all(pr.reduce_one(alns.size(), op_fast_add, 0), pr.reduce_one(num_bad_alns, op_fast_add, 0));
    auto fut_report = fut_reduce.then([=](auto all_num_alns, auto all_num_bad_alns) {
      if (all_num_bad_alns) SLOG_VERBOSE("Dropped ", perc_str(all_num_bad_alns, all_num_alns), " low quality alns\n");
    });
    Timings::set_pending(when_all(fut_progbar, fut_report));
  }

  void write_aln_depths(string fname, const vector<string> &read_group_names) {
    DBG(__FILEFUNC__, "\n");
    if (read_group_names.size() != num_read_groups)
      SDIE("Wong size if read_group_names.  Expecting ", num_read_groups, " got ", read_group_names.size());
    finish_all();
    shared_ptr<upcxx_utils::dist_ofstream> ctg_ofstream;
    if (fname != "") ctg_ofstream = make_shared<upcxx_utils::dist_ofstream>(fname);

    if (!upcxx::rank_me() && ctg_ofstream) {
      *ctg_ofstream << "contigName\tcontigLen\ttotalAvgDepth";
      for (auto rg_name : read_group_names) {
        string shortname = upcxx_utils::get_basename(rg_name);
        *ctg_ofstream << "\t" << shortname << "-avg_depth\t" << shortname << "-var_depth";
      }
      *ctg_ofstream << "\n";
    }

    // FIXME: the depths need to be in the same order as the contigs in the final_assembly.fasta file. This is an inefficient
    // way of ensuring that

    future<> fut_chain = make_future();
    for (auto it = ctgs.begin(); it != ctgs.end(); it++) {
      auto &ctg = *it;
      if ((int)ctg.seq.length() < min_ctg_len) continue;
      auto fut_rg_avg_vars = fut_get_depth(ctg.id);
      auto fut_ready = when_all(fut_chain, fut_rg_avg_vars);
      fut_chain =
          fut_ready.then([&ctg, ctg_ofstream, num_read_groups = this->num_read_groups](const vector<AvgVar<float>> &rg_avg_vars) {
            assert(rg_avg_vars.size() == num_read_groups);
            double tot_depth = 0.0;
            for (int rg = 0; rg < num_read_groups; rg++) {
              tot_depth += rg_avg_vars[rg].avg;
            }
            ctg.depth = tot_depth;
            if (ctg_ofstream) {
              *ctg_ofstream << "Contig" << ctg.id << "\t" << ctg.seq.length() << "\t" << tot_depth;
              for (int rg = 0; rg < num_read_groups; rg++) {
                *ctg_ofstream << "\t" << rg_avg_vars[rg].avg << "\t" << rg_avg_vars[rg].var;
              }
              *ctg_ofstream << "\n";
            }
          });
      limit_outstanding_futures(fut_chain).wait();
      upcxx::progress();
    }
    flush_outstanding_futures();
    fut_chain.wait();
    barrier();
    if (fname != "") {
      assert(ctg_ofstream);
      DBG("Prepared contig depths for '", fname, "\n");
      ctg_ofstream->close();  // sync and print stats
    }

  }  // write_aln_depths

};   // class CtgsDepths

// wrapper method for scaffolding to just calculate depths
// so squash all read_groups into the total
void compute_aln_depths_scaffolding(Contigs &ctgs, const Alns &alns, int max_kmer_len, int min_ctg_len,
                                    bool double_count_merged_region) {
  auto sh_ctgs_depths = CtgsDepths::build_ctgs_depths(ctgs, min_ctg_len, 1);
  sh_ctgs_depths->compute_aln_depths_by_read_group(alns, double_count_merged_region, -1);
  vector<string> names;
  names.push_back("");
  sh_ctgs_depths->write_aln_depths("", names);
}

// methods for post-asm
shared_ptr<CtgsDepths> build_ctgs_depths(Contigs &ctgs, int min_ctg_len, int num_read_groups) {
  return CtgsDepths::build_ctgs_depths(ctgs, min_ctg_len, num_read_groups);
}

void compute_aln_depths_post_asm(CtgsDepths &ctg_depths, const Alns &alns, bool double_count_merged_region, int read_group_id) {
  DBG("alns=", alns.size(), "\n");
  ctg_depths.compute_aln_depths_by_read_group(alns, double_count_merged_region, read_group_id);
}

void write_aln_depths(CtgsDepths &ctg_depths, string fname, const vector<string> &read_group_names) {
  ctg_depths.write_aln_depths(fname, read_group_names);
}