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

struct CtgBaseDepths {
  using base_count_t = vector<uint16_t>;
  using read_group_base_count_t = vector<base_count_t>;
  int64_t cid;
  read_group_base_count_t read_group_base_counts;
  CtgBaseDepths(int64_t cid, int num_read_groups, int clen)
      : cid(cid)
      , read_group_base_counts(num_read_groups) {
    for (auto &base_count : read_group_base_counts) base_count.resize(clen, 0);
  }

  // UPCXX_SERIALIZED_FIELDS(cid, read_group_base_counts);
};

template <typename T>
struct AvgVar {
  T avg, var;
  AvgVar()
      : avg(0)
      , var(0) {}
};

class CtgsDepths {
 private:
  using local_ctgs_depths_map_t = HASH_TABLE<int64_t, CtgBaseDepths>;
  using ctgs_depths_map_t = upcxx::dist_object<local_ctgs_depths_map_t>;
  ctgs_depths_map_t ctgs_depths;
  int edge_base_len, num_read_groups;
  HASH_TABLE<int64_t, CtgBaseDepths>::iterator ctgs_depths_iter;
  struct CtgAlnDepth {
    int64_t cid;
    int read_group_id, aln_start, aln_stop, aln_merge_start, aln_merge_stop;
  };
  ThreeTierAggrStore<CtgAlnDepth> ctg_aln_depth_store;

  size_t get_target_rank(int64_t cid) { return std::hash<int64_t>{}(cid) % upcxx::rank_n(); }

  void _update_ctg_aln_depth(int64_t cid, int read_group_id, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
    const auto it = ctgs_depths->find(cid);
    if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
    CtgBaseDepths::read_group_base_count_t &rg_ctg = it->second.read_group_base_counts;
    assert(read_group_id < rg_ctg.size());
    auto &ctg_base_counts = rg_ctg[read_group_id];
    // DBG_VERBOSE("cid=", cid, " counting aln_start=", aln_start, " aln_stop=", aln_stop, " read_group_id=", read_group_id,
    //    " contig.size()=", ctg_base_counts.size(), "\n");
    assert(aln_start >= 0 && "Align start >= 0");
    assert(aln_stop <= ctg_base_counts.size() && "Align stop <= contig.size()");
    for (int i = aln_start; i < aln_stop; i++) {
      ctg_base_counts[i]++;
    }
    // disabled by default -- counts are "better" if this does not happen
    // instead, count *insert coverage* not *read coverage*
    // insert coverage should be closer to a Poisson distribution

    // this is the merged region in a merged read - will be counted double
    if (aln_merge_start != -1 && aln_merge_stop != -1) {
      assert(aln_merge_start >= 0 && "merge start >= 0");
      assert(aln_merge_stop <= ctg_base_counts.size() && "merge_stop <= size");
      for (int i = aln_merge_start; i < aln_merge_stop; i++) {
        ctg_base_counts[i]++;
      }
    }
  }

 public:
  CtgsDepths(int edge_base_len, int num_read_groups)
      : ctgs_depths(local_ctgs_depths_map_t{})
      , edge_base_len(edge_base_len)
      , num_read_groups(num_read_groups)
      , ctg_aln_depth_store() {
    int64_t mem_to_use = 0.1 * get_free_mem() / local_team().rank_n();
    auto max_store_bytes = std::max(mem_to_use, (int64_t)sizeof(CtgAlnDepth) * 100);
    ctg_aln_depth_store.set_size("Ctg Aln Depths", max_store_bytes);
    ctg_aln_depth_store.set_update_func([&self = *this](CtgAlnDepth cad) {
      self._update_ctg_aln_depth(cad.cid, cad.read_group_id, cad.aln_start, cad.aln_stop, cad.aln_merge_start, cad.aln_merge_stop);
    });
  }
  void flush() {
    ctg_aln_depth_store.flush_updates();
    ctg_aln_depth_store.clear();
  }

  int64_t get_num_ctgs() { return reduce_one((int64_t)ctgs_depths->size(), upcxx::op_fast_add, 0).wait(); }

  void add_new_ctg(int64_t cid, int num_read_groups, int clen) {
    upcxx::rpc(
        get_target_rank(cid),
        [](ctgs_depths_map_t &ctgs_depths, int64_t cid, int num_read_groups, int clen) {
          CtgBaseDepths newctg(cid, num_read_groups, clen);
          ctgs_depths->insert({cid, std::move(newctg)});
        },
        ctgs_depths, cid, num_read_groups, clen)
        .wait();
  }

  void update_ctg_aln_depth(int64_t cid, int read_group_id, int aln_start, int aln_stop, int aln_merge_start, int aln_merge_stop) {
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
  future<vector<AvgVar<float>>> fut_get_depth(int64_t cid) {
    auto target_rank = get_target_rank(cid);
    // DBG_VERBOSE("Sending rpc to ", target_rank, " for cid=", cid, "\n");
    return upcxx::rpc(
        target_rank,
        [](ctgs_depths_map_t &ctgs_depths, int64_t cid, int edge_base_len) -> vector<AvgVar<float>> {
          const auto it = ctgs_depths->find(cid);
          if (it == ctgs_depths->end()) DIE("could not fetch vertex ", cid, "\n");
          CtgBaseDepths::read_group_base_count_t &rg_ctg = it->second.read_group_base_counts;
          int num_read_groups = rg_ctg.size();
          vector<AvgVar<double>> stats;  // calculate in double, send in float
          stats.resize(num_read_groups + 1);
          auto &avg_depth = stats[num_read_groups].avg;
          auto &variance = stats[num_read_groups].var;
          size_t clen = rg_ctg[0].size() - 2 * edge_base_len;
          if (clen <= 0) clen = 1;
          for (int rg = 0; rg < num_read_groups; rg++) {
            auto &ctg_base_depths = rg_ctg[rg];
            for (int i = edge_base_len; i < (int)ctg_base_depths.size() - edge_base_len; i++) {
              avg_depth += ctg_base_depths[i];
              stats[rg].avg += ctg_base_depths[i];
            }
            stats[rg].avg /= clen;
          }
          avg_depth /= clen;
          for (int rg = 0; rg < num_read_groups; rg++) {
            auto &ctg_base_depths = rg_ctg[rg];
            for (int i = edge_base_len; i < (int)ctg_base_depths.size() - edge_base_len; i++) {
              variance += pow((double)ctg_base_depths[i] - avg_depth, 2.0);
              stats[rg].var += pow((double)ctg_base_depths[i] - stats[rg].avg, 2.0);
            }
            stats[rg].var /= clen;
          }
          variance /= clen;
          // if (avg_depth < 2) avg_depth = 2;
          vector<AvgVar<float>> ret_vector;
          ret_vector.resize(num_read_groups + 1);
          for (int i = 0; i < num_read_groups + 1; i++) {
            ret_vector[i].avg = stats[i].avg;
            ret_vector[i].var = stats[i].var;
          }
          return ret_vector;
        },
        ctgs_depths, cid, edge_base_len);
  }
};

void compute_aln_depths(const string &fname, Contigs &ctgs, Alns &alns, int kmer_len, int min_ctg_len, vector<string> read_groups,
                        bool double_count_merged_region) {
  BarrierTimer timer(__FILEFUNC__);
  int edge_base_len = (min_ctg_len >= 75 ? 75 : 0);
  size_t bases = 0;
  for (auto &ctg : ctgs) {
    int clen = ctg.seq.length();
    if (clen < min_ctg_len) continue;
    bases += clen;
  }
  int num_read_groups = (read_groups.empty() ? 1 : read_groups.size());
  LOG("Locally processing ", ctgs.size(), " contigs with ", bases,
      " total bases. mem=", get_size_str(bases * num_read_groups * sizeof(CtgBaseDepths::base_count_t::value_type)), "\n");
  CtgsDepths ctgs_depths(edge_base_len, read_groups.size());
  LOG_MEM("Before allocating per_base ctgs_depths");
  SLOG_VERBOSE("Processing contigs, using an edge base length of ", edge_base_len, " and a min ctg len of ", min_ctg_len, "\n");
  for (auto &ctg : ctgs) {
    int clen = ctg.seq.length();
    if (clen < min_ctg_len) continue;
    ctgs_depths.add_new_ctg(ctg.id, num_read_groups, clen);
    upcxx::progress();
  }
  LOG_MEM("After allocating per_base ctgs_depths");
  barrier();
  auto unmerged_rlen = alns.calculate_unmerged_rlen();
  int64_t num_bad_alns = 0;
  auto num_ctgs = ctgs_depths.get_num_ctgs();
  SLOG_VERBOSE("Computing aln depths for ", num_ctgs, " ctgs\n");
  ProgressBar progbar(alns.size(), "Processing alignments");
  for (auto &aln : alns) {
    progbar.update();
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
    ctgs_depths.update_ctg_aln_depth(aln.cid, aln.read_group_id, adjusted_start, adjusted_stop, aln_cstart_merge, aln_cstop_merge);
    upcxx::progress();
  }
  auto fut_progbar = progbar.set_done();
  ctgs_depths.flush();
  fut_progbar.wait();
  barrier();
  auto all_num_alns = reduce_one(alns.size(), op_fast_add, 0).wait();
  auto all_num_bad_alns = reduce_one(num_bad_alns, op_fast_add, 0).wait();
  if (all_num_bad_alns) SLOG_VERBOSE("Dropped ", perc_str(all_num_bad_alns, all_num_alns), " low quality alns\n");
  // get string to dump
  shared_ptr<upcxx_utils::dist_ofstream> ctg_ofstream;
  if (fname != "") {
    if (read_groups.empty()) SDIE("No read groups passed in for file names - this is incorrect usage");
    ctg_ofstream = make_shared<upcxx_utils::dist_ofstream>(fname);
    if (!upcxx::rank_me()) {
      *ctg_ofstream << "contigName\tcontigLen\ttotalAvgDepth";
      for (auto rg_name : read_groups) {
        string shortname = upcxx_utils::get_basename(rg_name);
        *ctg_ofstream << "\t" << shortname << "-avg_depth\t" << shortname << "-var_depth";
      }
      *ctg_ofstream << "\n";
    }
  }
  // FIXME: the depths need to be in the same order as the contigs in the final_assembly.fasta file. This is an inefficient
  // way of ensuring that

  future<> fut_chain = make_future();
  for (auto it = ctgs.begin(); it != ctgs.end(); it++) {
    auto &ctg = *it;
    if ((int)ctg.seq.length() < min_ctg_len) continue;
    auto fut_rg_avg_vars = ctgs_depths.fut_get_depth(ctg.id);
    auto fut_ready = when_all(fut_chain, fut_rg_avg_vars);
    fut_chain = fut_ready.then([&ctg_ofstream, it = it, &num_read_groups](const vector<AvgVar<float>> &rg_avg_vars) {
      auto &ctg = *it;
      if (ctg_ofstream) {
        *ctg_ofstream << "Contig" << ctg.id << "\t" << ctg.seq.length() << "\t" << rg_avg_vars[num_read_groups].avg;
        for (int rg = 0; rg < num_read_groups; rg++) {
          *ctg_ofstream << "\t" << rg_avg_vars[rg].avg << "\t" << rg_avg_vars[rg].var;
        }
        *ctg_ofstream << "\n";
      }
      ctg.depth = rg_avg_vars[num_read_groups].avg;
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
}
