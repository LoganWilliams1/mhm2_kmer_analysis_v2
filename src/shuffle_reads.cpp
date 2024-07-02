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
#include <experimental/random>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/limit_outstanding.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/three_tier_aggr_store.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/reduce_prefix.hpp"
#include "utils.hpp"
#include "hash_funcs.h"
#include "kmer.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

// each read can map to multiple contigs, so first pass is to deterime exactly one cid per read (the highest score)
// generate a map of cids to read vectors from alignments
// determine location of every read using atomics and generate a read_id to location map
// process reads and put in correct locations

using kmer_t = Kmer<32>;

static intrank_t get_target_rank(int64_t val) {
  return MurmurHash3_x64_64(reinterpret_cast<const void *>(&val), sizeof(int64_t)) % rank_n();
}

static intrank_t get_kmer_target_rank(kmer_t &kmer) { return kmer.hash() % rank_n(); }
// static intrank_t get_kmer_target_rank(kmer_t &kmer) { return kmer.minimizer_hash_fast(15) % rank_n(); }

using cid_to_reads_map_t = HASH_TABLE<int64_t, vector<int64_t>>;
using read_to_target_map_t = HASH_TABLE<int64_t, int>;
using kmer_to_cid_map_t = HASH_TABLE<uint64_t, int64_t>;

static dist_object<kmer_to_cid_map_t> compute_kmer_to_cid_map(Contigs &ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  assert(SHUFFLE_KMER_LEN < 32);
  kmer_t::set_k(SHUFFLE_KMER_LEN);
  // create kmer-cid hash table - this will only be for k < 32
  dist_object<kmer_to_cid_map_t> kmer_to_cid_map(kmer_to_cid_map_t{});
  ThreeTierAggrStore<pair<uint64_t, int64_t>> kmer_cid_store;
  kmer_cid_store.set_update_func([&kmer_to_cid_map](pair<uint64_t, int64_t> &&kmer_cid_info) {
    auto it = kmer_to_cid_map->find(kmer_cid_info.first);
    if (it == kmer_to_cid_map->end()) kmer_to_cid_map->insert(kmer_cid_info);
    // else
    //  WARN("Found duplicate kmer in cids - this shouldn't happen!");
  });
  int est_update_size = sizeof(pair<uint64_t, int64_t>);
  auto local_n = local_team().rank_n();
  int64_t mem_to_use = 0.1 * get_free_mem(true) / local_n;
  auto max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  // estimate the number of updates
  uint64_t num_updates = 0;
  for (auto &ctg : ctgs) {
    num_updates += ctg.seq.length();
  }
  kmer_cid_store.set_size("kmer cid store", max_store_bytes, local_n * 2, num_updates);
  for (auto &ctg : ctgs) {
    vector<kmer_t> kmers;
    assert(ctg.seq.length() > 0);
    kmer_t::get_kmers(SHUFFLE_KMER_LEN, ctg.seq, kmers);
    // can skip kmers to make it more efficient
    for (int i = 0; i < kmers.size(); i += 1) {
      kmer_t kmer = kmers[i];
      auto kmer_rc = kmer.revcomp();
      if (kmer < kmer_rc) kmer = kmer_rc;
      kmer_cid_store.update(get_kmer_target_rank(kmer), {kmer.get_longs()[0], ctg.id});
    }
  }
  kmer_cid_store.flush_updates();
  barrier();
  auto tot_map_size = reduce_one(kmer_to_cid_map->size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Inserted ", tot_map_size, " kmers in the map of kmers to contig ids\n");
  return kmer_to_cid_map;
}

struct KmerReqBuf {
  vector<uint64_t> kmers;
  vector<int64_t> read_ids;
  void add(uint64_t kmer, int64_t read_id) {
    kmers.push_back(kmer);
    read_ids.push_back(read_id);
  }
  void clear() {
    kmers.clear();
    read_ids.clear();
  }
};

// hack to avoid agg store update() within progress context
struct UpdateCidReadsBuffer {
  int32_t target;
  int64_t cid, rid;
};
static void update_cid_reads_buffer(ThreeTierAggrStore<pair<int64_t, int64_t>> &cid_reads_store,
                                    vector<UpdateCidReadsBuffer> &update_buffer) {
  if (update_buffer.empty()) return;
  vector<UpdateCidReadsBuffer> tmp;
  tmp.swap(update_buffer);
  for (auto &ub : tmp) {
    assert(ub.cid >= 0);
    cid_reads_store.update(ub.target, {ub.cid, ub.rid});
  }
}
static future<> update_cid_reads(intrank_t target, KmerReqBuf &kmer_req_buf, dist_object<kmer_to_cid_map_t> &kmer_to_cid_map,
                                 ThreeTierAggrStore<pair<int64_t, int64_t>> &cid_reads_store,
                                 vector<UpdateCidReadsBuffer> &update_buffer) {
  auto fut = rpc(
                 target,
                 [](dist_object<kmer_to_cid_map_t> &kmer_to_cid_map, const vector<uint64_t> &kmers) -> vector<int64_t> {
                   vector<int64_t> cids;
                   cids.resize(kmers.size());
                   for (int i = 0; i < kmers.size(); i++) {
                     const auto it = kmer_to_cid_map->find(kmers[i]);
                     cids[i] = (it == kmer_to_cid_map->end() ? -1 : it->second);
                   }
                   return cids;
                 },
                 kmer_to_cid_map, kmer_req_buf.kmers)
                 .then([read_ids = kmer_req_buf.read_ids, &cid_reads_store, &update_buffer](const vector<int64_t> &cids) {
                   if (cids.size() != read_ids.size()) WARN("buff size is wrong, ", cids.size(), " != ", read_ids.size());
                   for (int i = 0; i < cids.size(); i++) {
                     // if (cids[i] != -1) cid_reads_store.update(get_target_rank(cids[i]), {cids[i], read_ids[i]});
                     UpdateCidReadsBuffer ub{.target = get_target_rank(cids[i]), .cid = cids[i], .rid = read_ids[i]};
                     if (cids[i] != -1) update_buffer.push_back(ub);
                   }
                 });
  kmer_req_buf.clear();
  return fut;
}

static dist_object<cid_to_reads_map_t> compute_cid_to_reads_map(PackedReadsList &packed_reads_list,
                                                                dist_object<kmer_to_cid_map_t> &kmer_to_cid_map, int64_t num_ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  dist_object<cid_to_reads_map_t> cid_to_reads_map(cid_to_reads_map_t{});
  auto tot_num_ctgs = reduce_one(num_ctgs, op_fast_add, 0).wait();
  cid_to_reads_map->reserve(num_ctgs);
  ThreeTierAggrStore<pair<int64_t, int64_t>> cid_reads_store;
  cid_reads_store.set_update_func([&cid_to_reads_map](pair<int64_t, int64_t> &&cid_reads_info) {
    auto &[cid, read_id] = cid_reads_info;
    auto it = cid_to_reads_map->find(cid);
    if (it == cid_to_reads_map->end())
      cid_to_reads_map->insert({cid, {read_id}});
    else
      it->second.push_back(read_id);
  });
  int est_update_size = sizeof(pair<int64_t, int64_t>);
  auto local_n = local_team().rank_n();
  int64_t mem_to_use = 0.1 * get_free_mem(true) / local_n;
  auto max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  // estimate the number of updates
  uint64_t num_updates = 0;
  for (auto packed_reads : packed_reads_list) {
    num_updates += packed_reads->get_local_num_reads() * (packed_reads->get_max_read_len() - SHUFFLE_KMER_LEN + 1);
  }
  cid_reads_store.set_size("Read cid store", max_store_bytes, local_n * 2, num_updates);
  string read_id_str, read_seq, read_quals;
  const int MAX_REQ_BUFF = 1000;
  vector<KmerReqBuf> kmer_req_bufs;
  kmer_req_bufs.resize(rank_n());
  vector<UpdateCidReadsBuffer> cid_update_buffer;
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    int64_t prev_read_id = -1;
    for (int i = 0; i < packed_reads->get_local_num_reads(); i += 2) {
      progress();
      auto &packed_read1 = (*packed_reads)[i];
      auto &packed_read2 = (*packed_reads)[i + 1];
      // first in pair is indicated by a negative index
      auto read_id = abs(packed_read1.get_id());
      if (i > 0 && prev_read_id == read_id)
        DIE("Duplicate read id across pairs, i ", i, " read id ", read_id, " ", packed_reads->get_full_read_id(i), " ",
            packed_reads->get_full_read_id(i - 2));
      prev_read_id = read_id;
      if (-packed_read1.get_id() != packed_read2.get_id())
        DIE("Paired reads id mismatch: ", -packed_read1.get_id(), " != ", packed_read2.get_id());
      for (auto &packed_read : {packed_read1, packed_read2}) {
        vector<kmer_t> kmers;
        packed_read.unpack(read_id_str, read_seq, read_quals, packed_reads->get_qual_offset());
        if (read_seq.length() < SHUFFLE_KMER_LEN) continue;
        kmer_t::get_kmers(SHUFFLE_KMER_LEN, read_seq, kmers);
        for (int j = 0; j < kmers.size(); j += 32) {
          kmer_t kmer = kmers[j];
          auto kmer_rc = kmer.revcomp();
          if (kmer < kmer_rc) kmer = kmer_rc;
          auto target = get_kmer_target_rank(kmer);
          kmer_req_bufs[target].add(kmer.get_longs()[0], read_id);
          if (kmer_req_bufs[target].kmers.size() == MAX_REQ_BUFF) {
            auto fut = update_cid_reads(target, kmer_req_bufs[target], kmer_to_cid_map, cid_reads_store, cid_update_buffer);
            limit_outstanding_futures(fut).wait();
            update_cid_reads_buffer(cid_reads_store, cid_update_buffer);
          }
        }
      }
    }
  }
  for (auto target : upcxx_utils::foreach_rank_by_node()) {  // stagger by rank_me, round robin by node
    if (!kmer_req_bufs[target].kmers.empty()) {
      auto fut = update_cid_reads(target, kmer_req_bufs[target], kmer_to_cid_map, cid_reads_store, cid_update_buffer);
      limit_outstanding_futures(fut).wait();
      update_cid_reads_buffer(cid_reads_store, cid_update_buffer);
    }
  }
  flush_outstanding_futures();
  update_cid_reads_buffer(cid_reads_store, cid_update_buffer);
  cid_reads_store.flush_updates();
  barrier();
  auto all_map_size = reduce_one(cid_to_reads_map->size(), op_fast_add, 0).wait();
  auto all_num_ctgs = reduce_one(num_ctgs, op_fast_add, 0).wait();
  SLOG_VERBOSE("Ctg IDs that map to reads: ", perc_str(all_map_size, all_num_ctgs), "\n");
  return cid_to_reads_map;
}

static dist_object<cid_to_reads_map_t> process_alns(PackedReadsList &packed_reads_list, Alns &alns, int64_t num_ctgs) {
  BarrierTimer timer(__FILEFUNC__);
  using read_to_cid_map_t = HASH_TABLE<int64_t, pair<int64_t, int>>;
  dist_object<read_to_cid_map_t> read_to_cid_map(read_to_cid_map_t{});
  ThreeTierAggrStore<tuple<int64_t, int64_t, int>> read_cid_store;
  read_cid_store.set_update_func([&read_to_cid_map](tuple<int64_t, int64_t, int> &&read_cid_info) {
    auto &[read_id, cid, score] = read_cid_info;
    auto it = read_to_cid_map->find(read_id);
    if (it == read_to_cid_map->end()) {
      read_to_cid_map->insert({read_id, {cid, score}});
    } else if (it->second.second < score) {
      it->second.first = cid;
      it->second.second = score;
    }
  });
  int est_update_size = sizeof(tuple<int64_t, int64_t, int>);
  auto local_n = local_team().rank_n();
  int64_t mem_to_use = 0.1 * get_free_mem(true) / local_n;
  auto max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  read_cid_store.set_size("Read cid store", max_store_bytes, local_n * 2, alns.size());

  for (auto &aln : alns) {
    progress();
    // use abs to ensure both reads in a pair map to the same ctg
    int64_t packed_read_id = abs(PackedRead::to_packed_id(aln.read_id));
    read_cid_store.update(get_target_rank(packed_read_id), {packed_read_id, aln.cid, aln.score1});
  }
  read_cid_store.flush_updates();
  barrier();

  dist_object<cid_to_reads_map_t> cid_to_reads_map(cid_to_reads_map_t{});
  cid_to_reads_map->reserve(num_ctgs);
  ThreeTierAggrStore<pair<int64_t, int64_t>> cid_reads_store;
  cid_reads_store.set_update_func([&cid_to_reads_map](pair<int64_t, int64_t> &&cid_reads_info) {
    auto &[cid, read_id] = cid_reads_info;
    auto it = cid_to_reads_map->find(cid);
    if (it == cid_to_reads_map->end())
      cid_to_reads_map->insert({cid, {read_id}});
    else
      it->second.push_back(read_id);
  });
  est_update_size = sizeof(pair<int64_t, int64_t>);
  mem_to_use = 0.1 * get_free_mem(true) / local_n;
  max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  cid_reads_store.set_size("Read cid store", max_store_bytes, local_n * 2, read_to_cid_map->size());
  for (auto &[read_id, cid_elem] : *read_to_cid_map) {
    progress();
    cid_reads_store.update(get_target_rank(cid_elem.first), {cid_elem.first, read_id});
  }
  cid_reads_store.flush_updates();
  barrier();
  return cid_to_reads_map;
}

struct ReadTarget {
  int64_t read_id;
  int64_t target;
};

static dist_object<read_to_target_map_t> compute_read_locations(dist_object<cid_to_reads_map_t> &cid_to_reads_map,
                                                                int64_t tot_num_reads) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_mapped_reads = 0;
  for (auto &[cid, read_ids] : *cid_to_reads_map) num_mapped_reads += read_ids.size();
  // counted read pairs
  num_mapped_reads *= 2;
  auto &pr = Timings::get_promise_reduce();
  auto fut_read_slot = upcxx_utils::reduce_prefix(num_mapped_reads, upcxx::op_fast_add);
  auto fut_all_num_mapped_reads = pr.reduce_all(num_mapped_reads, op_fast_add);
  auto fut_max_num_mapped_reads = pr.reduce_one(num_mapped_reads, op_fast_max, 0);
  // complete pending reductions
  pr.fulfill().wait();
  auto max_num_mapped_reads = fut_max_num_mapped_reads.wait();
  auto all_num_mapped_reads = fut_all_num_mapped_reads.wait();
  auto read_slot = fut_read_slot.wait() - num_mapped_reads;  // get my starting read slot
  barrier();
  auto avg_num_mapped_reads = all_num_mapped_reads / rank_n();
  SLOG_VERBOSE("Avg mapped reads per rank ", avg_num_mapped_reads, " max ", max_num_mapped_reads, " balance ",
               (double)avg_num_mapped_reads / max_num_mapped_reads, "\n");

  dist_object<read_to_target_map_t> read_to_target_map(read_to_target_map_t{});
  upcxx_utils::ThreeTierAggrStore<ReadTarget> read_target_store;
  auto local_n = local_team().rank_n();
  int64_t mem_to_use = 0.1 * get_free_mem(true) / local_n;
  auto max_store_bytes = max(mem_to_use, (int64_t)sizeof(ReadTarget) * 100);
  // estimate the number of updates
  uint64_t num_updates = 0;
  read_target_store.set_size("Read-Targets", max_store_bytes, local_n * 2, num_mapped_reads / 2);
  read_target_store.set_update_func([&read_to_target_map](ReadTarget rt) { read_to_target_map->insert({rt.read_id, rt.target}); });
  read_to_target_map->reserve(avg_num_mapped_reads);
  barrier();
  ProgressBar progbar(num_mapped_reads, "Updating read targets");
  int64_t block = (all_num_mapped_reads + rank_n() - 1) / rank_n();
  LOG("read_slot=", read_slot, " block ", block, " num_mapped_reads=", num_mapped_reads,
      " max_num_mapped_reads=", max_num_mapped_reads, " all_num_mapped_reads=", all_num_mapped_reads, "\n");
  int64_t num_reads_found = 0;
  for (auto &[cid, read_ids] : *cid_to_reads_map) {
    progress();
    if (read_ids.empty()) continue;
    for (auto read_id : read_ids) {
      num_reads_found++;
      ReadTarget rt{read_id, read_slot / block};
      assert(rt.target < rank_n() && "Target is on a rank's block");
      read_target_store.update(get_target_rank(read_id), rt);
      // each entry is a pair
      read_slot += 2;
    }
    progbar.update(read_ids.size() * 2);
  }
  assert(read_slot == fut_read_slot.wait() && "updated all Read-Targets");
  auto fut_progbar = progbar.set_done();
  read_target_store.flush_updates();
  barrier();
  read_target_store.clear();
  Timings::wait_pending();
  fut_progbar.wait();
  auto tot_reads_found = reduce_one(read_to_target_map->size(), op_fast_add, 0).wait();
  SLOG_VERBOSE("Number of read pairs mapping to contigs is ", perc_str(tot_reads_found, tot_num_reads / 2), "\n");
  return read_to_target_map;
}

// hack to avoid agg store update() within progress context
// avoid agg store update() within progress
using PairPackedRead = pair<PackedRead, PackedRead>;
struct MoveReadsToTargetUpdatesBuffer {
  int32_t target;
  PairPackedRead read_pair;
};
static void move_reads_to_target_update(ThreeTierAggrStore<PairPackedRead> &read_seq_store,
                                        vector<MoveReadsToTargetUpdatesBuffer> &updates_buffer) {
  if (!updates_buffer.empty()) {
    vector<MoveReadsToTargetUpdatesBuffer> tmp_ub;
    tmp_ub.swap(updates_buffer);
    for (auto &ub : tmp_ub) {
      read_seq_store.update(ub.target, std::move(ub.read_pair));
    }
  }
}
static void move_reads_to_targets(PackedReadsList &packed_reads_list, dist_object<read_to_target_map_t> &read_to_target_map,
                                  int64_t all_num_reads, dist_object<PackedReadsContainer> &new_packed_reads) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_not_found = 0;
  ThreeTierAggrStore<PairPackedRead> read_seq_store;
  // FIXME read_seq_store.set_restricted_updates();
  read_seq_store.set_update_func([&new_packed_reads](PairPackedRead &&read_pair_info) {
    new_packed_reads->emplace_back(std::move(read_pair_info.first));
    new_packed_reads->emplace_back(std::move(read_pair_info.second));
  });
  int est_update_size = 600;
  auto local_n = local_team().rank_n();
  int64_t mem_to_use = 0.1 * get_free_mem(true) / local_n;
  auto max_store_bytes = max(mem_to_use, (int64_t)est_update_size * 100);
  // estimate the number of updates and bases transferred
  uint64_t num_updates = 0, num_bases = 0;
  for (auto packed_reads : packed_reads_list) {
    num_updates += packed_reads->get_local_num_reads();
    num_bases += packed_reads->get_local_num_reads() * packed_reads->get_max_read_len();
  }
  // reduce max_store_bytes by the extra_size of this non-trivial updated data
  num_updates = reduce_all(num_updates, upcxx::op_fast_max).wait();
  auto non_trivial_size = sizeof(PairPackedRead) + ((num_bases + num_updates - 1) / num_updates);
  max_store_bytes = max_store_bytes * sizeof(PairPackedRead) / non_trivial_size;
  read_seq_store.set_size("Read seq store", max_store_bytes, local_n * 2, num_updates);
  // hack to avoid agg store update() within progress context
  // avoid agg store update() within progress
  vector<MoveReadsToTargetUpdatesBuffer> updates_buffer;
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    for (int i = 0; i < packed_reads->get_local_num_reads(); i += 2) {
      progress();
      auto &packed_read1 = (*packed_reads)[i];
      auto &packed_read2 = (*packed_reads)[i + 1];
      auto read_id = abs(packed_read1.get_id());
      auto fut = rpc(
                     get_target_rank(read_id),
                     [](dist_object<read_to_target_map_t> &read_to_target_map, int64_t read_id) -> int {
                       const auto it = read_to_target_map->find(read_id);
                       if (it == read_to_target_map->end()) return -1;
                       return it->second;
                     },
                     read_to_target_map, read_id)
                     .then([&updates_buffer, &num_not_found, &read_seq_store, &packed_read1, &packed_read2](int target) {
                       if (target == -1) {
                         num_not_found++;
                         target = std::experimental::randint(0, rank_n() - 1);
                       }
                       if (target < 0 || target >= rank_n()) DIE("target out of range ", target);
                       assert(target >= 0 && target < rank_n());
                       MoveReadsToTargetUpdatesBuffer ub{.target = target, {std::move(packed_read1), std::move(packed_read2)}};
                       updates_buffer.emplace_back(std::move(ub));
                     });
      limit_outstanding_futures(fut).wait();
      move_reads_to_target_update(read_seq_store, updates_buffer);
    }
  }
  flush_outstanding_futures();
  move_reads_to_target_update(read_seq_store, updates_buffer);
  read_seq_store.flush_updates();
  barrier();
  auto all_num_not_found = reduce_one(num_not_found, op_fast_add, 0).wait();
  SLOG_VERBOSE("Didn't find contig targets for ", perc_str(all_num_not_found, all_num_reads / 2), " pairs\n");
  return;
}

void shuffle_reads(int qual_offset, PackedReadsList &packed_reads_list, Contigs &ctgs) {
  BarrierTimer timer(__FILEFUNC__);

  int64_t num_reads = 0;
  for (auto packed_reads : packed_reads_list) num_reads += packed_reads->get_local_num_reads();
  auto msm_num_reads = min_sum_max_reduce_one(num_reads).wait();
  auto all_num_reads = broadcast(msm_num_reads.sum, 0).wait();

  auto kmer_to_cid_map = compute_kmer_to_cid_map(ctgs);
  auto cid_to_reads_map = compute_cid_to_reads_map(packed_reads_list, kmer_to_cid_map, ctgs.size());
  auto read_to_target_map = compute_read_locations(cid_to_reads_map, all_num_reads);

  dist_object<PackedReadsContainer> new_packed_reads(PackedReadsContainer{});
  new_packed_reads->reserve(num_reads * 1.3);
  move_reads_to_targets(packed_reads_list, read_to_target_map, all_num_reads, new_packed_reads);

  LOG("Had ", num_reads, " shuffled to ", new_packed_reads->size(), "\n");

  // now copy the new packed reads to the old
  for (auto packed_reads : packed_reads_list) delete packed_reads;
  packed_reads_list.clear();
  packed_reads_list.push_back(new PackedReads(qual_offset, *new_packed_reads, "combined-shuffled-reads"));
  packed_reads_list[0]->set_max_read_len();
  assert(packed_reads_list.size() == 1);
  uint64_t num_reads_received = packed_reads_list[0]->get_local_num_reads();
  assert(num_reads_received == new_packed_reads->size());
  new_packed_reads->clear();
  auto &pr = Timings::get_promise_reduce();
  auto fut_msm_num_reads_received = pr.msm_reduce_one(num_reads_received);
  uint64_t num_bases_received = packed_reads_list[0]->get_local_bases();
  auto fut_msm_bases_received = pr.msm_reduce_one(num_bases_received);

  auto fut = when_all(Timings::get_pending(), fut_msm_num_reads_received, fut_msm_bases_received)
                 .then([msm_num_reads, all_num_reads](const MinSumMax<uint64_t> &shuffled_msm_num_reads,
                                                      const MinSumMax<uint64_t> &shuffled_msm_num_bases) {
                   SLOG_VERBOSE("initial  num_reads: ", msm_num_reads.to_string(), "\n");
                   SLOG_VERBOSE("shuffled num_reads: ", shuffled_msm_num_reads.to_string(), "\n");
                   SLOG_VERBOSE("shuffled num_bases: ", shuffled_msm_num_bases.to_string(), "\n");
                   auto all_num_new_reads = shuffled_msm_num_reads.sum;
                   if (all_num_new_reads != all_num_reads)
                     SWARN("Not all reads shuffled, expected ", all_num_reads, " but only shuffled ", all_num_new_reads);
                 });
  // finish pending reductions
  pr.fulfill().wait();
  Timings::wait_pending();
  fut.wait();
  barrier();
}
