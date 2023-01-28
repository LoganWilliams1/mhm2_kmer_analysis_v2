#pragma once

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
#include <upcxx/upcxx.hpp>
#include <memory>

#include "alignments.hpp"
#include "contigs.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"

#include "localassm-gpu/localassm_struct.hpp"

namespace localassm_core {

struct CtgInfo {
  int64_t cid;
  char orient;
  char side;
  CtgInfo()
      : cid{}
      , orient{}
      , side{} {}
  CtgInfo(int64_t _cid, char _orient, char _side)
      : cid(_cid)
      , orient(_orient)
      , side(_side) {}
  UPCXX_SERIALIZED_FIELDS(cid, orient, side);
};

struct ReadCtgInfo {
  string read_id;
  CtgInfo ctg_info;
  UPCXX_SERIALIZED_FIELDS(read_id, ctg_info);
};

class ReadsToCtgsDHT {
 public:
  using reads_to_ctgs_map_t = HASH_TABLE<string, vector<CtgInfo>>;
  using dist_reads_to_ctgs_map_t = dist_object<reads_to_ctgs_map_t>;
  static size_t get_target_rank(const string &read_id);

 private:
  dist_reads_to_ctgs_map_t reads_to_ctgs_map;
  upcxx_utils::ThreeTierAggrStore<ReadCtgInfo> rtc_store;

 public:
  ReadsToCtgsDHT(int64_t initial_size, int64_t num_alns);

  void clear();

  void add(const string &read_id, int64_t cid, char orient, char side);

  void flush_updates();

  int64_t get_num_mappings();

  vector<CtgInfo> get_ctgs(string &read_id);

  future<vector<vector<CtgInfo>>> get_ctgs(intrank_t target_rank, vector<string> &read_ids);
};

struct CtgData {
  int64_t cid;
  string seq;
  double depth;
  UPCXX_SERIALIZED_FIELDS(cid, seq, depth);
};

struct CtgReadData {
  CtgReadData()
      : cid{}
      , side{}
      , read_seq{} {}
  CtgReadData(int64_t _cid, char _side, const ReadSeq _read_seq)
      : cid(_cid)
      , side(_side)
      , read_seq(_read_seq) {}
  int64_t cid;
  char side;
#if UPCXX_VERSION < 20210300L
  char pad[7];  // necessary in upcxx <= 2021.03 see upcxx Issue #427
  ReadSeq read_seq;

  UPCXX_SERIALIZED_FIELDS(cid, side, pad, read_seq);
#else
  ReadSeq read_seq;

  UPCXX_SERIALIZED_FIELDS(cid, side, read_seq);
#endif
};

class CtgsWithReadsDHT {
 public:
  using ctgs_map_t = HASH_TABLE<int64_t, CtgWithReads>;
  static size_t get_target_rank(int64_t cid);

 private:
  dist_object<ctgs_map_t> ctgs_map;
  ctgs_map_t::iterator ctgs_map_iter;
  upcxx_utils::ThreeTierAggrStore<CtgData> ctg_store;
  upcxx_utils::ThreeTierAggrStore<CtgReadData> ctg_read_store;

 public:
  CtgsWithReadsDHT(int64_t num_ctgs, int64_t num_ctg_bases);

  void add_ctg(Contig &ctg);

  void add_read(int64_t cid, char side, const ReadSeq read_seq);

  void add_read(const CtgReadData &ctg_read_data);

  void add_reads(vector<CtgReadData> &_ctg_read_datas);

  void flush_ctg_updates();

  void prep_ctg_read_store(int64_t num_reads, int64_t num_read_bases);

  void flush_read_updates();

  int64_t get_num_ctgs();

  int64_t get_local_num_ctgs();

  CtgWithReads *get_first_local_ctg();

  CtgWithReads *get_next_local_ctg();
};

struct WalkMetrics {
  int64_t num_walks = 0, sum_clen = 0, sum_ext = 0, max_walk_len = 0, num_reads = 0, num_sides = 0, max_num_reads = 0,
          excess_reads = 0;
  std::array<int64_t, 3> term_counts = {0};
  WalkMetrics &operator+(const WalkMetrics &wm) {
    num_walks += wm.num_walks;
    sum_clen += wm.sum_clen;
    max_walk_len += wm.max_walk_len;
    num_reads += wm.num_reads;
    num_sides += wm.num_sides;
    max_num_reads = std::max(max_num_reads, wm.max_num_reads);
    excess_reads += wm.excess_reads;
    term_counts[0] += wm.term_counts[0];
    term_counts[1] += wm.term_counts[1];
    term_counts[2] += wm.term_counts[2];
    return *this;
  }
};

void extend_ctg(CtgWithReads *ctg, WalkMetrics &wm, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                int qual_offset, int walk_len_limit, upcxx_utils::IntermittentTimer &count_mers_timer,
                upcxx_utils::IntermittentTimer &walk_mers_timer);

void add_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs);
void process_reads(unsigned kmer_len, PackedReadsList &packed_reads_list, ReadsToCtgsDHT &reads_to_ctgs,
                   CtgsWithReadsDHT &ctgs_dht);
void process_alns(const Alns &alns, ReadsToCtgsDHT &reads_to_ctgs, int insert_avg, int insert_stddev);

}  // namespace localassm_core