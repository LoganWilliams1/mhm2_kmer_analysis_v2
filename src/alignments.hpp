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

#include <upcxx/upcxx.hpp>

#include "contigs.hpp"
#include "version.h"
#include "upcxx_utils/ofstream.hpp"

using upcxx_utils::dist_ofstream;

struct Aln {
  // optimal packing of data fields (does not match constructor exactly)
  string read_id;
  int64_t cid;
  int cstart, cstop, clen;
  int rstart, rstop, rlen;  // TODO can these be int16_t (for short reads only)?
  int score1, score2;       // TODO can this be uint16_t (for short reads only)?
  int mismatches;           // TODO can this be uint16_t (for short reads only)?
  int identity;
  int mapq;
  string sam_string;
  string cigar;
  int16_t read_group_id;
  char orient;  // TODO can this be bool?

  Aln();
  Aln(const string &read_id, int64_t cid, int rlen, int cstart, int clen, char orient);
  Aln(const string &read_id, int64_t cid, int rstart, int rstop, int rlen, int cstart, int cstop, int clen, char orient, int score1,
      int score2, int mismatches, int read_group_id);

  void set(int ref_begin, int ref_end, int query_begin, int query_end, int top_score, int next_best_score, int aln_mismatches,
           int aln_read_group_id);
  void set_sam_string(string cigar);
  void add_cigar_pair_info(int64_t other_cid, int other_aln_cstart, char other_orient, int read_len);
  // writes out in the format meraligner uses
  string to_paf_string() const;
  string to_blast6_string() const;
  bool is_valid() const;
  std::pair<int, int> get_unaligned_overlaps() const;
  void set_identity();
  bool check_quality() const;
  static bool cmp(const Aln &aln1, const Aln &aln2);
  friend bool operator==(const Aln &aln1, const Aln &aln2);
  friend bool operator!=(const Aln &aln1, const Aln &aln2);
};  // class Aln

class Alns {
  using alns_t = std::vector<Aln>;
  alns_t alns;
  int64_t num_dups;
  int64_t num_bad;
  int read_len;

  bool set_pair_info(const string &read_id, vector<size_t> &read1_aln_indexes, vector<size_t> &read2_aln_indexes);
  upcxx::future<> write_sam_alignments(dist_ofstream &of, int min_contig_len) const;

 public:
  enum class Format { PAF, BLAST, SAM };

  Alns();
  Alns(int read_len);

  void clear();

  bool check_dup(Aln &aln);

  void add_aln(Aln &aln);

  void append(Alns &more_alns);

  const Aln &get_aln(int64_t i) const;

  Aln &get_aln(int64_t i);

  size_t size() const;

  bool empty() const { return alns.empty(); }

  void reserve(size_t capacity);

  void reset();

  int64_t get_num_dups();

  int64_t get_num_bad();

  inline auto begin() { return alns.begin(); }
  inline auto end() { return alns.end(); };
  inline auto begin() const { return alns.begin(); }
  inline auto end() const { return alns.end(); };

  template <typename OSTREAM>
  void dump_all(OSTREAM &os, Format fmt, int min_ctg_len = 0) const;
  void dump_single_file(const string fname, Format fmt) const;
  static upcxx::future<> write_sam_header(dist_ofstream &of, const vector<string> &read_group_names, const Contigs &ctgs,
                                          int min_ctg_len);
  void dump_sam_file(const string fname, const vector<string> &read_group_names, const Contigs &ctgs, int min_contig_len = 0) const;
  void dump_sam_file(const string fname, int min_ctg_len) const;
  void dump_rank_file(const string fname, Format fmt) const;

  int calculate_unmerged_rlen() const;

  void sort_alns();

  void compute_stats(size_t &num_reads_aligned, size_t &num_bases_aligned);

  void select_pairs(size_t &num_proper_pairs);
};  // class Alns
