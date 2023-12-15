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

#include <array>
#include <string>
#include <vector>

using std::string;
using std::vector;

#ifdef TNF_PATH_RESOLUTION
#define nTNF 136
using tnf_t = std::array<double, nTNF>;
#endif

using cid_t = int64_t;

struct Contig {
  cid_t id;
  string seq;
  double depth;
#ifdef TNF_PATH_RESOLUTION
  // histogram of teranucleotide frequencies for contig similarity check in cgraph walks
  tnf_t tnf;
#endif
  uint16_t get_uint16_t_depth() { return (depth > UINT16_MAX ? UINT16_MAX : depth); }
};

class Contigs {
  vector<Contig> contigs;

 public:
  Contigs() {}

  void clear();

  void set_capacity(int64_t sz);

  void add_contig(const Contig &contig);
  void add_contig(Contig &&contig);

  size_t size() const;

  auto begin() { return contigs.begin(); }

  auto end() { return contigs.end(); }

  auto begin() const { return contigs.begin(); }

  auto end() const { return contigs.end(); }

  void print_stats(unsigned min_ctg_len) const;

  void dump_contigs(const string &fname, unsigned min_ctg_len, const string &prefix);

  void load_contigs(const string &ctgs_fname);

  size_t get_num_ctg_kmers(int kmer_len) const;
};
