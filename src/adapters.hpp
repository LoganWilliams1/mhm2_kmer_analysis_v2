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

#include <string>
#include <vector>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/timers.hpp"
#include "kmer.hpp"
#include "ssw.hpp"
#include "utils.hpp"

using std::string;

#define MAX_ADAPTER_K 32
// kmer mapping to {index for adapter seq in adapter_seqs vector, offset within adapter seq}
using adapter_hash_table_t = HASH_TABLE<Kmer<MAX_ADAPTER_K>, vector<pair<int, int>>>;
using adapter_sequences_t = std::vector<string>;

class Adapters {
  adapter_sequences_t adapter_seqs;
  adapter_hash_table_t kmer_adapter_map;
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  int64_t bases_trimmed = 0;
  int64_t reads_removed = 0;
  int adapter_k;
  bool use_blastn_scores;
  upcxx_utils::BaseTimer trim_timer;
  upcxx_utils::BaseTimer trim_timer_ssw;

  void load_adapter_seqs(const string &fname);
  bool trim_seq(const string &id, string &seq, bool is_read_1);

 public:
  Adapters(int adapter_k, const string &fname, bool use_blastn_scores);
  void done(size_t all_bases_read, size_t all_num_pairs);
  bool trim(const string &id1, string &seq1, string &quals1, const string &id2, string &seq2, string &quals2);
};