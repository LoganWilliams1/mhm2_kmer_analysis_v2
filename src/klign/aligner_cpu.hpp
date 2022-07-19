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

#include <string_view>
#include <algorithm>
#include <iostream>

#include <upcxx/upcxx.hpp>

#include "upcxx_utils.hpp"
#include "kmer.hpp"
#include "ssw.hpp"
#include "utils.hpp"
#include "alignments.hpp"

using std::string;
using std::string_view;
using upcxx_utils::IntermittentTimer;

struct AlnScoring {
  int match, mismatch, gap_opening, gap_extending, ambiguity;

  string to_string() {
    std::ostringstream oss;
    oss << "match " << match << " mismatch " << mismatch << " gap open " << gap_opening << " gap extend " << gap_extending
        << " ambiguity " << ambiguity;
    return oss.str();
  }
};

// encapsulate the data for a kernel to run a block independently
struct AlignBlockData {
  vector<Aln> kernel_alns;
  vector<string> ctg_seqs;
  vector<string> read_seqs;
  shared_ptr<Alns> alns;
  int64_t max_clen;
  int64_t max_rlen;
  int read_group_id;
  AlnScoring aln_scoring;

  AlignBlockData(vector<Aln> &_kernel_alns, vector<string> &_ctg_seqs, vector<string> &_read_seqs, int64_t max_clen,
                 int64_t max_rlen, int read_group_id, AlnScoring &aln_scoring);
};

int get_cigar_length(const string &cigar);
void set_sam_string(Aln &aln, string_view read_seq, string cigar);

struct CPUAligner {
  // default aligner and filter
  StripedSmithWaterman::Aligner ssw_aligner;
  StripedSmithWaterman::Filter ssw_filter;
  AlnScoring aln_scoring;

  CPUAligner(bool compute_cigar);

  static void ssw_align_read(StripedSmithWaterman::Aligner &ssw_aligner, StripedSmithWaterman::Filter &ssw_filter, Alns *alns,
                             AlnScoring &aln_scoring, Aln &aln, const string_view &cseq, const string_view &rseq,
                             int read_group_id);

  void ssw_align_read(Alns *alns, Aln &aln, const string &cseq, const string &rseq, int read_group_id);

  upcxx::future<> ssw_align_block(shared_ptr<AlignBlockData> aln_block_data, Alns *alns, IntermittentTimer &aln_kernel_timer);
};
