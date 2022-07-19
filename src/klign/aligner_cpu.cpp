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

#include "aligner_cpu.hpp"

#ifdef __PPC64__  // FIXME remove after solving Issues #60 #35 #49
#define NO_KLIGN_CPU_WORK_STEAL
#endif

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

AlignBlockData::AlignBlockData(vector<Aln> &_kernel_alns, vector<string> &_ctg_seqs, vector<string> &_read_seqs, int64_t max_clen,
                               int64_t max_rlen, int read_group_id, AlnScoring &aln_scoring)
    : max_clen(max_clen)
    , max_rlen(max_rlen)
    , read_group_id(read_group_id)
    , aln_scoring(aln_scoring) {
  // copy/swap/reserve necessary data and configs
  size_t batch_sz = std::max(kernel_alns.size(), _kernel_alns.size());
  kernel_alns.swap(_kernel_alns);
  _kernel_alns.reserve(batch_sz);
  ctg_seqs.swap(_ctg_seqs);
  _ctg_seqs.reserve(batch_sz);
  read_seqs.swap(_read_seqs);
  _read_seqs.reserve(batch_sz);
  alns = make_shared<Alns>();
  alns->reserve(kernel_alns.size());
}

int get_cigar_length(const string &cigar) {
  // check that cigar string length is the same as the sequence, but only if the sequence is included
  int base_count = 0;
  string num = "";
  for (char c : cigar) {
    switch (c) {
      case 'M':
      case 'S':
      case '=':
      case 'X':
      case 'I':
        base_count += stoi(num);
        num = "";
        break;
      case 'D':
        // base_count -= stoi(num);
        num = "";
        break;
      default:
        if (!isdigit(c)) DIE("Invalid char detected in cigar: '", c, "'");
        num += c;
        break;
    }
  }
  return base_count;
}

CPUAligner::CPUAligner(bool compute_cigar)
    : ssw_aligner() {
  // default for normal alignments in the pipeline, but for final alignments, uses minimap2 defaults
  if (!compute_cigar)
    aln_scoring = {.match = ALN_MATCH_SCORE,
                   .mismatch = ALN_MISMATCH_COST,
                   .gap_opening = ALN_GAP_OPENING_COST,
                   .gap_extending = ALN_GAP_EXTENDING_COST,
                   .ambiguity = ALN_AMBIGUITY_COST};
  else
    // these are BLASTN defaults (https://www.arabidopsis.org/Blast/BLASToptions.jsp) - except match is actually 2
    // there is no BLAST value for ambiguity
    aln_scoring = {.match = 2, .mismatch = 3, .gap_opening = 5, .gap_extending = 2, .ambiguity = 1};
  SLOG_VERBOSE("Alignment scoring parameters: ", aln_scoring.to_string(), "\n");

  // aligner construction: SSW internal defaults are 2 2 3 1
  ssw_aligner.Clear();
  ssw_aligner.ReBuild(aln_scoring.match, aln_scoring.mismatch, aln_scoring.gap_opening, aln_scoring.gap_extending,
                      aln_scoring.ambiguity);
  ssw_filter.report_cigar = compute_cigar;
}

void CPUAligner::ssw_align_read(StripedSmithWaterman::Aligner &ssw_aligner, StripedSmithWaterman::Filter &ssw_filter, Alns *alns,
                                AlnScoring &aln_scoring, Aln &aln, const string_view &cseq, const string_view &rseq,
                                int read_group_id) {
  // debugging with these alignments
  // missing from mhm:
  // CP000510.1-101195/2	Contig6043	96.667	150	5	0	1	150	109	258	1.44e-66	249
  // probably too close to this one, which is found by mhm:
  // CP000510.1-101195/2	Contig6043	92.708	96	7	0	6	101	1	96	1.54e-34	142

  assert(aln.clen >= cseq.length() && "contig seq is contained within the greater contig");
  assert(aln.rlen >= rseq.length() && "read seq is contained with the greater read");

  StripedSmithWaterman::Alignment ssw_aln;

  // align query, ref, reflen
  ssw_aligner.Align(rseq.data(), rseq.length(), cseq.data(), cseq.length(), ssw_filter, &ssw_aln,
                    max((int)(rseq.length() / 2), 15));
  aln.set(ssw_aln.ref_begin, ssw_aln.ref_end, ssw_aln.query_begin, ssw_aln.query_end, ssw_aln.sw_score, ssw_aln.sw_score_next_best,
          ssw_aln.mismatches, read_group_id);
  if (ssw_filter.report_cigar) aln.set_sam_string(rseq, ssw_aln.cigar_string);
  alns->add_aln(aln);
  /*
  if (ssw_filter.report_cigar && (aln.read_id == "CP000510.1-101195/2")) {
    cout << KLGREEN << "aln: " << aln.to_string() << KNORM << endl;
    cout << KLGREEN << "rseq " << rseq << KNORM << endl;
    cout << KLGREEN << "cseq " << cseq << KNORM << endl;
    cout << KLGREEN << "ssw aln: ref begin " << ssw_aln.ref_begin << " ref_end " << ssw_aln.ref_end << " query begin "
         << ssw_aln.query_begin << " query end " << ssw_aln.query_end << KNORM << endl;
  }*/
}

void CPUAligner::ssw_align_read(Alns *alns, Aln &aln, const string &cseq, const string &rseq, int read_group_id) {
  ssw_align_read(ssw_aligner, ssw_filter, alns, aln_scoring, aln, cseq, rseq, read_group_id);
}

upcxx::future<> CPUAligner::ssw_align_block(shared_ptr<AlignBlockData> aln_block_data, Alns *alns, IntermittentTimer &aln_kernel_timer) {
  AsyncTimer t("ssw_align_block (thread)");
  future<> fut = upcxx_utils::execute_in_thread_pool(
      [&ssw_aligner = this->ssw_aligner, &ssw_filter = this->ssw_filter, &aln_scoring = this->aln_scoring, aln_block_data, t]() {
        t.start();
        assert(!aln_block_data->kernel_alns.empty());
        DBG_VERBOSE("Starting _ssw_align_block of ", aln_block_data->kernel_alns.size(), "\n");
        auto alns_ptr = aln_block_data->alns.get();
        for (int i = 0; i < aln_block_data->kernel_alns.size(); i++) {
          Aln &aln = aln_block_data->kernel_alns[i];
          string &cseq = aln_block_data->ctg_seqs[i];
          string &rseq = aln_block_data->read_seqs[i];
          DBG_VERBOSE("aligning ", i, " of ", aln_block_data->kernel_alns.size(), " ", aln.read_id, "\n");
          ssw_align_read(ssw_aligner, ssw_filter, alns_ptr, aln_scoring, aln, cseq, rseq, aln_block_data->read_group_id);
        }
        t.stop();
      });
  fut = fut.then([alns = alns, aln_block_data, t, &aln_kernel_timer]() {
    SLOG_VERBOSE("Finished CPU SSW aligning block of ", aln_block_data->kernel_alns.size(), " in ", t.get_elapsed(), " s (",
                 (t.get_elapsed() > 0 ? aln_block_data->kernel_alns.size() / t.get_elapsed() : 0.0), " aln/s)\n");
    DBG_VERBOSE("appending and returning ", aln_block_data->alns->size(), "\n");
    alns->append(*(aln_block_data->alns));
    aln_kernel_timer += t;
  });

  return fut;
}
