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

#include "alignments.hpp"

#include <fcntl.h>
#include <limits>
#include <sstream>
#include <string>
#include <upcxx/upcxx.hpp>
#include <unordered_set>

#include "contigs.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/thread_pool.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "version.h"
#include "zstr.hpp"

using namespace upcxx_utils;
using std::ostringstream;
using std::string;

using upcxx::future;

//
// class Aln
//
Aln::Aln()
    : read_id("")
    , cid(-1)
    , cstart(0)
    , cstop(0)
    , clen(0)
    , rstart(0)
    , rstop(0)
    , rlen(0)
    , score1(0)
    , score2(0)
    , mismatches(0)
    , sam_string({})
    , read_group_id(-1)
    , orient() {}

Aln::Aln(const string &read_id, int64_t cid, int rstart, int rstop, int rlen, int cstart, int cstop, int clen, char orient,
         int score1, int score2, int mismatches, int read_group_id)
    : read_id(read_id)
    , cid(cid)
    , cstart(cstart)
    , cstop(cstop)
    , clen(clen)
    , rstart(rstart)
    , rstop(rstop)
    , rlen(rlen)
    , score1(score1)
    , score2(score2)
    , mismatches(mismatches)
    , sam_string({})
    , read_group_id(read_group_id)
    , orient(orient) {
  // DBG_VERBOSE(read_id, " cid=", cid, " RG=", read_group_id, " mismatches=", mismatches, "\n");
}

void Aln::set(int ref_begin, int ref_end, int query_begin, int query_end, int top_score, int next_best_score, int aln_mismatches,
              int aln_read_group_id) {
  cstop = cstart + ref_end + 1;
  cstart += ref_begin;
  rstop = rstart + query_end + 1;
  rstart += query_begin;
  if (orient == '-') switch_orient(rstart, rstop, rlen);
  score1 = top_score;
  score2 = next_best_score;
  auto [unaligned_left, unaligned_right] = get_unaligned_overlaps();
  // FIXME: mismatches should include unaligned overlaps
  mismatches = aln_mismatches;  // + unaligned_left + unaligned_right;
  // FIXME: need to increase the start and stop values to include the unaligned overlap
  // cstart -= unaligned_left;
  // cstop += unaligned_right;
  // rstart -= unaligned_left;
  // rstop += unaligned_right;
  read_group_id = aln_read_group_id;
}

void Aln::set_sam_string(std::string_view read_seq, string cigar) {
  assert(is_valid());
  sam_string = read_id + "\t";
  string tmp;
  if (orient == '-') {
    sam_string += "16\t";
    if (read_seq != "*") {
      tmp = revcomp(string(read_seq.data(), read_seq.size()));
      read_seq = string_view(tmp.data(), tmp.size());
    }
    // reverse(read_quals.begin(), read_quals.end());
  } else {
    sam_string += "0\t";
  }
  sam_string += "Contig" + std::to_string(cid) + "\t" + std::to_string(cstart + 1) + "\t";
  uint32_t mapq;
  // for perfect match, set to same maximum as used by minimap or bwa
  if (score2 == 0) {
    mapq = 60;
  } else {
    mapq = -4.343 * log(1 - (double)abs(score1 - score2) / (double)score1);
    mapq = (uint32_t)(mapq + 4.99);
    mapq = mapq < 254 ? mapq : 254;
  }
  sam_string += std::to_string(mapq) + "\t";
  // FIXME: need to add the unaligned left and right to the cigar string as skipped bases
  // auto [unaligned_left, unaligned_right] = get_unaligned_overlaps();
  // cigar = to_string(unaligned_left) + 'S' + cigar + to_string(unaligned_right) + 'S';
  sam_string += cigar + "\t*\t0\t" + std::to_string(cstop - cstart + 1);
  // Don't output either the read sequence or quals - that causes the SAM file to bloat up hugely, and that info is already
  // available in the read files
  // aln.sam_string += read_subseq + "\t*\t";
  sam_string += "\t*\t*\t";
  sam_string +=
      "AS:i:" + std::to_string(score1) + "\tNM:i:" + std::to_string(mismatches) + "\tRG:Z:" + std::to_string(read_group_id);
  // for debugging
  // sam_string += " rstart " + to_string(rstart) + " rstop " + to_string(rstop) + " cstop " + to_string(cstop) +
  //                  " clen " + to_string(clen) + " alnlen " + to_string(rstop - rstart);
  /*
#ifdef DEBUG
  // only used if we actually include the read seq and quals in the SAM, which we don't
  int base_count = get_cigar_length(cigar);
  if (base_count != read_seq.length())
    DIE("number of bases in cigar != aln rlen, ", base_count, " != ", read_subseq.length(), "\nsam string ", aln.sam_string);
#endif
  */
}

// minimap2 PAF output format
string Aln::to_paf_string() const {
  ostringstream os;
  os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
     << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t" << (orient == '+' ? "Plus" : "Minus") << "\t"
     << score1 << "\t" << score2;
  return os.str();
}

// blast6 output format
string Aln::to_blast6_string() const {
  ostringstream os;
  // we don't track gap opens
  int gap_opens = 0;
  int aln_len = std::max(rstop - rstart, abs(cstop - cstart));
  double identity = 100.0 * (aln_len - mismatches) / aln_len;
  os << read_id << "\t"
     << "Contig" << cid << "\t" << std::fixed << std::setprecision(3) << identity << "\t" << aln_len << "\t" << mismatches << "\t"
     << gap_opens << "\t" << rstart + 1 << "\t" << rstop << "\t";
  // subject start and end reversed when orientation is minus
  if (orient == '+')
    os << cstart + 1 << "\t" << cstop;
  else
    os << cstop << "\t" << cstart + 1;
  // evalue and bitscore, which we don't have here
  os << "\t0\t0";
  return os.str();
}

bool Aln::is_valid() const {
  assert(rstart >= 0 && "start >= 0");
  assert(rstop <= rlen && "stop <= len");
  assert(cstart >= 0 && "cstart >= 0");
  assert(cstop <= clen && "cstop <= clen");

  return read_group_id >= 0 && (orient == '+' || orient == '-') && mismatches >= 0 && cid >= 0 && read_id.size() > 0;
}

std::pair<int, int> Aln::get_unaligned_overlaps() const {
  int fwd_cstart = cstart, fwd_cstop = cstop;
  if (orient == '-') switch_orient(fwd_cstart, fwd_cstop, clen);
  int unaligned_left = std::min(rstart, fwd_cstart);
  int unaligned_right = std::min(rlen - rstop, clen - fwd_cstop);
  return {unaligned_left, unaligned_right};
}

//
// class Alns
//

Alns::Alns()
    : num_dups(0)
    , num_bad(0) {}

void Alns::clear() {
  alns.clear();
  vector<Aln>().swap(alns);
}

void Alns::add_aln(Aln &aln) {
#ifdef DEBUG
  // check for duplicate alns to this read - do this backwards because only the most recent entries could be for this read
  for (auto it = alns.rbegin(); it != alns.rend(); ++it) {
    // we have no more entries for this read
    if (it->read_id != aln.read_id || it->cid != aln.cid) break;
    // now check for equality
    if (it->rstart == aln.rstart && it->rstop == aln.rstop && it->cstart == aln.cstart && it->cstop == aln.cstop) {
      num_dups++;
      return;
    }
  }
#endif
  if (!aln.is_valid()) DIE("Invalid alignment: ", aln.to_paf_string());
  assert(aln.is_valid());
  // FIXME: we'd like to require high value alignments, but can't do this because mismatch counts are not yet supported in ADEPT
  // if (aln.identity >= KLIGN_ALN_IDENTITY_CUTOFF)
  // Currently, we just filter based on excessive unaligned overlap
  // Only filter out if the SAM string is not set, i.e. we are using the alns internally rather than for post processing output
  auto [unaligned_left, unaligned_right] = aln.get_unaligned_overlaps();
  auto unaligned = unaligned_left + unaligned_right;
  int aln_len = std::max(aln.rstop - aln.rstart + unaligned, abs(aln.cstop - aln.cstart + unaligned));
  double identity = 100.0 * (aln_len - aln.mismatches - unaligned) / aln_len;
  if (!aln.sam_string.empty() ||
      (unaligned_left <= KLIGN_UNALIGNED_THRES && unaligned_right <= KLIGN_UNALIGNED_THRES))  // && identity >= 95))
    alns.push_back(aln);
  else
    num_bad++;
}

void Alns::append(Alns &more_alns) {
  alns.insert(alns.end(), more_alns.alns.begin(), more_alns.alns.end());
  num_dups += more_alns.num_dups;
  more_alns.clear();
}

const Aln &Alns::get_aln(int64_t i) const { return alns[i]; }

Aln &Alns::get_aln(int64_t i) { return alns[i]; }

size_t Alns::size() const { return alns.size(); }

void Alns::reserve(size_t capacity) { alns.reserve(capacity); }

void Alns::reset() { alns.clear(); }

int64_t Alns::get_num_dups() { return upcxx::reduce_one(num_dups, upcxx::op_fast_add, 0).wait(); }

int64_t Alns::get_num_bad() { return upcxx::reduce_one(num_bad, upcxx::op_fast_add, 0).wait(); }

void Alns::dump_rank_file(string fname) const {
  get_rank_path(fname, rank_me());
  zstr::ofstream f(fname);
  dump_all(f, false);
  f.close();
  upcxx::barrier();
}

void Alns::dump_single_file(const string fname) const {
  dist_ofstream of(fname);
  dump_all(of, false);
  of.close();
  upcxx::barrier();
}

void Alns::dump_sam_file(const string fname, const vector<string> &read_group_names, const Contigs &ctgs, int min_ctg_len) const {
  BarrierTimer timer(__FILEFUNC__);

  string out_str = "";

  dist_ofstream of(fname);
  future<> all_done = make_future();

  // First all ranks dump Sequence tags - @SQ	SN:Contig0	LN:887
  for (const auto &ctg : ctgs) {
    if (ctg.seq.length() < min_ctg_len) continue;
    assert(ctg.id >= 0);
    of << "@SQ\tSN:Contig" << std::to_string(ctg.id) << "\tLN:" << std::to_string(ctg.seq.length()) << "\n";
  }
  // all @SQ headers aggregated to the top of the file
  all_done = of.flush_collective();

  // rank 0 continues with header
  if (!upcxx::rank_me()) {
    // add ReadGroup tags - @RG ID:[0-n] DS:filename
    for (int i = 0; i < read_group_names.size(); i++) {
      string basefilename = upcxx_utils::get_basename(read_group_names[i]);
      of << "@RG\tID:" << std::to_string(i) << "\tDS:" << basefilename << "\n";
    }
    // add program information
    of << "@PG\tID:MHM2\tPN:MHM2\tVN:" << string(MHM2_VERSION) << "\n";
  }

  // next alignments.  rank0 will be first with the remaining header fields
  dump_all(of, true, min_ctg_len);

  all_done = when_all(all_done, of.close_async());
  all_done.wait();
  of.close_and_report_timings().wait();
}

int Alns::calculate_unmerged_rlen() {
  BarrierTimer timer(__FILEFUNC__);
  // get the unmerged read length - most common read length
  HASH_TABLE<int, int64_t> rlens;
  int64_t sum_rlens = 0;
  for (auto &aln : alns) {
    rlens[aln.rlen]++;
    sum_rlens += aln.rlen;
  }
  auto all_sum_rlens = upcxx::reduce_all(sum_rlens, op_fast_add).wait();
  auto all_nalns = upcxx::reduce_all(alns.size(), op_fast_add).wait();
  auto avg_rlen = all_sum_rlens / all_nalns;
  int most_common_rlen = avg_rlen;
  int64_t max_count = 0;
  for (auto &rlen : rlens) {
    if (rlen.second > max_count) {
      max_count = rlen.second;
      most_common_rlen = rlen.first;
    }
  }
  SLOG_VERBOSE("Computed unmerged read length as ", most_common_rlen, " with a count of ", max_count, " and average of ", avg_rlen,
               "\n");
  return most_common_rlen;
}

future<> Alns::sort_alns() {
  AsyncTimer timer(__FILEFUNC__);
  // execute this in a separate thread so master can continue to communicate freely
  auto fut = execute_in_thread_pool([&alns = this->alns, timer]() {
    timer.start();
    // sort the alns by name and then for the read from best score to worst - this is needed in later stages
    std::sort(alns.begin(), alns.end(), [](const Aln &elem1, const Aln &elem2) {
      if (elem1.read_id == elem2.read_id) {
        // sort by score, then contig len then last by cid to get a deterministic ordering
        if (elem1.score1 == elem2.score1) {
          if (elem1.clen == elem2.clen) return elem1.cid > elem2.cid;
          return elem1.clen > elem2.clen;
        }
        return elem1.score1 > elem2.score1;
      }
      if (elem1.read_id.length() == elem2.read_id.length()) {
        auto rlen = elem1.read_id.length();
        auto cmp = elem1.read_id.compare(0, rlen - 2, elem2.read_id, 0, rlen - 2);
        if (cmp == 0) return (elem1.read_id[rlen - 1] == '1');
        return cmp > 0;
      }
      return elem1.read_id > elem2.read_id;
    });
    timer.stop();
    return timer;
  });
  return fut.then([](AsyncTimer timer) {
    // TODO record timer and initate reports after waiting...
  });
}
