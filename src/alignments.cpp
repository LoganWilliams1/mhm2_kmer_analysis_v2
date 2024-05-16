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
#include <bitset>

#include "contigs.hpp"
#include "utils.hpp"
#include "version.h"
#include "zstr.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/thread_pool.hpp"
#include "upcxx_utils/timers.hpp"
#include "upcxx_utils/mem_profile.hpp"

using namespace upcxx_utils;
using std::bitset;
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
    , identity(0)
    , sam_string({})
    , cigar({})
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
    , identity(0)
    , sam_string({})
    , cigar({})
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
  // auto [unaligned_left, unaligned_right] = get_unaligned_overlaps();
  //  FIXME: mismatches should include unaligned overlaps
  mismatches = aln_mismatches;  // + unaligned_left + unaligned_right;
  // FIXME: need to increase the start and stop values to include the unaligned overlap
  // cstart -= unaligned_left;
  // cstop += unaligned_right;
  // rstart -= unaligned_left;
  // rstop += unaligned_right;
  read_group_id = aln_read_group_id;
  identity = calc_identity();
}

static string construct_sam_string(const string &read_id, char orient, int64_t cid, int cstart, int mapq, const string &cigar,
                                   string pair_read_id, int pair_pos, int tlen, int score1, int mismatches, int read_group_id,
                                   int identity) {
  string sam_string = "";
  // query name
  sam_string = read_id + "\t";
  // flags: https://www.samformat.info/sam-format-flag
  if (orient == '-')
    sam_string += "16\t";
  else
    sam_string += "0\t";
  // reference name
  sam_string += "scaffold_" + std::to_string(cid) + "\t";
  // first mapping position in reference
  sam_string += std::to_string(cstart + 1) + "\t";
  // mapping quality score
  sam_string += std::to_string(mapq) + "\t";
  // cigar
  sam_string += cigar + "\t";
  // reference name that mate read maps to
  sam_string += pair_read_id + "\t";
  // position of mate read in reference
  sam_string += std::to_string(pair_pos) + "\t";
  // template length (should be aln length)
  sam_string += std::to_string(tlen) + "\t";
  // read sequence - empty to save space
  sam_string += "*\t";
  // read quals - empty to save space
  sam_string += "*\t";
  // tag - alignment score
  sam_string += "AS:i:" + std::to_string(score1) + "\t";
  // tag - number of mismatches
  sam_string += "NM:i:" + std::to_string(mismatches) + "\t";
  // tag - read group (sample)
  sam_string += "RG:Z:" + std::to_string(read_group_id) + "\t";
  // tag - identity
  sam_string += "YI:f:" + std::to_string(identity);
  return sam_string;
}

void Aln::set_sam_string(string cigar) {
  assert(is_valid());
  mismatches = 0;
  // need to count the number of mismatches because these are not available from the GPU computation
  string count_str = "";
  int aln_len = 0;
  int first_trim = 0, last_trim = cigar.length();
  for (int i = 0; i < cigar.length(); i++) {
    if (isdigit(cigar[i])) {
      count_str += cigar[i];
    } else {
      if (cigar[i] == 'X' || cigar[i] == 'I' || cigar[i] == 'D') {
        mismatches += std::stoi(count_str);
      } else if (cigar[i] == 'S') {
        if (i == cigar.length() - 1)
          last_trim = i - count_str.length();
        else
          first_trim = i + 1;
      } else {
        aln_len++;
      }
      count_str = "";
    }
  }
  this->cigar = cigar.substr(first_trim, last_trim - first_trim);
  identity = calc_identity();
  // this is not being accurately computed since we don't have all the possible secondary aln scores from SSW or Adept
  uint32_t mapq = 255;
  /*
  // for perfect match, set to same maximum as appears to be generated by bbmap
  if (score2 == 0) {
    mapq = 45;
  } else {
    mapq = -4.343 * log(1 - (double)abs(score1 - score2) / (double)score1);
    mapq = (uint32_t)(mapq + 4.99);
    mapq = std::min(mapq, (uint32_t)45);
  }
  */

  sam_string = construct_sam_string(read_id, orient, cid, cstart, mapq, cigar, "*", 0, cstop - cstart + 1, score1, mismatches,
                                    read_group_id, identity);
}

void Aln::add_cigar_pair_info(int64_t other_cid, int other_aln_cstart) {
  string other_cid_str = "*";
  if (other_cid == cid)
    other_cid_str = "=";
  else if (other_cid > -1)
    other_cid_str = "scaffold_" + std::to_string(other_cid);
  sam_string = construct_sam_string(read_id, orient, cid, cstart, 255, cigar, other_cid_str, other_aln_cstart, cstop - cstart + 1,
                                    score1, mismatches, read_group_id, identity);
}

// minimap2 PAF output format
string Aln::to_paf_string() const {
  ostringstream os;
  os << read_id << "\t" << rstart + 1 << "\t" << rstop << "\t" << rlen << "\t"
     << "Contig" << cid << "\t" << cstart + 1 << "\t" << cstop << "\t" << clen << "\t" << (orient == '+' ? "Plus" : "Minus") << "\t"
     << score1 << "\t" << score2;
  //<< "0";
  return os.str();
}

// blast6 output format
string Aln::to_blast6_string() const {
  ostringstream os;
  // we don't track gap opens
  int gap_opens = 0;
  int aln_len = std::max(rstop - rstart, abs(cstop - cstart));
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

double Aln::calc_identity() const {
  int aln_len = std::max(rstop - rstart, abs(cstop - cstart));
  int num_matches = aln_len - mismatches;
  auto [unaligned_left, unaligned_right] = get_unaligned_overlaps();
  aln_len += unaligned_left + unaligned_right;
  return 100.0 * (double)num_matches / aln_len;
}

bool Aln::check_quality() const {
  int aln_len = std::max(rstop - rstart, abs(cstop - cstart));
  double perc_id = 100.0 * (aln_len - mismatches) / aln_len;
  int cigar_aln_len = 0;
  int cigar_mismatches = 0;
  string num_str = "";
  for (int i = 0; i < cigar.length(); i++) {
    if (isdigit(cigar[i])) {
      num_str += cigar[i];
      continue;
    }
    int count = stoi(num_str);
    if (cigar[i] != 'S') cigar_aln_len += count;
    num_str = "";
    switch (cigar[i]) {
      case 'X': cigar_mismatches++; break;
      case 'I':
        cigar_mismatches++;
        cigar_aln_len++;
        break;
      case 'D':
        cigar_mismatches++;
        cigar_aln_len--;
        break;
      case '=':
      case 'M':
      case 'S': break;
      default: WARN("unexpected type in cigar: '", cigar[i], "'");
    };
  }
  if (aln_len != cigar_aln_len) DIE(cigar, " ", aln_len, " ", mismatches, " [", cigar_aln_len, " ", cigar_mismatches, "]");
  return true;
}

bool Aln::cmp(const Aln &aln1, const Aln &aln2) {
  if (aln1.read_id == aln2.read_id) {
    // sort by score, then cstart, then rstart, to get the same ordering across runs (can't use cid since that can change
    // between runs
    if (aln1.score1 != aln2.score1) return aln1.score1 > aln2.score1;
    if (aln1.clen != aln2.clen) return aln1.clen > aln2.clen;
    if (aln1.cstart != aln2.cstart) return aln1.cstart < aln2.cstart;
    if (aln1.rstart != aln2.rstart) return aln1.rstart < aln2.rstart;
    if (aln1.cstop != aln2.cstop) return aln1.cstop < aln2.cstop;
    if (aln1.rstop != aln2.rstop) return aln1.rstop < aln2.rstop;
    if (aln1.orient != aln2.orient) return aln1.orient < aln2.orient;
    if (aln1.cid != aln2.cid) return aln1.cid < aln2.cid;
    // WARN("Duplicate alns: ", aln1.to_paf_string(), " ", aln2.to_paf_string());
    // if (aln1 != aln2) WARN("BUT NOT EQUAL!");
  }
  if (aln1.read_id.length() == aln2.read_id.length()) {
    auto id_len = aln1.read_id.length();
    auto id_cmp = aln1.read_id.compare(0, id_len - 2, aln2.read_id, 0, id_len - 2);
    if (id_cmp == 0) return (aln1.read_id[id_len - 1] == '1');
    return id_cmp > 0;
  }
  return aln1.read_id > aln2.read_id;
}

bool operator==(const Aln &aln1, const Aln &aln2) {
  if (aln1.read_id != aln2.read_id) return false;
  if (aln1.cid != aln2.cid) return false;
  if (aln1.cstart != aln2.cstart) return false;
  if (aln1.cstop != aln2.cstop) return false;
  if (aln1.clen != aln2.clen) return false;
  if (aln1.rstart != aln2.rstart) return false;
  if (aln1.rstop != aln2.rstop) return false;
  if (aln1.rlen != aln2.rlen) return false;
  if (aln1.score1 != aln2.score1) return false;
  if (aln1.orient != aln2.orient) return false;
  return true;
}

bool operator!=(const Aln &aln1, const Aln &aln2) { return (!(aln1 == aln2)); }

//
// class Alns
//

Alns::Alns()
    : num_dups(0)
    , num_bad(0) {}

void Alns::clear() {
  alns.clear();
  alns_t().swap(alns);
}

void Alns::add_aln(Aln &aln) {
  auto new_identity = aln.calc_identity();
  // This is not done in bbmap - poorer alns are kept. This filtering is done when computing aln depths
  // if (new_identity < 97) {
  //   num_bad++;
  //   return;
  // }
  //  check for multiple read-ctg alns. Check backwards from most recent entry, since all alns for a read are grouped
  for (auto it = alns.rbegin(); it != alns.rend();) {
    // we have no more entries for this read
    if (it->read_id != aln.read_id) break;
    if (it->cid == aln.cid) {
      num_dups++;
      auto old_identity = it->identity;
      auto new_identity = aln.calc_identity();
      // SLOG("multi aln: ", it->read_id, " ", it->cid, " ", it->score1, " ", aln.score1, " ", old_identity, " ", new_identity, "
      // ", num_dups, "\n");
      auto old_aln_len = it->rstop - it->rstart;
      it++;
      if ((new_identity > old_identity) || (new_identity == old_identity && (aln.rstop - aln.rstart > old_aln_len))) {
        // new one is better - erase the old one
        auto fit = it.base();
        if (fit != alns.end()) alns.erase(fit);
        // can only happen once because previous add_aln calls will have ensured there is only the best single aln for that cid
        break;
      } else {
        // new one is no better - don't add
        return;
      }
    } else {
      it++;
    }
  }
  if (!aln.is_valid()) DIE("Invalid alignment: ", aln.to_paf_string());
  assert(aln.is_valid());
  // FIXME: we'd like to require high value alignments, but can't do this because mismatch counts are not yet supported in ADEPT
  // if (aln.identity >= KLIGN_ALN_IDENTITY_CUTOFF)
  // Currently, we just filter based on excessive unaligned overlap
  // Only filter out if the SAM string is not set, i.e. we are using the alns internally rather than for post processing output
  auto [unaligned_left, unaligned_right] = aln.get_unaligned_overlaps();
  auto unaligned = unaligned_left + unaligned_right;
  // int aln_len = std::max(aln.rstop - aln.rstart + unaligned, abs(aln.cstop - aln.cstart + unaligned));
  if (!aln.sam_string.empty() || (unaligned_left <= KLIGN_UNALIGNED_THRES && unaligned_right <= KLIGN_UNALIGNED_THRES))
    alns.push_back(aln);
  else
    num_bad++;
}

void Alns::append(Alns &more_alns) {
  DBG("Appending ", more_alns.size(), " alignments to ", alns.size(), "\n");
  reserve(alns.size() + more_alns.alns.size());
  for (auto &a : more_alns.alns) {
    alns.emplace_back(std::move(a));
  }
  num_dups += more_alns.num_dups;
  num_bad += more_alns.num_bad;
  more_alns.clear();
}

const Aln &Alns::get_aln(int64_t i) const { return alns[i]; }

Aln &Alns::get_aln(int64_t i) { return alns[i]; }

size_t Alns::size() const { return alns.size(); }

void Alns::reserve(size_t capacity) { alns.reserve(capacity); }

void Alns::reset() { alns.clear(); }

int64_t Alns::get_num_dups() { return num_dups; }

int64_t Alns::get_num_bad() { return num_bad; }

template <typename OSTREAM>
void Alns::dump_all(OSTREAM &os, Format fmt, int min_ctg_len) const {
  // all ranks dump their valid alignments
  for (const auto &aln : alns) {
    if (aln.cid == -1) continue;
    // DBG(aln.to_paf_string(), "\n");
    if (aln.clen < min_ctg_len) continue;
    switch (fmt) {
      case Format::PAF: os << aln.to_paf_string() << "\n"; break;
      case Format::BLAST: os << aln.to_blast6_string() << "\n"; break;
      case Format::SAM: os << aln.sam_string << "\n"; break;
    }
  }
}

void Alns::dump_rank_file(string fname, Format fmt) const {
  get_rank_path(fname, rank_me());
  zstr::ofstream f(fname);
  dump_all(f, fmt);
  f.close();
  upcxx::barrier();
}

void Alns::dump_single_file(const string fname, Format fmt) const {
  dist_ofstream of(fname);
  dump_all(of, fmt);
  of.close();
  upcxx::barrier();
}

upcxx::future<> Alns::write_sam_header(dist_ofstream &of, const vector<string> &read_group_names, const Contigs &ctgs,
                                       int min_ctg_len) {
  // First all ranks dump Sequence tags - @SQ	SN:Contig0	LN:887
  for (const auto &ctg : ctgs) {
    if (ctg.seq.length() < min_ctg_len) continue;
    assert(ctg.id >= 0);
    of << "@SQ\tSN:scaffold_" << std::to_string(ctg.id) << "\tLN:" << std::to_string(ctg.seq.length()) << "\n";
  }
  // all @SQ headers aggregated to the top of the file
  auto all_done = of.flush_collective();

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
  return all_done;
}

upcxx::future<> Alns::write_sam_alignments(dist_ofstream &of, int min_ctg_len) const {
  // next alignments.  rank0 will be first with the remaining header fields
  dump_all(of, Format::SAM, min_ctg_len);
  return of.flush_collective();
}

void Alns::dump_sam_file(const string fname, const vector<string> &read_group_names, const Contigs &ctgs, int min_ctg_len) const {
  BarrierTimer timer(__FILEFUNC__);
  string out_str = "";
  dist_ofstream of(fname);
  auto fut = write_sam_header(of, read_group_names, ctgs, min_ctg_len);
  auto fut2 = write_sam_alignments(of, min_ctg_len);
  when_all(fut, fut2, of.close_async()).wait();
  of.close_and_report_timings().wait();
}

int Alns::calculate_unmerged_rlen() const {
  BarrierTimer timer(__FILEFUNC__);
  // get the unmerged read length - most common read length
  DBG("calculate_unmerged_rlen alns.size=", alns.size(), "\n");
  HASH_TABLE<int, int64_t> rlens;
  int64_t sum_rlens = 0;
  for (auto &aln : alns) {
    if (aln.cid == -1) continue;
    // DBG("aln.rlen=", aln.rlen, " ", aln.to_paf_string(), "\n");
    rlens[aln.rlen]++;
    sum_rlens += aln.rlen;
  }
  auto all_sum_rlens = upcxx::reduce_all(sum_rlens, op_fast_add).wait();
  auto all_nalns = upcxx::reduce_all(alns.size(), op_fast_add).wait();
  auto avg_rlen = all_nalns > 0 ? all_sum_rlens / all_nalns : 0;
  int most_common_rlen = avg_rlen;
  int64_t max_count = 0;
  for (auto &rlen : rlens) {
    if (rlen.second > max_count) {
      max_count = rlen.second;
      most_common_rlen = rlen.first;
    }
  }
  SLOG_VERBOSE("Computed unmerged read length as ", most_common_rlen, " with a count of ", max_count, " and average of ", avg_rlen,
               " over ", all_nalns, " alignments\n");
  return most_common_rlen;
}

void Alns::sort_alns() {
  BaseTimer timer(__FILEFUNC__);
  timer.start();
  // sort the alns by name and then for the read from best score to worst - this is needed in later stages
  std::sort(alns.begin(), alns.end(), Aln::cmp);
  // now purge any duplicates
  auto start_size = alns.size();
  auto ip = unique(alns.begin(), alns.end());
  alns.resize(distance(alns.begin(), ip));
  timer.stop();
  auto num_dups = start_size - alns.size();
  LOG("Sorted alns and removed ", num_dups, " duplicates in ", timer.get_elapsed(), " s\n");
}

void Alns::compute_stats(size_t &num_reads_mapped, size_t &num_bases_mapped, size_t &num_proper_pairs) {
  auto get_read_id = [](const string &read_id) {
    string pair_id_str = read_id.substr(read_id.length() - 2, 2);
    int pair_id = 0;
    if (pair_id_str == "/1") pair_id = 1;
    if (pair_id_str == "/2") pair_id = 2;
    if (pair_id == 0) return make_pair(read_id, 0);
    return make_pair(read_id.substr(0, read_id.length() - 2), pair_id);
  };

  BaseTimer timer(__FILEFUNC__);
  // bitset needs to be defined at compile time and the max paired read length we expect is 250. This means we use 2x the
  // space needed for a read length of 101, with 250 taking 32 and 101 taking 16
  HASH_TABLE<string, bitset<256>> mapped_reads(alns.size());
  for (auto &aln : alns) {
    if (aln.cid == -1) continue;
    auto elem = mapped_reads.find(aln.read_id);
    if (elem == mapped_reads.end()) elem = mapped_reads.insert({aln.read_id, bitset<256>()}).first;
    for (int i = aln.rstart; i < aln.rstop; i++) elem->second[i] = 1;
  }
  num_reads_mapped = mapped_reads.size();
  num_bases_mapped = 0;
  for (auto &mapped_read : mapped_reads) num_bases_mapped += mapped_read.second.count();

  num_proper_pairs = 0;
  // count proper pairs: both sides of the pair map to the same contig in the correct orientation, less than 32kbp apart
  for (size_t i = 1; i < alns.size(); i++) {
    if (alns[i].cid == -1) continue;
    auto [read_id, pair_id] = get_read_id(alns[i].read_id);
    if (pair_id != 2) continue;
    // now at second read pair - check for previous one
    auto [prev_read_id, prev_pair_id] = get_read_id(alns[i - 1].read_id);
    if (prev_read_id == read_id && prev_pair_id == 1 && alns[i - 1].cid == alns[i].cid && alns[i - 1].orient != alns[i].orient) {
      int d = (alns[i - 1].cstart < alns[i].cstart ? alns[i].cstart - alns[i - 1].cstop : alns[i - 1].cstart - alns[i].cstop);
      if (d < 32000) num_proper_pairs++;
    }
  }
  LOG_MEM("After alns.compute_stats");
}

#define PRINT_READ_ID "CP000510.1-10145 X"

void Alns::select_pairs() {
  auto get_highest_ri = [&alns = this->alns](vector<size_t> read_aln_indexes, const string &read_id) {
    if (read_aln_indexes.empty()) return (int64_t)-1;
    int64_t highest_ri = read_aln_indexes[0];
    for (int i = 1; i < read_aln_indexes.size(); i++) {
      int64_t ri = read_aln_indexes[i];
      if (alns[ri].score1 > alns[highest_ri].score1)
        highest_ri = ri;
      else if (alns[ri].score1 == alns[highest_ri].score1 && alns[ri].cid < alns[highest_ri].cid)
        highest_ri = ri;
    }
    if (read_id == PRINT_READ_ID)
      WARN("found highest aln for read ", alns[highest_ri].read_id, " cid ", alns[highest_ri].cid, " score ",
           alns[highest_ri].score1);
    return highest_ri;
  };

  auto clear_other_alns = [&alns = this->alns](int64_t highest_ri, vector<size_t> read_aln_indexes, size_t &num_cleared) {
    for (auto ri : read_aln_indexes) {
      if (ri != highest_ri) alns[ri].cid = -1;
    }
    num_cleared += read_aln_indexes.size() - 1;
    if (highest_ri == -1) num_cleared++;
  };

  vector<size_t> read1_aln_indexes;
  vector<size_t> read2_aln_indexes;
  string curr_read_id = "";
  size_t num_alns_cleared = 0;
  // chose only one pair of alns per read pair
  for (size_t i = 0; i < alns.size(); i++) {
    Aln &aln = alns[i];
    if (aln.identity < 50) {
      aln.cid = -1;
      continue;
    }
    string read_id = aln.read_id.substr(0, aln.read_id.length() - 2);
    if (i == 0) curr_read_id = read_id;
    if (read_id == PRINT_READ_ID) WARN(aln.sam_string);
    if (read_id != curr_read_id) {
      size_t highest_r1 = get_highest_ri(read1_aln_indexes, curr_read_id);
      size_t highest_r2 = get_highest_ri(read2_aln_indexes, curr_read_id);
      if (highest_r1 != -1 && highest_r2 != -1) {
        if (alns[highest_r1].cid != alns[highest_r2].cid) {
          // find highest scoring pair
          int max_score = 0;
          for (auto aln_r1 : read1_aln_indexes) {
            for (auto aln_r2 : read2_aln_indexes) {
              // int score = alns[aln_r1].score1 + alns[aln_r2].score1;
              int score = alns[aln_r1].score1 * alns[aln_r1].identity + alns[aln_r2].score1 * alns[aln_r2].identity;
              if (score > max_score) {
                highest_r1 = aln_r1;
                highest_r2 = aln_r2;
                max_score = score;
              }
            }
          }
        }
      }
      if (highest_r1 != -1) {
        if (highest_r2 != -1)
          alns[highest_r1].add_cigar_pair_info((alns[highest_r2].cid), alns[highest_r2].cstart + 1);
        else
          alns[highest_r1].add_cigar_pair_info(-1, 0);
      }
      if (highest_r2 != -1) {
        if (highest_r1 != -1)
          alns[highest_r2].add_cigar_pair_info((alns[highest_r1].cid), alns[highest_r1].cstart + 1);
        else
          alns[highest_r2].add_cigar_pair_info(-1, 0);
      }
      clear_other_alns(highest_r1, read1_aln_indexes, num_alns_cleared);
      clear_other_alns(highest_r2, read2_aln_indexes, num_alns_cleared);
      read1_aln_indexes.clear();
      read2_aln_indexes.clear();
      i--;
      curr_read_id = read_id;
      continue;
    }
    char pair = aln.read_id[aln.read_id.length() - 1];
    switch (pair) {
      case '1': read1_aln_indexes.push_back(i); break;
      case '2': read2_aln_indexes.push_back(i); break;
      default: DIE("Incorrect pair information for read ", aln.read_id, " found '", pair, "'");
    }
  }
  auto all_num_alns_cleared = reduce_one(num_alns_cleared, op_fast_add, 0).wait();
  SLOG(KLGREEN, "Number other alns dropped ", all_num_alns_cleared, KNORM, "\n");
}