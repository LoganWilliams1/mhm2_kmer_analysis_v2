#include <sstream>
#include <string>
#include <vector>

#include "ssw.hpp"
#include "gtest/gtest.h"

#include "klign/klign.hpp"
#include "klign/aligner_cpu.hpp"
#include "klign/adept-sw/driver.hpp"

#ifdef ENABLE_GPUS
#include "gpu-utils/gpu_utils.hpp"
#endif

using std::max;
using std::min;
using std::string;
using std::vector;

using namespace StripedSmithWaterman;

void translate_adept_to_ssw(Alignment &aln, const adept_sw::AlignmentResults &aln_results, int idx) {
  aln.sw_score = aln_results.top_scores[idx];
  aln.sw_score_next_best = 0;
#ifndef REF_IS_QUERY
  aln.ref_begin = aln_results.ref_begin[idx];
  aln.ref_end = aln_results.ref_end[idx];
  aln.query_begin = aln_results.query_begin[idx];
  aln.query_end = aln_results.query_end[idx];
#else
  aln.query_begin = aln_results.ref_begin[idx];
  aln.query_end = aln_results.ref_end[idx];
  aln.ref_begin = aln_results.query_begin[idx];
  aln.ref_end = aln_results.query_end[idx];
#endif
  aln.ref_end_next_best = 0;
  aln.mismatches = 0;
  aln.cigar_string.clear();
  aln.cigar.clear();
}

#ifdef ENABLE_GPUS
void test_aligns_gpu(vector<Alignment> &alns, vector<string> query, vector<string> ref, adept_sw::GPUDriver &gpu_driver) {
  alns.reserve(query.size());
  unsigned max_q_len = 0, max_ref_len = 0;
  for (int i = 0; i < query.size(); i++) {
    if (max_q_len < query[i].size()) max_q_len = query[i].size();
    if (max_ref_len < ref[i].size()) max_ref_len = ref[i].size();
  }
  gpu_driver.run_kernel_forwards(query, ref, max_q_len, max_ref_len);
  gpu_driver.kernel_block();
  gpu_driver.run_kernel_backwards(query, ref, max_q_len, max_ref_len);
  gpu_driver.kernel_block();

  auto aln_results = gpu_driver.get_aln_results();

  for (int i = 0; i < query.size(); i++) {
    alns.push_back({});
    Alignment &alignment = alns[i];
    translate_adept_to_ssw(alignment, aln_results, i);
  }
}

void check_alns_gpu(vector<Alignment> &alns, vector<int> qstart, vector<int> qend, vector<int> rstart, vector<int> rend) {
  int i = 0;
  for (Alignment &aln : alns) {
    if (i == 15) {  // mismatch test
      EXPECT_TRUE(aln.ref_end - aln.ref_begin <= 3) << "adept.ref_begin:" << aln.ref_begin << "\tadept.ref_end:" << aln.ref_end;
      EXPECT_TRUE(aln.query_end - aln.query_begin <= 3)
          << "\tadept.query_begin:" << aln.query_begin << "\tadept.query_end:" << aln.query_end;
      EXPECT_TRUE(aln.sw_score <= 4);
      EXPECT_TRUE(aln.sw_score_next_best == 0);
    } else {
      EXPECT_EQ(aln.ref_begin, rstart[i]) << "adept.ref_begin:" << aln.ref_begin << "\t"
                                          << "correct ref_begin:" << rstart[i];
      EXPECT_EQ(aln.ref_end, rend[i]) << "\tadept.ref_end:" << aln.ref_end << "\tcorrect ref_end:" << rend[i];
      EXPECT_EQ(aln.query_begin, qstart[i]) << "\tadept.query_begin:" << aln.query_begin << "\tcorrect query_begin:" << qstart[i];
      EXPECT_EQ(aln.query_end, qend[i]) << "\tadept.query_end:" << aln.query_end << "\tcorrect query end:" << qend[i];
    }
    i++;
  }
}
#endif

string aln2string(Alignment &aln) {
  std::stringstream ss;
  ss << "score=" << aln.sw_score << " score2=" << aln.sw_score_next_best;
  ss << " rbegin=" << aln.ref_begin << " rend=" << aln.ref_end;
  ss << " qbegin=" << aln.query_begin << " qend=" << aln.query_end;
  ss << " rend2=" << aln.ref_end_next_best << " mismatches=" << aln.mismatches;
  ss << " cigarstr=" << aln.cigar_string;
  return ss.str();
}

AlnScoring aln_scoring = {.match = ALN_MATCH_SCORE,
                          .mismatch = ALN_MISMATCH_COST,
                          .gap_opening = ALN_GAP_OPENING_COST,
                          .gap_extending = ALN_GAP_EXTENDING_COST,
                          .ambiguity = ALN_AMBIGUITY_COST};
AlnScoring cigar_aln_scoring = {.match = 2, .mismatch = 4, .gap_opening = 4, .gap_extending = 2, .ambiguity = 1};

Aligner ssw_aligner;
Aligner ssw_aligner_mhm2(aln_scoring.match, aln_scoring.mismatch, aln_scoring.gap_opening, aln_scoring.gap_extending,
                         aln_scoring.ambiguity);
Aligner ssw_aligner_cigar(cigar_aln_scoring.match, cigar_aln_scoring.mismatch, cigar_aln_scoring.gap_opening,
                          cigar_aln_scoring.gap_extending, cigar_aln_scoring.ambiguity);

Filter ssw_filter(true, false, 0, 32767), ssw_filter_cigar(true, true, 0, 32767);

void test_aligns(vector<Alignment> &alns, string query, string ref) {
  alns.resize(6);
  auto reflen = ref.size();
  auto qlen = query.size();
  auto masklen = max((int)min(reflen, qlen) / 2, 15);
  ssw_aligner.Align(query.c_str(), ref.c_str(), reflen, ssw_filter, &alns[0], masklen);
  ssw_aligner.Align(query.c_str(), ref.c_str(), reflen, ssw_filter_cigar, &alns[1], masklen);

  ssw_aligner_mhm2.Align(query.c_str(), ref.c_str(), reflen, ssw_filter, &alns[2], masklen);
  ssw_aligner_mhm2.Align(query.c_str(), ref.c_str(), reflen, ssw_filter_cigar, &alns[3], masklen);

  ssw_aligner_cigar.Align(query.c_str(), ref.c_str(), reflen, ssw_filter, &alns[4], masklen);
  ssw_aligner_cigar.Align(query.c_str(), ref.c_str(), reflen, ssw_filter_cigar, &alns[5], masklen);
}

void check_alns(vector<Alignment> &alns, int qstart, int qend, int rstart, int rend, int mismatches, string query = "",
                string ref = "", string cigar = "") {
  for (Alignment &aln : alns) {
    EXPECT_EQ(aln.ref_begin, rstart) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.ref_end, rend) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.query_begin, qstart) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_EQ(aln.query_end, qend) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    if (!aln.cigar_string.empty()) {  // mismatches should be recorded...
      EXPECT_EQ(aln.mismatches, mismatches) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
      if (!cigar.empty())
        EXPECT_STREQ(aln.cigar_string.c_str(), cigar.c_str()) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    }
  }
}
void check_not_alns(vector<Alignment> &alns, string query = "", string ref = "") {
  for (Alignment &aln : alns) {
    EXPECT_TRUE(aln.ref_end - aln.ref_begin <= 2) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_TRUE(aln.query_end - aln.query_begin <= 2) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_TRUE(aln.sw_score <= 4) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
    EXPECT_TRUE(aln.sw_score_next_best == 0) << "query=" << query << " ref=" << ref << " aln=" << aln2string(aln);
  }
}

TEST(MHMTest, ssw) {
  // arrange
  // act
  // assert

  EXPECT_EQ(ssw_filter.report_cigar, false);
  EXPECT_EQ(ssw_filter_cigar.report_cigar, true);

  vector<Alignment> alns;
  string ref = "ACGT";
  string query = ref;
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 0, 3, 0, query, ref, "4=");
  ref = "AACGT";
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 1, 4, 0, query, ref, "4=");
  ref = "ACGTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 0, 3, 0, query, ref, "4=");

  ref = "ACGT";

  query = "TACGT";
  test_aligns(alns, query, ref);
  check_alns(alns, 1, 4, 0, 3, 0, query, ref, "1S4=");
  query = "TTACGT";
  test_aligns(alns, query, ref);
  check_alns(alns, 2, 5, 0, 3, 0, query, ref, "2S4=");
  query = "ACGTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 0, 3, 0, query, ref, "4=1S");
  query = "ACGTTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 0, 3, 0, 3, 0, query, ref, "4=2S");

  query = "TACGTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 1, 4, 0, 3, 0, query, ref, "1S4=1S");
  query = "TTACGTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 2, 5, 0, 3, 0, query, ref, "2S4=1S");
  query = "TACGTTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 1, 4, 0, 3, 0, query, ref, "1S4=2S");
  query = "TTACGTTT";
  test_aligns(alns, query, ref);
  check_alns(alns, 2, 5, 0, 3, 0, query, ref, "2S4=2S");

  string r = "AAAATTTTCCCCGGGG";
  string q = "AAAATTTTCCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 15, 0, 15, 0, q, r, "16=");

  // 1 subst
  q = "AAAATTTTACCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 15, 0, 15, 1, q, r, "8=1X7=");

  // 1 insert
  q = "AAAATTTTACCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 16, 0, 15, 1, q, r, "8=1I8=");

  // 1 del
  q = "AAAATTTCCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 14, 0, 15, 1, q, r, "4=1D11=");

  // no match
  q = "GCTAGCTAGCTAGCTA";
  test_aligns(alns, q, r);
  check_not_alns(alns, q, r);

  // soft clip start
  q = "GCTAAAATTTTCCCCGGGG";
  test_aligns(alns, q, r);
  check_alns(alns, 3, 18, 0, 15, 0, q, r, "3S16=");

  // soft clip end
  q = "AAAATTTTCCCCGGGGACT";
  test_aligns(alns, q, r);
  check_alns(alns, 0, 15, 0, 15, 0, q, r, "16=3S");
}

TEST(MHMTest, AdeptSW) {
  // arrange
  // act
  // assert

  double time_to_initialize;
  int device_count;
  size_t total_mem;
#ifdef ENABLE_GPUS
  gpu_utils::initialize_gpu(time_to_initialize, 0);
  //  if (device_count > 0) {
  //    EXPECT_TRUE(total_mem > 32 * 1024 * 1024);  // >32 MB
  //  }

  double init_time = 0;
  adept_sw::GPUDriver gpu_driver(0, 1, (short)aln_scoring.match, (short)-aln_scoring.mismatch, (short)-aln_scoring.gap_opening,
                                 (short)-aln_scoring.gap_extending, 300, init_time);
  std::cout << "Initialized gpu in " << time_to_initialize << "s and " << init_time << "s\n";
#endif

  vector<Alignment> alns;
  vector<string> refs, queries;
  vector<int> qstarts, qends, rstarts, rends;
  // first test
  string ref = "ACGT";
  string query = ref;
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(3);
  rstarts.push_back(0);
  rends.push_back(3);
  // second test
  ref = "AACGT";
  query = "ACGT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(3);
  rstarts.push_back(1);
  rends.push_back(4);
  // third test
  ref = "ACGTT";
  query = "ACGT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(3);
  rstarts.push_back(0);
  rends.push_back(3);
  // fourth test
  ref = "ACGT";
  query = "TACGT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(1);
  qends.push_back(4);
  rstarts.push_back(0);
  rends.push_back(3);
  // fifth test
  ref = "ACGT";
  query = "TTACGT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(2);
  qends.push_back(5);
  rstarts.push_back(0);
  rends.push_back(3);
  // sixth test
  ref = "ACGT";
  query = "ACGTT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(3);
  rstarts.push_back(0);
  rends.push_back(3);
  // seventh test
  ref = "ACGT";
  query = "ACGTTT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(3);
  rstarts.push_back(0);
  rends.push_back(3);
  // eighth test
  ref = "ACGT";
  query = "TACGTT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(1);
  qends.push_back(4);
  rstarts.push_back(0);
  rends.push_back(3);
  // ninth test
  ref = "ACGT";
  query = "TTACGTT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(2);
  qends.push_back(5);
  rstarts.push_back(0);
  rends.push_back(3);
  // tenth test
  ref = "ACGT";
  query = "TACGTTT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(1);
  qends.push_back(4);
  rstarts.push_back(0);
  rends.push_back(3);
  // eleventh test
  ref = "ACGT";
  query = "TTACGTTT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(2);
  qends.push_back(5);
  rstarts.push_back(0);
  rends.push_back(3);
  // twelvth test
  ref = "AAAATTTTCCCCGGGG";
  query = "AAAATTTTCCCCGGGG";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(15);
  rstarts.push_back(0);
  rends.push_back(15);
  // thirteenth test
  ref = "AAAATTTTCCCCGGGG";
  query = "AAAATTTTACCCGGGG";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(15);
  rstarts.push_back(0);
  rends.push_back(15);
  // 1 insert // fourteenth test
  ref = "AAAATTTTCCCCGGGG";
  query = "AAAATTTTACCCCGGGG";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(16);
  rstarts.push_back(0);
  rends.push_back(15);
  // 1 del // fifteenth test
  ref = "AAAATTTTCCCCGGGG";
  query = "AAAATTTCCCCGGGG";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(14);
  rstarts.push_back(0);
  rends.push_back(15);
  // no match // sixteenth
  ref = "AAAATTTTCCCCGGGG";
  query = "GCTAGCTAGCTAGCTA";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(14);
  rstarts.push_back(0);
  rends.push_back(15);

  // soft clip start // seventeenth test
  ref = "AAAATTTTCCCCGGGG";
  query = "GCTAAAATTTTCCCCGGGG";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(3);
  qends.push_back(18);
  rstarts.push_back(0);
  rends.push_back(15);
  // soft clip end // eighteenth test
  ref = "AAAATTTTCCCCGGGG";
  query = "AAAATTTTCCCCGGGGACT";
  queries.push_back(query);
  refs.push_back(ref);
  qstarts.push_back(0);
  qends.push_back(15);
  rstarts.push_back(0);
  rends.push_back(15);

#ifdef ENABLE_GPUS
  // run kernel
  test_aligns_gpu(alns, queries, refs, gpu_driver);
  // verify results
  check_alns_gpu(alns, qstarts, qends, rstarts, rends);
  // cuda tear down happens in driver destructor
#endif
}

TEST(MHMTest, Issue118) {
  /* From stacktrace:
  [14] #9  0x00000000107a864c in banded_sw (ref=0x200f6c71926c "\001\001\002", read=0x200f6c71a39c "\001\001\002", refLen=138,
  readLen=138, score=248, weight_gapO=4, weight_gapE=2, band_width2=1 , mat=0x4a93aba0
  "\002\374\374\374\377\374\002\374\374\377\374\374\002\374\377\374\374\374\002\377\377\377\377\377\377", n=5) at
  /ccs/home/rsegan/workspace/mhm2/src/ssw/ssw_core.cpp:943 [14] #10 0x00000000107a9bcc in ssw_align (prof=0x200f6c361ab0,
  ref=0x200f6c719260 "", refLen=138, weight_gapO=4 '\004', weight_gapE=2 '\002', flag=15 '\017', filters=0, filterd=32767, maskLen=
  75) at /ccs/home/rsegan/workspace/mhm2/src/ssw/ssw_core.cpp:1212
  [14] #11 0x00000000107a2b98 in StripedSmithWaterman::Aligner::Align (this=0x7fffe6dccb60, query=0x4b9156c0
  "AGCGGTGAATCGCCGATCGAGACGGTGCCGCCCGACGCCCGCACCACGCCTGCGACGATCGCGAGCCCCAGGCCGCTGCCCCCG
  GCATCCCTCGCCCGTCCTTCGTCGAGGCGAACGAAGCGTTCGAACACCCGCTCTCGCTCCGAGGCC", query_len=@0x200f68cddc6c: 150, ref=0x4b915760
  "AGCGGTCACTATCCGATCGAGACGGTGCCGCGCGACGCCCGCACCACGCCTGCGACGTTCGCGAGCCCCAGGCCG
  CTGCCCCCGGCATCCCTCGACCGTCCTTAGTCGAGGCTAACGAAGCGTTCGAACACCCGCTCTCGCTCCGAGGCC", ref_len=@0x200f68cddc68: 150, filter=...,
  alignment=0x200f68cddc70, maskLen=75) at /ccs/home/rsegan/workspace/mhm2 /src/ssw/ssw.cpp:457 [14] #12 0x0000000010c054f4 in
  CPUAligner::ssw_align_read (ssw_aligner=..., ssw_filter=..., alns=0x6b92bdd0, aln_scoring=..., aln=..., cseq=..., rseq=...,
  read_group_id=0) at /ccs/home/rsegan/ workspace/mhm2/src/klign/aligner_cpu.cpp:169
  */

  Alns alns;
  Aln aln[10] = {};
  string query, ref, exp_sam, exp_aln;

  query = string("AGCGGTGAATCGCCGATCGAGACGGTGCCGCCCGACGCCCGCACCACGCCTGCGACGATCGCGAGCCCCAGGCCGCTGCCCCCGGCATCCCTCGCCCGTCCTTCGTCGAGGCG"
                 "AACGAAGCGTTCGAACACCCGCTCTCGCTCCGAGGCC");
  ref = string("AGCGGTCACTATCCGATCGAGACGGTGCCGCGCGACGCCCGCACCACGCCTGCGACGTTCGCGAGCCCCAGGCCGCTGCCCCCGGCATCCCTCGACCGTCCTTAGTCGAGGCTAA"
               "CGAAGCGTTCGAACACCCGCTCTCGCTCCGAGGCC");

  aln[0] = Aln("a", 0, 0, query.size() - 1, query.size(), 0, ref.size() - 1, ref.size(), '+', 0, 0, 0, -1);
  exp_aln = "a\t1\t150\t150\tContig0\t1\t150\t150\tPlus\t248\t114";
  exp_sam = "a\t0\tContig0\t1\t7\t6=1X1=2I2=2D19=1X25=1X36=1X8=1X8=1X37=\t*\t0\t151\t*\t*\tAS:i:248\tNM:i:10\tRG:Z:0";
  // minimap cigar: 12S19=1X25=1X36=1X8=1X8=1X37=

  CPUAligner::ssw_align_read(ssw_aligner_cigar, ssw_filter_cigar, &alns, cigar_aln_scoring, aln[0], string_view(query),
                             string_view(ref), 0);
  EXPECT_EQ(alns.size(), 1) << "did not align query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(0).to_paf_string(), exp_aln) << "did align correctly query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(0).sam_string, exp_sam) << "did align correctly cigar query=" << query << " ref=" << ref;

  aln[1] = Aln("b", 0, 0, ref.size() - 1, ref.size(), 0, query.size() - 1, query.size(), '+', 0, 0, 0, -1);
  exp_aln = "b\t1\t150\t150\tContig0\t1\t150\t150\tPlus\t248\t114";
  exp_sam = "b\t0\tContig0\t1\t7\t6=1X1=2D2=2I19=1X25=1X36=1X8=1X8=1X37=\t*\t0\t151\t*\t*\tAS:i:248\tNM:i:10\tRG:Z:0";

  CPUAligner::ssw_align_read(ssw_aligner_cigar, ssw_filter_cigar, &alns, cigar_aln_scoring, aln[1], string_view(ref),
                             string_view(query), 0);
  EXPECT_EQ(alns.size(), 2) << "did not align query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(1).to_paf_string(), exp_aln) << "did align correctly query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(1).sam_string, exp_sam) << "did align correctly cigar query=" << query << " ref=" << ref;

  /*
  [20] #8  <signal handler called>
  [20] #9  0x00000000107a12a8 in (anonymous namespace)::ConvertAlignment (s_al=..., query_len=@0x200f68cddc5c: 150,
  al=0x200f68cddc60) at /ccs/home/rsegan/workspace/mhm2/src/ssw/ssw.cpp:60 [20] #10 0x00000000107a2bc0 in
  StripedSmithWaterman::Aligner::Align (this=0x7fffd8932ac0, query=0x56a4cf70
  "CAAGTCAACAAACGTAAAGTGATGGGTATGTTCTGACTCTTTGATTTTAAATTTCGAAATCTGAGCTTTTTGGGGGATGTGGCG
  TGTGAAAGCAGCTAAATCATTCCTCCCACTCAAATTTCAGGCAACGCCATTGAGTACAGGTTGTGA", query_len=@0x200f68cddc5c: 150, ref=0x555f4830
  "CCAATCAACACGCATAAAGTGATGGAGCGGTTCTGATCCCTTGGTTTAAAATTTCGAAATCTGAGCTTTTCGAGG
  GATGTGACTTGCGAAAGCAGCTAAATAATTCCTCCCGCTCAAATTTCAGGCAACGCCATTGAGTACAGGTTGTGA", ref_len=@0x200f68cddc58: 150, filter=...,
  alignment=0x200f68cddc60, maskLen=75) at /ccs/home/rsegan/workspace/mhm2 /src/ssw/ssw.cpp:463 [20] #11 0x0000000010c05394 in
  CPUAligner::ssw_align_read (ssw_aligner=..., ssw_filter=..., alns=0x76e32000, aln_scoring=..., aln=..., cseq=..., rseq=...,
  read_group_id=0) at /ccs/home/rsegan/ workspace/mhm2/src/klign/aligner_cpu.cpp:169 [20] #12 0x0000000010c05a54 in
  CPUAligner::<lambda()>::operator()(void) const (__closure=0x7569a2f8) at
  /ccs/home/rsegan/workspace/mhm2/src/klign/aligner_cpu.cpp:207
  */

  query = string("CAAGTCAACAAACGTAAAGTGATGGGTATGTTCTGACTCTTTGATTTTAAATTTCGAAATCTGAGCTTTTTGGGGGATGTGGCGTGTGAAAGCAGCTAAATCATTCCTCCCAC"
                 "TCAAATTTCAGGCAACGCCATTGAGTACAGGTTGTGA");
  ref = string("CCAATCAACACGCATAAAGTGATGGAGCGGTTCTGATCCCTTGGTTTAAAATTTCGAAATCTGAGCTTTTCGAGGGATGTGACTTGCGAAAGCAGCTAAATAATTCCTCCCGCTC"
               "AAATTTCAGGCAACGCCATTGAGTACAGGTTGTGA");

  aln[2] = Aln("c", 0, 0, query.size() - 1, query.size(), 0, ref.size() - 1, ref.size(), '+', 0, 0, 0, -1);
  exp_aln = "c\t2\t150\t150\tContig0\t1\t150\t150\tPlus\t186\t68";
  exp_sam = "c\t0\tContig0\t1\t8\t1S3=1D5=2D3=2I11=1I1=1D2X7=1I1=1D1=1X3=1X3=1X22=1X1=1X8=1X1=1X2=1X14=1X9=1X38=\t*\t0\t151\t*\t*"
            "\tAS:i:186\tNM:i:21\tRG:Z:0";

  CPUAligner::ssw_align_read(ssw_aligner_cigar, ssw_filter_cigar, &alns, cigar_aln_scoring, aln[2], string_view(query),
                             string_view(ref), 0);
  EXPECT_EQ(alns.size(), 3) << "did not align query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(2).to_paf_string(), exp_aln) << "did align correctly query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(2).sam_string, exp_sam) << "did align correctly cigar query=" << query << " ref=" << ref;

  aln[3] = Aln("d", 0, 0, ref.size() - 1, ref.size(), 0, query.size() - 1, query.size(), '+', 0, 0, 0, -1);
  exp_aln = "d\t1\t150\t150\tContig0\t2\t150\t150\tPlus\t186\t68";
  exp_sam =
      "d\t0\tContig0\t2\t8\t3=1I5=2I3=2D11=1D1=1I2X7=1I1=1D1=1X3=1X3=1X22=1X1=1X8=1X1=1X2=1X14=1X9=1X38=\t*\t0\t150\t*\t*\tAS:"
      "i:186\tNM:i:21\tRG:Z:0";

  CPUAligner::ssw_align_read(ssw_aligner_cigar, ssw_filter_cigar, &alns, cigar_aln_scoring, aln[3], string_view(ref),
                             string_view(query), 0);
  EXPECT_EQ(alns.size(), 4) << "did not align query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(3).to_paf_string(), exp_aln) << "did align correctly query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(3).sam_string, exp_sam) << "did align correctly cigar query=" << query << " ref=" << ref;

  /*
   * cseq=CTTTTAGTTTCCCATTTCTAATTTCTCATTTGTCATTTCTCATTTTAAAGTCTTTCTCATTTTTTTAAAGCTTTTCCTCTCAAAAATGCCTTGATTGTTGGGAATAATCACTAAACAAACTTAAAAATTACCGATAAAATCA
   * rseq=AGGGATGGGACTCTTTTTTTTTAGTTTTCCATTTCTCATTTGTCATTTGTCATTTCTCATTTTTTTAAAGTTTTTCCTCTCAAAAATGCCTTGATTGTTGGGAATAATCACTAAACAAACTTAAAAATTACCGTTAAAATCA
   * */
  query = string("AGGGATGGGACTCTTTTTTTTTAGTTTTCCATTTCTCATTTGTCATTTGTCATTTCTCATTTTTTTAAAGTTTTTCCTCTCAAAAATGCCTTGATTGTTGGGAATAATCACTA"
                 "AACAAACTTAAAAATTACCGTTAAAATCA");
  ref = string("CTTTTAGTTTCCCATTTCTAATTTCTCATTTGTCATTTCTCATTTTAAAGTCTTTCTCATTTTTTTAAAGCTTTTCCTCTCAAAAATGCCTTGATTGTTGGGAATAATCACTAAA"
               "CAAACTTAAAAATTACCGATAAAATCA");

  aln[4] = Aln("e", 0, 0, query.size() - 1, query.size(), 0, ref.size() - 1, ref.size(), '+', 0, 0, 0, -1);
  exp_aln = "e\t22\t142\t142\tContig0\t26\t142\t142\tPlus\t190\t72";
  exp_sam = "e\t0\tContig0\t26\t8\t21S4=1I5=1X6=1X6=4I3=1D18=1X62=1X8=\t*\t0\t118\t*\t*\tAS:i:190\tNM:i:10\tRG:Z:0";

  CPUAligner::ssw_align_read(ssw_aligner_cigar, ssw_filter_cigar, &alns, cigar_aln_scoring, aln[4], string_view(query),
                             string_view(ref), 0);
  EXPECT_EQ(alns.size(), 5) << "did not align query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(4).to_paf_string(), exp_aln) << "did align correctly query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(4).sam_string, exp_sam) << "did align correctly cigar query=" << query << " ref=" << ref;

  aln[5] = Aln("f", 0, 0, ref.size() - 1, ref.size(), 0, query.size() - 1, query.size(), '+', 0, 0, 0, -1);
  exp_aln = "f\t26\t142\t142\tContig0\t22\t142\t142\tPlus\t190\t78";
  exp_sam = "f\t0\tContig0\t22\t7\t25S4=1D5=1X6=1X6=4D3=1I18=1X62=1X8=\t*\t0\t122\t*\t*\tAS:i:190\tNM:i:10\tRG:Z:0";
  // minimap cigar                 25S4=1D5=1X6=1X6=4D3=1I18=1X62=1X8=

  CPUAligner::ssw_align_read(ssw_aligner_cigar, ssw_filter_cigar, &alns, cigar_aln_scoring, aln[5], string_view(ref),
                             string_view(query), 0);
  EXPECT_EQ(alns.size(), 6) << "did not align query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(5).to_paf_string(), exp_aln) << "did align correctly query=" << query << " ref=" << ref;
  EXPECT_EQ(alns.get_aln(5).sam_string, exp_sam) << "did align correctly cigar query=" << query << " ref=" << ref;
}
