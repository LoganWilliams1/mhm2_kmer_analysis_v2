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

#include <math.h>
#include <stdarg.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#ifdef __x86_64__
#include <emmintrin.h>
#include <immintrin.h>
#include <x86intrin.h>
#endif
#include <functional>
#include <string_view>
#include <upcxx/upcxx.hpp>
#include <utility>

using namespace std;
using namespace upcxx;

#include "fastq.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "utils.hpp"
#include "zstr.hpp"
#include "kmer.hpp"
#include "ssw.hpp"

using namespace upcxx_utils;

static const double Q2Perror[] = {
    1.0,       0.7943,    0.6309,    0.5012,    0.3981,    0.3162,    0.2512,    0.1995,    0.1585,    0.1259,     0.1,
    0.07943,   0.06310,   0.05012,   0.03981,   0.03162,   0.02512,   0.01995,   0.01585,   0.01259,   0.01,       0.007943,
    0.006310,  0.005012,  0.003981,  0.003162,  0.002512,  0.001995,  0.001585,  0.001259,  0.001,     0.0007943,  0.0006310,
    0.0005012, 0.0003981, 0.0003162, 0.0002512, 0.0001995, 0.0001585, 0.0001259, 0.0001,    7.943e-05, 6.310e-05,  5.012e-05,
    3.981e-05, 3.162e-05, 2.512e-05, 1.995e-05, 1.585e-05, 1.259e-05, 1e-05,     7.943e-06, 6.310e-06, 5.012e-06,  3.981e-06,
    3.162e-06, 2.512e-06, 1.995e-06, 1.585e-06, 1.259e-06, 1e-06,     7.943e-07, 6.310e-07, 5.012e-07, 3.981e-07,  3.1622e-07,
    2.512e-07, 1.995e-07, 1.585e-07, 1.259e-07, 1e-07,     7.943e-08, 6.310e-08, 5.012e-08, 3.981e-08, 3.1622e-08, 2.512e-08,
    1.995e-08, 1.585e-08, 1.259e-08, 1e-08};

static pair<uint64_t, int> estimate_num_reads(vector<string> &reads_fname_list, int ranks_per_file = 7) {
  // estimate reads in this rank's section of all the files
  future<> progress_fut = make_future();

  BarrierTimer timer(__FILEFUNC__);
  size_t total_size = FastqReaders::open_all(reads_fname_list);

  // Issue #61 - reduce the # of reading ranks to fix excessively long estimates on poor filesystems
  // only a handful of ranks per file are needed to perform the estimate (i.e. 7)
  SLOG_VERBOSE("Estimating with about ", ranks_per_file, " ranks for each of the ", reads_fname_list.size(), " file(s)\n");
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t total_records_processed = 0;
  int64_t total_file_bytes_read = 0;
  string id, seq, quals;
  int max_read_len = 0;
  int read_file_idx = 0;
  int64_t my_total_file_size = 0;
  int64_t my_estimated_total_records = 0;
  int64_t estimated_total_records = 0;
  std::vector<int> file_bytes_per_record(reads_fname_list.size(), 0);

  ProgressBar progbar(total_size / rank_n(), "Scanning reads file to estimate number of reads");
  for (auto const &reads_fname : reads_fname_list) {
    discharge();
    // assume this rank will not read this file
    size_t my_file_size = 0;
    bool will_read = false;
    FastqReader &fqr = FastqReaders::get(reads_fname);
    if (FastqReaders::is_open(reads_fname) && fqr.my_file_size() > 0) {
      // test if my block crosses one of the evenly spaced boundaries on this file
      auto division_blocks = fqr.get_file_size().wait() / (ranks_per_file + 1);
      auto my_start = fqr.tellg();
      my_file_size = fqr.my_file_size();
      auto my_end = my_start + my_file_size;
      for (int i = 0; i < ranks_per_file; i++) {
        auto block_boundary = division_blocks * (i + 1);
        if (my_start <= block_boundary && my_end >= block_boundary) {
          DBG("Will read block ", i, " for the estimate of file#", read_file_idx, ":", fqr.get_fname(), " my_start=", my_start,
              " <= boundary=", block_boundary, " <= my_end=", my_end, "\n");
          will_read = true;
        }
        if (will_read) break;
      }
    }

    size_t tot_bytes_read = 0, file_bytes_read = 0;
    int64_t records_processed = 0, total_records = 0;
    my_total_file_size += my_file_size;
    if (will_read) {
      auto pos = fqr.tellg();
      while (true) {
        size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
        if (!bytes_read) break;
        num_lines += 4;
        num_reads++;
        tot_bytes_read += bytes_read;
        progbar.update(tot_bytes_read);
        records_processed++;
        // do not read the entire data set for just an estimate
        if (records_processed > 50000) break;
      }
      file_bytes_read = fqr.tellg() - pos;
      total_file_bytes_read += file_bytes_read;
      DBG("processed ", records_processed, " records in ", fqr.get_fname(), " over ", file_bytes_read, " file bytes (",
          tot_bytes_read, " stream bytes)\n");
      if (tot_bytes_read < file_bytes_read) file_bytes_read = tot_bytes_read;  // use the minimum of the two measures
      fqr.reset();  // rewind this file for the next reading as this was only an estimation
    }
    auto file_size = fqr.get_file_size().wait();
    auto fname = fqr.get_fname();
    total_records_processed += records_processed;
    int bytes_per = (int)(records_processed > 0 ? (file_bytes_read / records_processed) : std::numeric_limits<int>::max());
    auto fut_reduce =
        reduce_all(bytes_per, op_fast_min)
            .then([file_size, my_file_size, fname, &my_estimated_total_records, &estimated_total_records](int min_bytes_per) {
              DBG("Found min_bytes_per=", min_bytes_per, " for file ", fname, " my_size=", my_file_size, "\n");
              assert(min_bytes_per < std::numeric_limits<int>::max());
              assert(min_bytes_per > 0);
              my_estimated_total_records += my_file_size / min_bytes_per;
              estimated_total_records += file_size / min_bytes_per;
            });
    progress_fut = when_all(progress_fut, fut_reduce);
    max_read_len = max(fqr.get_max_read_len(), max_read_len);
    read_file_idx++;
  }
  progress_fut = when_all(progress_fut, progbar.set_done());
  auto fut_max_read_len = reduce_all(max_read_len, op_fast_max);
  DBG("This rank processed ", num_lines, " lines (", num_reads, " reads) with max_read_len=", max_read_len,
      " my_est=", my_estimated_total_records, "\n");
  progress_fut.wait();
  max_read_len = fut_max_read_len.wait();

  timer.initate_exit_barrier();  // barrier ensures all have completed for next reduction
  SLOG_VERBOSE("Found maximum read length of ", max_read_len, " and (rank 0's) max estimated total ", estimated_total_records,
               " per rank\n");
  return {my_estimated_total_records, max_read_len};
}

// returns the number of mismatches if it is <= max or a number greater than max (but no the actual count)
int16_t fast_count_mismatches(const char *a, const char *b, int len, int16_t max) {
  assert(len < 32768);
  int16_t mismatches = 0;
  int16_t jumpSize, jumpLen;

#if defined(__APPLE__) && defined(__MACH__)
#else
#if defined(__x86_64__)
  // 128-bit SIMD
  if (len >= 16) {
    jumpSize = sizeof(__m128i);
    jumpLen = len / jumpSize;
    for (int16_t i = 0; i < jumpLen; i++) {
      __m128i aa = _mm_loadu_si128((const __m128i *)a);     // load 16 bytes from a
      __m128i bb = _mm_loadu_si128((const __m128i *)b);     // load 16 bytes from b
      __m128i matched = _mm_cmpeq_epi8(aa, bb);             // bytes that are equal are now 0xFF, not equal are 0x00
      uint32_t myMaskMatched = _mm_movemask_epi8(matched);  // mask of most significant bit for each byte
      // count mismatches
      mismatches += _popcnt32((~myMaskMatched) & 0xffff);  // over 16 bits
      if (mismatches > max) break;
      a += jumpSize;
      b += jumpSize;
    }
    len -= jumpLen * jumpSize;
  }
#endif
#endif
  // CPU version and fall through 8 bytes at a time
  if (mismatches <= max) {
    assert(len >= 0);
    jumpSize = sizeof(int64_t);
    jumpLen = len / jumpSize;
    for (int16_t i = 0; i < jumpLen; i++) {
      int64_t *aa = (int64_t *)a, *bb = (int64_t *)b;
      if (*aa != *bb) {  // likely
        for (int j = 0; j < jumpSize; j++) {
          if (a[j] != b[j]) mismatches++;
        }
        if (mismatches > max) break;
      }  // else it matched
      a += jumpSize;
      b += jumpSize;
    }
    len -= jumpLen * jumpSize;
  }
  // do the remaining bytes, if needed
  if (mismatches <= max) {
    assert(len >= 0);
    for (int j = 0; j < len; j++) {
      mismatches += ((a[j] == b[j]) ? 0 : 1);
    }
  }
  return mismatches;
}

#define MAX_ADAPTER_K 32

// using string_view so we don't store the string again for every kmer
using adapter_hash_table_t = HASH_TABLE<Kmer<MAX_ADAPTER_K>, vector<string_view>>;
using adapter_sequences_t = vector<string>;

static void load_adapter_seqs(const string &fname, adapter_sequences_t &adapter_seqs, adapter_hash_table_t &adapters,
                              int adapter_k) {
  adapters.clear();

  // avoid every rank reading this small file
  adapter_sequences_t new_seqs;
  vector<size_t> sizes;
  if (!rank_me()) {
    ifstream f(fname);
    if (!f.is_open()) DIE("Could not open adapters file '", fname, "': ", strerror(errno));
    string line;
    string name;
    int num = 0;
    while (getline(f, line)) {
      if (line[0] == '>') {
        name = line;
        continue;
      }
      num++;
      if (line.length() < adapter_k) {
        SWARN("adapter seq for ", name, " is too short ", line.length(), " < ", adapter_k);
        continue;
      }
      new_seqs.push_back(line);
      sizes.push_back(line.size());
    }
  }
  //
  // broadcast the new sequences
  // non trivial broadcasts are not allowed (yet), so send a series of fixed-size broadcasts
  //

  // broadcast the number of sequences and allocate
  auto num_seqs = upcxx::broadcast(sizes.size(), 0).wait();
  sizes.resize(num_seqs);
  new_seqs.resize(num_seqs);
  adapter_seqs.reserve(adapter_seqs.size() + num_seqs * 2);

  // broadcast the sequence lengths
  upcxx::broadcast(sizes.data(), num_seqs, 0).wait();

  // broadcast the concatenated new sequences
  string concat_seqs;
  size_t total_seq_size = 0;
  for (int i = 0; i < num_seqs; i++) {
    concat_seqs += new_seqs[i];
    total_seq_size += sizes[i];
  }
  if (rank_me() == 0) assert(concat_seqs.size() == total_seq_size);
  concat_seqs.resize(total_seq_size);
  upcxx::broadcast(concat_seqs.data(), total_seq_size, 0).wait();

  // partition back into separate sequences
  size_t cursor = 0;
  for (int i = 0; i < num_seqs; i++) {
    if (rank_me() == 0) assert(new_seqs[i].size() == sizes[i]);
    new_seqs[i].resize(sizes[i]);
    if (rank_me() == 0) assert(new_seqs[i].compare(concat_seqs.substr(cursor, sizes[i])) == 0);
    new_seqs[i] = concat_seqs.substr(cursor, sizes[i]);
    cursor += sizes[i];
    adapter_seqs.push_back(new_seqs[i]);
    // revcomped adapters are very rare, so we don't bother with them
    // insert both kmer and kmer revcomp so we don't have to revcomp kmers in reads, which takes more time and since this
    // is such a small table storing both kmer and kmer_rc is fine
    adapter_seqs.push_back(revcomp(new_seqs[i]));
  }
  concat_seqs.clear();

  for (auto const &seq : adapter_seqs) {
    vector<Kmer<MAX_ADAPTER_K>> kmers;
    Kmer<MAX_ADAPTER_K>::set_k(adapter_k);
    Kmer<MAX_ADAPTER_K>::get_kmers(adapter_k, seq, kmers, false);
    for (int i = 0; i < kmers.size(); i++) {
      auto kmer = kmers[i];
      auto it = adapters.find(kmer);

      if (it == adapters.end())
        adapters.insert({kmer, {string_view(seq)}});
      else
        it->second.push_back(string_view(seq));
    }
  }
  SLOG_VERBOSE("Loaded ", adapter_seqs.size() / 2, " adapters, with a total of ", adapters.size(), " kmers\n");
  /*
  #ifdef DEBUG
  barrier();
  if (!rank_me()) {
    for (auto [kmer, seqs] : adapters) {
      for (auto seq : seqs) {
        DBG("adapter: ", kmer, " ", seq.first, " ", seq.second, "\n");
      }
    }
  }
  barrier();
  #endif
  */
}

static bool trim_adapters(StripedSmithWaterman::Aligner &ssw_aligner, StripedSmithWaterman::Filter &ssw_filter,
                          adapter_hash_table_t &adapters, const string &rname, string &seq, bool is_read_1, int adapter_k,
                          int64_t &bases_trimmed, int64_t &reads_removed, BaseTimer &time_overhead, BaseTimer &time_ssw) {
  time_overhead.start();
  vector<Kmer<MAX_ADAPTER_K>> kmers;
  // Kmer<MAX_ADAPTER_K>::set_k(adapter_k);
  Kmer<MAX_ADAPTER_K>::get_kmers(adapter_k, seq, kmers, false);
  double best_identity = 0;
  int best_trim_pos = seq.length();
  string best_adapter_seq;
  HASH_TABLE<string_view, bool> adapters_matching;
  for (auto &kmer : kmers) {
    auto it = adapters.find(kmer);
    if (it != adapters.end()) {
      for (auto &adapter_seq : it->second) {
        if (adapters_matching.find(adapter_seq) != adapters_matching.end()) continue;
        time_ssw.start();
        adapters_matching[adapter_seq] = true;
        StripedSmithWaterman::Alignment ssw_aln;
        ssw_aligner.Align(adapter_seq.data(), adapter_seq.length(), seq.data(), seq.length(), ssw_filter, &ssw_aln,
                          max((int)(seq.length() / 2), 15));
        int max_match_len = min(adapter_seq.length(), seq.length() - ssw_aln.ref_begin);
        double identity = (double)ssw_aln.sw_score / (double)ALN_MATCH_SCORE / (double)(max_match_len);
        if (identity >= best_identity) {
          best_identity = identity;
          best_trim_pos = ssw_aln.ref_begin;
          best_adapter_seq = adapter_seq;
        }
        time_ssw.stop();
        break;
      }
    }
  }
  time_overhead.stop();
  if (best_identity >= 0.5) {
    if (best_trim_pos < 12) best_trim_pos = 0;
    //DBG("Read ", rname, " is trimmed at ", best_trim_pos, " best identity ", best_identity, "\n", best_adapter_seq, "\n", seq, "\n");
    if (!best_trim_pos) reads_removed++;
    bases_trimmed += seq.length() - best_trim_pos;
    seq.resize(best_trim_pos);
    return true;
  }
  return false;
}

void merge_reads(vector<string> reads_fname_list, int qual_offset, double &elapsed_write_io_t,
                 vector<PackedReads *> &packed_reads_list, bool checkpoint, const string &adapter_fname, int min_kmer_len,
                 int subsample_pct) {
  assert(subsample_pct > 0 && subsample_pct <= 100);
  BarrierTimer timer(__FILEFUNC__);
  Timer merge_time(__FILEFUNC__ + " merging all");

  string fake_qual;
  fake_qual += (char)qual_offset;

  adapter_sequences_t adapter_seqs;
  adapter_hash_table_t adapters;
  StripedSmithWaterman::Aligner ssw_aligner(ALN_MATCH_SCORE, ALN_MISMATCH_COST, ALN_GAP_OPENING_COST, ALN_GAP_EXTENDING_COST,
                                            ALN_AMBIGUITY_COST);
  StripedSmithWaterman::Filter ssw_filter;

  ssw_filter.report_cigar = false;
  if (!adapter_fname.empty()) {
    load_adapter_seqs(adapter_fname, adapter_seqs, adapters, min_kmer_len);
  }

  size_t total_size = FastqReaders::open_all(reads_fname_list, subsample_pct);
  vector<string> merged_reads_fname_list;

  using shared_of = shared_ptr<upcxx_utils::dist_ofstream>;
  std::vector<shared_of> all_outputs;

  int64_t tot_bytes_read = 0;
  int64_t tot_num_ambiguous = 0;
  int64_t tot_num_merged = 0;
  int tot_max_read_len = 0;
  int64_t tot_bases = 0;
  // for unique read id need to estimate number of reads in our sections of all files
  auto [my_num_reads_estimate, read_len] = estimate_num_reads(reads_fname_list);
  auto max_num_reads = reduce_all(my_num_reads_estimate, op_fast_max).wait();
  auto tot_num_reads = reduce_all(my_num_reads_estimate, op_fast_add).wait();
  SLOG_VERBOSE("Estimated total number of reads as ", tot_num_reads, ", and max for any rank ", max_num_reads, "\n");
  // 2 reads per pair, 5x the block size estimate to be sure that we have no overlap. The read ids do not have to be contiguous
  auto read_id_block = (max_num_reads + 10000) * 2 * 5;
  uint64_t read_id = rank_me() * read_id_block;
  uint64_t start_read_id = read_id;
  DBG("starting read_id=", start_read_id, " max_num_reads=", max_num_reads, " read_id_block=", read_id_block, "\n");
  IntermittentTimer dump_reads_t("dump_reads");
  future<> wrote_all_files_fut = make_future();
  promise<> summary_promise;
  future<> fut_summary = summary_promise.get_future();
  int ri = 0;
  ProgressBar progbar(total_size / rank_n(), "Merging reads");
  for (auto const &reads_fname : reads_fname_list) {
    Timer merge_file_timer("merging " + get_basename(reads_fname));
    BaseTimer trim_timer("Adapter Trim overhead");
    BaseTimer trim_timer_ssw("Adapter Trim SSW");
    merge_file_timer.initiate_entrance_reduction();

    string out_fname = get_merged_reads_fname(reads_fname);
    if (file_exists(out_fname)) SWARN("File ", out_fname, " already exists, will overwrite...");

    FastqReader &fqr = FastqReaders::get(reads_fname);
    fqr.advise(true);
    auto my_file_size = fqr.my_file_size();

    shared_of sh_out_file;
    if (checkpoint) {
      auto merged_name = get_merged_reads_fname(reads_fname);
      sh_out_file = make_shared<upcxx_utils::dist_ofstream>(merged_name);
      all_outputs.push_back(sh_out_file);
      merged_reads_fname_list.push_back(merged_name);
    }
    int max_read_len = 0;
    int64_t overlap_len = 0;
    int64_t merged_len = 0;

    const int16_t MIN_OVERLAP = 12;
    const int16_t EXTRA_TEST_OVERLAP = 2;
    const int16_t MAX_MISMATCHES = 3;  // allow up to 3 mismatches, with MAX_PERROR
    const int Q2PerrorSize = sizeof(Q2Perror) / sizeof(*Q2Perror);
    assert(qual_offset == 33 || qual_offset == 64);

    // illumina reads generally accumulate errors at the end, so allow more mismatches in the overlap as long as differential
    // quality indicates a clear winner
    const double MAX_PERROR = 0.025;  // max 2.5% accumulated mismatch prob of error within overlap by differential quality score
    const int16_t EXTRA_MISMATCHES_PER_1000 = (int)150;  // allow addtl mismatches per 1000 bases overlap before aborting test
    const uint8_t MAX_MATCH_QUAL = 41 + qual_offset;

    string id1, seq1, quals1, id2, seq2, quals2, tmp_id, tmp_seq, tmp_quals;
    int64_t num_pairs = 0;
    int64_t bytes_read = 0;
    int64_t num_ambiguous = 0;
    int64_t num_merged = 0;
    int64_t num_reads = 0;
    DBG("Starting merge on ", fqr.get_fname(), " read_id=", read_id, " tell=", (int64_t)(fqr.my_file_size() > 0 ? fqr.tellg() : -1),
        " sz=", fqr.my_file_size(), "\n");
    int64_t bases_trimmed = 0;
    int64_t reads_removed = 0;
    int64_t bases_read = 0;
    int64_t missing_read1 = 0;
    int64_t missing_read2 = 0;

    bool skip_read1 = false, skip_read2 = false;
    for (;; num_pairs++) {
      discharge();
      //DBG_VERBOSE("Merging num_pair=", num_pairs, " read_id=", read_id, "\n");
      if (!fqr.is_paired()) {
        // unpaired reads get dummy read2 just like merged reads
        int64_t bytes_read1 = fqr.get_next_fq_record(id1, seq1, quals1);
        if (!bytes_read1) {
          DBG("Found end on ", fqr.get_fname(), " after read_id=", read_id, "\n");
          break;
        };
        bytes_read += bytes_read1;
        progbar.update(bytes_read);
        packed_reads_list[ri]->add_read("r" + to_string(read_id) + "/1", seq1, quals1);
        packed_reads_list[ri]->add_read("r" + to_string(read_id) + "/2", "N", fake_qual);
        read_id += 2;
        if (checkpoint) {
          *sh_out_file << "@r" << read_id << "/1\n" << seq1 << "\n+\n" << quals1 << "\n";
          *sh_out_file << "@r" << read_id << "/2\nN\n+\n" << fake_qual << "\n";
        }
        continue;
      }
      if (skip_read1 && skip_read2) break;
      if (!skip_read1) {
        int64_t bytes_read1 = fqr.get_next_fq_record(id1, seq1, quals1);
        if (!bytes_read1) break;  // end of file
        //DBG("Read1: ", id1, " ", seq1.length(), "\n");

        bytes_read += bytes_read1;
        bases_read += seq1.length();
      } else {
        // use the last read as read1
        assert(!tmp_id.empty());
        //DBG("Using deferred Read1: ", tmp_id, " ", tmp_seq.length(), "\n");
        id1 = tmp_id;
        seq1 = tmp_seq;
        quals1 = tmp_quals;
        tmp_id.clear();
        skip_read1 = false;
      }
      id2.clear();

      if (id1.length() > 2 && id1[id1.length() - 1] == '1') {
        skip_read2 = false;
      } else {
        assert(id1.empty() || id1[id1.length() - 1] == '2');
        //DBG("Missing read1, faking it\n");
        // got read 2: missing read 1 of expected pair! (Issue 117 to be robust to this missing read)
        missing_read1++;
        // set read2
        id2 = id1;
        seq2 = seq1;
        quals2 = quals1;
        // generate a fake read1
        id1[id1.length() - 1] = '1';
        seq1 = "N";
        quals1 = fake_qual;
        skip_read2 = true;
      }

      if (!skip_read2) {
        int64_t bytes_read2 = fqr.get_next_fq_record(id2, seq2, quals2);
        if (!bytes_read2) {
          // record missing read2
          id2.clear();
          skip_read1 = skip_read2 = true;
        }
        //DBG("Read2: ", id2, " ", seq2.length(), "\n");
        bytes_read += bytes_read2;
        bases_read += seq2.length();
      } else {
        skip_read2 = false;
      }

      if (id2.length() > 2 && id2[id2.length() - 1] == '2' && id1.compare(0, id1.length() - 1, id2, 0, id2.length() - 1) == 0) {
        skip_read1 = false;
      } else {
        if (skip_read1 && (skip_read2 || id2.empty())) break;  // end of file
        assert(id2.empty() || id2[id2.length() - 1] == '1' ||
               id2[id2.length() - 1] == '2');  // can miss both this read2 and the next read1, getting the next read2
        //DBG("Missing read2, faking it\n");
        // got read1 : missing read2 of expected pair! (Issue 117 to be robust to this missing read)
        missing_read2++;
        // preserve this as the *next* read1 (may actually be a read2)
        assert(tmp_id.empty());
        tmp_id = id2;
        tmp_seq = seq2;
        tmp_quals = quals2;
        // generate a fake read2
        if (id2.empty()) skip_read2 = true;  // end of file
        id2 = id1;
        id2[id2.length() - 1] = '2';
        seq2 = "N";
        quals2 = fake_qual;
        skip_read1 = true;
      }

      progbar.update(bytes_read);

      if (id1.compare(0, id1.length() - 2, id2, 0, id2.length() - 2) != 0) DIE("Mismatched pairs ", id1, " ", id2);
      if (id1[id1.length() - 1] != '1' || id2[id2.length() - 1] != '2') DIE("Mismatched pair numbers ", id1, " ", id2);

      if (!adapters.empty()) {
        bool trim1 = trim_adapters(ssw_aligner, ssw_filter, adapters, id1, seq1, true, min_kmer_len, bases_trimmed, reads_removed,
                                   trim_timer, trim_timer_ssw);
        bool trim2 = trim_adapters(ssw_aligner, ssw_filter, adapters, id2, seq2, false, min_kmer_len, bases_trimmed, reads_removed,
                                   trim_timer, trim_timer_ssw);
        // trim to same length - like the tpe option in bbduk
        if ((trim1 || trim2) && seq1.length() > 1 && seq2.length() > 1) {
          auto min_seq_len = min(seq1.length(), seq2.length());
          seq1.resize(min_seq_len);
          seq2.resize(min_seq_len);
          quals1.resize(min_seq_len);
          quals2.resize(min_seq_len);
        }
        // it's possible that really short reads could be merged, but unlikely and they'd still be short, so drop all below min
        // kmer length
        if (seq1.length() < min_kmer_len && seq2.length() < min_kmer_len) continue;
      }

      bool is_merged = 0;
      int8_t abort_merge = 0;

      // revcomp the second mate pair and reverse the second quals
      string rc_seq2 = revcomp(seq2);
      string rev_quals2 = quals2;
      reverse(rev_quals2.begin(), rev_quals2.end());

      // use start_i to offset inequal lengths which can be very different but still overlap near the end.  250 vs 178..
      int16_t len = (rc_seq2.length() < seq1.length()) ? rc_seq2.length() : seq1.length();
      int16_t start_i = ((len == (int16_t)seq1.length()) ? 0 : seq1.length() - len);
      int16_t found_i = -1;
      int16_t best_i = -1;
      int16_t best_mm = len;
      double best_perror = -1.0;

      // slide along seq1
      for (int16_t i = 0; i < len - MIN_OVERLAP + EXTRA_TEST_OVERLAP; i++) {  // test less overlap than MIN_OVERLAP
        if (abort_merge) break;
        int16_t overlap = len - i;
        int16_t this_max_mismatch = MAX_MISMATCHES + (EXTRA_MISMATCHES_PER_1000 * overlap / 1000);
        int16_t error_max_mismatch = this_max_mismatch * 4 / 3 + 1;  // 33% higher
        if (fast_count_mismatches(seq1.c_str() + start_i + i, rc_seq2.c_str(), overlap, error_max_mismatch) > error_max_mismatch)
          continue;
        int16_t matches = 0, mismatches = 0, bothNs = 0, Ncount = 0;
        int16_t overlapChecked = 0;
        double perror = 0.0;
        for (int16_t j = 0; j < overlap; j++) {
          overlapChecked++;
          char ps = seq1[start_i + i + j];
          char rs = rc_seq2[j];
          if (ps == rs) {
            matches++;
            if (ps == 'N') {
              Ncount += 2;
              if (bothNs++) {
                abort_merge++;
                num_ambiguous++;
                break;  // do not match multiple Ns in the same position -- 1 is okay
              }
            }
          } else {
            mismatches++;
            if (ps == 'N') {
              mismatches++;  // N still counts as a mismatch
              Ncount++;
              quals1[start_i + i + j] = qual_offset;
              assert(rev_quals2[j] - qual_offset < Q2PerrorSize);
              assert(rev_quals2[j] - qual_offset >= 0);
              perror += Q2Perror[rev_quals2[j] - qual_offset];
            } else if (rs == 'N') {
              Ncount++;
              mismatches++;  // N still counts as a mismatch
              rev_quals2[j] = qual_offset;
              assert(quals1[start_i + i + j] - qual_offset < Q2PerrorSize);
              assert(quals1[start_i + i + j] - qual_offset >= 0);
              perror += Q2Perror[quals1[start_i + i + j] - qual_offset];
            }
            if (MAX_PERROR > 0.0) {
              assert(quals1[start_i + i + j] >= qual_offset);
              assert(rev_quals2[j] >= qual_offset);
              uint8_t q1 = quals1[start_i + i + j] - qual_offset;
              uint8_t q2 = rev_quals2[j] - qual_offset;
              if (q1 < 0 || q2 < 0 || q1 >= Q2PerrorSize || q2 >= Q2PerrorSize)
                DIE("Invalid quality score for read ", id1, " '", quals1[start_i + i + j], "' ", id2, " '", rev_quals2[j],
                    "' assuming common qual_offset of ", qual_offset,
                    ". Check the data and make sure it follows a single consistent quality scoring model ",
                    "(phred+64 vs. phred+33)");

              // sum perror as the difference in q score perrors
              uint8_t diffq = (q1 > q2) ? q1 - q2 : q2 - q1;
              if (diffq <= 2) {
                perror += 0.5;  // cap at flipping a coin when both quality scores are close
              } else {
                assert(diffq < Q2PerrorSize);
                perror += Q2Perror[diffq];
              }
            }
          }
          if (Ncount > 3) {
            abort_merge++;
            num_ambiguous++;
            break;  // do not match reads with many Ns
          }
          if (mismatches > error_max_mismatch) break;
        }
        int16_t match_thres = overlap - this_max_mismatch;
        if (match_thres < MIN_OVERLAP) match_thres = MIN_OVERLAP;
        if (matches >= match_thres && overlapChecked == overlap && mismatches <= this_max_mismatch &&
            perror / overlap <= MAX_PERROR) {
          if (best_i < 0 && found_i < 0) {
            best_i = i;
            best_mm = mismatches;
            best_perror = perror;
          } else {
            // another good or ambiguous overlap detected
            num_ambiguous++;
            best_i = -1;
            best_mm = len;
            best_perror = -1.0;
            break;
          }
        } else if (overlapChecked == overlap && mismatches <= error_max_mismatch && perror / overlap <= MAX_PERROR * 4 / 3) {
          // lower threshold for detection of an ambigious overlap
          found_i = i;
          if (best_i >= 0) {
            // ambiguous mapping found after a good one was
            num_ambiguous++;
            best_i = -1;
            best_mm = len;
            best_perror = -1.0;
            break;
          }
        }
      }

      if (best_i >= 0 && !abort_merge) {
        int16_t i = best_i;
        int16_t overlap = len - i;
        // pick the base with the highest quality score for the overlapped region
        for (int16_t j = 0; j < overlap; j++) {
          if (seq1[start_i + i + j] == rc_seq2[j]) {
            // match boost quality up to the limit
            uint16_t newQual = quals1[start_i + i + j] + rev_quals2[j] - qual_offset;
            quals1[start_i + i + j] = ((newQual > MAX_MATCH_QUAL) ? MAX_MATCH_QUAL : newQual);
            assert(quals1[start_i + i + j] >= quals1[start_i + i + j]);
            // FIXME: this fails for a CAMISIM generated dataset. I don't even know what this is checking...
            // assert(quals1[start_i + i + j] >= rev_quals2[j]);
          } else {
            uint8_t newQual;
            if (quals1[start_i + i + j] < rev_quals2[j]) {
              // use rev base and discount quality
              newQual = rev_quals2[j] - quals1[start_i + i + j] + qual_offset;
              seq1[start_i + i + j] = rc_seq2[j];
            } else {
              // keep prev base, but still discount quality
              newQual = quals1[start_i + i + j] - rev_quals2[j] + qual_offset;
            }
            // a bit better than random chance here
            quals1[start_i + i + j] = ((newQual > (2 + qual_offset)) ? newQual : (2 + qual_offset));
          }
          assert(quals1[start_i + i + j] >= qual_offset);
        }

        // include the remainder of the rc_seq2 and quals
        seq1 = seq1.substr(0, start_i + i + overlap) + rc_seq2.substr(overlap);
        quals1 = quals1.substr(0, start_i + i + overlap) + rev_quals2.substr(overlap);

        is_merged = true;
        num_merged++;

        int read_len = seq1.length();  // caculate new merged length
        if (max_read_len < read_len) max_read_len = read_len;
        merged_len += read_len;
        overlap_len += overlap;

        packed_reads_list[ri]->add_read("r" + to_string(read_id) + "/1", seq1, quals1);
        packed_reads_list[ri]->add_read("r" + to_string(read_id) + "/2", "N", fake_qual);
        if (checkpoint) {
          *sh_out_file << "@r" << read_id << "/1\n" << seq1 << "\n+\n" << quals1 << "\n";
          *sh_out_file << "@r" << read_id << "/2\nN\n+\n" << fake_qual << "\n";
        }
      }
      if (!is_merged) {
        // write without the revcomp
        packed_reads_list[ri]->add_read("r" + to_string(read_id) + "/1", seq1, quals1);
        packed_reads_list[ri]->add_read("r" + to_string(read_id) + "/2", seq2, quals2);
        if (checkpoint) {
          *sh_out_file << "@r" << read_id << "/1\n" << seq1 << "\n+\n" << quals1 << "\n";
          *sh_out_file << "@r" << read_id << "/2\n" << seq2 << "\n+\n" << quals2 << "\n";
        }
      }
      // inc by 2 so that we can use a later optimization of treating the even as /1 and the odd as /2
      read_id += 2;
    }
    DBG("Merged my set of reads. num_merged=", num_merged, " num_ambig=", num_ambiguous, " bytes_read=", bytes_read, "\n");
    fqr.advise(false);  // free kernel memory

    if (checkpoint) {
      // close this file, but do not wait for it yet
      dump_reads_t.start();
      wrote_all_files_fut = when_all(wrote_all_files_fut, sh_out_file->close_async());
      dump_reads_t.stop();
    }

    tot_num_merged += num_merged;
    tot_num_ambiguous += num_ambiguous;
    tot_max_read_len = std::max(tot_max_read_len, max_read_len);
    tot_bytes_read += bytes_read;
    tot_bases += bases_read;

    // start the collective reductions
    // delay the summary output for when they complete
    auto fut_sh_ssw_timings = trim_timer_ssw.reduce_timings();
    auto fut_sh_trim_overhead = trim_timer.reduce_timings();
    auto fut_reductions = when_all(
        reduce_one(num_pairs, op_fast_add, 0), reduce_one(num_merged, op_fast_add, 0), reduce_one(num_ambiguous, op_fast_add, 0),
        reduce_one(merged_len, op_fast_add, 0), reduce_one(overlap_len, op_fast_add, 0), reduce_one(max_read_len, op_fast_max, 0),
        reduce_one(bases_trimmed, op_fast_add, 0), reduce_one(reads_removed, op_fast_add, 0),
        reduce_one(bases_read, op_fast_add, 0), reduce_one(bytes_read, op_fast_add, 0), reduce_one(missing_read1, op_fast_add, 0),
        reduce_one(missing_read2, op_fast_add, 0), fut_sh_ssw_timings, fut_sh_trim_overhead);

    fut_summary = when_all(fut_summary, fut_reductions)
                      .then([reads_fname, bytes_read, &adapters](
                                int64_t all_num_pairs, int64_t all_num_merged, int64_t all_num_ambiguous, int64_t all_merged_len,
                                int64_t all_overlap_len, int all_max_read_len, int64_t all_bases_trimmed, int64_t all_reads_removed,
                                int64_t all_bases_read, int64_t all_bytes_read, int64_t all_missing_read1,
                                int64_t all_missing_read2, ShTimings sh_ssw_timings, ShTimings sh_trim_overhead) {
                        SLOG_VERBOSE("Merged reads in file ", reads_fname, ":\n");
                        SLOG_VERBOSE("  merged ", perc_str(all_num_merged, all_num_pairs), " of ", all_num_pairs, " pairs\n");
                        SLOG_VERBOSE("  ambiguous ", perc_str(all_num_ambiguous, all_num_pairs), " ambiguous pairs\n");
                        SLOG_VERBOSE("  missing pair1 ", all_missing_read1, " pair2 ", all_missing_read2, "\n");
                        SLOG_VERBOSE("  average merged length ", (double)all_merged_len / all_num_merged, "\n");
                        SLOG_VERBOSE("  average overlap length ", (double)all_overlap_len / all_num_merged, "\n");
                        if (!adapters.empty()) {
                          SLOG_VERBOSE("  adapter bases trimmed ", perc_str(all_bases_trimmed, all_bases_read), "\n");
                          SLOG_VERBOSE("  adapter reads removed ", perc_str(all_reads_removed, all_num_pairs * 2), "\n");
                          SLOG_VERBOSE("  adapter SSW timings: ", sh_ssw_timings->to_string(), "\n");
                          SLOG_VERBOSE("  adapter trim total overhead: ", sh_trim_overhead->to_string(), "\n");
                        }
                        SLOG_VERBOSE("  max read length ", all_max_read_len, "\n");
                        SLOG_VERBOSE("Rank0 bytes read ", bytes_read, " of ", all_bytes_read, "\n");
                      });

    num_reads += num_pairs * 2;
    ri++;
    FastqReaders::close(reads_fname);
  }
  auto prog_done = progbar.set_done();
  wrote_all_files_fut = when_all(wrote_all_files_fut, prog_done);

  DBG("last read_id=", read_id, " last should not be > ", start_read_id + read_id_block, "\n");
  if (read_id >= start_read_id + read_id_block) WARN("Invalid read_id=", read_id, " start_read_id=", start_read_id, " read_id_block=", read_id_block, "\n");
  assert(read_id < start_read_id + read_id_block);
  merge_time.initiate_exit_reduction();

  //#ifdef DEBUG
  // ensure there is no overlap in read_ids which will cause a crash later
  using SSPair = std::pair<uint64_t, uint64_t>;
  SSPair start_stop(start_read_id, read_id);
  dist_object<SSPair> dist_ss(world(), start_stop);
  future<> rpc_tests = make_future();
  // check next rank
  assert(dist_ss->first <= dist_ss->second);
  if (rank_me() < rank_n() - 1) {
    auto fut = rpc(
        rank_me() + 1,
        [](dist_object<pair<uint64_t, uint64_t>> &dist_ss, SSPair ss) {
          if (!(ss.first < dist_ss->first && ss.second < dist_ss->first))
            DIE("Invalid read ids from previous rank: ", rank_me(), "=", dist_ss->first, "-", dist_ss->second,
                " prev rank=", ss.first, "-", ss.second, "\n");
        },
        dist_ss, *dist_ss);
    rpc_tests = when_all(rpc_tests, fut);
  }
  if (rank_me() > 0) {
    auto fut = rpc(
        rank_me() - 1,
        [](dist_object<pair<uint64_t, uint64_t>> &dist_ss, SSPair ss) {
          if (!(ss.first > dist_ss->second && ss.second > dist_ss->second))
            DIE("Invalid read ids from next rank: ", rank_me(), "=", dist_ss->first, "-", dist_ss->second, " next rank=", ss.first,
                "-", ss.second, "\n");
        },
        dist_ss, *dist_ss);
    rpc_tests = when_all(rpc_tests, fut);
  }
  rpc_tests.wait();
  //#endif

  // finish all file writing and report
  dump_reads_t.start();
  wrote_all_files_fut.wait();
  for (auto sh_of : all_outputs) {
    wrote_all_files_fut = when_all(wrote_all_files_fut, sh_of->report_timings());
  }
  wrote_all_files_fut.wait();

  dump_reads_t.stop();
  elapsed_write_io_t = dump_reads_t.get_elapsed();
  dump_reads_t.done();

  summary_promise.fulfill_anonymous(1);
  fut_summary.wait();

  timer.initate_exit_barrier();
}
