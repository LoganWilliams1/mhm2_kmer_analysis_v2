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

#include "upcxx_utils/log.hpp"
#include "adapters.hpp"

using namespace std;

void Adapters::load_adapter_seqs(const string &fname) {
  kmer_adapter_map.clear();

  // avoid every rank reading this small file
  adapter_sequences_t new_seqs;
  vector<uint64_t> sizes;
  if (!rank_me()) {
    ifstream f(fname);
    if (!f.is_open()) DIE("Could not open adapters file '", fname, "': ", strerror(errno));
    string line;
    string name;
    int num = 0;
    int num_short = 0;
    while (getline(f, line)) {
      if (line[0] == '>') {
        name = line;
        continue;
      }
      num++;
      if (line.length() < adapter_k) {
        num_short++;
        continue;
      }
      new_seqs.push_back(line);
      sizes.push_back(line.size());
    }
    if (num_short) SLOG_VERBOSE("Ignoring ", num_short, " adapters of length less than ", adapter_k, "\n");
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
  uint64_t total_seq_size = 0;
  for (int i = 0; i < num_seqs; i++) {
    concat_seqs += new_seqs[i];
    total_seq_size += sizes[i];
  }
  if (rank_me() == 0) assert(concat_seqs.size() == total_seq_size);
  concat_seqs.resize(total_seq_size);
  upcxx::broadcast(concat_seqs.data(), total_seq_size, 0).wait();

  // partition back into separate sequences
  uint64_t cursor = 0;
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

  for (int i = 0; i < adapter_seqs.size(); i++) {
    auto &seq = adapter_seqs[i];
    vector<Kmer<MAX_ADAPTER_K>> kmers;
    Kmer<MAX_ADAPTER_K>::set_k(adapter_k);
    Kmer<MAX_ADAPTER_K>::get_kmers(adapter_k, seq, kmers, false);
    for (int j = 0; j < kmers.size(); j++) {
      auto kmer = kmers[j];
      auto it = kmer_adapter_map.find(kmer);
      if (it == kmer_adapter_map.end())
        kmer_adapter_map.insert({kmer, {{i, j}}});
      else
        it->second.push_back({i, j});
    }
  }
  SLOG_VERBOSE("Loaded ", adapter_seqs.size() / 2, " adapters, with a total of ", kmer_adapter_map.size(), " kmers\n");
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

bool Adapters::trim_seq(const string &id, string &seq, bool is_read_1) {
  trim_timer.start();
  vector<Kmer<MAX_ADAPTER_K>> kmers;
  Kmer<MAX_ADAPTER_K>::get_kmers(adapter_k, seq, kmers, false);
  double best_identity = 0;
  int best_match_len = 0;
  int best_trim_pos = seq.length();
  string best_adapter_seq = "";

  vector<bool> adapters_matching(adapter_seqs.size(), false);
  bool found = false;
#ifdef MERGE_READS_TRIM_WITH_SSW
  const int STEP = 4;
#else
  const int STEP = 1;
#endif
  for (int i = 0; i < kmers.size(); i += STEP) {
    auto &kmer = kmers[i];
    auto it = kmer_adapter_map.find(kmer);
    if (it != kmer_adapter_map.end()) {
      for (auto adapter_record : it->second) {
        int adapter_index = adapter_record.first;
        int kmer_offset = adapter_record.second;
        if (adapters_matching[adapter_index]) continue;
        auto &adapter_seq = adapter_seqs[adapter_index];
        trim_timer_ssw.start();
        adapters_matching[adapter_index] = true;
#ifdef MERGE_READS_TRIM_WITH_SSW
        StripedSmithWaterman::Alignment ssw_aln;

        int adapter_seq_start = max(0, kmer_offset - i - 2);
        int adapter_seq_match_len = min(adapter_seq_start + seq.length() + 2, adapter_seq.length());
        auto adapter_subseq = adapter_seq.substr(adapter_seq_start, adapter_seq_match_len);

        // FIXME: use the kmer location to align only the necessary subsequence in the adapter seq
        ssw_aligner.Align(adapter_subseq.data(), adapter_subseq.length(), seq.data(), seq.length(), ssw_filter, &ssw_aln,
                          max((int)(seq.length() / 2), 15));

        int max_match_len = min(adapter_seq.length(), seq.length() - ssw_aln.ref_begin);
        double identity = (double)ssw_aln.sw_score / (double)ssw_aligner.get_match_score() / (double)(max_match_len);
        if (identity >= best_identity) {
          best_identity = identity;
          best_trim_pos = ssw_aln.ref_begin;
          best_adapter_seq = adapter_seq;
          if (identity > 0.97) found = true;
        }
#else
        int num_mismatches = 0;
        for (int j = 0;; j++) {
          int seq_pos = adapter_k + i + j;
          int adapter_pos = adapter_k + kmer_offset + j;
          if (seq_pos >= seq.length() || adapter_pos >= adapter_seq.length()) break;
          if (adapter_seq[adapter_pos] != seq[seq_pos]) {
            num_mismatches++;
            if (num_mismatches > 1) {
              int match_len = adapter_k + j;
              if (match_len > best_match_len) {
                best_identity = (double)match_len / (double)adapter_seq.length();
                best_trim_pos = i;
                best_adapter_seq = adapter_seq;
                best_match_len = match_len;
                if (match_len >= adapter_seq.length() - 1) found = true;
              }
              break;
            }
          }
        }
#endif
        trim_timer_ssw.stop();
        break;
      }
    }
    if (found) break;
  }
  trim_timer.stop();

  if (best_identity >= 0.5) {
    if (best_trim_pos < 12) best_trim_pos = 0;
    // DBG("Read ", rname, " is trimmed at ", best_trim_pos, " best identity ", best_identity, "\n", best_adapter_seq, "\n", seq,
    // "\n");
    if (!best_trim_pos) reads_removed++;
    bases_trimmed += seq.length() - best_trim_pos;
    seq.resize(best_trim_pos);
    return true;
  }
  return false;
}

Adapters::Adapters(int adapter_k, const string &fname, bool use_blastn_scores)
    : adapter_k(adapter_k)
    , use_blastn_scores(use_blastn_scores) {
  ssw_aligner.Clear();
  if (!ssw_aligner.ReBuild(to_string(use_blastn_scores ? BLASTN_ALN_SCORES : ALTERNATE_ALN_SCORES)))
    SDIE("Failed to set aln scores");
  ssw_filter.report_cigar = false;
  if (!fname.empty()) load_adapter_seqs(fname);
}

void Adapters::done(size_t all_bases_read, size_t all_num_pairs) {
  if (adapter_seqs.empty()) return;
  auto ssw_elapsed_time = trim_timer_ssw.get_elapsed();
  auto trim_elapsed_time = trim_timer.get_elapsed();
  auto all_bases_trimmed = reduce_one(bases_trimmed, op_fast_add, 0).wait();
  auto all_reads_removed = reduce_one(reads_removed, op_fast_add, 0).wait();
  SLOG_VERBOSE("Adapter trimming:\n");
  SLOG_VERBOSE("  bases trimmed ", upcxx_utils::perc_str(all_bases_trimmed, all_bases_read), "\n");
  SLOG_VERBOSE("  reads removed ", upcxx_utils::perc_str(all_reads_removed, all_num_pairs * 2), "\n");
  SLOG_VERBOSE("  SSW time ", std::fixed, std::setprecision(3), ssw_elapsed_time, " s\n");
  SLOG_VERBOSE("  total time ", std::fixed, std::setprecision(3), trim_elapsed_time, " s\n");
}

bool Adapters::trim(const string &id1, string &seq1, string &quals1, const string &id2, string &seq2, string &quals2) {
  if (!adapter_seqs.empty()) {
    bool trim1 = trim_seq(id1, seq1, true);
    bool trim2 = trim_seq(id2, seq2, false);
    // trim to same length - like the tpe option in bbduk
    if ((trim1 || trim2) && seq1.length() > 1 && seq2.length() > 1) {
      auto min_seq_len = min(seq1.length(), seq2.length());
      seq1.resize(min_seq_len);
      seq2.resize(min_seq_len);
      quals1.resize(min_seq_len);
      quals2.resize(min_seq_len);
    }
    return trim1 || trim2;
  }
  return false;
}