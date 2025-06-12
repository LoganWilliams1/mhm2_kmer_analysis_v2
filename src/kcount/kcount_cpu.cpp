/*
 HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
 through Lawrence Berkeley National Laboratory (subject to receipt of any required
 approvals from the U.S. Dept. of Energy).  All rights reserved.

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

#include "upcxx_utils.hpp"
#include "kcount.hpp"
#include "kmer_dht.hpp"
#include "prime.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

// #define SLOG_CPU_HT(...) SLOG(KLMAGENTA, __VA_ARGS__, KNORM)
#define SLOG_CPU_HT(...) SLOG_VERBOSE(__VA_ARGS__)

template <int MAX_K>
struct SeqBlockInserter<MAX_K>::SeqBlockInserterState {
  int64_t bytes_kmers_sent = 0;
  int64_t bytes_supermers_sent = 0;
  int64_t num_kmers = 0;
  vector<Kmer<MAX_K>> kmers;
};

template <int MAX_K>
SeqBlockInserter<MAX_K>::SeqBlockInserter(int qual_offset, int minimizer_len) {
  state = new SeqBlockInserterState();
}

template <int MAX_K>
SeqBlockInserter<MAX_K>::~SeqBlockInserter() {
  if (state) delete state;
}

template <int MAX_K>
void SeqBlockInserter<MAX_K>::process_seq(string &seq, kmer_count_t depth, dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  if (!depth) depth = 1;
  auto kmer_len = Kmer<MAX_K>::get_k();
  Kmer<MAX_K>::get_kmers(kmer_len, seq, state->kmers);
  for (size_t i = 0; i < state->kmers.size(); i++) {
    state->bytes_kmers_sent += sizeof(KmerAndExt<MAX_K>);
    Kmer<MAX_K> kmer_rc = state->kmers[i].revcomp();
    if (kmer_rc < state->kmers[i]) state->kmers[i] = kmer_rc;
  }

  Supermer supermer{.seq = seq.substr(0, kmer_len + 1), .count = (kmer_count_t)depth};
  auto prev_target_rank = kmer_dht->get_kmer_target_rank(state->kmers[1]);
  for (int i = 1; i < (int)(seq.length() - kmer_len); i++) {
    auto &kmer = state->kmers[i];
    auto target_rank = kmer_dht->get_kmer_target_rank(kmer);
    if (target_rank == prev_target_rank) {
      supermer.seq += seq[i + kmer_len];
    } else {
      state->bytes_supermers_sent += supermer.get_bytes();
      kmer_dht->add_supermer(supermer, prev_target_rank);
      supermer.seq = seq.substr(i - 1, kmer_len + 2);
      prev_target_rank = target_rank;
    }
  }
  if (supermer.seq.length() >= kmer_len + 2) {
    state->bytes_supermers_sent += supermer.get_bytes();
    kmer_dht->add_supermer(supermer, prev_target_rank);
  }
  state->num_kmers += seq.length() - 2 - kmer_len;
}

template <int MAX_K>
void SeqBlockInserter<MAX_K>::done_processing(dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  auto tot_supermers_bytes_sent = reduce_one(state->bytes_supermers_sent, op_fast_add, 0).wait();
  auto tot_kmers_bytes_sent = reduce_one(state->bytes_kmers_sent, op_fast_add, 0).wait();
  SLOG_CPU_HT("Total bytes sent in compressed supermers ", get_size_str(tot_supermers_bytes_sent), " (compression is ", fixed,
              setprecision(3), (double)tot_kmers_bytes_sent / tot_supermers_bytes_sent, " over kmers)\n");
  auto all_num_kmers = reduce_one(state->num_kmers, op_fast_add, 0).wait();
  SLOG_CPU_HT("Processed a total of ", all_num_kmers, " kmers\n");
}

struct ExtCounts {
  kmer_count_t count_A;
  kmer_count_t count_C;
  kmer_count_t count_G;
  kmer_count_t count_T;

  void set(uint16_t *counts) {
    count_A = counts[0];
    count_C = counts[1];
    count_G = counts[2];
    count_T = counts[3];
  }

  void set(uint32_t *counts) {
    count_A = static_cast<uint16_t>(counts[0]);
    count_C = static_cast<uint16_t>(counts[1]);
    count_G = static_cast<uint16_t>(counts[2]);
    count_T = static_cast<uint16_t>(counts[3]);
  }

  std::array<std::pair<char, int>, 4> get_sorted() {
    std::array<std::pair<char, int>, 4> counts = {std::make_pair('A', (int)count_A), std::make_pair('C', (int)count_C),
                                                  std::make_pair('G', (int)count_G), std::make_pair('T', (int)count_T)};
    std::sort(std::begin(counts), std::end(counts), [](const auto &elem1, const auto &elem2) {
      if (elem1.second == elem2.second)
        return elem1.first > elem2.first;
      else
        return elem1.second > elem2.second;
    });
    return counts;
  }

  bool is_zero() {
    if (count_A + count_C + count_G + count_T == 0) return true;
    return false;
  }

  kmer_count_t inc_with_limit(int count1, int count2) {
    count1 += count2;
    return std::min(count1, (int)std::numeric_limits<kmer_count_t>::max());
  }

  void inc(char ext, int count) {
    switch (ext) {
      case 'A': count_A = inc_with_limit(count_A, count); break;
      case 'C': count_C = inc_with_limit(count_C, count); break;
      case 'G': count_G = inc_with_limit(count_G, count); break;
      case 'T': count_T = inc_with_limit(count_T, count); break;
    }
  }

  void add(ExtCounts &ext_counts) {
    count_A = inc_with_limit(count_A, ext_counts.count_A);
    count_C = inc_with_limit(count_C, ext_counts.count_C);
    count_G = inc_with_limit(count_G, ext_counts.count_G);
    count_T = inc_with_limit(count_T, ext_counts.count_T);
  }

  char get_ext(kmer_count_t count) {
    auto sorted_counts = get_sorted();
    int top_count = sorted_counts[0].second;
    int runner_up_count = sorted_counts[1].second;
    // set dynamic_min_depth to 1.0 for single depth data (non-metagenomes)
    int dmin_dyn = std::max((int)((1.0 - DYN_MIN_DEPTH) * count), _dmin_thres);
    if (top_count < dmin_dyn) return 'X';
    if (runner_up_count >= dmin_dyn) return 'F';
    return sorted_counts[0].first;
  }

  string to_string() {
    ostringstream os;
    os << count_A << "," << count_C << "," << count_G << "," << count_T;
    return os.str();
  }
};

struct KmerExtsCounts {
  ExtCounts left_exts;
  ExtCounts right_exts;
  kmer_count_t count;
  bool from_ctg;

  char get_left_ext() { return left_exts.get_ext(count); }

  char get_right_ext() { return right_exts.get_ext(count); }
};

// template <int MAX_K>
// using KmerMapExts = HASH_TABLE<Kmer<MAX_K>, KmerExtsCounts>;

template <int MAX_K>
class KmerMapExts {
  size_t capacity = 0;
  size_t num_elems = 0;
  size_t num_dropped = 0;
  size_t num_singleton_overrides = 0;
  vector<Kmer<MAX_K>> keys;
  size_t sum_probe_lens = 0;
  size_t max_probe_len = 0;
  vector<KmerExtsCounts> counts;
  size_t iter_pos = 0;
  const int N_LONGS = Kmer<MAX_K>::get_N_LONGS();
  const uint64_t KEY_EMPTY = 0xffffffffffffffff;

 public:
  void reserve(size_t max_elems) {
    primes::Prime prime;
    prime.set(max_elems, true);
    capacity = prime.get();
    num_elems = 0;
    keys.resize(capacity);
    memset((void *)keys.data(), 0xff, sizeof(Kmer<MAX_K>) * capacity);
    counts.resize(capacity, {0});
    SLOG_CPU_HT("Capacity is set to ", capacity, " for ", max_elems, " max elements, size ",
                get_size_str(capacity * (sizeof(Kmer<MAX_K>) + sizeof(KmerExtsCounts))), "\n");
  }

  pair<KmerExtsCounts *, bool> insert(const Kmer<MAX_K> &kmer, bool override_singletons) {
    size_t slot = kmer.hash() % capacity;
    size_t start_slot = slot;
    const size_t MAX_PROBE = (capacity < KCOUNT_HT_MAX_PROBE ? capacity : KCOUNT_HT_MAX_PROBE);
    for (size_t i = 1; i <= MAX_PROBE; i++) {
      if (keys[slot].get_longs()[N_LONGS - 1] == KEY_EMPTY) {
        keys[slot] = kmer;
        sum_probe_lens += i;
        if (i > max_probe_len) max_probe_len = i;
        num_elems++;
        return {&(counts[slot]), true};
      } else if (kmer == keys[slot]) {
        return {&(counts[slot]), false};
      }
      slot = (slot + 1) % capacity;
    }
    // if we reach here we have not found the kmer and there is no space to insert it.
    // if overriding singletons, we have to repeat the insertion, but this time replacing the first singleton
    // this approach is more computationally expensive, but it allows us to reuse slots that would be purged for ctg kmers
    if (override_singletons) {
      // reset variables for search
      slot = start_slot;
      for (size_t i = 1; i <= MAX_PROBE; i++) {
        assert(kmer != keys[slot]);  // FIXME? probe_lens[slot] != 0
        if (counts[slot].count == 1) {
          num_singleton_overrides++;
          keys[slot] = kmer;
          if (i > max_probe_len) max_probe_len = i;
          sum_probe_lens += i;
          return {&(counts[slot]), true};
        }
        slot = (slot + 1) % capacity;
      }
    }
    num_dropped++;
    return {nullptr, false};
  }

  size_t size() { return num_elems; }

  double load_factor() { return (double)num_elems / capacity; }

  size_t get_num_dropped() { return num_dropped; }

  void clear_num_dropped() { num_dropped = 0; }

  size_t get_num_singleton_overrides() { return num_singleton_overrides; }

  size_t get_sum_probe_lens() { return sum_probe_lens; }

  size_t get_max_probe_len() { return max_probe_len; }

  void begin_iterate() { iter_pos = 0; }

  pair<Kmer<MAX_K> *, KmerExtsCounts *> get_next() {
    for (; iter_pos < capacity; iter_pos++) {
      if (keys[iter_pos].get_longs()[N_LONGS - 1] != KEY_EMPTY) {
        iter_pos++;
        return {&keys[iter_pos - 1], &counts[iter_pos - 1]};
      }
    }
    return {nullptr, nullptr};
  }
};

// Another idea is that we basically have two arrays, one which contains keys + one ext each side (packed into 1 byte), and the
// other which contains the full KmerExtCounts. Also, make the keys array say 2x bigger than the values array. Then when we
// insert for the first time, don't insert a full kmerextcounts, but rather just the kmer and its packed exts byte. On the
// subsequent inserts of the same kmer, then we insert the full extcounts. And given that at least 50% of all kmers are
// singletons, we should save a lot of space, e.g. for k=21, the kmer is 8 bytes, so the first array will require 9 bytes, whereas
// the second array will require 19 bytes. So for singletons, we use 9 bytes and others we use 9+19=28 bytes. Thus our memory
// requirement for 50% singletons is now 0.5*9+0.5*28=19, instead of 28, so a savings of ~44%. For kmer of length 2, we'd have
// 0.5*(17+36)=27 compared to 36, so a savings of 25%. Still pretty worthwhile. The only catch is that if we size the arrays so
// that the keys is 2x the non-keys, then we cannot save more than those fractions, even for higher singleton counts. eg if we
// have 80% singletons we'll still get 44% and 25% savings.

template <int MAX_K>
static void get_kmers_and_exts(Supermer &supermer, vector<KmerAndExt<MAX_K>> &kmers_and_exts) {
  vector<bool> quals;
  quals.resize(supermer.seq.length());
  for (int i = 0; i < supermer.seq.length(); i++) {
    quals[i] = isupper(supermer.seq[i]);
    if (supermer.seq[i] >= 'a' && supermer.seq[i] <= 'z') supermer.seq[i] += ('A' - 'a');
  }
  auto kmer_len = Kmer<MAX_K>::get_k();
  vector<Kmer<MAX_K>> kmers;
  Kmer<MAX_K>::get_kmers(kmer_len, supermer.seq, kmers);
  kmers_and_exts.clear();
  for (int i = 1; i < (int)(supermer.seq.length() - kmer_len); i++) {
    Kmer<MAX_K> kmer = kmers[i];
    char left_ext = supermer.seq[i - 1];
    if (!quals[i - 1]) left_ext = '0';
    char right_ext = supermer.seq[i + kmer_len];
    if (!quals[i + kmer_len]) right_ext = '0';
    // get the lexicographically smallest
    Kmer<MAX_K> kmer_rc = kmer.revcomp();
    if (kmer_rc < kmer) {
      kmer = kmer_rc;
      swap(left_ext, right_ext);
      left_ext = comp_nucleotide(left_ext);
      right_ext = comp_nucleotide(right_ext);
    };
    kmers_and_exts.push_back({.kmer = kmer, .count = supermer.count, .left = left_ext, .right = right_ext});
  }
}

template <int MAX_K>
static void insert_supermer_from_read(Supermer &supermer, dist_object<KmerMapExts<MAX_K>> &kmers) {
  auto kmer_len = Kmer<MAX_K>::get_k();
  vector<KmerAndExt<MAX_K>> kmers_and_exts;
  kmers_and_exts.reserve(supermer.seq.length() - kmer_len);
  get_kmers_and_exts(supermer, kmers_and_exts);
  for (auto &kmer_and_ext : kmers_and_exts) {
    // find it - if it isn't found then insert it - this doen't set or change the value
    auto [exts_counts, is_new] = kmers->insert(kmer_and_ext.kmer, false);
    // no space - had to drop it
    if (!exts_counts) continue;
    int count = exts_counts->count + kmer_and_ext.count;
    if (count > numeric_limits<kmer_count_t>::max()) count = numeric_limits<kmer_count_t>::max();
    exts_counts->count = count;
    exts_counts->left_exts.inc(kmer_and_ext.left, kmer_and_ext.count);
    exts_counts->right_exts.inc(kmer_and_ext.right, kmer_and_ext.count);
  }
}

template <int MAX_K>
static void insert_supermer_from_ctg(Supermer &supermer, dist_object<KmerMapExts<MAX_K>> &kmers) {
  auto kmer_len = Kmer<MAX_K>::get_k();
  vector<KmerAndExt<MAX_K>> kmers_and_exts;
  kmers_and_exts.reserve(supermer.seq.length() - kmer_len);
  get_kmers_and_exts(supermer, kmers_and_exts);
  for (auto &kmer_and_ext : kmers_and_exts) {
    // insert a new kmer derived from the previous round's contigs
    auto [exts_counts, is_new] = kmers->insert(kmer_and_ext.kmer, true);
    // no space - had to drop it
    if (!exts_counts) continue;
    bool insert_it = false;
    if (is_new) {
      insert_it = true;
    } else if (!exts_counts->from_ctg) {
      // existing entry is from a read
      if (exts_counts->count == 1) {
        // singleton read kmer, replace - this will just be purged anyway
        insert_it = true;
      } else {
        char left_ext = exts_counts->get_left_ext();
        char right_ext = exts_counts->get_right_ext();
        // non-UU, replace
        if (left_ext == 'X' || left_ext == 'F' || right_ext == 'X' || right_ext == 'F') insert_it = true;
      }
    } else {
      // existing entry from contig
      if (exts_counts->count) {
        // will always insert, although it may get purged later for a conflict
        insert_it = true;
        char left_ext = exts_counts->get_left_ext();
        char right_ext = exts_counts->get_right_ext();
        if (left_ext != kmer_and_ext.left || right_ext != kmer_and_ext.right) {
          // if the two contig kmers disagree on extensions, set up to purge by setting the count to 0
          kmer_and_ext.count = 0;
        } else {
          // multiple occurrences of the same kmer derived from different contigs or parts of contigs
          // The only way this kmer could have been already found in the contigs only is if it came from a localassm
          // extension. In which case, all such kmers should not be counted again for each contig, because each
          // contig can use the same reads independently, and the depth will be oversampled.
          kmer_and_ext.count = min(kmer_and_ext.count, exts_counts->count);
        }
      }
    }
    if (insert_it) {
      *exts_counts = {.left_exts = {0}, .right_exts = {0}, .count = kmer_and_ext.count, .from_ctg = true};
      exts_counts->left_exts.inc(kmer_and_ext.left, kmer_and_ext.count);
      exts_counts->right_exts.inc(kmer_and_ext.right, kmer_and_ext.count);
    }
  }
}

template <int MAX_K>
struct HashTableInserter<MAX_K>::HashTableInserterState {
  bool using_ctg_kmers = false;
  dist_object<KmerMapExts<MAX_K>> kmers;
  upcxx_utils::BaseTimer insert_timer, kernel_timer;

  HashTableInserterState()
      : kmers(KmerMapExts<MAX_K>{}) {}
};

template <int MAX_K>
HashTableInserter<MAX_K>::HashTableInserter() {}

template <int MAX_K>
HashTableInserter<MAX_K>::~HashTableInserter() {
  if (state) delete state;
}

template <int MAX_K>
void HashTableInserter<MAX_K>::init(size_t max_elems, size_t max_ctg_elems, size_t num_errors, bool use_qf) {
  state = new HashTableInserterState();
  state->using_ctg_kmers = false;
  double free_mem = get_free_mem(true);
  SLOG_CPU_HT("There is ", get_size_str(free_mem), " free memory for hash table allocations\n");
  double avail_mem = free_mem / local_team().rank_n();

  size_t elem_size = sizeof(Kmer<MAX_K>) + sizeof(KmerExtsCounts);
  // expected size of compact hash table
  size_t compact_elem_size = sizeof(Kmer<MAX_K>) + sizeof(KmerCounts);
  double elem_size_ratio = (double)compact_elem_size / (double)elem_size;
  SLOG_CPU_HT("Element size for main HT ", elem_size, " and for compact HT ", compact_elem_size, " (ratio ", fixed, setprecision(3),
              elem_size_ratio, ")\n");

  double target_load_factor = 0.66;
  double load_multiplier = 1.0 / target_load_factor;
  size_t max_read_kmers = load_multiplier * (max_elems + max_ctg_elems + num_errors);
  size_t read_kmers_size = max_read_kmers * elem_size;
  size_t max_compact_kmers = load_multiplier * (max_elems + max_ctg_elems);
  size_t compact_kmers_size = max_compact_kmers * compact_elem_size;

  SLOG_CPU_HT("Element counts for a target load factor of ", fixed, setprecision(3), target_load_factor, ": read kmers ",
              max_read_kmers, ", compact ht ", max_compact_kmers, "\n");
  size_t tot_size = read_kmers_size + compact_kmers_size;
  SLOG_CPU_HT("Hash table sizes: read kmers ", read_kmers_size, ", compact ht ", compact_kmers_size, ", total ", tot_size, "\n");

  // keep some in reserve for all the various other requirements, e.g. rpcs, local allocs, etc
  double mem_ratio = (double)(0.8 * avail_mem) / tot_size;
  if (mem_ratio > 3.0) mem_ratio = 3.0;
  if (mem_ratio < 0.9)
    SWARN("Insufficent memory for ", fixed, setprecision(3), target_load_factor,
          " load factor across all data structures; reducing to ", (mem_ratio * target_load_factor),
          ". This may result in an OOM or dropped kmers\n");
  max_read_kmers *= mem_ratio;
  max_compact_kmers *= mem_ratio;
  SLOG_CPU_HT("Adjusted element counts by ", fixed, setprecision(3), mem_ratio, ": read kmers ", max_read_kmers, " compact ht ",
              max_compact_kmers, "\n");
  state->kmers->reserve(max_read_kmers);
  auto free_mem_after = get_free_mem(true);
  double used_mem = free_mem - free_mem_after;
  SLOG_CPU_HT("Memory available after hash table allocation: ", get_size_str(free_mem_after), ", used ", get_size_str(used_mem),
              "\n");
}

template <int MAX_K>
void HashTableInserter<MAX_K>::init_ctg_kmers(size_t max_elems) {
  state->using_ctg_kmers = true;
}

template <int MAX_K>
void HashTableInserter<MAX_K>::insert_supermer(const std::string &supermer_seq, kmer_count_t supermer_count) {
  Supermer supermer = {.seq = supermer_seq, .count = supermer_count};
  state->kernel_timer.start();
  for (int i = 0; i < supermer.seq.length(); i++) {
    char base = supermer.seq[i];
    if (base >= 'a' && base <= 'z') base += ('A' - 'a');
    if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N')
      DIE("bad char '", supermer.seq[i], "' in supermer seq int val ", (int)supermer.seq[i], " length ", supermer.seq.length(),
          " supermer ", supermer.seq);
  }
  if (!state->using_ctg_kmers)
    insert_supermer_from_read(supermer, state->kmers);
  else
    insert_supermer_from_ctg(supermer, state->kmers);
  state->kernel_timer.stop();
}

template <int MAX_K>
void HashTableInserter<MAX_K>::flush_inserts() {
  int64_t tot_num_kmers = reduce_one(state->kmers->size(), op_fast_add, 0).wait();
  SLOG_CPU_HT("After inserting ", state->using_ctg_kmers ? "ctg kmers" : "read kmers", " in hash table:\n");
  SLOG_CPU_HT("  Number of elements: ", tot_num_kmers, "\n");
  auto avg_load_factor = reduce_one(state->kmers->load_factor(), op_fast_add, 0).wait() / upcxx::rank_n();
  auto max_load_factor = reduce_one(state->kmers->load_factor(), op_fast_max, 0).wait();
  SLOG_CPU_HT("  load factor: ", avg_load_factor, " avg, ", max_load_factor, " max, load balance ", fixed, setprecision(3),
              (double)avg_load_factor / max_load_factor, "\n");
  auto avg_probe_len = (double)reduce_one(state->kmers->get_sum_probe_lens(), op_fast_add, 0).wait() / tot_num_kmers;
  auto max_probe_len = (double)reduce_one(state->kmers->get_max_probe_len(), op_fast_max, 0).wait();
  SLOG_CPU_HT("  probe lengths: ", avg_probe_len, " avg, ", max_probe_len, " max\n");
  auto tot_num_dropped = reduce_one(state->kmers->get_num_dropped(), op_fast_add, 0).wait();
  state->kmers->clear_num_dropped();
  auto tot_kmers = tot_num_kmers + tot_num_dropped;
  if (tot_num_dropped) SLOG_CPU_HT("  Number dropped ", perc_str(tot_num_dropped, tot_kmers), "\n");
  auto tot_num_overrides = reduce_one(state->kmers->get_num_singleton_overrides(), op_fast_add, 0).wait();
  if (tot_num_overrides) SLOG_CPU_HT("  Number singleton overrides ", perc_str(tot_num_overrides, tot_kmers), "\n");
  SLOG_VERBOSE("CPU kcount found total of ", tot_num_kmers + tot_num_dropped, " unique kmers including singletons and dropped\n");
  if (100.0 * tot_num_dropped / tot_kmers > 0.1)
    SWARN("Lack of memory caused ", perc_str(tot_num_dropped, tot_kmers), " kmers to be dropped (singleton overrides ",
          perc_str(tot_num_overrides, tot_kmers), ")\n");
  barrier();
  auto avg_kmers_processed = reduce_one(state->kmers->size(), op_fast_add, 0).wait() / rank_n();
  auto max_kmers_processed = reduce_one(state->kmers->size(), op_fast_max, 0).wait();
  SLOG_CPU_HT("  Avg kmers per rank ", avg_kmers_processed, " (balance ", (double)avg_kmers_processed / max_kmers_processed, ")\n");
}

template <int MAX_K>
double HashTableInserter<MAX_K>::insert_into_local_hashtable(dist_object<KmerMap<MAX_K>> &local_kmers) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_good_kmers = state->kmers->size();
  int64_t max_kmer_counts = 0;
  state->insert_timer.start();
  state->kmers->begin_iterate();
  while (true) {
    auto [kmer, kmer_ext_counts] = state->kmers->get_next();
    if (!kmer) break;
    if (kmer_ext_counts->count > max_kmer_counts) max_kmer_counts = kmer_ext_counts->count;
    if ((kmer_ext_counts->count < 2) || (kmer_ext_counts->left_exts.is_zero() && kmer_ext_counts->right_exts.is_zero()))
      num_good_kmers--;
  }
  auto fut_msm_max_kmer_count = upcxx_utils::min_sum_max_reduce_all(max_kmer_counts);
  LOG_MEM("Before inserting into local hashtable");
  SLOG_CPU_HT("Reserving compact hash table for ", num_good_kmers, " elements, requires ",
              get_size_str(num_good_kmers * (sizeof(Kmer<MAX_K>) + sizeof(KmerCounts))), "\n");
  auto free_mem = get_free_mem(true);
  local_kmers->reserve(num_good_kmers);
  auto free_mem_after = get_free_mem(true);
  double used_mem = free_mem - free_mem_after;
  SLOG_CPU_HT("Memory available after compact hash table allocation: ", get_size_str(free_mem_after), ", used ",
              get_size_str(used_mem), "\n");
  LOG_MEM("After reserving for local hashtable");
  int64_t num_purged = 0, num_inserted = 0;
  auto msm_max_kmer_count = fut_msm_max_kmer_count.wait();
  if (!rank_me()) LOG("High count (max) for kmers: ", msm_max_kmer_count.to_string(), "\n");
  int64_t high_count_threshold = msm_max_kmer_count.avg;
  state->kmers->begin_iterate();
  uint64_t sum_kmer_counts = 0;
  while (true) {
    auto [kmer, kmer_ext_counts] = state->kmers->get_next();
    if (!kmer) break;
    if (kmer_ext_counts->count < 2) {
      num_purged++;
      continue;
    }
    if (kmer_ext_counts->count >= high_count_threshold) {
      NET_LOG("High count kmer: k = ", Kmer<MAX_K>::get_k(), " count = ", kmer_ext_counts->count, " kmer = ", kmer->to_string(),
              "\n");
    }
    KmerCounts kmer_counts = {.uutig_frag = nullptr,
                              .count = kmer_ext_counts->count,
                              .left = kmer_ext_counts->get_left_ext(),
                              .right = kmer_ext_counts->get_right_ext()};
    if (kmer_counts.left == 'X' || kmer_counts.right == 'X' || kmer_counts.left == 'F' || kmer_counts.right == 'F') {
      // these are never used in the dbjg traversal, and are overwritten by ctg kmers if those have proper extensions
      num_purged++;
      continue;
    }
    const auto it = local_kmers->find(*kmer);
    if (it != local_kmers->end())
      WARN("Found a duplicate kmer ", kmer->to_string(), " - shouldn't happen: existing count ", it->second.count, " new count ",
           kmer_counts.count);
    local_kmers->insert({*kmer, kmer_counts});
    num_inserted++;
    sum_kmer_counts += kmer_counts.count;
  }
  state->insert_timer.stop();
  if (num_inserted > num_good_kmers) WARN("Inserted ", num_inserted, " but was expecting only ", num_good_kmers, "\n");
  barrier();
  LOG_MEM("After inserting into local hashtable");
  auto tot_num_purged = reduce_one(num_purged, op_fast_add, 0).wait();
  auto tot_num_kmers = reduce_one(state->kmers->size(), op_fast_add, 0).wait();
  SLOG_CPU_HT("CPU Hashtable: purged ", perc_str(tot_num_purged, tot_num_kmers), " singleton kmers out of ", tot_num_kmers, "\n");
  auto final_load_factor = (double)reduce_one(local_kmers->load_factor(), op_fast_add, 0).wait() / rank_n();
  auto max_final_load_factor = (double)reduce_one(local_kmers->load_factor(), op_fast_max, 0).wait();
  auto all_kmers_inserted = reduce_all(num_inserted, op_fast_add).wait();
  SLOG_CPU_HT("Compact hash table has ", num_inserted, " kmers and load factor ", fixed, setprecision(3), final_load_factor,
              " avg ", max_final_load_factor, " max\n");
  auto all_sum_kmer_counts = reduce_all(sum_kmer_counts, op_fast_add).wait();
  double avg_kmer_count = (double)all_sum_kmer_counts / all_kmers_inserted;
  SLOG_CPU_HT("For ", all_kmers_inserted, " kmers, average kmer count (depth): ", fixed, setprecision(2), avg_kmer_count, "\n");

  SLOG_CPU_HT("Total kmer count sum: ", all_sum_kmer_counts, "\n");

  return avg_kmer_count;
}

template <int MAX_K>
void HashTableInserter<MAX_K>::get_elapsed_time(double &insert_time, double &kernel_time) {
  insert_time = state->insert_timer.get_elapsed();
  kernel_time = state->kernel_timer.get_elapsed();
}

#define SEQ_BLOCK_INSERTER_K(KMER_LEN) template struct SeqBlockInserter<KMER_LEN>;
#define HASH_TABLE_INSERTER_K(KMER_LEN) template class HashTableInserter<KMER_LEN>;

SEQ_BLOCK_INSERTER_K(32);
HASH_TABLE_INSERTER_K(32);
#if MAX_BUILD_KMER >= 64
SEQ_BLOCK_INSERTER_K(64);
HASH_TABLE_INSERTER_K(64);
#endif
#if MAX_BUILD_KMER >= 96
SEQ_BLOCK_INSERTER_K(96);
HASH_TABLE_INSERTER_K(96);
#endif
#if MAX_BUILD_KMER >= 128
SEQ_BLOCK_INSERTER_K(128);
HASH_TABLE_INSERTER_K(128);
#endif
#if MAX_BUILD_KMER >= 160
SEQ_BLOCK_INSERTER_K(160);
HASH_TABLE_INSERTER_K(160);
#endif
#undef SEQ_BLOCK_INSERTER_K
#undef HASH_TABLE_INSERTER_K
