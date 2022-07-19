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

#include <stdarg.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "zstr.hpp"

#include "stage_timers.hpp"
#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

//#define DBG_INS_CTG_KMER DBG
#define DBG_INS_CTG_KMER(...)
//#define DBG_INSERT_KMER DBG
#define DBG_INSERT_KMER(...)

void Supermer::pack(const string &unpacked_seq) {
  // each position in the sequence is an upper or lower case nucleotide, not including Ns
  seq = string(unpacked_seq.length() / 2 + unpacked_seq.length() % 2, 0);
  for (int i = 0; i < unpacked_seq.length(); i++) {
    char packed_val = 0;
    switch (unpacked_seq[i]) {
      case 'a': packed_val = 1; break;
      case 'c': packed_val = 2; break;
      case 'g': packed_val = 3; break;
      case 't': packed_val = 4; break;
      case 'A': packed_val = 5; break;
      case 'C': packed_val = 6; break;
      case 'G': packed_val = 7; break;
      case 'T': packed_val = 8; break;
      case 'N': packed_val = 9; break;
      default: DIE("Invalid value encountered when packing '", unpacked_seq[i], "' ", (int)unpacked_seq[i]);
    };
    seq[i / 2] |= (!(i % 2) ? (packed_val << 4) : packed_val);
    if (seq[i / 2] == '_') DIE("packed byte is same as sentinel _");
  }
}

void Supermer::unpack() {
  static const char to_base[] = {0, 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T', 'N'};
  string unpacked_seq;
  for (int i = 0; i < seq.length(); i++) {
    unpacked_seq += to_base[(seq[i] & 240) >> 4];
    int right_ext = seq[i] & 15;
    if (right_ext) unpacked_seq += to_base[right_ext];
  }
  seq = unpacked_seq;
}

int Supermer::get_bytes() { return seq.length() + sizeof(kmer_count_t); }

template <int MAX_K>
KmerDHT<MAX_K>::KmerDHT(uint64_t my_num_kmers, size_t max_kmer_store_bytes, int max_rpcs_in_flight, bool useHHSS, bool use_qf)
    : local_kmers({})
    , ht_inserter({})
    , kmer_store()
    , max_kmer_store_bytes(max_kmer_store_bytes)
    , my_num_kmers(my_num_kmers)
    , max_rpcs_in_flight(max_rpcs_in_flight)
    , num_supermer_inserts(0) {
  // minimizer len depends on k
  minimizer_len = Kmer<MAX_K>::get_k() * 2 / 3 + 1;
  if (minimizer_len < 15) minimizer_len = 15;
  if (minimizer_len > 27) minimizer_len = 27;
  SLOG_VERBOSE("Using a minimizer length of ", minimizer_len, "\n");
  // main purpose of the timer here is to track memory usage
  BarrierTimer timer(__FILEFUNC__);
  auto node0_cores = upcxx::local_team().rank_n();
  // check if we have enough memory to run - conservative because we don't want to run out of memory
  double adjustment_factor = 0.2;
  auto my_adjusted_num_kmers = my_num_kmers * adjustment_factor;
  double required_space = estimate_hashtable_memory(my_adjusted_num_kmers, sizeof(Kmer<MAX_K>) + sizeof(KmerCounts)) * node0_cores;
  auto max_reqd_space = upcxx::reduce_all(required_space, upcxx::op_fast_max).wait();
  auto free_mem = get_free_mem();
  auto lowest_free_mem = upcxx::reduce_all(free_mem, upcxx::op_fast_min).wait();
  auto highest_free_mem = upcxx::reduce_all(free_mem, upcxx::op_fast_max).wait();
  SLOG_VERBOSE("With adjustment factor of ", adjustment_factor, " require ", get_size_str(max_reqd_space), " per node (",
               my_adjusted_num_kmers, " kmers per rank), and there is ", get_size_str(lowest_free_mem), " to ",
               get_size_str(highest_free_mem), " available on the nodes\n");
  if (lowest_free_mem * 0.80 < max_reqd_space)
    SWARN("Insufficient memory available: this could crash with OOM (lowest=", get_size_str(lowest_free_mem),
          " vs reqd=", get_size_str(max_reqd_space), ")");

  kmer_store.set_size("kmers", max_kmer_store_bytes, max_rpcs_in_flight, useHHSS);

  barrier();
  // in this case we have to roughly estimate the hash table size because the space is reserved now
  // err on the side of excess because the whole point of doing this is speed and we don't want a
  // hash table resize
  // Unfortunately, this estimate depends on the depth of the sample - high depth means more wasted memory,
  // but low depth means potentially resizing the hash table, which is very expensive
  double kmers_space_reserved = my_adjusted_num_kmers * (sizeof(Kmer<MAX_K>) + sizeof(KmerCounts));
  SLOG_VERBOSE("Reserving at least ", get_size_str(node0_cores * kmers_space_reserved), " for kmer hash tables with ",
               node0_cores * my_adjusted_num_kmers, " entries on node 0\n");
  double init_free_mem = get_free_mem();
  if (my_adjusted_num_kmers <= 0) DIE("no kmers to reserve space for");
  kmer_store.set_update_func(
      [&ht_inserter = this->ht_inserter, &num_supermer_inserts = this->num_supermer_inserts](Supermer supermer) {
        num_supermer_inserts++;
        ht_inserter->insert_supermer(supermer.seq, supermer.count);
      });
  // this is conservative, actually varies from around 1/2 to 1/5
  if (use_qf) my_num_kmers *= 0.6;
  ht_inserter->init(my_num_kmers, use_qf);
  barrier();
}

template <int MAX_K>
void KmerDHT<MAX_K>::clear_stores() {
  kmer_store.clear();
}

template <int MAX_K>
KmerDHT<MAX_K>::~KmerDHT() {
  local_kmers->clear();
  KmerMap<MAX_K>().swap(*local_kmers);
  clear_stores();
}

template <int MAX_K>
void KmerDHT<MAX_K>::init_ctg_kmers(int64_t max_elems) {
  using_ctg_kmers = true;
  ht_inserter->init_ctg_kmers(max_elems);
}

template <int MAX_K>
int KmerDHT<MAX_K>::get_minimizer_len() {
  return minimizer_len;
}

template <int MAX_K>
uint64_t KmerDHT<MAX_K>::get_num_kmers(bool all) {
  if (!all)
    return reduce_one((uint64_t)local_kmers->size(), op_fast_add, 0).wait();
  else
    return reduce_all((uint64_t)local_kmers->size(), op_fast_add).wait();
}

template <int MAX_K>
int64_t KmerDHT<MAX_K>::get_local_num_kmers(void) {
  return local_kmers->size();
}

template <int MAX_K>
upcxx::intrank_t KmerDHT<MAX_K>::get_kmer_target_rank(const Kmer<MAX_K> &kmer, const Kmer<MAX_K> *kmer_rc) const {
  assert(&kmer != kmer_rc && "Can be a palindrome, cannot be the same Kmer instance");
  return kmer.minimizer_hash_fast(minimizer_len, kmer_rc) % rank_n();
}

template <int MAX_K>
KmerCounts *KmerDHT<MAX_K>::get_local_kmer_counts(Kmer<MAX_K> &kmer) {
  const auto it = local_kmers->find(kmer);
  if (it == local_kmers->end()) return nullptr;
  return &it->second;
}

template <int MAX_K>
int64_t KmerDHT<MAX_K>::get_num_supermer_inserts() {
  return num_supermer_inserts;
}

template <int MAX_K>
bool KmerDHT<MAX_K>::kmer_exists(Kmer<MAX_K> kmer_fw) {
  const Kmer<MAX_K> kmer_rc = kmer_fw.revcomp();
  const Kmer<MAX_K> *kmer = (kmer_rc < kmer_fw) ? &kmer_rc : &kmer_fw;

  return rpc(
             get_kmer_target_rank(kmer_fw, &kmer_rc),
             [](Kmer<MAX_K> kmer, dist_object<KmerMap<MAX_K>> &local_kmers) -> bool {
               const auto it = local_kmers->find(kmer);
               if (it == local_kmers->end()) return false;
               return true;
             },
             *kmer, local_kmers)
      .wait();
}

template <int MAX_K>
void KmerDHT<MAX_K>::add_supermer(Supermer &supermer, int target_rank) {
  kmer_store.update(target_rank, supermer);
}

template <int MAX_K>
void KmerDHT<MAX_K>::flush_updates() {
  BarrierTimer timer(__FILEFUNC__);
  kmer_store.flush_updates();
  barrier();
  ht_inserter->flush_inserts();
}

template <int MAX_K>
void KmerDHT<MAX_K>::finish_updates() {
  ht_inserter->insert_into_local_hashtable(local_kmers);
  double insert_time, kernel_time;
  ht_inserter->get_elapsed_time(insert_time, kernel_time);
  stage_timers.kernel_kmer_analysis->inc_elapsed(kernel_time);
  LOG("Total time for local hashtables: insert_time=", insert_time, "kernel_time=", kernel_time, "\n");
}

// one line per kmer, format:
// KMERCHARS LR N
// where L is left extension and R is right extension, one char, either X, F or A, C, G, T
// where N is the count of the kmer frequency
template <int MAX_K>
void KmerDHT<MAX_K>::dump_kmers() {
  BarrierTimer timer(__FILEFUNC__);
  int k = Kmer<MAX_K>::get_k();
  string dump_fname = "kmers-" + to_string(k) + ".txt.gz";
  get_rank_path(dump_fname, rank_me());
  zstr::ofstream dump_file(dump_fname);
  ostringstream out_buf;
  ProgressBar progbar(local_kmers->size(), "Dumping kmers to " + dump_fname);
  int64_t i = 0;
  for (auto &elem : *local_kmers) {
    out_buf << elem.first << " " << elem.second.count << " " << elem.second.left << " " << elem.second.right;
    out_buf << endl;
    i++;
    if (!(i % 1000)) {
      dump_file << out_buf.str();
      out_buf = ostringstream();
    }
    progbar.update();
  }
  if (!out_buf.str().empty()) dump_file << out_buf.str();
  dump_file.close();
  progbar.done();
  SLOG_VERBOSE("Dumped ", this->get_num_kmers(), " kmers\n");
}

template <int MAX_K>
typename KmerMap<MAX_K>::iterator KmerDHT<MAX_K>::local_kmers_begin() {
  return local_kmers->begin();
}

template <int MAX_K>
typename KmerMap<MAX_K>::iterator KmerDHT<MAX_K>::local_kmers_end() {
  return local_kmers->end();
}

template <int MAX_K>
int32_t KmerDHT<MAX_K>::get_time_offset_us() {
  std::chrono::duration<double> t_elapsed = CLOCK_NOW() - start_t;
  return std::chrono::duration_cast<std::chrono::microseconds>(t_elapsed).count();
}

#define KMER_DHT_K(KMER_LEN) template class KmerDHT<KMER_LEN>

KMER_DHT_K(32);
#if MAX_BUILD_KMER >= 64
KMER_DHT_K(64);
#endif
#if MAX_BUILD_KMER >= 96
KMER_DHT_K(96);
#endif
#if MAX_BUILD_KMER >= 128
KMER_DHT_K(128);
#endif
#if MAX_BUILD_KMER >= 160
KMER_DHT_K(160);
#endif

#undef KMER_DHT_K
