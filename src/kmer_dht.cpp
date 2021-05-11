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

#include "kmer_dht.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

//#define DBG_INS_CTG_KMER DBG
#define DBG_INS_CTG_KMER(...)
//#define DBG_INSERT_KMER DBG
#define DBG_INSERT_KMER(...)

// global variables to avoid passing dist objs to rpcs
static int64_t _num_kmers_counted = 0;
static int64_t _num_kmers_counted_locally = 0;

template <int MAX_K>
auto KmerDHT<MAX_K>::update_bloom_set_func() {
  return [&bloom_filter1 = this->bloom_filter1, &bloom_filter2 = this->bloom_filter2](const Kmer<MAX_K> &kmer) {
    // look for it in the first bloom filter - if not found, add it just to the first bloom filter
    // if found, add it to the second bloom filter
    if (!bloom_filter1->possibly_contains(kmer.get_bytes()))
      bloom_filter1->add(kmer.get_bytes());
    else
      bloom_filter2->add(kmer.get_bytes());
  };
}

template <int MAX_K>
auto KmerDHT<MAX_K>::update_ctg_bloom_set_func() {
  return [&bloom_filter2 = this->bloom_filter2](Kmer<MAX_K> kmer) {
    // only add to bloom_filter2
    bloom_filter2->add(kmer.get_bytes());
  };
}

template <int MAX_K>
auto KmerDHT<MAX_K>::update_bloom_count_func() {
  // second pass, actually count if in second bloom filter
  return [&bloom_filter = this->bloom_filter2, &kmers = this->kmers](KmerAndExt kmer_and_ext) {
    if (!bloom_filter->possibly_contains(kmer_and_ext.kmer.get_bytes())) return;
    update_count(kmer_and_ext, kmers);
  };
}

template <int MAX_K>
auto KmerDHT<MAX_K>::update_count_func() {
  return [&kmers = this->kmers](KmerAndExt kmer_and_ext) { update_count(kmer_and_ext, kmers); };
}

template <int MAX_K>
void KmerDHT<MAX_K>::update_count(KmerAndExt kmer_and_ext, dist_object<KmerMap> &kmers) {
#if defined(ENABLE_GPUS) & defined(KCOUNT_ENABLE_GPUS)
  // FIXME: buffer hash table entries, and when full, copy across to
  // gpu_driver->insert_kmer(kmer_and_ext.kmer.get_longs(), kmer_and_ext.count, kmer_and_ext.left, kmer_and_ext.right);
#endif

  //#else
  // find it - if it isn't found then insert it, otherwise increment the counts
  const auto it = kmers->find(kmer_and_ext.kmer);
  if (it == kmers->end()) {
    KmerCounts kmer_counts = {.left_exts = {0},
                              .right_exts = {0},
                              .uutig_frag = nullptr,
                              .count = kmer_and_ext.count,
                              .left = 'X',
                              .right = 'X',
                              .from_ctg = false};
    kmer_counts.left_exts.inc(kmer_and_ext.left, kmer_and_ext.count);
    kmer_counts.right_exts.inc(kmer_and_ext.right, kmer_and_ext.count);
    auto prev_bucket_count = kmers->bucket_count();
    kmers->insert({kmer_and_ext.kmer, kmer_counts});
    // since sizes are an estimate this could happen, but it will impact performance
    if (prev_bucket_count < kmers->bucket_count())
      SWARN("Hash table on rank 0 was resized from ", prev_bucket_count, " to ", kmers->bucket_count());
    DBG_INSERT_KMER("inserted kmer ", kmer_and_ext.kmer.to_string(), " with count ", kmer_counts.count, "\n");
  } else {
    auto kmer_count = &it->second;
    int count = kmer_count->count + kmer_and_ext.count;
    if (count > numeric_limits<kmer_count_t>::max()) count = numeric_limits<kmer_count_t>::max();
    kmer_count->count = count;
    kmer_count->left_exts.inc(kmer_and_ext.left, kmer_and_ext.count);
    kmer_count->right_exts.inc(kmer_and_ext.right, kmer_and_ext.count);
  }
  //#endif
}

template <int MAX_K>
auto KmerDHT<MAX_K>::update_ctg_kmers_count_func() {
  return [&kmers = this->kmers](KmerAndExt kmer_and_ext) {
    // insert a new kmer derived from the previous round's contigs
    const auto it = kmers->find(kmer_and_ext.kmer);
    bool insert = false;
    if (it == kmers->end()) {
      // if it isn't found then insert it
      insert = true;
    } else {
      auto kmer_counts = &it->second;
      if (!kmer_counts->from_ctg) {
        // existing kmer is from a read, only replace with new contig kmer if the existing kmer is not UU
        char left_ext = kmer_counts->get_left_ext();
        char right_ext = kmer_counts->get_right_ext();
        if (left_ext == 'X' || left_ext == 'F' || right_ext == 'X' || right_ext == 'F') {
          // non-UU, replace
          insert = true;
          // but keep the count from the read kmer
          // or could sum the depths
          DBG_INS_CTG_KMER("replace non-UU read kmer\n");
        }
      } else {
        // existing kmer from previous round's contigs
        // update kmer counts
        if (!kmer_counts->count) {
          // previously must have been a conflict and set to zero, so don't do anything
          DBG_INS_CTG_KMER("skip conflicted kmer, depth 0\n");
        } else {
          // will always insert, although it may get purged later for a conflict
          insert = true;
          char left_ext = kmer_counts->get_left_ext();
          char right_ext = kmer_counts->get_right_ext();
          if (left_ext != kmer_and_ext.left || right_ext != kmer_and_ext.right) {
            // if the two contig kmers disagree on extensions, set up to purge by setting the count to 0
            kmer_and_ext.count = 0;
            DBG_INS_CTG_KMER("set to purge conflict: prev ", left_ext, ", ", right_ext, " new ", kmer_and_ext.left, ", ",
                             kmer_and_ext.right, "\n");
          } else {
            // multiple occurrences of the same kmer derived from different contigs or parts of contigs
            // The only way this kmer could have been already found in the contigs only is if it came from a localassm
            // extension. In which case, all such kmers should not be counted again for each contig, because each
            // contig can use the same reads independently, and the depth will be oversampled.
            kmer_and_ext.count = min(kmer_and_ext.count, kmer_counts->count);
            DBG_INS_CTG_KMER("increase count of existing ctg kmer from ", kmer_counts->count, " to ", kmer_and_ext.count, "\n");
          }
        }
      }
    }
    if (insert) {
      kmer_count_t count = kmer_and_ext.count;
      KmerCounts kmer_counts = {
          .left_exts = {0}, .right_exts = {0}, .uutig_frag = nullptr, .count = count, .left = 'X', .right = 'X', .from_ctg = true};
      kmer_counts.left_exts.inc(kmer_and_ext.left, count);
      kmer_counts.right_exts.inc(kmer_and_ext.right, count);
      (*kmers)[kmer_and_ext.kmer] = kmer_counts;
    }
  };
}

template <int MAX_K>
KmerDHT<MAX_K>::KmerDHT(uint64_t my_num_kmers, int max_kmer_store_bytes, int max_rpcs_in_flight, bool force_bloom, bool useHHSS)
    : kmers({})
    , bloom_filter1({})
    , bloom_filter2({})
    , kmer_store_bloom()
    , kmer_store()
    , max_kmer_store_bytes(max_kmer_store_bytes)
    , initial_kmer_dht_reservation(0)
    , my_num_kmers(my_num_kmers)
    , max_rpcs_in_flight(max_rpcs_in_flight)
    , bloom1_cardinality(0)
    , estimated_error_rate(0.0) {
  // minimizer len depends on k
  minimizer_len = Kmer<MAX_K>::get_k() * 2 / 3 + 1;
  if (minimizer_len < 15) minimizer_len = 15;
  if (minimizer_len > 27) minimizer_len = 27;
  SLOG_VERBOSE("Using a minimizer length of ", minimizer_len, "\n");
  // main purpose of the timer here is to track memory usage
  BarrierTimer timer(__FILEFUNC__);
  auto node0_cores = upcxx::local_team().rank_n();
  // check if we have enough memory to run without bloom - require 2x the estimate for non-bloom - conservative because we don't
  // want to run out of memory
  double adjustment_factor = KCOUNT_NO_BLOOM_ADJUSTMENT_FACTOR;
  auto my_adjusted_num_kmers = my_num_kmers * adjustment_factor;
  double required_space = estimate_hashtable_memory(my_adjusted_num_kmers, sizeof(Kmer<MAX_K>) + sizeof(KmerCounts)) * node0_cores;
  auto max_reqd_space = upcxx::reduce_all(required_space, upcxx::op_fast_max).wait();
  auto free_mem = get_free_mem();
  auto lowest_free_mem = upcxx::reduce_all(free_mem, upcxx::op_fast_min).wait();
  auto highest_free_mem = upcxx::reduce_all(free_mem, upcxx::op_fast_max).wait();
  SLOG_VERBOSE("Without bloom filters and adjustment factor of ", adjustment_factor, " require ", get_size_str(max_reqd_space),
               " per node (", my_adjusted_num_kmers, " kmers per rank), and there is ", get_size_str(lowest_free_mem), " to ",
               get_size_str(highest_free_mem), " available on the nodes\n");
  if (force_bloom) {
    use_bloom = true;
    SLOG_VERBOSE("Using bloom (--force-bloom set)\n");
  } else {
    if (lowest_free_mem * 0.80 < max_reqd_space) {
      use_bloom = true;
      SLOG("Insufficient memory available: enabling bloom filters assuming 80% of free mem is available for hashtables\n");
    } else {
      use_bloom = false;
      SLOG_VERBOSE("Sufficient memory available; not using bloom filters\n");
    }
  }

  if (use_bloom)
    kmer_store_bloom.set_size("bloom", max_kmer_store_bytes, max_rpcs_in_flight, useHHSS);
  else
    kmer_store.set_size("kmers", max_kmer_store_bytes, max_rpcs_in_flight, useHHSS);

  if (use_bloom) {
#if defined(ENABLE_GPUS) & defined(KCOUNT_ENABLE_GPUS)
    SDIE("Cannot use bloom with GPU kmer counting");
#endif

    // in this case we get an accurate estimate of the hash table size after the first bloom round, so the hash table space
    // is reserved then
    double init_mem_free = get_free_mem();
    bloom_filter1->init(my_num_kmers, KCOUNT_BLOOM_FP);
    // second bloom has far fewer kmers
    bloom_filter2->init(my_num_kmers * adjustment_factor, KCOUNT_BLOOM_FP);
    SLOG_VERBOSE("Bloom filters used ", get_size_str(init_mem_free - get_free_mem()), " memory on node 0\n");
  } else {
    barrier();
    // in this case we have to roughly estimate the hash table size because the space is reserved now
    // err on the side of excess because the whole point of doing this is speed and we don't want a
    // hash table resize
    // Unfortunately, this estimate depends on the depth of the sample - high depth means more wasted memory,
    // but low depth means potentially resizing the hash table, which is very expensive
    initial_kmer_dht_reservation = my_adjusted_num_kmers;
    double kmers_space_reserved = my_adjusted_num_kmers * (sizeof(Kmer<MAX_K>) + sizeof(KmerCounts));
    SLOG_VERBOSE("Reserving at least ", get_size_str(node0_cores * kmers_space_reserved), " for kmer hash tables with ",
                 node0_cores * my_adjusted_num_kmers, " entries on node 0\n");
    double init_free_mem = get_free_mem();
    if (my_adjusted_num_kmers <= 0) DIE("no kmers to reserve space for");
    kmers->reserve(my_adjusted_num_kmers);
#if defined(ENABLE_GPUS) & defined(KCOUNT_ENABLE_GPUS)
    // vector<GPUKmerCounts> gpu_kmer_counts;
    // gpu_kmer_counts.reserve(KCOUNT_GPU_MAX_HT_ENTRIES);
    if (gpu_utils::get_num_node_gpus() <= 0) {
      DIE("GPUs are enabled but no GPU could be configured for kmer counting");
    } else {
      // calculate total slots for hash table. Reserve space for parse and pack
      int bytes_for_pnp = KCOUNT_GPU_SEQ_BLOCK_SIZE * (2 + Kmer<MAX_K>::get_N_LONGS() * sizeof(uint64_t) + sizeof(int));
      int max_dev_id = reduce_one(gpu_utils::get_gpu_device_pci_id(), op_fast_max, 0).wait();
      auto gpu_avail_mem = gpu_utils::get_free_gpu_mem() * max_dev_id / upcxx::local_team().rank_n() - bytes_for_pnp;
      auto gpu_tot_mem = gpu_utils::get_tot_gpu_mem() * max_dev_id / upcxx::local_team().rank_n() - bytes_for_pnp;
      SLOG_VERBOSE(KLMAGENTA, "Available GPU memory per rank for kcount hash table is ", get_size_str(gpu_avail_mem),
                   " out of a max of ", get_size_str(gpu_tot_mem), KNORM, "\n");
      // don't use up all the memory
      // gpu_avail_mem *= 0.9;
      double init_time;
      gpu_driver = new kcount_gpu::HashTableGPUDriver(rank_me(), rank_n(), Kmer<MAX_K>::get_k(), Kmer<MAX_K>::get_N_LONGS(),
                                                      gpu_avail_mem, init_time);
      SLOG_VERBOSE(KLMAGENTA, "Initialized hash table GPU driver in ", std::fixed, std::setprecision(3), init_time, " s", KNORM,
                   "\n");
      auto num_ht_slots = gpu_driver->get_num_ht_slots();
      SLOG_VERBOSE(KLMAGENTA, "GPU hash table has ", num_ht_slots, " slots/rank", KNORM, "\n");
      auto gpu_free_mem = gpu_utils::get_free_gpu_mem() * max_dev_id / upcxx::local_team().rank_n();
      SLOG_VERBOSE(KLMAGENTA, "After initializing GPU hash table, there is ", get_size_str(gpu_free_mem),
                   " memory available per rank, with ", get_size_str(bytes_for_pnp), " reserved for parse and pack" KNORM, "\n");
    }
#endif
    barrier();
  }
  start_t = CLOCK_NOW();
}

template <int MAX_K>
void KmerDHT<MAX_K>::clear() {
  kmers->clear();
  KmerMap().swap(*kmers);
  clear_stores();
}

template <int MAX_K>
void KmerDHT<MAX_K>::clear_stores() {
  kmer_store_bloom.clear();
  kmer_store.clear();
}

template <int MAX_K>
KmerDHT<MAX_K>::~KmerDHT() {
  clear();
#if defined(ENABLE_GPUS) & defined(KCOUNT_ENABLE_GPUS)
  delete gpu_driver;
#endif
}

template <int MAX_K>
void KmerDHT<MAX_K>::reserve_space_and_clear_bloom1() {
  BarrierTimer timer(__FILEFUNC__);
  // at this point we're done with generating the bloom filters, so we can drop the first bloom filter and
  // allocate the kmer hash table

  // purge the kmer store and prep the kmer + count
  kmer_store_bloom.clear();
  kmer_store.set_size("kmers", max_kmer_store_bytes, max_rpcs_in_flight);

  int64_t cardinality1 = bloom_filter1->estimate_num_items();
  int64_t cardinality2 = bloom_filter2->estimate_num_items();
  bloom1_cardinality = cardinality1;
  SLOG_VERBOSE("Rank 0: first bloom filter size estimate is ", cardinality1, " and second size estimate is ", cardinality2,
               " ratio is ", (double)cardinality2 / cardinality1, "\n");
  bloom_filter1->clear();  // no longer need it

  barrier();
  // two bloom false positive rates applied
  initial_kmer_dht_reservation = (int64_t)(cardinality2 * (1 + KCOUNT_BLOOM_FP) * (1 + KCOUNT_BLOOM_FP) + 20000);
  auto node0_cores = upcxx::local_team().rank_n();
  double kmers_space_reserved = initial_kmer_dht_reservation * (sizeof(Kmer<MAX_K>) + sizeof(KmerCounts));
  SLOG_VERBOSE("Reserving at least ", get_size_str(node0_cores * kmers_space_reserved), " for kmer hash tables with ",
               node0_cores * initial_kmer_dht_reservation, " entries on node 0\n");
  barrier(local_team());
  double init_free_mem = get_free_mem();
  barrier(local_team());
  kmers->reserve(initial_kmer_dht_reservation);
  barrier(local_team());
  SLOG_VERBOSE("Actually used ", get_size_str(init_free_mem - get_free_mem()), "\n");
  barrier();
}

template <int MAX_K>
bool KmerDHT<MAX_K>::get_use_bloom() {
  return use_bloom;
}

template <int MAX_K>
pair<int64_t, int64_t> KmerDHT<MAX_K>::get_bytes_sent() {
  auto all_bytes_sent = reduce_one(bytes_sent, op_fast_add, 0).wait();
  auto max_bytes_sent = reduce_one(bytes_sent, op_fast_max, 0).wait();
  return {all_bytes_sent, max_bytes_sent};
}

template <int MAX_K>
void KmerDHT<MAX_K>::set_pass(PASS_TYPE pass_type) {
  _num_kmers_counted = 0;
  _num_kmers_counted_locally = 0;
  this->pass_type = pass_type;
  switch (pass_type) {
    case BLOOM_SET_PASS: kmer_store_bloom.set_update_func(update_bloom_set_func()); break;
    case CTG_BLOOM_SET_PASS: kmer_store_bloom.set_update_func(update_ctg_bloom_set_func()); break;
    case BLOOM_COUNT_PASS: kmer_store.set_update_func(update_bloom_count_func()); break;
    case NO_BLOOM_PASS: kmer_store.set_update_func(update_count_func()); break;
    case CTG_KMERS_PASS: kmer_store.set_update_func(update_ctg_kmers_count_func()); break;
  };
}

template <int MAX_K>
int KmerDHT<MAX_K>::get_minimizer_len() {
  return minimizer_len;
}

template <int MAX_K>
int64_t KmerDHT<MAX_K>::get_num_kmers(bool all) {
  if (!all)
    return reduce_one(kmers->size(), op_fast_add, 0).wait();
  else
    return reduce_all(kmers->size(), op_fast_add).wait();
}

template <int MAX_K>
float KmerDHT<MAX_K>::max_load_factor() {
  return reduce_one(kmers->max_load_factor(), op_fast_max, 0).wait();
}

template <int MAX_K>
void KmerDHT<MAX_K>::print_load_factor() {
  int64_t num_kmers_est = initial_kmer_dht_reservation * rank_n();
  int64_t num_kmers = get_num_kmers();
  SLOG_VERBOSE("Originally reserved ", num_kmers_est, " and now have ", num_kmers, " elements\n");
  auto avg_load_factor = reduce_one(kmers->load_factor(), op_fast_add, 0).wait() / upcxx::rank_n();
  SLOG_VERBOSE("kmer DHT load factor: ", avg_load_factor, "\n");
}

template <int MAX_K>
int64_t KmerDHT<MAX_K>::get_local_num_kmers(void) {
  return kmers->size();
}

template <int MAX_K>
double KmerDHT<MAX_K>::get_estimated_error_rate() {
  return estimated_error_rate;
}

template <int MAX_K>
upcxx::intrank_t KmerDHT<MAX_K>::get_kmer_target_rank(const Kmer<MAX_K> &kmer, const Kmer<MAX_K> *kmer_rc) const {
  assert(&kmer != kmer_rc && "Can be a palindrome, cannot be the same Kmer instance");
  return kmer.minimizer_hash_fast(minimizer_len, kmer_rc) % rank_n();
}

template <int MAX_K>
KmerCounts *KmerDHT<MAX_K>::get_local_kmer_counts(Kmer<MAX_K> &kmer) {
  const auto it = kmers->find(kmer);
  if (it == kmers->end()) return nullptr;
  return &it->second;
}

#ifdef DEBUG
template <int MAX_K>
bool KmerDHT<MAX_K>::kmer_exists(Kmer<MAX_K> kmer_fw) {
  const Kmer<MAX_K> kmer_rc = kmer_fw.revcomp();
  const Kmer<MAX_K> *kmer = (kmer_rc < kmer_fw) ? &kmer_rc : &kmer_fw;

  return rpc(get_kmer_target_rank(kmer_fw, &kmer_rc),
             [](Kmer<MAX_K> kmer, dist_object<KmerMap> &kmers) -> bool {
               const auto it = kmers->find(kmer);
               if (it == kmers->end()) return false;
               return true;
             },
             *kmer, kmers)
      .wait();
}
#endif

template <int MAX_K>
void KmerDHT<MAX_K>::add_kmer(Kmer<MAX_K> kmer, char left_ext, char right_ext, kmer_count_t count, int target_rank) {
  _num_kmers_counted++;
  if (!count) count = 1;
  if (target_rank == -1) {
    // get the lexicographically smallest
    Kmer<MAX_K> kmer_rc = kmer.revcomp();
    if (kmer_rc < kmer) {
      kmer.swap(kmer_rc);
      swap(left_ext, right_ext);
      left_ext = comp_nucleotide(left_ext);
      right_ext = comp_nucleotide(right_ext);
    }
    target_rank = get_kmer_target_rank(kmer, &kmer_rc);
  }
  if (pass_type == BLOOM_SET_PASS || pass_type == CTG_BLOOM_SET_PASS) {
    if (count) {
      kmer_store_bloom.update(target_rank, kmer);
      bytes_sent += sizeof(kmer);
    }
  } else {
    KmerAndExt kmer_and_ext = {.kmer = kmer, .count = count, .left = left_ext, .right = right_ext};
    kmer_store.update(target_rank, kmer_and_ext);
    bytes_sent += sizeof(kmer_and_ext);
  }
}

template <int MAX_K>
void KmerDHT<MAX_K>::flush_updates() {
  BarrierTimer timer(__FILEFUNC__);
  if (pass_type == BLOOM_SET_PASS || pass_type == CTG_BLOOM_SET_PASS)
    kmer_store_bloom.flush_updates();
  else
    kmer_store.flush_updates();
#if defined(ENABLE_GPUS) & defined(KCOUNT_ENABLE_GPUS)
    // FIXME: call gpus to process any outstanding buffered local updates
    // FIXME: copy from gpu to cpu and insert in cpu unordered_map
#endif
  auto avg_kmers_processed = reduce_one(_num_kmers_counted, op_fast_add, 0).wait() / rank_n();
  auto max_kmers_processed = reduce_one(_num_kmers_counted, op_fast_max, 0).wait();
  auto avg_local_kmers = reduce_one(_num_kmers_counted_locally, op_fast_add, 0).wait() / rank_n();
  auto max_local_kmers = reduce_one(_num_kmers_counted_locally, op_fast_max, 0).wait();
  SLOG_VERBOSE("Avg kmers processed per rank ", avg_kmers_processed, " (balance ",
               (double)avg_kmers_processed / max_kmers_processed, ")\n");
  SLOG_VERBOSE("Avg local kmers processed per rank ", perc_str(avg_local_kmers, avg_kmers_processed), " (balance ",
               (double)avg_local_kmers / max_local_kmers, ")\n");
}

template <int MAX_K>
void KmerDHT<MAX_K>::purge_kmers(int threshold) {
  BarrierTimer timer(__FILEFUNC__);
  auto num_prior_kmers = get_num_kmers();
  int64_t num_purged = 0;
  for (auto it = kmers->begin(); it != kmers->end();) {
    auto kmer_counts = make_shared<KmerCounts>(it->second);
    if ((kmer_counts->count < threshold) || (kmer_counts->left_exts.is_zero() && kmer_counts->right_exts.is_zero())) {
      num_purged++;
      it = kmers->erase(it);
    } else {
      ++it;
    }
  }
  auto all_num_purged = reduce_one(num_purged, op_fast_add, 0).wait();
  SLOG_VERBOSE("Purged ", perc_str(all_num_purged, num_prior_kmers), " kmers below frequency threshold of ", threshold, "\n");
  estimated_error_rate = 1.0 - pow(1.0 - (double)all_num_purged / (double)num_prior_kmers, 1.0 / (double)Kmer<MAX_K>::get_k());
  SLOG_VERBOSE("Estimated per-base error rate from purge: ", estimated_error_rate, "\n");
}

template <int MAX_K>
void KmerDHT<MAX_K>::compute_kmer_exts() {
  BarrierTimer timer(__FILEFUNC__);
  for (auto &elem : *kmers) {
    auto kmer_counts = &elem.second;
    kmer_counts->left = kmer_counts->get_left_ext();
    kmer_counts->right = kmer_counts->get_right_ext();
  }
}

// one line per kmer, format:
// KMERCHARS LR N
// where L is left extension and R is right extension, one char, either X, F or A, C, G, T
// where N is the count of the kmer frequency
template <int MAX_K>
void KmerDHT<MAX_K>::dump_kmers(int k) {
  BarrierTimer timer(__FILEFUNC__);
  string dump_fname = "kmers-" + to_string(k) + ".txt.gz";
  get_rank_path(dump_fname, rank_me());
  zstr::bgzf_ofstream dump_file(dump_fname);
  ostringstream out_buf;
  ProgressBar progbar(kmers->size(), "Dumping kmers to " + dump_fname);
  int64_t i = 0;
  for (auto &elem : *kmers) {
    out_buf << elem.first << " " << elem.second.count << " " << elem.second.left << " " << elem.second.right << endl;
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
typename KmerDHT<MAX_K>::KmerMap::iterator KmerDHT<MAX_K>::local_kmers_begin() {
  return kmers->begin();
}

template <int MAX_K>
typename KmerDHT<MAX_K>::KmerMap::iterator KmerDHT<MAX_K>::local_kmers_end() {
  return kmers->end();
}

template <int MAX_K>
int32_t KmerDHT<MAX_K>::get_time_offset_us() {
  std::chrono::duration<double> t_elapsed = CLOCK_NOW() - start_t;
  return std::chrono::duration_cast<std::chrono::microseconds>(t_elapsed).count();
}
