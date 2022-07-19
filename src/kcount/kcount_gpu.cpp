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

#include "upcxx_utils.hpp"
#include "kcount.hpp"
#include "kmer_dht.hpp"
#include "devices_gpu.hpp"

#include "gpu-utils/gpu_utils.hpp"
#include "kcount-gpu/parse_and_pack.hpp"
#include "kcount-gpu/gpu_hash_table.hpp"

//#define SLOG_GPU(...) SLOG(KLMAGENTA, __VA_ARGS__, KNORM)
#define SLOG_GPU SLOG_VERBOSE

using namespace std;
using namespace upcxx_utils;
using namespace upcxx;
using namespace kcount_gpu;

template <int MAX_K>
struct SeqBlockInserter<MAX_K>::SeqBlockInserterState {
  ParseAndPackGPUDriver *pnp_gpu_driver;
  int64_t num_pnp_gpu_waits = 0;
  int num_block_calls = 0;
  int64_t num_kmers = 0;
  int64_t bytes_kmers_sent = 0;
  int64_t bytes_supermers_sent = 0;
  string seq_block;
  vector<kmer_count_t> depth_block;

  SeqBlockInserterState()
      : pnp_gpu_driver(nullptr) {
    seq_block.reserve(KCOUNT_SEQ_BLOCK_SIZE);
    depth_block.reserve(KCOUNT_SEQ_BLOCK_SIZE);
  }
};

template <int MAX_K>
SeqBlockInserter<MAX_K>::SeqBlockInserter(int qual_offset, int minimizer_len) {
  double init_time;
  state = new SeqBlockInserterState();
  state->pnp_gpu_driver = new ParseAndPackGPUDriver(rank_me(), rank_n(), qual_offset, Kmer<MAX_K>::get_k(),
                                                    Kmer<MAX_K>::get_N_LONGS(), minimizer_len, init_time);
  SLOG_GPU("Initialized PnP GPU driver in ", fixed, setprecision(3), init_time, " s\n");
}

template <int MAX_K>
SeqBlockInserter<MAX_K>::~SeqBlockInserter() {
  if (state->pnp_gpu_driver) {
    delete state->pnp_gpu_driver;
    state->pnp_gpu_driver = nullptr;
  }
  if (state) delete state;
}

template <int MAX_K>
static void process_block(SeqBlockInserter<MAX_K> *seq_block_inserter, dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  unsigned int num_valid_kmers = 0;
  auto state = seq_block_inserter->state;
  bool from_ctgs = !state->depth_block.empty();
  state->num_block_calls++;
  future<bool> fut = execute_in_thread_pool(
      [&state, &num_valid_kmers] { return state->pnp_gpu_driver->process_seq_block(state->seq_block, num_valid_kmers); });
  while (!fut.ready()) {
    state->num_pnp_gpu_waits++;
    progress();
  }
  bool success = fut.wait();
  if (!success) DIE("seq length is too high, ", state->seq_block.length(), " >= ", KCOUNT_SEQ_BLOCK_SIZE);
  state->bytes_kmers_sent += sizeof(KmerAndExt<MAX_K>) * num_valid_kmers;
  future<> fut_pnp = execute_in_thread_pool([&state] { state->pnp_gpu_driver->pack_seq_block(state->seq_block); });
  while (!fut_pnp.ready()) {
    state->num_pnp_gpu_waits++;
    progress();
  }
  fut_pnp.wait();
  int num_targets = (int)state->pnp_gpu_driver->supermers.size();
  for (int i = 0; i < num_targets; i++) {
    auto target = state->pnp_gpu_driver->supermers[i].target;
    auto offset = state->pnp_gpu_driver->supermers[i].offset;
    auto len = state->pnp_gpu_driver->supermers[i].len;
    Supermer supermer;
    int packed_len = len / 2;
    if (offset % 2 || len % 2) packed_len++;
    supermer.seq = state->pnp_gpu_driver->packed_seqs.substr(offset / 2, packed_len);
    if (offset % 2) supermer.seq[0] &= 15;
    if ((offset + len) % 2) supermer.seq[supermer.seq.length() - 1] &= 240;
    supermer.count = (from_ctgs ? state->depth_block[offset + 1] : (kmer_count_t)1);
    state->bytes_supermers_sent += supermer.get_bytes();
    kmer_dht->add_supermer(supermer, target);
    state->num_kmers += (2 * supermer.seq.length() - Kmer<MAX_K>::get_k());
    progress();
  }
}

template <int MAX_K>
void SeqBlockInserter<MAX_K>::process_seq(string &seq, kmer_count_t depth, dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  if (seq.length() >= KCOUNT_SEQ_BLOCK_SIZE)
    DIE("Oh dear, my laziness is revealed: the ctg seq is too long ", seq.length(), " for this GPU implementation ",
        KCOUNT_SEQ_BLOCK_SIZE);
  if (state->seq_block.length() + 1 + seq.length() >= KCOUNT_SEQ_BLOCK_SIZE) {
    process_block(this, kmer_dht);
    state->seq_block.clear();
    state->depth_block.clear();
  }
  state->seq_block += seq;
  state->seq_block += "_";
  if (depth) state->depth_block.insert(state->depth_block.end(), seq.length() + 1, depth);
}

template <int MAX_K>
void SeqBlockInserter<MAX_K>::done_processing(dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  if (kmer_dht->using_ctg_kmers) {
    DBG("in done_processing with seq_block empty? ", state->seq_block.empty(), "\n");
    // WARN("in done_processing with seq_block empty? ", state->seq_block.empty(), "\n");
  }
  if (!state->seq_block.empty()) process_block(this, kmer_dht);
  if (kmer_dht->using_ctg_kmers) {
    DBG("done_processing\n");
    // WARN("done_processing\n");
  }
  barrier();
  SLOG_GPU("GPU Parse-N-Pack stats:\n");
  SLOG_GPU("  number of calls to progress (on rank 0): ", state->num_pnp_gpu_waits, "\n");
  SLOG_GPU("  number of calls to PnP GPU kernel (on rank 0): ", state->num_block_calls, "\n");
  auto [gpu_time_tot, gpu_time_kernel] = state->pnp_gpu_driver->get_elapsed_times();
  SLOG_GPU("  elapsed times: ", fixed, setprecision(3), " total ", gpu_time_tot, ", kernel ", gpu_time_kernel, "\n");
  auto tot_supermers_bytes_sent = reduce_one(state->bytes_supermers_sent, op_fast_add, 0).wait();
  auto tot_kmers_bytes_sent = reduce_one(state->bytes_kmers_sent, op_fast_add, 0).wait();
  SLOG_VERBOSE("Total bytes sent in compressed supermers ", get_size_str(tot_supermers_bytes_sent), " (compression is ", fixed,
               setprecision(3), (double)tot_kmers_bytes_sent / tot_supermers_bytes_sent, " over kmers)\n");
  auto all_num_kmers = reduce_one(state->num_kmers, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", all_num_kmers, " kmers\n");
  barrier();
}

template <int MAX_K>
struct HashTableInserter<MAX_K>::HashTableInserterState {
  HashTableGPUDriver<MAX_K> ht_gpu_driver;

  HashTableInserterState()
      : ht_gpu_driver({}) {}
};

template <int MAX_K>
HashTableInserter<MAX_K>::HashTableInserter() {}

template <int MAX_K>
HashTableInserter<MAX_K>::~HashTableInserter() {
  if (state != nullptr) delete state;
}

template <int MAX_K>
void HashTableInserter<MAX_K>::init(size_t max_elems, bool use_qf) {
  this->use_qf = use_qf;
  state = new HashTableInserterState();
  double init_time;
  // calculate total slots for hash table. Reserve space for parse and pack
  size_t bytes_for_pnp = KCOUNT_SEQ_BLOCK_SIZE * (2 + Kmer<MAX_K>::get_N_LONGS() * sizeof(uint64_t) + sizeof(int));
  size_t gpu_bytes_reqd = 0, ht_bytes_used = 0, qf_bytes_used = 0;
  auto init_gpu_mem = gpu_utils::get_gpu_avail_mem();
  auto gpu_avail_mem_per_rank = (get_avail_gpu_mem_per_rank() - bytes_for_pnp) * 0.9;
  SLOG_GPU("Available GPU memory per rank for kmers hash table is ", get_size_str(gpu_avail_mem_per_rank), "\n");
  assert(state != nullptr);
  state->ht_gpu_driver.init(rank_me(), rank_n(), Kmer<MAX_K>::get_k(), max_elems, gpu_avail_mem_per_rank, init_time, gpu_bytes_reqd,
                            ht_bytes_used, qf_bytes_used, use_qf);
  auto capacity = state->ht_gpu_driver.get_capacity();
  SLOG_GPU("GPU read kmers hash table has capacity per rank of ", capacity, " for ", (int64_t)max_elems, " elements\n");
  SLOG_GPU("Using ", get_size_str(ht_bytes_used), " for the GPU hash table and ", get_size_str(qf_bytes_used), " for the QF\n");
  if (capacity < max_elems)
    SLOG_VERBOSE("GPU read kmers hash table has less than requested capacity: ", perc_str(capacity, max_elems),
                 "; full capacity requires ", get_size_str(gpu_bytes_reqd), " memory on GPU but only have ",
                 get_size_str(gpu_avail_mem_per_rank), "\n");
  SLOG_GPU("Initialized hash table GPU driver in ", fixed, setprecision(3), init_time, " s\n");
  auto gpu_used_mem = init_gpu_mem - gpu_utils::get_gpu_avail_mem();
  SLOG_GPU("GPU read kmers hash table used ", get_size_str(gpu_used_mem), " memory on GPU out of ",
           get_size_str(gpu_utils::get_gpu_tot_mem()), "\n");
}

template <int MAX_K>
void HashTableInserter<MAX_K>::init_ctg_kmers(size_t max_elems) {
  assert(state != nullptr);
  auto init_gpu_mem = gpu_utils::get_gpu_avail_mem();
  // we don't need to reserve space for either pnp or the read kmers because those have already reduced the gpu_avail_mem
  auto gpu_avail_mem_per_rank = get_avail_gpu_mem_per_rank();
  SLOG_GPU("Available GPU memory per rank for ctg kmers hash table is ", get_size_str(gpu_avail_mem_per_rank), "\n");
  state->ht_gpu_driver.init_ctg_kmers(max_elems, gpu_avail_mem_per_rank);
  SLOG_GPU("GPU ctg kmers hash table has capacity per rank of ", state->ht_gpu_driver.get_capacity(), " for ", fixed, max_elems,
           " elements\n");
  auto gpu_used_mem = init_gpu_mem - gpu_utils::get_gpu_avail_mem();
  SLOG_GPU("GPU ctg kmers hash table used ", get_size_str(gpu_used_mem), " memory on GPU out of ",
           get_size_str(gpu_utils::get_gpu_tot_mem()), "\n");
}

template <int MAX_K>
void HashTableInserter<MAX_K>::insert_supermer(const std::string &supermer_seq, kmer_count_t supermer_count) {
  assert(state != nullptr);
  state->ht_gpu_driver.insert_supermer(supermer_seq, supermer_count);
}

template <int MAX_K>
void HashTableInserter<MAX_K>::flush_inserts() {
  state->ht_gpu_driver.flush_inserts();
  auto avg_num_gpu_calls = reduce_one(state->ht_gpu_driver.get_num_gpu_calls(), op_fast_add, 0).wait() / rank_n();
  auto max_num_gpu_calls = reduce_one(state->ht_gpu_driver.get_num_gpu_calls(), op_fast_max, 0).wait();
  if (state->ht_gpu_driver.pass_type == kcount_gpu::READ_KMERS_PASS)
    SLOG_GPU("GPU hash table stats for read kmers pass:\n");
  else
    SLOG_GPU("GPU hash table stats for ctg kmers pass:\n");
  SLOG_GPU("  number of calls to hash table GPU driver: ", avg_num_gpu_calls, " avg, ", max_num_gpu_calls, " max\n");
  // a bunch of stats about the hash table on the GPU
  auto insert_stats = state->ht_gpu_driver.get_stats();
  uint64_t num_dropped_elems = reduce_one((uint64_t)insert_stats.dropped, op_fast_add, 0).wait();
  uint64_t num_attempted_inserts = reduce_one((uint64_t)insert_stats.attempted, op_fast_add, 0).wait();
  uint64_t num_inserts = reduce_one((uint64_t)insert_stats.new_inserts, op_fast_add, 0).wait();
  uint64_t capacity = state->ht_gpu_driver.get_capacity();
  uint64_t all_capacity = reduce_one(capacity, op_fast_add, 0).wait();
  if (num_dropped_elems) {
    if (num_dropped_elems > num_attempted_inserts / 10000)
      SWARN("GPU hash table: failed to insert ", perc_str(num_dropped_elems, num_attempted_inserts), " elements; capacity ",
            all_capacity);
    else
      SLOG_GPU("  failed to insert ", perc_str(num_dropped_elems, num_attempted_inserts), " elements; capacity ", all_capacity,
               "\n");
  }
  if (use_qf && state->ht_gpu_driver.pass_type == kcount_gpu::READ_KMERS_PASS) {
    uint64_t num_unique_qf = reduce_one((uint64_t)insert_stats.num_unique_qf, op_fast_add, 0).wait();
    // SLOG_GPU("  QF found ", perc_str(num_unique_qf, num_inserts), " unique kmers ", num_inserts, "\n");
    SLOG_GPU("  QF filtered out ", perc_str(num_unique_qf - num_inserts, num_unique_qf), " singletons\n");
    auto qf_max_load = reduce_one(state->ht_gpu_driver.get_qf_load_factor(), op_fast_max, 0).wait();
    auto qf_tot_load = reduce_one(state->ht_gpu_driver.get_qf_load_factor(), op_fast_add, 0).wait();
    double qf_avg_load = (double)qf_tot_load / rank_n();
    SLOG_GPU("  QF load factor ", fixed, setprecision(2), qf_avg_load, " avg ", qf_max_load, " max ", qf_avg_load / qf_max_load,
             " balance\n");
    uint64_t qf_failures = (uint64_t)state->ht_gpu_driver.get_qf_failures();
    //if (qf_failures) WARN("GQF failed to insert ", qf_failures, " items, load factor ", state->ht_gpu_driver.get_qf_load_factor());
    auto all_qf_failures = reduce_one(qf_failures, op_fast_add, 0).wait();
    if (all_qf_failures) SWARN("GQF failed to insert ", all_qf_failures, " items");
  }
  double load = (double)(insert_stats.new_inserts) / capacity;
  double avg_load_factor = reduce_one(load, op_fast_add, 0).wait() / rank_n();
  double max_load_factor = reduce_one(load, op_fast_max, 0).wait();
  SLOG_GPU("  load factor ", fixed, setprecision(2), avg_load_factor, " avg, ", max_load_factor, " max\n");
  SLOG_GPU("  final size per rank is ", insert_stats.new_inserts, " entries\n");
}

template <int MAX_K>
void HashTableInserter<MAX_K>::get_elapsed_time(double &insert_time, double &kernel_time) {
  state->ht_gpu_driver.get_elapsed_time(insert_time, kernel_time);
}

template <int MAX_K>
void HashTableInserter<MAX_K>::insert_into_local_hashtable(dist_object<KmerMap<MAX_K>> &local_kmers) {
  barrier();
  IntermittentTimer insert_timer("gpu insert to cpu timer");
  insert_timer.start();
  if (state->ht_gpu_driver.pass_type == CTG_KMERS_PASS) {
    int attempted_inserts = 0, dropped_inserts = 0, new_inserts = 0;
    state->ht_gpu_driver.done_ctg_kmer_inserts(attempted_inserts, dropped_inserts, new_inserts);
    barrier();
    auto num_dropped_elems = reduce_one((uint64_t)dropped_inserts, op_fast_add, 0).wait();
    auto num_attempted_inserts = reduce_one((uint64_t)attempted_inserts, op_fast_add, 0).wait();
    auto num_new_inserts = reduce_one((uint64_t)new_inserts, op_fast_add, 0).wait();
    SLOG_GPU("GPU ctg kmers hash table: inserted ", num_new_inserts, " new elements into read kmers hash table\n");
    auto all_capacity = reduce_one((uint64_t)state->ht_gpu_driver.get_capacity(), op_fast_add, 0).wait();
    if (num_dropped_elems) {
      if (num_dropped_elems > num_attempted_inserts / 10000)
        SWARN("GPU read kmers hash table: failed to insert ", perc_str(num_dropped_elems, num_attempted_inserts),
              " ctg kmers; total capacity ", all_capacity);
      else
        SLOG_VERBOSE("GPU read kmers hash table: failed to insert ", perc_str(num_dropped_elems, num_attempted_inserts),
                     " ctg kmers; total capacity ", all_capacity, "\n");
    }
  }
  barrier();
  int num_dropped = 0, num_entries = 0, num_purged = 0;
  state->ht_gpu_driver.done_all_inserts(num_dropped, num_entries, num_purged);
  barrier();
  if (num_dropped)
    WARN("GPU dropped ", num_dropped, " entries out of ", num_entries, " when compacting to output hash table" KNORM "\n");

  auto all_capacity = reduce_one((uint64_t)state->ht_gpu_driver.get_final_capacity(), op_fast_add, 0).wait();
  auto all_num_purged = reduce_one((uint64_t)num_purged, op_fast_add, 0).wait();
  auto all_num_entries = reduce_one((uint64_t)num_entries, op_fast_add, 0).wait();
  auto prepurge_num_entries = all_num_entries + all_num_purged;
  SLOG_GPU("GPU hash table: purged ", perc_str(all_num_purged, prepurge_num_entries), " singleton kmers out of ",
           prepurge_num_entries, "\n");
  SLOG_GPU("GPU hash table final size is ", all_num_entries, " entries and final load factor is ",
           ((double)all_num_entries / all_capacity), "\n");
  barrier();

  // add some space for the ctg kmers
  local_kmers->reserve(num_entries * 1.5);
  uint64_t invalid = 0;
  while (true) {
    auto [kmer_array, count_exts] = state->ht_gpu_driver.get_next_entry();
    if (!kmer_array) break;
    // empty slot
    if (!count_exts->count) continue;
    if ((char)count_exts->left == 'X' && (char)count_exts->right == 'X') {
      // these are eliminated during purging in CPU version
      invalid++;
      continue;
    }
    if ((count_exts->count < 2)) {
      WARN("Found a kmer that should have been purged, count is ", count_exts->count);
      invalid++;
      continue;
    }
    KmerCounts kmer_counts = {.uutig_frag = nullptr,
                              .count = static_cast<kmer_count_t>(min(count_exts->count, static_cast<count_t>(UINT16_MAX))),
                              .left = (char)count_exts->left,
                              .right = (char)count_exts->right};
    Kmer<MAX_K> kmer(reinterpret_cast<const uint64_t *>(kmer_array->longs));
    const auto it = local_kmers->find(kmer);
    if (it != local_kmers->end())
      WARN("Found a duplicate kmer ", kmer.to_string(), " - shouldn't happen: existing count ", it->second.count, " new count ",
           kmer_counts.count);
    local_kmers->insert({kmer, kmer_counts});
  }
  insert_timer.stop();

  auto all_avg_elapsed_time = reduce_one(insert_timer.get_elapsed(), op_fast_add, 0).wait() / rank_n();
  auto all_max_elapsed_time = reduce_one(insert_timer.get_elapsed(), op_fast_max, 0).wait();
  SLOG_GPU("Inserting kmers from GPU to cpu hash table took ", all_avg_elapsed_time, " avg, ", all_max_elapsed_time, " max\n");
  auto all_kmers_size = reduce_one((uint64_t)local_kmers->size(), op_fast_add, 0).wait();
  if (local_kmers->size() != (num_entries - invalid))
    WARN("kmers->size() is ", local_kmers->size(), " != ", (num_entries - invalid), " num_entries");
  auto all_invalid = reduce_one((uint64_t)invalid, op_fast_add, 0).wait();
  if (all_kmers_size != all_num_entries - all_invalid)
    SWARN("CPU kmer counts not equal to gpu kmer counts: ", all_kmers_size, " != ", (all_num_entries - all_invalid),
          " all_num_entries: ", all_num_entries, " all_invalid: ", all_invalid);
  double gpu_insert_time = 0, gpu_kernel_time = 0;
  state->ht_gpu_driver.get_elapsed_time(gpu_insert_time, gpu_kernel_time);
  auto avg_gpu_insert_time = reduce_one(gpu_insert_time, op_fast_add, 0).wait() / rank_n();
  auto max_gpu_insert_time = reduce_one(gpu_insert_time, op_fast_max, 0).wait();
  auto avg_gpu_kernel_time = reduce_one(gpu_kernel_time, op_fast_add, 0).wait() / rank_n();
  auto max_gpu_kernel_time = reduce_one(gpu_kernel_time, op_fast_max, 0).wait();
  SLOG_GPU("Elapsed GPU time for kmer hash tables:\n");
  double load_balance = (double)avg_gpu_insert_time / max_gpu_insert_time;
  SLOG_GPU("  insert: ", fixed, setprecision(3), (load_balance < 0.5 ? KLRED : ""), avg_gpu_insert_time, " avg, ",
           max_gpu_insert_time, " max, load balance ", load_balance, KNORM, "\n");
  SLOG_GPU("  kernel: ", fixed, setprecision(3), avg_gpu_kernel_time, " avg, ", max_gpu_kernel_time, " max\n");
  barrier();
}

#define seq_block_inserter_K(KMER_LEN) template struct SeqBlockInserter<KMER_LEN>;
#define HASH_TABLE_INSERTER_K(KMER_LEN) template class HashTableInserter<KMER_LEN>;

seq_block_inserter_K(32);
HASH_TABLE_INSERTER_K(32);
#if MAX_BUILD_KMER >= 64
seq_block_inserter_K(64);
HASH_TABLE_INSERTER_K(64);
#endif
#if MAX_BUILD_KMER >= 96
seq_block_inserter_K(96);
HASH_TABLE_INSERTER_K(96);
#endif
#if MAX_BUILD_KMER >= 128
seq_block_inserter_K(128);
HASH_TABLE_INSERTER_K(128);
#endif
#if MAX_BUILD_KMER >= 160
seq_block_inserter_K(160);
HASH_TABLE_INSERTER_K(160);
#endif
#undef seq_block_inserter_K
#undef HASH_TABLE_INSERTER_K
