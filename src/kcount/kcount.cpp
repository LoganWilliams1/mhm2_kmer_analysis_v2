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
#include "utils.hpp"
#include "kcount.hpp"

//#define DBG_ADD_KMER DBG
#define DBG_ADD_KMER(...)

using namespace std;
using namespace upcxx_utils;
using namespace upcxx;

template <int MAX_K>
static void count_kmers(unsigned kmer_len, int qual_offset, vector<PackedReads *> &packed_reads_list,
                        dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_reads = 0;
  int64_t num_lines = 0;
  int64_t num_bad_quals = 0;
  int64_t tot_read_len = 0;

  string progbar_prefix = "";
  IntermittentTimer t_pp(__FILENAME__ + string(":kmer parse and pack"));
  barrier();
  SeqBlockInserter<MAX_K> seq_block_inserter(qual_offset, kmer_dht->get_minimizer_len());
  int64_t total_local_num_reads = PackedReads::get_total_local_num_reads(packed_reads_list);
  int64_t total_local_raw_kmers = 0;
  ProgressBar progbar(total_local_num_reads, "Processing reads to count kmers");

  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    string id, seq, quals;
    while (true) {
      if (!packed_reads->get_next_read(id, seq, quals)) break;
      num_reads++;
      progbar.update();
      if (seq.length() < kmer_len) continue;
      tot_read_len += seq.length();
      for (int i = 0; i < seq.length(); i++) {
        if (quals[i] < qual_offset + KCOUNT_QUAL_CUTOFF) {
          seq[i] = tolower(seq[i]);
          num_bad_quals++;
        }
      }
      total_local_raw_kmers += seq.length() - kmer_len + 1;
      seq_block_inserter.process_seq(seq, 0, kmer_dht);
      progress();
    }
  }
  seq_block_inserter.done_processing(kmer_dht);
  progbar.done();
  kmer_dht->flush_updates();
  auto all_num_reads = reduce_one(num_reads, op_fast_add, 0).wait();
  auto all_raw_kmers = reduce_one(total_local_raw_kmers, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", all_num_reads, " reads ", all_raw_kmers, " raw kmers\n");
  auto avg_supermer_inserts = reduce_one(kmer_dht->get_num_supermer_inserts(), op_fast_add, 0).wait() / rank_n();
  auto max_supermer_inserts = reduce_one(kmer_dht->get_num_supermer_inserts(), op_fast_max, 0).wait();
  SLOG_VERBOSE("Avg supermer inserts ", avg_supermer_inserts, " max ", max_supermer_inserts, " load ", setprecision(3), fixed,
       (double)avg_supermer_inserts / max_supermer_inserts, "\n");
  auto all_num_bad_quals = reduce_one(num_bad_quals, op_fast_add, 0).wait();
  auto all_tot_read_len = reduce_one(tot_read_len, op_fast_add, 0).wait();
  if (all_num_bad_quals) SLOG_VERBOSE("Found ", perc_str(all_num_bad_quals, all_tot_read_len), " bad quality positions\n");
};

template <int MAX_K>
static void add_ctg_kmers(unsigned kmer_len, unsigned prev_kmer_len, Contigs &ctgs, dist_object<KmerDHT<MAX_K>> &kmer_dht) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_prev_kmers = kmer_dht->get_num_kmers();
  int64_t total_local_raw_kmers = 0;

  ProgressBar progbar(ctgs.size(), "Adding extra contig kmers from kmer length " + to_string(prev_kmer_len));
  auto start_local_num_kmers = kmer_dht->get_local_num_kmers();

  SeqBlockInserter<MAX_K> seq_block_inserter(0, kmer_dht->get_minimizer_len());
  barrier();
  DBG("After seq_block_inserter constructor, with ", ctgs.size(), " ctgs\n");
  // WARN("After seq_block_inserter constructor, with ", ctgs.size(), " ctgs\n");
  //  estimate number of kmers from ctgs
  int64_t max_kmers = 0;
  for (auto &ctg : ctgs) {
    if (ctg.seq.length() > kmer_len) max_kmers += ctg.seq.length() - kmer_len + 1;
  }
  int64_t all_max_kmers = reduce_all(max_kmers, op_fast_add).wait();
  // increase max kmers to allow for load factor 0.67
  kmer_dht->init_ctg_kmers(1.5 * all_max_kmers / rank_n());
  barrier();
  DBG("after kmer_dht->init_ctg_kmers\n");
  DBG("looping over ", ctgs.size(), " ctgs\n");
  // WARN("after kmer_dht->init_ctg_kmers\n");
  // WARN("looping over ", ctgs.size(), " ctgs\n");
  for (auto it = ctgs.begin(); it != ctgs.end(); ++it) {
    auto ctg = it;
    progbar.update();
    if (ctg->seq.length() < kmer_len + 2) continue;
    total_local_raw_kmers += ctg->seq.length() - kmer_len + 1;
    seq_block_inserter.process_seq(ctg->seq, ctg->get_uint16_t_depth(), kmer_dht);
  }
  DBG("after ctgs loop\n");
  // WARN("after ctgs loop\n");
  seq_block_inserter.done_processing(kmer_dht);
  progbar.done();
  kmer_dht->flush_updates();
  auto all_num_ctgs = reduce_one(ctgs.size(), op_fast_add, 0).wait();
  auto all_raw_kmers = reduce_one(total_local_raw_kmers, op_fast_add, 0).wait();
  SLOG_VERBOSE("Processed a total of ", all_num_ctgs, " contigs and ", total_local_raw_kmers, " raw kmers\n");
};

template <int MAX_K>
void analyze_kmers(unsigned kmer_len, unsigned prev_kmer_len, int qual_offset, vector<PackedReads *> &packed_reads_list,
                   int dmin_thres, Contigs &ctgs, dist_object<KmerDHT<MAX_K>> &kmer_dht, bool dump_kmers) {
  BarrierTimer timer(__FILEFUNC__);
  auto fut_has_contigs = upcxx::reduce_all(ctgs.size(), upcxx::op_fast_max).then([](size_t max_ctgs) { return max_ctgs > 0; });
  _dmin_thres = dmin_thres;

  count_kmers(kmer_len, qual_offset, packed_reads_list, kmer_dht);
  barrier();
  if (fut_has_contigs.wait()) {
    add_ctg_kmers(kmer_len, prev_kmer_len, ctgs, kmer_dht);
    barrier();
  }
  kmer_dht->finish_updates();
  if (dump_kmers) kmer_dht->dump_kmers();
  barrier();
  kmer_dht->clear_stores();
};
