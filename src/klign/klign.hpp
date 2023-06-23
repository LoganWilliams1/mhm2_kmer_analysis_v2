#pragma once

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

#include <memory>
#include <upcxx/upcxx.hpp>

#include "alignments.hpp"
#include "contigs.hpp"
#include "kmer.hpp"
#include "packed_reads.hpp"
#include "utils.hpp"

#include "upcxx_utils/timers.hpp"

using std::shared_ptr;
using upcxx::global_ptr;

#include "upcxx_utils/timers.hpp"

struct KlignTimers {
  upcxx_utils::IntermittentTimer fetch_ctg_maps, compute_alns, rget_ctg_seqs, aln_kernel, aln_kernel_mem, aln_kernel_block,
      aln_kernel_launch, sort_t;

  KlignTimers()
      : fetch_ctg_maps("klign: fetch ctg maps")
      , compute_alns("klign: compute alns")
      , rget_ctg_seqs("klign: rget ctg seqs")
      , aln_kernel("klign: aln kernel")
      , aln_kernel_mem("klign: aln kernel mem")
      , aln_kernel_block("klign: aln kernel block")
      , aln_kernel_launch("klign: aln kernel launch")
      , sort_t("klign: sort alns") {}

  void done_all() {
    fetch_ctg_maps.done_all_async();
    compute_alns.done_all_async();
    rget_ctg_seqs.done_all_async();
    aln_kernel.done_all_async();
    aln_kernel_mem.done_all_async();
    aln_kernel_block.done_all_async();
    aln_kernel_launch.done_all_async();
    sort_t.done_all_async();
  }

  void clear() {
    fetch_ctg_maps.clear();
    compute_alns.clear();
    rget_ctg_seqs.clear();
    aln_kernel.clear();
    aln_kernel_mem.clear();
    aln_kernel_block.clear();
    aln_kernel_launch.clear();
    sort_t.clear();
  }
};  // struct KlignTimers

struct CtgLoc {
  cid_t cid;
  global_ptr<char> seq_gptr;
  int32_t clen;
  float depth;
  int32_t pos : 31;    // pack int 31 bits
  uint32_t is_rc : 1;  // bool
};

struct CtgAndReadLoc {
  CtgLoc ctg_loc;
  int32_t cstart;
  int32_t pos_in_read : 31;  // pack into 31 bits
  uint32_t read_is_rc : 1;   // bool
};

template <int MAX_K>
struct KmerAndCtgLoc {
  Kmer<MAX_K> kmer;
  CtgLoc ctg_loc;
  UPCXX_SERIALIZED_FIELDS(kmer, ctg_loc);
};

struct CtgLocAndKmerIdx {
  CtgLoc ctg_loc;
  int kmer_i;
};

using CtgAndReadLocsMap = HASH_TABLE<cid_t, vector<CtgAndReadLoc>>;

struct ReadRecord {
  int64_t index : 40;
  int64_t rlen : 24;

  CtgAndReadLocsMap aligned_ctgs_map;

  ReadRecord()
      : index(-1)
      , rlen(0) {}

  ReadRecord(int index, int rlen)
      : index(index)
      , rlen(rlen)
      , aligned_ctgs_map{} {}

  bool is_valid() const { return index >= 0 && rlen > 0; }
};  // struct ReadRecord

struct ReadRecordPtr {
  ReadRecord *read_record;
  int read_offset;
  bool is_rc;
};  // struct ReadRecordPtr

template <int MAX_K>
class KmerCtgDHT;

template <int MAX_K>
double find_alignments(unsigned kmer_len, PackedReadsList &packed_reads_list, int max_store_size, int max_rpcs_in_flight,
                       Contigs &ctgs, Alns &alns, int seed_space, int rlen_limit, bool report_cigar, bool use_blastn_scores,
                       int min_ctg_len);

template <int MAX_K>
shared_ptr<KmerCtgDHT<MAX_K>> build_kmer_ctg_dht(unsigned, int, int, Contigs &, int, bool);

template <int MAX_K>
void compute_alns(PackedReads *, vector<ReadRecord> &, Alns &, int, int, bool, bool, int64_t, KlignTimers &);

template <int MAX_K>
void fetch_ctg_maps(KmerCtgDHT<MAX_K> &, PackedReads *, vector<ReadRecord> &, int, KlignTimers &);

// Reduce compile time by instantiating templates of common types
// extern template declarations are in kmer.hpp
// template instantiations each happen in src/CMakeLists via kmer-extern-template.in.cpp

#define __MACRO_KLIGN__(KMER_LEN, MODIFIER)                                                                                        \
  MODIFIER double find_alignments<KMER_LEN>(unsigned, PackedReadsList &, int, int, Contigs &, Alns &, int, int, bool, bool, int);  \
  MODIFIER shared_ptr<KmerCtgDHT<KMER_LEN>> build_kmer_ctg_dht<KMER_LEN>(unsigned, int, int, Contigs &, int, bool);                \
  MODIFIER void compute_alns<KMER_LEN>(PackedReads *, vector<ReadRecord> &, Alns &, int, int, bool, bool, int64_t, KlignTimers &); \
  MODIFIER void fetch_ctg_maps<KMER_LEN>(KmerCtgDHT<KMER_LEN> &, PackedReads *, vector<ReadRecord> &, int, KlignTimers &);

__MACRO_KLIGN__(32, extern template);
#if MAX_BUILD_KMER >= 64
__MACRO_KLIGN__(64, extern template);
#endif
#if MAX_BUILD_KMER >= 96
__MACRO_KLIGN__(96, extern template);
#endif
#if MAX_BUILD_KMER >= 128
__MACRO_KLIGN__(128, extern template);
#endif
#if MAX_BUILD_KMER >= 160
__MACRO_KLIGN__(160, extern template);
#endif
