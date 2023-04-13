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

#include <fstream>
#include <iostream>
#include <regex>
#include <upcxx/upcxx.hpp>
#include <memory>

#include "alignments.hpp"
#include "contigs.hpp"
#include "kmer_dht.hpp"
#include "packed_reads.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"
#include "localassm_core.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;
using namespace localassm_core;

void extend_ctgs(CtgsWithReadsDHT &ctgs_dht, Contigs &ctgs, int insert_avg, int insert_stddev, int max_kmer_len, int kmer_len,
                 int qual_offset, unsigned max_read_size);

void localassm(int max_kmer_len, int kmer_len, PackedReadsList &packed_reads_list, int insert_avg, int insert_stddev,
               int qual_offset, Contigs &ctgs, const Alns &alns) {
  BarrierTimer timer(__FILEFUNC__);
  // get the contig count for ctgs_dht size hints
  int64_t num_ctg_bases = 0;
  for (auto &ctg : ctgs) num_ctg_bases += ctg.seq.length();
  // get the read and base counts for reads_to_ctgs size hints
  int64_t num_reads = 0;
  int64_t num_read_bases = 0;
  unsigned max_read_len = 0;
  for (auto &packed_reads : packed_reads_list) {
    auto l_reads = packed_reads->get_local_num_reads();
    max_read_len = std::max(max_read_len, (unsigned) packed_reads->get_max_read_len());
    num_reads += l_reads;
  }
  num_read_bases = num_reads * max_read_len;
  LOG("num_reads=", num_reads, " num_read_bases=", num_read_bases, " max_read_len=", max_read_len, "\n");
  CtgsWithReadsDHT ctgs_dht(ctgs.size(), num_ctg_bases);
  add_ctgs(ctgs_dht, ctgs);
  ctgs_dht.prep_ctg_read_store(num_reads, num_read_bases);
  ReadsToCtgsDHT reads_to_ctgs(100, alns.size());
  // extract read id to ctg id mappings from alignments
  process_alns(alns, reads_to_ctgs, insert_avg, insert_stddev);
  // extract read seqs and add to ctgs
  process_reads(max_kmer_len, packed_reads_list, reads_to_ctgs, ctgs_dht);
  // free the reads_to_contigs map
  reads_to_ctgs.clear();
  // clear out the local contigs
  ctgs.clear();
  ctgs.set_capacity(ctgs_dht.get_local_num_ctgs());
  // extend contigs using locally mapped reads
  extend_ctgs(ctgs_dht, ctgs, insert_avg, insert_stddev, max_kmer_len, kmer_len, qual_offset, max_read_len);
}
