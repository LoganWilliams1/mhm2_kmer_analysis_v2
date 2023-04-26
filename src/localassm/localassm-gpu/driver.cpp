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

#include <numeric>
#include <cassert>
#include <memory>

#include "gpu-utils/gpu_compatibility.hpp"
#include "gpu-utils/gpu_utils.hpp"
#include "gpu-utils/gpu_common.hpp"

#include "driver.hpp"
#include "kernel.hpp"

using namespace std;

static void revcomp(char *str, char *str_rc, unsigned int size) {
  if (size == FULL || size == EMPTY) {
    printf("Bad revcomp\n");
    return;
  }
  int size_rc = 0;
  for (int i = size - 1; i >= 0; i--) {
    switch (str[i]) {
      case 'A':
      case 'a': str_rc[size_rc] = 'T'; break;
      case 'C':
      case 'c': str_rc[size_rc] = 'G'; break;
      case 'G':
      case 'g': str_rc[size_rc] = 'C'; break;
      case 'T':
      case 't': str_rc[size_rc] = 'A'; break;
      case 'N':
      case 'n': str_rc[size_rc] = 'N'; break;
      case 'U':
      case 'R':
      case 'Y':
      case 'K':
      case 'M':
      case 'S':
      case 'W':
      case 'B':
      case 'D':
      case 'H':
      case 'V': str_rc[size_rc] = 'N'; break;
      default:
        std::cerr << "Illegal revcomp char " << ((str[i] >= 32 && str[i] <= 126) ? str[i] : ' ') << " int=" << (int)str[i] << "\n";
        break;
    }
    size_rc++;
  }
}

static std::string revcomp(std::string instr) {
  std::string str_rc;
  if (instr.size() == FULL || instr.size() == EMPTY) {
    printf("Bad revcomp2\n");
    return str_rc;
  }
  for (int i = instr.size() - 1; i >= 0; i--) {
    switch (instr[i]) {
      case 'A':
      case 'a': str_rc += 'T'; break;
      case 'C':
      case 'c': str_rc += 'G'; break;
      case 'G':
      case 'g': str_rc += 'C'; break;
      case 'T':
      case 't': str_rc += 'A'; break;
      case 'N':
      case 'n': str_rc += 'N'; break;
      case 'U':
      case 'R':
      case 'Y':
      case 'K':
      case 'M':
      case 'S':
      case 'W':
      case 'B':
      case 'D':
      case 'H':
      case 'V': str_rc += 'N'; break;
      default:
        std::cerr << "Illegal revcomp char:" << ((instr[i] >= 32 && instr[i] <= 126) ? instr[i] : ' ') << " int=" << (int)instr[i]
                  << "\n";
        break;
    }
  }
  return str_rc;
}

void localassm_driver::localassm_driver(vector<CtgWithReads *> &data_in, uint32_t max_ctg_size, uint32_t max_read_size,
                                        uint32_t max_r_count, uint32_t max_l_count, int mer_len, int max_kmer_len,
                                        accum_data &sizes_vecs, int walk_len_limit, int qual_offset, int my_rank,
                                        size_t gpu_mem_avail, int debug_line) {
  gpu_utils::set_gpu_device(my_rank);
  size_t TOO_BIG = 2147483647L;
  int max_mer_len = max_kmer_len;  // mer_len;// max_mer_len needs to come from macro (121) and mer_len is the mer_len for current
                                   // go
  uint64_t tot_extensions = data_in.size();
  uint32_t max_read_count = max_r_count > max_l_count ? max_r_count : max_l_count;
  int max_walk_len = walk_len_limit;
  uint64_t ht_tot_size = accumulate(sizes_vecs.ht_sizes.begin(), sizes_vecs.ht_sizes.end(), (uint64_t) 0);
  uint64_t total_r_reads = accumulate(sizes_vecs.r_reads_count.begin(), sizes_vecs.r_reads_count.end(), (uint64_t) 0);
  uint64_t total_l_reads = accumulate(sizes_vecs.l_reads_count.begin(), sizes_vecs.l_reads_count.end(), (uint64_t) 0);
  uint64_t total_ctg_len = accumulate(sizes_vecs.ctg_sizes.begin(), sizes_vecs.ctg_sizes.end(), (uint64_t) 0);

  size_t gpu_mem_req = sizeof(int32_t) * tot_extensions * 2                            // prefix_ht_size_d, ctg_seq_offsets_d
                       + sizeof(int32_t) * tot_extensions * 2                          // rds_l_cnt_offset_d, rds_r_cnt_offset_d
                       + sizeof(int32_t) * tot_extensions                              // final_walk_lens_d
                       + sizeof(int64_t) * tot_extensions * 2                          // cid_d
                       + sizeof(int32_t) * (total_l_reads + total_r_reads)             // reads_l_offset_d, reads_r_offset_d
                       + sizeof(char) * total_ctg_len                                  // ctg_seqs_d (contigs' sequence)
                       + sizeof(char) * total_l_reads * max_read_size * 2              // reads_left_d, quals_left_d
                       + sizeof(char) * total_r_reads * max_read_size * 2              // reads_right_d, quals_right_d
                       + sizeof(double) * tot_extensions                               // depth_h
                       + sizeof(int32_t) * 3                                           // term_counts_d
                       + sizeof(loc_ht) * (ht_tot_size + 1)                            // d_ht // changed to try the new method
                       + sizeof(loc_ht_bool) * (tot_extensions * max_walk_len + 1)     // d_ht_bool
                       + sizeof(char) * tot_extensions * (max_walk_len + max_mer_len)  // mer_walk_temp_d
                       + sizeof(char) * tot_extensions * max_walk_len                  // longest_walks_d
      ;

  assert(gpu_mem_avail > 0);
  // factor is to buffer for the extra mem that is used when allocating once and using again
  float factor = 0.80;
  uint64_t iterations = ceil(((double)gpu_mem_req) / ((double)gpu_mem_avail * factor));
  iterations = std::max(iterations, (max_ctg_size * tot_extensions + TOO_BIG - 1) / TOO_BIG);
  iterations = std::max(iterations, (ht_tot_size + TOO_BIG - 1) / TOO_BIG);
  iterations = std::max(iterations, (total_r_reads * max_read_size + TOO_BIG - 1) / TOO_BIG);
  iterations = std::max(iterations, (total_l_reads * max_read_size + TOO_BIG - 1) / TOO_BIG);
  iterations = std::max(iterations, (total_ctg_len + TOO_BIG - 1) / TOO_BIG);
  assert(iterations > 0);

  // run more iterations to ensure no int32 overflows
  if (max_ctg_size * tot_extensions / iterations > TOO_BIG || ht_tot_size / iterations > TOO_BIG ||
      total_r_reads * max_read_size / iterations > TOO_BIG || total_l_reads * max_read_size / iterations > TOO_BIG ||
      total_ctg_len / iterations > TOO_BIG)
    fprintf(stderr,
            "WARN potential overflow rank=%d tot_extensions=%lu iterations=%lu max_ctg_size=%u ht_tot_size=%lu total_r_reads=%lu "
            "total_l_reads=%lu max_read_size=%u total_ctg_len=%lu\n",
            my_rank, tot_extensions, iterations, max_ctg_size, ht_tot_size, total_r_reads, total_l_reads, max_read_size,
            total_ctg_len);

  const uint64_t max_slice_size = (tot_extensions + iterations - 1) / iterations;
  assert(max_slice_size < TOO_BIG);
  assert(max_slice_size > 0);

  // to get the largest ht size for any iteration and allocate GPU memory for that (once)
  uint64_t max_ht = 0, max_r_rds_its = 0, max_l_rds_its = 0, max_ctg_len_its = 0, test_sum = 0;
  for (unsigned i = 0; i < iterations; i++) {
    uint64_t extensions_offset = max_slice_size * i;
    uint64_t num_extensions = tot_extensions - extensions_offset;
    if (num_extensions > max_slice_size) num_extensions = max_slice_size;
    assert(num_extensions > 0 && num_extensions <= max_slice_size);
    uint64_t temp_max_ht = 0, temp_max_r_rds = 0, temp_max_l_rds = 0, temp_max_ctg_len = 0;

    temp_max_ht = accumulate(sizes_vecs.ht_sizes.begin() + extensions_offset,
                             sizes_vecs.ht_sizes.begin() + extensions_offset + num_extensions, (uint64_t) 0);
    temp_max_r_rds = accumulate(sizes_vecs.r_reads_count.begin() + extensions_offset,
                                sizes_vecs.r_reads_count.begin() + extensions_offset + num_extensions, (uint64_t) 0);
    temp_max_l_rds = accumulate(sizes_vecs.l_reads_count.begin() + extensions_offset,
                                sizes_vecs.l_reads_count.begin() + extensions_offset + num_extensions, (uint64_t) 0);
    temp_max_ctg_len = accumulate(sizes_vecs.ctg_sizes.begin() + extensions_offset,
                                  sizes_vecs.ctg_sizes.begin() + extensions_offset + num_extensions, (uint64_t) 0);

    if (temp_max_ht > max_ht) max_ht = temp_max_ht;
    if (temp_max_r_rds > max_r_rds_its) max_r_rds_its = temp_max_r_rds;
    if (temp_max_l_rds > max_l_rds_its) max_l_rds_its = temp_max_l_rds;
    if (temp_max_ctg_len > max_ctg_len_its) max_ctg_len_its = temp_max_ctg_len;
    test_sum += temp_max_ht;
  }
  if (max_ht > TOO_BIG || max_r_rds_its > TOO_BIG || max_l_rds_its > TOO_BIG || max_ctg_len_its > TOO_BIG)
    fprintf(stderr,
            "overflow with iterations=%lu max_slice_size=%lu for tot_extensions=%lu max_ht=%lu max_r_rds_its=%lu max_l_rds_its=%lu "
            "max_ctg_len_its=%lu\n",
            iterations, max_slice_size, tot_extensions, max_ht, max_r_rds_its, max_l_rds_its, max_ctg_len_its);

  // allocating maximum possible memory for a single iteration
  uint64_t all_walk_size = tot_extensions * max_walk_len;
#ifdef DEBUG
  if (!my_rank)
    fprintf(stderr,
            "Allocating memory for iterations=%lu max_slice_size=%lu and tot_extensions=%lu of max_walk_len=%u and 2x walks %lu\n",
            iterations, max_slice_size, tot_extensions, max_walk_len, all_walk_size);
#endif

  uint64_t max_ctg_seqs_h = max_ctg_size * max_slice_size;
  unique_ptr<char[]> ctg_seqs_h{new char[max_ctg_seqs_h]};
  unique_ptr<uint64_t[]> cid_h{new uint64_t[max_slice_size]};
  unique_ptr<char[]> ctgs_seqs_rc_h{
      new char[max_ctg_seqs_h]};  // revcomps not requried on GPU, ctg space will be re-used on GPU, but if we want
                                  // to do right left extensions in parallel, then we need separate space on GPU
  unique_ptr<uint32_t[]> ctg_seq_offsets_h{new uint32_t[max_slice_size]};
  unique_ptr<double[]> depth_h{new double[max_slice_size]};
  int32_t max_reads_left_h = max_read_size * max_l_rds_its;
  int32_t max_reads_right_h = max_read_size * max_r_rds_its;
  unique_ptr<char[]> reads_left_h{new char[max_reads_left_h]};
  unique_ptr<char[]> reads_right_h{new char[max_reads_right_h]};
  unique_ptr<char[]> quals_right_h{new char[max_reads_right_h]};
  unique_ptr<char[]> quals_left_h{new char[max_reads_left_h]};
  unique_ptr<uint32_t[]> reads_l_offset_h{new uint32_t[max_l_rds_its]};
  unique_ptr<uint32_t[]> reads_r_offset_h{new uint32_t[max_r_rds_its]};
  unique_ptr<uint32_t[]> rds_l_cnt_offset_h{new uint32_t[max_slice_size]};
  unique_ptr<uint32_t[]> rds_r_cnt_offset_h{new uint32_t[max_slice_size]};
  unique_ptr<uint32_t[]> term_counts_h{new uint32_t[3]};
  unique_ptr<char[]> longest_walks_r_h{new char[all_walk_size]};  // reserve memory for all the walks
  unique_ptr<char[]> longest_walks_l_h{new char[all_walk_size]};  // not needed on device, will re-use right walk memory
  unique_ptr<uint32_t[]> final_walk_lens_r_h{new uint32_t[max_slice_size * iterations]};  // reserve memory for all the walks.
  unique_ptr<uint32_t[]> final_walk_lens_l_h{new uint32_t[max_slice_size * iterations]};  // not needed on device, will re use right walk memory
  unique_ptr<uint32_t[]> prefix_ht_size_h{new uint32_t[max_slice_size]};

  gpu_mem_req = sizeof(int32_t) * max_slice_size * 2                            // prefix_ht_size_d, ctg_seq_offsets_d
                + sizeof(int32_t) * max_slice_size * 2                          // rds_l_cnt_offset_d, rds_r_cnt_offset_d
                + sizeof(int32_t) * max_slice_size                              // final_walk_lens_d
                + sizeof(int64_t) * max_slice_size * 2                          // cid_d
                + sizeof(int32_t) * (max_r_rds_its + max_l_rds_its)             // reads_l_offset_d, reads_r_offset_d
                + sizeof(char) * max_ctg_len_its                                // ctg_seqs_d (contigs' sequence)
                + sizeof(char) * max_l_rds_its * max_read_size * 2              // reads_left_d, quals_left_d
                + sizeof(char) * max_r_rds_its * max_read_size * 2              // reads_right_d, quals_right_d
                + sizeof(double) * max_slice_size                               // depth_h
                + sizeof(int32_t) * 3                                           // term_counts_d
                + sizeof(loc_ht) * (max_ht + 1)                                 // d_ht // changed to try the new method
                + sizeof(loc_ht_bool) * (max_slice_size * max_walk_len + 1)     // d_ht_bool
                + sizeof(char) * max_slice_size * (max_walk_len + max_mer_len)  // mer_walk_temp_d
                + sizeof(char) * max_slice_size * max_walk_len                  // longest_walks_d
      ;

  uint32_t *ctg_seq_offsets_d, *reads_l_offset_d, *reads_r_offset_d;
  uint64_t *cid_d;
  uint32_t *rds_l_cnt_offset_d, *rds_r_cnt_offset_d, *prefix_ht_size_d;

  char *ctg_seqs_d, *reads_left_d, *reads_right_d, *quals_left_d, *quals_right_d;
  char *longest_walks_d, *mer_walk_temp_d;
  double *depth_d;
  uint32_t *term_counts_d;
  loc_ht *d_ht;
  loc_ht_bool *d_ht_bool;
  uint32_t *final_walk_lens_d;
  // allocate GPU  memory
  ERROR_CHECK(Malloc(&prefix_ht_size_d, sizeof(uint32_t) * max_slice_size));
  ERROR_CHECK(Malloc(&cid_d, sizeof(uint64_t) * max_slice_size));
  ERROR_CHECK(Malloc(&ctg_seq_offsets_d, sizeof(uint32_t) * max_slice_size));
  ERROR_CHECK(Malloc(&reads_l_offset_d, sizeof(uint32_t) * max_l_rds_its));
  ERROR_CHECK(Malloc(&reads_r_offset_d, sizeof(uint32_t) * max_r_rds_its));
  ERROR_CHECK(Malloc(&rds_l_cnt_offset_d, sizeof(uint32_t) * max_slice_size));
  ERROR_CHECK(Malloc(&rds_r_cnt_offset_d, sizeof(uint32_t) * max_slice_size));
  ERROR_CHECK(Malloc(&ctg_seqs_d, sizeof(char) * max_ctg_len_its));
  ERROR_CHECK(Malloc(&reads_left_d, sizeof(char) * max_read_size * max_l_rds_its));
  ERROR_CHECK(Malloc(&reads_right_d, sizeof(char) * max_read_size * max_r_rds_its));
  ERROR_CHECK(Malloc(&depth_d, sizeof(double) * max_slice_size));
  ERROR_CHECK(Malloc(&quals_right_d, sizeof(char) * max_read_size * max_r_rds_its));
  ERROR_CHECK(Malloc(&quals_left_d, sizeof(char) * max_read_size * max_l_rds_its));
  ERROR_CHECK(Malloc(&term_counts_d, sizeof(uint32_t) * 3));
  // one local hashtable for each thread, so total hash_tables equal to vec_size i.e. total contigs
  ERROR_CHECK(Malloc(&d_ht, sizeof(loc_ht) * (max_ht + 1)));  // one more for FULL last entry
  ERROR_CHECK(Memset(d_ht + max_ht, FULL, sizeof(loc_ht)));   // set the FULL record
  ERROR_CHECK(Malloc(&longest_walks_d, sizeof(char) * max_slice_size * max_walk_len));
  ERROR_CHECK(Malloc(&mer_walk_temp_d, (max_mer_len + max_walk_len) * sizeof(char) * max_slice_size));
  ERROR_CHECK(Malloc(&d_ht_bool, sizeof(loc_ht_bool) * (max_slice_size * max_walk_len + 1)));   // one more for FULL last entry
  ERROR_CHECK(Memset(d_ht_bool + (max_slice_size * max_walk_len), FULL, sizeof(loc_ht_bool)));  // set the FULL record
  ERROR_CHECK(Malloc(&final_walk_lens_d, sizeof(uint32_t) * max_slice_size));

  for (unsigned slice = 0; slice < iterations; slice++) {
    auto extensions_offset = max_slice_size * slice;
    auto num_extensions = tot_extensions - extensions_offset;
    if (num_extensions > max_slice_size) num_extensions = max_slice_size;
    assert(num_extensions > 0 && num_extensions <= max_slice_size);

    vector<CtgWithReads *>::const_iterator slice_iter = data_in.begin() + extensions_offset;
    auto this_slice_size = num_extensions;
    uint32_t vec_size = this_slice_size;
    uint64_t ctgs_offset_sum = 0;
    uint64_t prefix_ht_sum = 0;
    uint64_t reads_r_offset_sum = 0;
    uint64_t reads_l_offset_sum = 0;
    uint64_t read_l_index = 0, read_r_index = 0;
    int num_bad = 0;
    for (unsigned i = 0; i < this_slice_size; i++) {
      const CtgWithReads &temp_data = *slice_iter[i];  // slice_data[i];
      cid_h[i] = temp_data.cid;
      depth_h[i] = temp_data.depth;
      // convert string to c-string
      auto ctg_seq_size = temp_data.seq.size();
      if (ctg_seq_size > TOO_BIG || ctgs_offset_sum + ctg_seq_size > max_ctg_seqs_h)
        printf("Invalid ctg_seq_size=%ld my_rank=%d i=%d of %lu cid=%ld ctgs_offset_sum=%lu max_ctg_seqs_h=%lu\n", ctg_seq_size,
               my_rank, i, this_slice_size, temp_data.cid, ctgs_offset_sum, max_ctg_seqs_h);
      char *ctgs_ptr = ctg_seqs_h.get() + ctgs_offset_sum;
      memcpy(ctgs_ptr, temp_data.seq.c_str(), ctg_seq_size);
      ctgs_offset_sum += ctg_seq_size;
      ctg_seq_offsets_h[i] = ctgs_offset_sum;
      prefix_ht_sum += temp_data.max_reads * max_read_size;
      prefix_ht_size_h[i] = prefix_ht_sum;

      for (unsigned j = 0; j < temp_data.reads_left.size(); j++) {
        auto read_size = temp_data.reads_left[j].seq.size();
        auto qual_size = temp_data.reads_left[j].quals.size();
        if (read_size != qual_size || read_size > max_read_size || reads_l_offset_sum + max_read_size > max_reads_left_h) {
          if (!num_bad)
            fprintf(stderr,
                    "WARN: myrank=%d: Invalid reads_left[%u of %lu].seq read_size=%lu qual_size=%lu i=%u of %lu, cid=%lu "
                    "reads_l_offset_sum=%lu max_reads_left_h=%u max_l_rds_its=%lu slice=%u of %lu\n",
                    my_rank, j, temp_data.reads_left.size(), read_size, qual_size, i, this_slice_size, temp_data.cid,
                    reads_l_offset_sum, max_reads_left_h, max_l_rds_its, slice, iterations);
          num_bad++;
          continue;
        }
        char *reads_l_ptr = reads_left_h.get() + reads_l_offset_sum;
        char *quals_l_ptr = quals_left_h.get() + reads_l_offset_sum;

        memcpy(reads_l_ptr, temp_data.reads_left[j].seq.c_str(), read_size);
        // quals offsets will be same as reads offset because quals and reads have same length
        memcpy(quals_l_ptr, temp_data.reads_left[j].quals.c_str(), qual_size);
        reads_l_offset_sum += read_size;
        reads_l_offset_h[read_l_index] = reads_l_offset_sum;
        read_l_index++;
      }
      rds_l_cnt_offset_h[i] = read_l_index;  // running sum of left reads count
      if (read_l_index > max_l_rds_its || reads_l_offset_sum > TOO_BIG)
        fprintf(stderr, "WARN: my_rank=%d read_l_index=%lu > max_l_rds_its=%lu reads_l_offset_sum=%lu\n", my_rank, read_l_index,
                max_l_rds_its, reads_l_offset_sum);

      for (unsigned j = 0; j < temp_data.reads_right.size(); j++) {
        auto read_size = temp_data.reads_right[j].seq.size();
        auto qual_size = temp_data.reads_right[j].quals.size();
        if (read_size != qual_size || read_size > max_read_size || reads_r_offset_sum + max_read_size > max_reads_right_h) {
          if (!num_bad)
            fprintf(stderr,
                    "WARN: myrank=%d: Invalid reads_right[%u of %lu].seq read_size=%lu quals_size=%lu i=%u of %lu, cid=%lu "
                    "reads_r_offset_sum=%lu max_reads_right_h=%u max_r_rds_its=%lu slice=%u of %lu\n",
                    my_rank, j, temp_data.reads_right.size(), read_size, qual_size, i, this_slice_size, temp_data.cid,
                    reads_r_offset_sum, max_reads_right_h, max_r_rds_its, slice, iterations);
          num_bad++;
          continue;
        }
        char *reads_r_ptr = reads_right_h.get() + reads_r_offset_sum;
        char *quals_r_ptr = quals_right_h.get() + reads_r_offset_sum;
        memcpy(reads_r_ptr, temp_data.reads_right[j].seq.c_str(), read_size);
        // quals offsets will be same as reads offset because quals and reads have same length
        memcpy(quals_r_ptr, temp_data.reads_right[j].quals.c_str(), qual_size);
        reads_r_offset_sum += read_size;
        reads_r_offset_h[read_r_index] = reads_r_offset_sum;
        read_r_index++;
      }
      rds_r_cnt_offset_h[i] = read_r_index;  // running sum of right reads count
      if (read_r_index > max_r_rds_its || reads_r_offset_sum > TOO_BIG)
        fprintf(stderr, "WARN: my_rank=%d read_r_index=%lu > max_r_rds_its=%lu reads_r_offset_sum=%lu\n", my_rank, read_r_index,
                max_r_rds_its, reads_r_offset_sum);

    }                                        // data packing for loop ends
    if (num_bad) fprintf(stderr, "WARN: myrank=%d found %d excess sized reads from line %d\n", my_rank, num_bad, debug_line);

    uint32_t total_r_reads_slice = read_r_index;
    uint32_t total_l_reads_slice = read_l_index;
    for (int i = 0; i < 3; i++) {
      term_counts_h[i] = 0;
    }

    ERROR_CHECK(Memcpy(prefix_ht_size_d, prefix_ht_size_h.get(), sizeof(uint32_t) * vec_size, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(cid_d, cid_h.get(), sizeof(uint64_t) * vec_size, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(ctg_seq_offsets_d, ctg_seq_offsets_h.get(), sizeof(uint32_t) * vec_size, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(reads_l_offset_d, reads_l_offset_h.get(), sizeof(uint32_t) * total_l_reads_slice, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(reads_r_offset_d, reads_r_offset_h.get(), sizeof(uint32_t) * total_r_reads_slice, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(rds_l_cnt_offset_d, rds_l_cnt_offset_h.get(), sizeof(uint32_t) * vec_size, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(rds_r_cnt_offset_d, rds_r_cnt_offset_h.get(), sizeof(uint32_t) * vec_size, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(ctg_seqs_d, ctg_seqs_h.get(), sizeof(char) * ctgs_offset_sum, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(depth_d, depth_h.get(), sizeof(double) * vec_size, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(term_counts_d, term_counts_h.get(), sizeof(uint32_t) * 3, MemcpyHostToDevice));

    // copy right reads
    ERROR_CHECK(Memcpy(reads_right_d, reads_right_h.get(), sizeof(char) * reads_r_offset_sum, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(quals_right_d, quals_right_h.get(), sizeof(char) * reads_r_offset_sum, MemcpyHostToDevice));

    // call kernel here, one thread per contig
    unsigned total_threads = vec_size * WARP_SIZE;  // we need one warp (32/64 threads) per extension, vec_size = extensions
    unsigned thread_per_blk = 512;
    unsigned blocks = (total_threads + thread_per_blk) / thread_per_blk;

    int64_t sum_ext = 0, num_walks = 0;
    uint32_t qual_offset_ = qual_offset;
    LaunchKernel(iterative_walks_kernel, blocks, thread_per_blk, cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_right_d, quals_right_d,
                 reads_r_offset_d, rds_r_cnt_offset_d, depth_d, d_ht, prefix_ht_size_d, d_ht_bool, mer_len, max_mer_len,
                 term_counts_d, num_walks, max_walk_len, sum_ext, max_read_size, max_read_count, qual_offset_, longest_walks_d,
                 mer_walk_temp_d, final_walk_lens_d, vec_size);

    // perform revcomp of contig sequences and launch kernel with left reads,

    for (unsigned j = 0; j < vec_size; j++) {
      int size_lst;
      char *curr_seq;
      char *curr_seq_rc;
      if (j == 0) {
        size_lst = ctg_seq_offsets_h[j];
        curr_seq = ctg_seqs_h.get();
        curr_seq_rc = ctgs_seqs_rc_h.get();
      } else {
        size_lst = ctg_seq_offsets_h[j] - ctg_seq_offsets_h[j - 1];
        curr_seq = ctg_seqs_h.get() + ctg_seq_offsets_h[j - 1];
        curr_seq_rc = ctgs_seqs_rc_h.get() + ctg_seq_offsets_h[j - 1];
      }
      revcomp(curr_seq, curr_seq_rc, size_lst);
    }
    ERROR_CHECK(Memcpy(longest_walks_r_h.get() + slice * max_walk_len * max_slice_size, longest_walks_d,
                       sizeof(char) * vec_size * max_walk_len, MemcpyDeviceToHost));
    ERROR_CHECK(Memcpy(final_walk_lens_r_h.get() + slice * max_slice_size, final_walk_lens_d, sizeof(int32_t) * vec_size,
                       MemcpyDeviceToHost));

    // cpying rev comped ctgs to device on same memory as previous ctgs
    ERROR_CHECK(Memcpy(ctg_seqs_d, ctgs_seqs_rc_h.get(), sizeof(char) * ctgs_offset_sum, MemcpyHostToDevice));

    // copy left reads
    ERROR_CHECK(Memcpy(reads_left_d, reads_left_h.get(), sizeof(char) * reads_l_offset_sum, MemcpyHostToDevice));
    ERROR_CHECK(Memcpy(quals_left_d, quals_left_h.get(), sizeof(char) * reads_l_offset_sum, MemcpyHostToDevice));

    LaunchKernel(iterative_walks_kernel, blocks, thread_per_blk, cid_d, ctg_seq_offsets_d, ctg_seqs_d, reads_left_d, quals_left_d,
                 reads_l_offset_d, rds_l_cnt_offset_d, depth_d, d_ht, prefix_ht_size_d, d_ht_bool, mer_len, max_mer_len,
                 term_counts_d, num_walks, max_walk_len, sum_ext, max_read_size, max_read_count, qual_offset_, longest_walks_d,
                 mer_walk_temp_d, final_walk_lens_d, vec_size);

    ERROR_CHECK(Memcpy(longest_walks_l_h.get() + slice * max_walk_len * max_slice_size, longest_walks_d,
                       sizeof(char) * vec_size * max_walk_len, MemcpyDeviceToHost));  // copy back left walks
    ERROR_CHECK(Memcpy(final_walk_lens_l_h.get() + slice * max_slice_size, final_walk_lens_d, sizeof(int32_t) * vec_size,
                       MemcpyDeviceToHost));
  }  // the for loop over all slices ends here

  // once all the alignments are on cpu, then go through them and stitch them with contigs in front and back.

  for (unsigned j = 0; j < iterations; j++) {
    auto extensions_offset = max_slice_size * j;
    auto num_extensions = tot_extensions - extensions_offset;
    if (num_extensions > max_slice_size) num_extensions = max_slice_size;
    assert(num_extensions > 0 && num_extensions <= max_slice_size);

    int loc_size = num_extensions;

    for (int i = 0; i < loc_size; i++) {
      auto &left_len = final_walk_lens_l_h[j * max_slice_size + i];
      if (left_len > TOO_BIG) {
        fprintf(stderr, "WARN: myrank=%d found TOO_BIG %d left_len line j=%d i=%d of %d\n", my_rank, left_len, j, i, loc_size);
        continue;
      }
      if (left_len > 0) {
        string left(longest_walks_l_h.get() + j * max_slice_size * max_walk_len + max_walk_len * i, left_len);
        string left_rc = revcomp(left);
        data_in[j * max_slice_size + i]->seq.insert(0, left_rc);
      }
      auto &right_len = final_walk_lens_r_h[j * max_slice_size + i];
      if (right_len > TOO_BIG) {
        fprintf(stderr, "WARN: myrank=%d found TOO_BIG %d right_len line j=%d i=%d of %d\n", my_rank, right_len, j, i, loc_size);
        continue;
      }
      if (right_len > 0) {
        string right(longest_walks_r_h.get() + j * max_slice_size * max_walk_len + max_walk_len * i,
                     final_walk_lens_r_h[j * max_slice_size + i]);
        data_in[j * max_slice_size + i]->seq += right;
      }
    }
  }

  ERROR_CHECK(Free(prefix_ht_size_d));
  ERROR_CHECK(Free(term_counts_d));
  ERROR_CHECK(Free(cid_d));
  ERROR_CHECK(Free(ctg_seq_offsets_d));
  ERROR_CHECK(Free(reads_l_offset_d));
  ERROR_CHECK(Free(reads_r_offset_d));
  ERROR_CHECK(Free(rds_l_cnt_offset_d));
  ERROR_CHECK(Free(rds_r_cnt_offset_d));
  ERROR_CHECK(Free(ctg_seqs_d));
  ERROR_CHECK(Free(reads_left_d));
  ERROR_CHECK(Free(reads_right_d));
  ERROR_CHECK(Free(depth_d));
  ERROR_CHECK(Free(quals_right_d));
  ERROR_CHECK(Free(quals_left_d));
  ERROR_CHECK(Free(d_ht));
  ERROR_CHECK(Free(longest_walks_d));
  ERROR_CHECK(Free(mer_walk_temp_d));
  ERROR_CHECK(Free(d_ht_bool));
  ERROR_CHECK(Free(final_walk_lens_d));
}
