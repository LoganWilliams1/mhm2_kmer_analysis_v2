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

#include <array>
#include <deque>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <upcxx/upcxx.hpp>

using std::deque;
using std::max;
using std::string;
using std::string_view;
using std::to_string;
using std::unique_ptr;
using std::vector;

#ifndef USE_PACKED_READS_LINEAR_ALLOCATOR
#define USE_PACKED_READS_LINEAR_ALLOCATOR 1
#endif

#include "linear_allocator_pool.hpp"

class PackedReads;
class PackedRead {
  static inline const std::array<char, 5> nucleotide_map = {'A', 'C', 'G', 'T', 'N'};
  static const uint64_t MAX_READ_ID = 1ULL << 46;  // 70 trillion read pairs.
  // This class is packed into 16 bytes, a pointer + bitfield packed into 64-bits: is_allocated, read_id and read_len
  uint64_t is_allocated : 1;
  // the pair number is indicated in the read id - negative means pair 1, positive means pair 2
  int64_t read_id : 47;
  // the read is not going to be larger than 65536 in length, but possibly larger than 256
  uint64_t read_len : 16;
  // each cached read packs the nucleotide into 3 bits (ACGTN), and the quality score into 5 bits
  unsigned char *bytes;

  // overall, we expect the compression to be around 50%. E.g. a read of 150bp would be
  // 16+150=166 vs 13+300=313

 public:
  // default, move and (deep) copy constructors to support containers without unique_ptr overhead
  PackedRead();
  PackedRead(const string &id_str, string_view seq, string_view quals, int qual_offset, PackedReads * = nullptr);
  PackedRead(const PackedRead &copy, PackedReads * = nullptr);
  PackedRead(const PackedRead &reference, bool no_new_allocation);
  PackedRead(PackedRead &&Move, PackedReads * = nullptr);
  PackedRead &operator=(const PackedRead &copy);
  PackedRead &operator=(PackedRead &&move);

  ~PackedRead();

  void unpack(string &read_id_str, string &seq, string &quals, int qual_offset) const;
  void unpack_seq(string &seq) const;
  void clear();

  int64_t get_id();
  string get_str_id();
  static int64_t to_packed_id(const string &id_str);

  uint16_t get_read_len() const;
  unsigned char *get_raw_bytes();

  struct upcxx_serialization {
    template <typename Writer>
    static void serialize(Writer &writer, PackedRead const &packed_read) {
      int64_t bitpacked_fields;
      memcpy(&bitpacked_fields, &packed_read, sizeof(int64_t));
      writer.write(bitpacked_fields);
      for (int i = 0; i < packed_read.read_len; i++) {
        writer.write(packed_read.bytes[i]);
      }
    }

    template <typename Reader>
    static PackedRead *deserialize(Reader &reader, void *storage) {
      PackedRead *packed_read = new (storage) PackedRead();
      int64_t bitpacked_fields = reader.template read<int64_t>();
      memcpy((void *)packed_read, &bitpacked_fields, sizeof(int64_t));
      packed_read->is_allocated = 0;
      packed_read->bytes = new unsigned char[packed_read->read_len];
      for (int i = 0; i < packed_read->read_len; i++) {
        packed_read->bytes[i] = reader.template read<unsigned char>();
      }
      return packed_read;
    }
  };
};

using PackedReadsList = deque<PackedReads *>;
using PackedReadsContainer = vector<PackedRead>;
class PackedReads {
  const static size_t ALLOCATION_BLOCK_SIZE = 4 * 1024 * 1024;
  PackedReadsContainer packed_reads;
  LinearAllocatorPool allocator;
  // this is only used when we need to know the actual name of the original reads
  deque<string> read_id_idx_to_str;
  unsigned max_read_len = 0;
  uint64_t _index = 0;
  uint64_t bases = 0;
  uint64_t name_bytes = 0;
  int qual_offset;
  string fname;
  bool str_ids;
  bool reads_are_paired;

 public:
  PackedReads(int qual_offset, const string &fname, bool str_ids = false);
  PackedReads(int qual_offset, PackedReadsContainer &packed_reads, const string &fname);
  ~PackedReads();

  bool get_next_read(string &id, string &seq, string &quals);
  uint64_t get_read_index() const;
  void get_read(uint64_t i, string &id, string &seq, string &quals) const;
  string get_full_read_id(uint64_t i);
  void get_read_seq(uint64_t i, string &seq) const;

  PackedRead &operator[](int i);

  void reset();

  void clear();

  void reserve(uint64_t num_reads, uint64_t num_bases);

  string get_fname() const;

  unsigned get_max_read_len() const;

  void set_max_read_len();

  int64_t get_local_num_reads() const;

  static int64_t get_total_local_num_reads(const PackedReadsList &packed_reads_list);

  void add_read(const string &read_id, const string &seq, const string &quals, const string &orig_id = "");

  upcxx::future<> load_reads_nb(const string &adapter_fname);

  void load_reads(const string &adapter_fname);

  static void load_reads_list(PackedReadsList &, const string &adapter_fname);

  void report_size();

  bool is_paired();

  int64_t get_local_bases() const;

  upcxx::future<uint64_t> fut_get_bases() const;
  uint64_t get_bases() const;

  int get_qual_offset();

  static uint64_t estimate_num_kmers(unsigned kmer_len, PackedReadsList &packed_reads_list);

  unsigned char *allocate_read(uint16_t read_len);
};
