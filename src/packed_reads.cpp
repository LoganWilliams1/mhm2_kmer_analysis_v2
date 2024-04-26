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
#include "packed_reads.hpp"

#include <iostream>
// Not available in gcc <= 7
// #include <charconv>
#include <fcntl.h>
#include <unistd.h>

#include <upcxx/upcxx.hpp>

#include "fastq.hpp"
#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "upcxx_utils/promise_collectives.hpp"

using std::string;
using std::string_view;
using std::to_string;

using upcxx::rank_me;
using upcxx::rank_n;

using namespace upcxx_utils;

using std::max;

PackedRead::PackedRead()
    : is_allocated(0)
    , read_id(0)
    , read_len(0)
    , bytes(nullptr) {
  // DBG("Constructed empty ", (void *)this, "\n");
}

PackedRead::PackedRead(const string &id_str, string_view seq, string_view quals, int qual_offset, PackedReads *packed_reads)
    : is_allocated(USE_PACKED_READS_LINEAR_ALLOCATOR && packed_reads) {
  read_id = strtol(id_str.c_str() + 1, nullptr, 10) + 1;
  assert(labs(read_id) < MAX_READ_ID);
  if (id_str[id_str.length() - 1] == '1') read_id *= -1;
  // read_id = strtol(id_str.c_str() + 1, nullptr, 10);
  // this uses from_chars because it's the fastest option out there
  // auto res = std::from_chars(id_str.data() + 2, id_str.data() + id_str.size() - 2, read_id);
  // if (res.ec != std::errc()) DIE("Failed to convert string to int64_t, ", res.ec);
  // negative if first of the pair
  // if (id_str[id_str.length() - 1] == '1') read_id *= -1;
  // packed is same length as sequence. Set first 3 bits to represent A,C,G,T,N
  // set next five bits to represent quality (from 0 to 32). This doesn't cover the full quality range (only up to 32)
  // but it's all we need since once quality is greater than the qual_thres (20), we treat the base as high quality
  if (is_allocated) {
    bytes = packed_reads->allocate_read(seq.length());
  } else {
    bytes = new unsigned char[seq.length()];
  }

  for (unsigned i = 0, len = seq.length(); i < len; i++) {
    switch (seq[i]) {
      case 'A':
      case 'a': bytes[i] = 0; break;
      case 'C':
      case 'c': bytes[i] = 1; break;
      case 'G':
      case 'g': bytes[i] = 2; break;
      case 'T':
      case 't': bytes[i] = 3; break;
      case 'N':
      case 'n': bytes[i] = 4; break;
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
      case 'V': bytes[i] = 4; break;
      default:
        DIE("Illegal char in comp nucleotide of '", ((seq[i] >= 32 && seq[i] <= 126) ? seq[i] : ' '), "' int=", (int)seq[i], "\n");
    }
    bytes[i] |= ((unsigned char)std::min(quals[i] - qual_offset, 31) << 3);
  }
  read_len = (uint16_t)seq.length();
  // DBG("Constructed ", (void *)this, " is_allocated=", is_allocated, " len=", read_len, " id=", read_id, "\n");
}

PackedRead::PackedRead(const PackedRead &copy, PackedReads *packed_reads)
    : is_allocated((USE_PACKED_READS_LINEAR_ALLOCATOR && packed_reads) ? 1 : 0)
    , read_id(copy.read_id)
    , read_len(copy.read_len)
    , bytes(nullptr) {
  if (is_allocated) {
    bytes = packed_reads->allocate_read(read_len);
  } else {
    bytes = new unsigned char[read_len];
  }
  // DBG("Constructed copy ", (void *)this, " from ", (void *)&copy, " allocate=", is_allocated, " was_allocated=",
  // copy.is_allocated,  " len=", read_len, " id=", read_id, "\n");
  memcpy(bytes, copy.bytes, read_len);
}

PackedRead::PackedRead(const PackedRead &copy, bool no_new_allocation_ignored)
    : is_allocated(1)
    , read_id(copy.read_id)
    , read_len(copy.read_len)
    , bytes(copy.bytes) {
  // DBG("Constructed non-allocated ", (void *)this, " from ", (void *)&copy, " allocate=", is_allocated,   " was_allocated=",
  // copy.is_allocated, " len=", read_len, " id=", read_id, "\n");
}

PackedRead::PackedRead(PackedRead &&move, PackedReads *packed_reads)
    : is_allocated((USE_PACKED_READS_LINEAR_ALLOCATOR && packed_reads) ? 1 : 0)
    , read_id(move.read_id)
    , read_len(move.read_len)
    , bytes(nullptr) {
  if (is_allocated || move.is_allocated) {
    // copy
    bytes = is_allocated ? packed_reads->allocate_read(read_len) : new unsigned char[read_len];
    memcpy(bytes, move.bytes, read_len);
  } else {
    bytes = move.bytes;
  }
  // DBG("Constructed move ", (void *)this, " from ", (void *)&move, " allocate=", is_allocated, " was_allocated=",
  // move.is_allocated,   " len=", read_len, " id=", read_id, "\n");
  move.bytes = nullptr;
  move.clear();
}

PackedRead &PackedRead::operator=(const PackedRead &copy) {
  // DBG("copy assigned ", (void *)this, " from ", (void *)&copy, "\n");
  PackedRead pr(copy);
  std::swap(*this, pr);
  return *this;
}

PackedRead &PackedRead::operator=(PackedRead &&move) {
  // DBG("move assigned ", (void *)this, " from ", (void *)&move, "\n");
  PackedRead pr(std::move(move));
  std::swap(*this, pr);
  return *this;
}

PackedRead::~PackedRead() { clear(); }

void PackedRead::clear() {
  if (bytes && !is_allocated) delete[] bytes;
  bytes = nullptr;
  read_len = 0;
  read_id = 0;
}

void PackedRead::unpack(string &read_id_str, string &seq, string &quals, int qual_offset) const {
  assert(bytes != nullptr);
  char pair_id = (read_id < 0 ? '1' : '2');
  read_id_str = "@r" + to_string(labs(read_id)) + '/' + pair_id;
  seq.resize(read_len);
  quals.resize(read_len);
  for (int i = 0; i < read_len; i++) {
    seq[i] = nucleotide_map[bytes[i] & 7];
    quals[i] = qual_offset + (bytes[i] >> 3);
  }
  assert(seq.length() == read_len);
  assert(quals.length() == read_len);
}

void PackedRead::unpack_seq(string &seq) const {
  assert(bytes != nullptr);
  seq.resize(read_len);
  for (int i = 0; i < read_len; i++) {
    seq[i] = nucleotide_map[bytes[i] & 7];
  }
  assert(seq.length() == read_len);
}

int64_t PackedRead::get_id() { return read_id; }

string PackedRead::get_str_id() {
  assert(labs(read_id) < MAX_READ_ID);
  char pair_id = (read_id < 0 ? '1' : '2');
  return "@r" + to_string(labs(read_id)) + '/' + pair_id;
}

int64_t PackedRead::to_packed_id(const string &id_str) {
  assert(id_str[0] == '@');
  int64_t read_id = strtol(id_str.c_str() + 2, nullptr, 10);
  if (id_str[id_str.length() - 1] == '1') read_id *= -1;
  assert(labs(read_id) < MAX_READ_ID);
  return read_id;
}

uint16_t PackedRead::get_read_len() const { return read_len; }

unsigned char *PackedRead::get_raw_bytes() { return bytes; }

PackedReads::PackedReads(int qual_offset, const string &fname, bool str_ids)
    : allocator(ALLOCATION_BLOCK_SIZE)
    , qual_offset(qual_offset)
    , fname(fname)
    , str_ids(str_ids) {}

PackedReads::PackedReads(int qual_offset, PackedReadsContainer &new_packed_reads, const string &fname)
    : allocator(ALLOCATION_BLOCK_SIZE)
    , packed_reads{}
    , index(0)
    , qual_offset(qual_offset)
    , fname(fname)
    , str_ids(false) {
  max_read_len = 0;
  // assert(!packed_reads.size());
  LOG("Constructed PackedReads ", (void *)this, ". Transferring ", new_packed_reads.size(), " to allocated storage fname=", fname,
      "\n");
  uint64_t num_bases = 0;
  for (auto &packed_read : new_packed_reads) num_bases += packed_read.get_read_len();
  reserve(new_packed_reads.size(), num_bases);
  for (const auto &packed_read : new_packed_reads) {
    packed_reads.emplace_back(packed_read, this);  // copy to this PackedReads using the allocator
    max_read_len = max((unsigned)packed_read.get_read_len(), max_read_len);
  }
}

PackedReads::~PackedReads() { clear(); }

bool PackedReads::get_next_read(string &id, string &seq, string &quals) {
  assert(qual_offset == 33 || qual_offset == 64);
  if (index == packed_reads.size()) return false;
  packed_reads[index].unpack(id, seq, quals, qual_offset);
  if (str_ids) id = read_id_idx_to_str[index];
  index++;
  return true;
}

unsigned char *PackedReads::allocate_read(uint16_t read_len) { return (unsigned char *)allocator.Allocate(read_len); }

uint64_t PackedReads::get_read_index() const { return index; }

void PackedReads::get_read(uint64_t index, string &id, string &seq, string &quals) const {
  if (index >= packed_reads.size()) DIE("Invalid get_read(", index, ") - size=", packed_reads.size());
  packed_reads[index].unpack(id, seq, quals, qual_offset);
  if (str_ids) id = read_id_idx_to_str[index];
}

void PackedReads::get_read_seq(uint64_t index, string &seq) const {
  if (index >= packed_reads.size()) DIE("Invalid get_read(", index, ") - size=", packed_reads.size());
  packed_reads[index].unpack_seq(seq);
}

void PackedReads::reset() { index = 0; }

void PackedReads::clear() {
  index = 0;
  fname.clear();
  PackedReadsContainer().swap(packed_reads);
  allocator.Reset();
  if (str_ids) deque<string>().swap(read_id_idx_to_str);
}

void PackedReads::reserve(uint64_t num_reads, uint64_t num_bases) {
  LOG("Reserving ", num_reads, " reads and ", num_bases, " bases: ", get_size_str(sizeof(PackedRead) * num_reads + num_bases),
      " for ", get_basename(fname), "\n");
  packed_reads.reserve(num_reads);
  allocator.Reserve(num_bases);
}

string PackedReads::get_fname() const { return fname; }

unsigned PackedReads::get_max_read_len() const { return max_read_len; }

void PackedReads::set_max_read_len() {
  max_read_len = 0;
  for (auto &packed_read : packed_reads) {
    max_read_len = max((unsigned)packed_read.get_read_len(), max_read_len);
  }
}

int64_t PackedReads::get_local_num_reads() const { return packed_reads.size(); }

int64_t PackedReads::get_total_local_num_reads(const PackedReadsList &packed_reads_list) {
  int64_t total_local_num_reads = 0;
  for (const PackedReads *pr : packed_reads_list) {
    total_local_num_reads += pr->get_local_num_reads();
  }
  return total_local_num_reads;
}

int PackedReads::get_qual_offset() { return qual_offset; }

void PackedReads::add_read(const string &read_id, const string &seq, const string &quals) {
  packed_reads.emplace_back(read_id, seq, quals, qual_offset, this);
  if (str_ids) {
    read_id_idx_to_str.push_back(read_id);
    name_bytes += sizeof(string) + read_id.size();
  }
  max_read_len = max(max_read_len, (unsigned)seq.length());
  bases += seq.length();
}

upcxx::future<> PackedReads::load_reads_nb(const string &adapter_fname) {
  // first estimate the number of records
  size_t tot_bytes_read = 0;
  int64_t num_records = 0;
  FastqReader &fqr = FastqReaders::get(fname);
  fqr.advise(true);
  string id, seq, quals;
  for (num_records = 0; num_records < 20000; num_records++) {
    size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
    if (!bytes_read) break;
    tot_bytes_read += bytes_read;
  }
  int64_t reserve_records = 0;
  int64_t estimated_records = 0;
  if (num_records > 0) {
    int64_t bytes_per_record = tot_bytes_read / num_records;
    estimated_records = fqr.my_file_size() / bytes_per_record;
    reserve_records = estimated_records * 1.10 + 10000;  // reserve more so there is not a big reallocation if it is under
  }
  fqr.reset();
  ProgressBar progbar(fqr.my_file_size(), "Loading reads from " + fname + " " + get_size_str(fqr.my_file_size()));
  tot_bytes_read = 0;
  int lines = 0;
  while (true) {
    size_t bytes_read = fqr.get_next_fq_record(id, seq, quals);
    if (!bytes_read) break;
    tot_bytes_read += bytes_read;
    progbar.update(tot_bytes_read);
    add_read(id, seq, quals);
  }
  DBG("Done loading reads in ", fname, "\n");
  fqr.advise(false);
  FastqReaders::close(fname);
  auto fut = progbar.set_done();
  int64_t underestimate = estimated_records - packed_reads.size();
  if (underestimate < 0 && reserve_records < packed_reads.size())
    LOG("NOTICE Underestimated by ", -underestimate, " estimated ", estimated_records, " found ", packed_reads.size(), "\n");
  auto &pr = Timings::get_promise_reduce();
  auto all_under_estimated_fut = pr.reduce_one(underestimate < 0 ? 1 : 0, upcxx::op_fast_add, 0);
  auto all_estimated_records_fut = pr.reduce_one(estimated_records, upcxx::op_fast_add, 0);
  auto all_num_records_fut = pr.reduce_one(packed_reads.size(), upcxx::op_fast_add, 0);
  auto all_num_bases_fut = pr.reduce_one(bases, upcxx::op_fast_add, 0);
  auto fut_report =
      when_all(fut, all_under_estimated_fut, all_estimated_records_fut, all_num_records_fut, all_num_bases_fut)
          .then([max_read_len = this->max_read_len](int64_t all_under_estimated, int64_t all_estimated_records,
                                                    int64_t all_num_records, int64_t all_num_bases) {
            SLOG_VERBOSE("Loaded ", all_num_records, " reads (estimated ", all_estimated_records, " with ", all_under_estimated,
                         " ranks underestimated) max_read=", max_read_len, " tot_bases=", all_num_bases, "\n");
          });
  Timings::set_pending(fut_report);
  return fut_report;
}

void PackedReads::load_reads_list(PackedReadsList &packed_reads_list, const string &adapter_fname) {
  BarrierTimer timer(__FILEFUNC__);
  upcxx::future<> all_done = upcxx::make_future();
  vector<string> file_names;
  for (auto pr : packed_reads_list) {
    file_names.push_back(pr->fname);
  }
  FastqReaders::open_all(file_names);
  for (auto pr : packed_reads_list) {
    upcxx::discharge();
    upcxx::progress();
    auto fut = pr->load_reads_nb(adapter_fname);
    all_done = when_all(all_done, fut);
  }
  FastqReaders::close_all();
  all_done.wait();
}

void PackedReads::load_reads(const string &adapter_fname) {
  BarrierTimer timer(__FILEFUNC__);
  auto fut = load_reads_nb(adapter_fname);
  Timings::wait_pending();
  fut.wait();
}

void PackedReads::report_size() {
  LOG("PackedReads num_reads=", get_local_num_reads(), " num_bases=", bases, " name_bytes=", name_bytes, ": ",
      get_size_str(sizeof(PackedRead) * get_local_num_reads() + bases + name_bytes), " for ", upcxx_utils::get_basename(fname),
      "\n");
  auto &pr = Timings::get_promise_reduce();
  auto all_num_records_fut = pr.reduce_one(packed_reads.size(), upcxx::op_fast_add, 0);
  auto all_num_bases_fut = pr.reduce_one(bases, upcxx::op_fast_add, 0);
  auto all_num_names_fut = pr.reduce_one(name_bytes, upcxx::op_fast_add, 0);
  auto min_rank_fut = pr.reduce_one(get_local_num_reads() > 0 ? rank_me() : rank_n(), upcxx::op_fast_min, 0);
  auto max_rank_fut = pr.reduce_one(get_local_num_reads() > 0 ? rank_me() : -1, upcxx::op_fast_max, 0);
  auto fut =
      when_all(Timings::get_pending(), all_num_records_fut, all_num_bases_fut, all_num_names_fut, min_rank_fut, max_rank_fut)
          .then([&self = *this](size_t all_num_records, size_t all_num_bases, size_t all_num_names, int min_rank, int max_rank) {
            auto num_ranks = max_rank - min_rank + 1;
            auto sz = all_num_records * sizeof(PackedRead) + all_num_bases + all_num_names;
            SLOG_VERBOSE("Loaded ", upcxx_utils::get_basename(self.fname), ": reads=", all_num_records,
                         " tot_bases=", all_num_bases, " names=", get_size_str(all_num_names), " PackedReads: ", get_size_str(sz),
                         " active_ranks=", num_ranks, " (", min_rank, "-", max_rank, ") avg=", get_size_str(sz / num_ranks), "\n");
          });
  Timings::set_pending(fut);  // include in reports
}

int64_t PackedReads::get_local_bases() const { return bases; }
upcxx::future<uint64_t> PackedReads::fut_get_bases() const {
  return upcxx_utils::Timings::get_promise_reduce().reduce_one(bases, upcxx::op_fast_add, 0);
}
uint64_t PackedReads::get_bases() const {
  auto fut = fut_get_bases();
  Timings::wait_pending();
  return fut.wait();
}

PackedRead &PackedReads::operator[](int index) {
  if (index >= packed_reads.size()) DIE("Array index out of bound ", index, " >= ", packed_reads.size());
  return packed_reads[index];
}

uint64_t PackedReads::estimate_num_kmers(unsigned kmer_len, PackedReadsList &packed_reads_list) {
  BarrierTimer timer(__FILEFUNC__);
  int64_t num_kmers = 0;
  int64_t num_reads = 0;
  int64_t tot_num_reads = PackedReads::get_total_local_num_reads(packed_reads_list);
  ProgressBar progbar(tot_num_reads, "Scanning reads to estimate number of kmers");
  for (auto packed_reads : packed_reads_list) {
    packed_reads->reset();
    string id, seq, quals;
    for (int i = 0; i < 50000; i++) {
      if (!packed_reads->get_next_read(id, seq, quals)) break;
      progbar.update();
      // do not read the entire data set for just an estimate
      if (seq.length() < kmer_len) continue;
      num_kmers += seq.length() - kmer_len + 1;
      num_reads++;
    }
  }
  auto fut = progbar.set_done();
  DBG("This rank processed ", num_reads, " reads, and found ", num_kmers, " kmers\n");
  auto &pr = Timings::get_promise_reduce();
  auto fut_all_num_reads = pr.reduce_one(num_reads, op_fast_add, 0);
  auto fut_all_tot_num_reads = pr.reduce_one(tot_num_reads, op_fast_add, 0);
  auto fut_all_num_kmers = pr.reduce_all(num_kmers, op_fast_add);

  auto fut_log = when_all(Timings::get_pending(), fut_all_num_reads, fut_all_tot_num_reads, fut_all_num_kmers, fut)
                     .then([](int64_t all_num_reads, int64_t all_tot_num_reads, int64_t all_num_kmers) {
                       SLOG_VERBOSE("Processed ", perc_str(all_num_reads, all_tot_num_reads), " reads, and estimated a maximum of ",
                                    (all_num_reads > 0 ? all_num_kmers * (all_tot_num_reads / all_num_reads) : 0), " kmers\n");
                     });
  Timings::set_pending(fut_log);
  Timings::wait_pending();
  fut_log.wait();
  return num_reads > 0 ? num_kmers * tot_num_reads / num_reads : 0;
}
