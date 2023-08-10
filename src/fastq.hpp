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

#include <map>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/timers.hpp"
#include "zstr.hpp"

using std::shared_ptr;
using std::string;
using std::to_string;

using upcxx::dist_object;
using upcxx::promise;
using upcxx::rank_me;
using upcxx::rank_n;

using upcxx_utils::AsyncTimer;
using upcxx_utils::BarrierTimer;
using upcxx_utils::get_basename;
using upcxx_utils::IntermittentTimer;
using upcxx_utils::Timer;

#define INT_CEIL(numerator, denominator) (((numerator)-1) / (denominator) + 1)
#define BUF_SIZE 2047

class FastqReader {
  string fname;
  // FIXME BGZF std::unique_ptr<zstr::base_ifstream> in;
  std::unique_ptr<ifstream> in;
  promise<> know_file_size;
  int64_t file_size;                         // may be bgzf_virtual_file_pointer int
  int64_t block_start = 0, block_size = -1;  // raw offsets. block_start <= start_read;  end_read<= block_start + block_size
  int64_t start_read;                        // may be bgzf_virtual_file_pointer int
  int64_t end_read;                          // may be bgzf_virtual_file_pointer int
  unsigned max_read_len;
  unsigned avg_bytes_per_read;
  uint64_t num_reads, num_pairs, num_bases;
  int subsample_pct = 100;
  string buf;
  int qual_offset;
  int read_count = 0;  // used in subsample
  shared_ptr<FastqReader> fqr2;
  bool first_file;
  bool _is_paired;        // file was declared as paired by the user
  bool _is_interleaved;   // 1 file that also _is_paired
  bool _fix_paired_name;  // Issue124 - file contains identical names of paired reads, so fix each paired read name to be unique
                          // on-the-fly
  bool _first_pair;       // alternate for pair1 (first_pair) and pair2 (!first_pair)
  bool _trim_comment;
  bool _is_bgzf;
  bool _first_read = true;

  IntermittentTimer io_t, read_io_t;
  struct PromStartStop {
    promise<int64_t> start_prom, stop_prom;
    upcxx::future<> set(FastqReader &fqr) {
      auto set_start = start_prom.get_future().then([&fqr](int64_t start) {
        DBG("Set start_read at ", start, " on ", fqr.fname, "\n");
        fqr.start_read = start;
        return start;
      });
      auto set_end = stop_prom.get_future().then([&fqr](int64_t stop) {
        DBG("Set end_read at ", stop, " on ", fqr.fname, "\n");
        fqr.end_read = stop;
        return stop;
      });
      return when_all(set_start, set_end, fqr.get_file_size(false)).then([&fqr](int64_t start, int64_t stop, int64_t sz) {
        LOG("Reading ", fqr.get_fname(), " from ", start, " to ", stop, "=", fqr.my_file_size(false), " of ", sz, "\n");
      });
    }
    upcxx::future<> get_future() {
      return when_all(start_prom.get_future(), stop_prom.get_future()).then([](int64_t start, int64_t stop) {});
    }
  };
  dist_object<PromStartStop> dist_prom;
  upcxx::future<> open_fut;

  void seekg(int64_t pos);

  inline static double overall_io_t = 0;

  static void rtrim(string &s);

  bool get_fq_name(string &header);

  bool is_sep(const char &sep) const;

  int64_t get_fptr_for_next_record(int64_t offset);

 public:
  FastqReader() = delete;                                    // no default constructor
  FastqReader(const string &_fname, upcxx::future<> first_wait = make_future(), bool is_second_file = false);
  FastqReader(const FastqReader &copy) = delete;             // no copy
  FastqReader(FastqReader &&move) = default;
  FastqReader &operator=(const FastqReader &copy) = delete;  // no copy
  FastqReader &operator=(FastqReader &&move) = delete;

  void set_subsample_pct(int pct) {
    assert(subsample_pct > 0 && subsample_pct <= 100);
    subsample_pct = pct;
    if (fqr2) fqr2->subsample_pct = pct;
  }

  // this happens within a separate thread
  upcxx::future<> continue_open();
  upcxx::future<> continue_open_default_per_rank_boundaries();

  ~FastqReader();
  void close();

  string get_fname() const;

  string get_ifstream_state() const;

  size_t my_file_size(bool include_file_2 = true) const;

  upcxx::future<int64_t> get_file_size(bool include_file2 = false) const;

  void advise(bool will_need);

  void set_block(int64_t start, int64_t size);

  size_t get_next_fq_record(string &id, string &seq, string &quals, bool wait_open = true);

  void set_max_read_len(int len);
  int get_max_read_len();

  double static get_io_time();

  void reset();

  upcxx::future<> get_open_fut() const { return open_fut; }

  bool is_paired() const { return _is_paired; }

  bool is_two_files() const {
    if (fqr2) assert(is_paired() && !is_interleaved());
    return (is_paired() && !is_interleaved());
  }

  bool is_interleaved() const {
    if (_is_interleaved) {
      assert(!fqr2);
      assert(_is_paired);
    }
    return _is_interleaved;
  }

  static upcxx::future<> set_matching_pair(FastqReader &fqr1, FastqReader &fqr2, dist_object<PromStartStop> &dist_start_stop1,
                                           dist_object<PromStartStop> &dist_start_stop2);
  void seek_start();
  int64_t tellg();
  bool is_open();
  bool is_open(bool wait_for_open);
  FastqReader &get_fqr2();
  void set_avg_bytes_per_read(unsigned bytes);
  unsigned get_avg_bytes_per_read() const;
  uint64_t get_num_bases() const { return num_bases + (fqr2 ? fqr2->get_num_bases() : 0); }
  uint64_t get_num_reads() const { return num_reads + (fqr2 ? fqr2->get_num_reads() : 0); }
  uint64_t get_num_pairs() const { return num_pairs; }
  bool is_first_file() const { return (_is_paired && fqr2) || _is_interleaved || !_is_paired; }

};  // FastqReader

class FastqReaders {
  // singleton class to hold as set of fastq readers open and re-usable
  using ShFastqReader = shared_ptr<FastqReader>;
  std::unordered_map<string, ShFastqReader> readers;

  FastqReaders();
  ~FastqReaders();

 public:
  static FastqReaders &getInstance();

  static bool is_open(const string fname, bool wait_for_open = true);

  static size_t get_open_file_size(const string fname, bool include_file_2 = false);

  static FastqReader &open(const string fname, int subsample_pct = 100, upcxx::future<> first_wait = make_future());

  // returns the total global size
  template <typename Container>
  static size_t open_all(Container &fnames, int subsample_pct = 100) {
    DBG("Open all ", fnames.size(), " files\n");
#ifdef NO_FASTQ_GLOBAL_BLOCKING
    return open_all_file_blocking(fnames, subsample_pct);
#else
    return open_all_global_blocking(fnames, subsample_pct);
#endif
  }

  // returns the total global size
  template <typename Container>
  static size_t open_all_file_blocking(Container &fnames, int subsample_pct = 100) {
    // every rank opens a partition of every file
    assert(subsample_pct > 0 && subsample_pct <= 100);
    auto fut_chain = make_future();
    size_t total_size = 0;
    AsyncTimer t_all("open_all_file_blocking all files");
    t_all.start();
    for (string &fname : fnames) {
      AsyncTimer t_finish_seek("open_all_file_blocking seeks on file:" + get_basename(fname));
      t_finish_seek.start();
      FastqReader &fqr = open(fname, subsample_pct);
      auto fut = fqr.get_open_fut().then([t_finish_seek]() { t_finish_seek.stop(); });
      fut_chain = when_all(fut_chain, fut, fqr.get_file_size(true).then([&total_size](int64_t sz) { total_size += sz; }));
    }
    fut_chain = fut_chain.then([t_all]() { t_all.stop(); });
    fut_chain.wait();
    DBG("All ", fnames.size(), " file (pairs) open.  total_size=", total_size, "\n");
    return total_size;
  }

  // returns the total global size
  template <typename Container>
  static size_t open_all_global_blocking(Container &fnames, int subsample_pct = 100) {
    BarrierTimer t(__FILEFUNC__);
    assert(!upcxx::in_progress());
    // opens only some files and reads a partition, as if the entire set of files is one single concatented file

    // all files need to be opened together.  Verify either all or none are open
    bool needs_blocking = false;
    size_t total_size = 0;
    for (string &fname : fnames) {
      if (!is_open(fname, true)) {
        DBG(fname, " is not yet open\n");
        needs_blocking = true;
      } else {
        total_size += get_open_file_size(fname, true);
      }
    }
    if (!needs_blocking) return total_size;  // all are open
    total_size = 0;                          // some are not open yet

    std::vector<promise<>> know_blocks(fnames.size());
    std::vector<int64_t> file_sizes(fnames.size());
    std::vector<int64_t> file_sizes1(fnames.size());
    assert(subsample_pct > 0 && subsample_pct <= 100);
    int filenum = 0;
    upcxx::future<> chain_fut = make_future();
    upcxx::future<> chain_open_fut = make_future();

    for (string &fname : fnames) {
      Timer open_t("Opening " + fname);
      if (is_open(fname)) {
        DBG("Closing open file ", fname, "\n");
        close(fname);
      }
      auto &fqr = open(fname, subsample_pct, know_blocks[filenum].get_future());
      if (!fqr.is_open()) {
        needs_blocking = true;
        AsyncTimer t_file_sizes("Get file sizes of " + get_basename(fqr.get_fname()));
        t_file_sizes.start();
        upcxx::future<> fut = when_all(fqr.get_file_size(true), fqr.get_file_size(false))
                                  .then([t_file_sizes, &fqr, &total_size, &file_size = file_sizes[filenum],
                                         &file_size1 = file_sizes1[filenum]](int64_t sz, int64_t sz1) {
                                    file_size = sz;
                                    file_size1 = sz1;
                                    total_size += sz;
                                    t_file_sizes.stop();
                                  });
        chain_fut = when_all(chain_fut, fut);
        chain_open_fut = when_all(chain_open_fut, fqr.get_open_fut());
      }
      filenum++;
    }

    AsyncTimer t_finish_seek("Finish global block seeks on all files");
    chain_fut = chain_fut.then([t_finish_seek, &fnames, &file_sizes, &file_sizes1, &total_size, &know_blocks]() {
      t_finish_seek.start();
      assert(fnames.size() == file_sizes.size());
      assert(file_sizes.size() == know_blocks.size());
      auto global_block = INT_CEIL(total_size, rank_n());
      auto my_global_start = global_block * rank_me();
      auto my_global_stop = my_global_start + global_block;
      int64_t my_read_offset = 0;
      int64_t my_read_remaining = global_block;
      DBG("Determining global blocks global_block=", global_block, " my_start=", my_global_start, " my_stop=", my_global_stop,
          "\n");
      for (int i = 0; i < file_sizes.size(); i++) {
        // Default is to not read this file
        auto file_start = file_sizes[i];
        auto file_stop = file_sizes[i];
        auto global_file_start = my_read_offset;
        auto global_file_stop = global_file_start + file_sizes[i];

        DBG("file[", i, "] size=", file_sizes[i], " read_offset/file_start=", my_read_offset, " file_stop=", global_file_stop,
            " remaining=", my_read_remaining, "\n");
        if (global_file_stop <= my_global_start) {
          DBG("file[", i, "] is completely before my block\n");
          assert(global_file_start <= my_global_start);
          assert(global_file_stop <= my_global_start);
          assert(global_file_start <= my_global_stop);
          assert(global_file_stop <= my_global_stop);
        } else if (my_global_stop <= global_file_start) {
          DBG("file[", i, "] is completely after my block. global_file_start=", global_file_start,
              " global_file_stop=", global_file_stop, "\n");
          assert(my_read_remaining == 0);
          assert(global_file_start >= my_global_start);
          assert(global_file_stop >= my_global_start);
          assert(global_file_start >= my_global_stop);
          assert(global_file_stop >= my_global_stop);
        } else {
          // file is at least partially within my block
          if (my_global_start <= global_file_start && global_file_stop <= my_global_stop) {
            // file is completely contained within my block
            DBG("file[", i, "] is completely contained within my block. global_file_start=", global_file_start,
                " global_file_stop=", global_file_stop, "\n");
            assert(global_file_start <= my_global_stop);
            file_start = 0;
            assert(file_stop == file_sizes[i]);
          } else if (global_file_start <= my_global_start && my_global_stop <= global_file_stop) {
            // my block is completely contained with this file
            DBG("file[", i, "] completely contains my block. ", global_file_start, " global_file_stop=", global_file_stop, "\n");
            file_start = my_global_start - global_file_start;
            file_stop = file_start + my_read_remaining;
          } else if (global_file_start <= my_global_start && my_global_start < global_file_stop &&
                     global_file_stop < my_global_stop) {
            // file crosses my block start and ends within my block
            DBG("file[", i, "] crosses just my block start. global_file_start=", global_file_start,
                " global_file_stop=", global_file_stop, "\n");
            assert(global_file_stop <= my_global_stop);
            file_start = my_global_start - global_file_start;
            assert(file_start < file_stop);
            assert(file_stop == file_sizes[i]);
          } else if (my_global_start < global_file_start && global_file_start < my_global_stop &&
                     my_global_stop <= global_file_stop) {
            // file starts within my block and crosses my block end
            DBG("file[", i, "] crosses just my block end. global_file_start=", global_file_start,
                " global_file_stop=", global_file_stop, "\n");
            assert(my_global_start <= global_file_start);
            assert(global_file_start + my_read_remaining <= global_file_stop);
            file_start = 0;
            file_stop = my_read_remaining;
            assert(file_stop <= file_sizes[i]);
          } else {
            DIE("Should not get here!");
          }
        }
        DBG("file_start=", file_start, " file_stop=", file_stop, " my_read_remaining=", my_read_remaining, "\n");
        assert(file_start <= file_stop);
        assert(file_start <= file_sizes[i]);
        assert(file_stop <= file_sizes[i]);
        int64_t file_read_len = file_stop - file_start;
        assert(file_read_len <= my_read_remaining);
        my_read_remaining -= file_read_len;
        assert(my_read_remaining >= 0);

        // set the partially open FQR block for this file
        auto it = getInstance().readers.find(fnames[i]);
        assert(it != getInstance().readers.end());
        if (it->second->is_two_files()) {
          // adjust blocks within both files
          int64_t fs1 = file_sizes1[i];
          int64_t fs2 = file_sizes[i] - fs1;
          int64_t file_start_1 = fs1, file_read_len_1 = 0, file_start_2 = fs2, file_read_len_2 = 0;
          if (file_start < file_sizes[i]) {
            double fs1_frac = 1. * fs1 / file_sizes[i];
            double fs2_frac = 1. * fs2 / file_sizes[i];
            file_start_1 = file_start * fs1_frac;
            file_read_len_1 = ((int64_t)(file_read_len * fs1_frac)) + 1;
            file_start_2 = file_start * fs2_frac;
            file_read_len_2 = ((int64_t)(file_read_len * fs2_frac)) + 1;
          }
          DBG("file_name=", fnames[i], " file_start=", file_start, " read_len=", file_read_len, " tot_file_size=", file_sizes[i],
              " fs1=", fs1, " fs2=", fs2, " start1=", file_start_1, " len1=", file_read_len_1, " start2=", file_start_2,
              " len2=", file_read_len_2, "\n");
          it->second->set_block(file_start_1, file_read_len_1);
          it->second->get_fqr2().set_block(file_start_2, file_read_len_2);
        } else {
          it->second->set_block(file_start, file_read_len);
          DBG("block for i=", i, " file=", fnames[i], " file_start=", file_start, " file_stop=", file_stop, " len=", file_read_len,
              " size=", file_sizes[i], "\n");
        }

        // notify FQR::continue_open to continue to open
        know_blocks[i].fulfill_anonymous(1);
        my_read_offset = global_file_stop;
      }
      t_finish_seek.stop();
    });

    // block and wait with variables still in scope
    chain_fut.wait();
    chain_open_fut.wait();
    t_finish_seek.initiate_start_reduction();
    t_finish_seek.initiate_stop_reduction();

    DBG("All ", fnames.size(), " file (pairs) open.  total_size=", total_size, "\n");
    return total_size;
  }

  static FastqReader &get(const string fname);

  static void close(const string fname);

  static void close_all();
};  // FastqReaders
