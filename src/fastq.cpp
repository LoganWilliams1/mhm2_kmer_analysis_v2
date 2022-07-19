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

#include "fastq.hpp"

#include <fcntl.h>
#include <unistd.h>

#include <iostream>
#include <memory>
#include <string>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/promise_collectives.hpp"
#include "upcxx_utils/thread_pool.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

#include "zstr.hpp"

using namespace upcxx_utils;
using std::make_shared;
using std::ostringstream;
using std::shared_ptr;
using std::string;

using upcxx::future;
using upcxx::local_team;
using upcxx::rank_me;

void FastqReader::rtrim(string &s) {
  auto pos = s.length() - 1;
  while (pos >= 0 && std::isspace(static_cast<unsigned char>(s[pos]))) pos--;
  s.resize(pos + 1);
}

bool FastqReader::get_fq_name(string &header) {
  if (header[0] != '@') {
    // WARN("unknown format ", header, " missing @\n");
    return false;
  }
  // trim off '@'
  header.erase(0, 1);
  // strip trailing spaces
  // header = std::regex_replace(header, std::regex("\\s+$"), std::string(""));
  rtrim(header);
  // convert if new illumina 2 format  or HudsonAlpha format
  unsigned len = header.length();
  if (len >= 3 && header[len - 2] != '/') {
    if (header[len - 2] == 'R') {
      // HudsonAlpha format  (@pair-R1,  @pair-R2)
      // replace with @pair/1 or @pair/2 */
      char rNum = header[len - 1];
      header[len - 3] = '/';
      header[len - 2] = rNum;
      header.resize(len - 1);
      return true;
    } else {
      // Latest Illumina header format
      auto end_pos = header.find('\t');
      if (end_pos == string::npos) {
        end_pos = header.find(' ');
        // no comment, just return without modification
        if (end_pos == string::npos) return true;
      }
      // check for @pair/1 @/pair2 vs @pair 1:Y:... @pair 2:Y:....
      if (end_pos > 3 && header[end_pos - 2] == '/' && (header[end_pos - 1] == '1' || header[end_pos - 1] == '2')) {
        // @pair/1 or @pair/2
        // truncate name at comment
        header.resize(end_pos);
        return true;
      }
      if ((len < end_pos + 7) || header[end_pos + 2] != ':' || header[end_pos + 4] != ':' || header[end_pos + 6] != ':' ||
          (header[end_pos + 1] != '1' && header[end_pos + 1] != '2')) {
        // unknown pairing format
        SWARN("unknown format ", header, " end pos ", end_pos, "\n");
        return false;
      }
      // @pair 1:Y:.:.: or @pair 2:Y:.:.:
      // replace with @pair/1 or @pair/2
      header[end_pos] = '/';
      header.resize(end_pos + 2);
    }
  }
  return true;
}

int64_t FastqReader::tellg() {
  return in->tellg();

}

void FastqReader::seekg(int64_t pos) {
  in->seekg(pos);
  assert(tellg() == pos);
}

int64_t FastqReader::get_fptr_for_next_record(int64_t offset) {
  // first record is the first record, include it.  Every other partition will be at least 1 full record after offset.
  // but read the first few lines anyway

  // eof - do not read anything
  if (offset >= file_size) return file_size;

  io_t.start();

  in->seekg(offset);
  if (!in->good()) DIE("Could not fseek in ", fname, " to ", offset, ": ", strerror(errno));
  
  if (offset != 0) {
    // skip first (likely partial) line after this offset to ensure we start at the beginning of a line
    std::getline(*in, buf);
  }
  if (buf.empty() && (in->eof() || in->fail())) {
    io_t.stop();
    DBG("Got eof, fail or empty getline at ", tellg(), " eof=", in->eof(), " fail=", in->fail(), " buf=", buf.size(), "\n");
    return file_size;
  }
  int64_t last_tell = tellg();

  // read up to 20 lines and store tellg() at each line, then process them below...
  std::vector<string> lines;
  lines.reserve(20);
  std::vector<int64_t> tells;
  tells.reserve(20);
  for (int i = 0; i < 20; i++) {
    tells.push_back(tellg());
    std::getline(*in, buf);
    if (in->eof() || in->fail() || buf.empty()) {
      DBG("Got eof, fail or empty getline at ", tellg(), " eof=", in->eof(), " fail=", in->fail(), " buf=", buf.size(), "\n");
      break;
    }
    assert(buf.back() != '\n');
    lines.push_back(buf);
  }

  char last_pair = '\0', this_pair = '\1';
  int i;
  string possible_header, last_header;
  for (i = 0; i < lines.size(); i++) {
    int64_t this_tell = tells[i];
    string &tmp = lines[i];
    if (tmp[0] == '@') {
      string header(tmp);
      possible_header = header;
      if (possible_header.compare(last_header) == 0) {
        // test and possibly fix identical read pair names without corresponding header field - Issue124
        if (_is_paired) {
          _fix_paired_name = true;
          LOG("Detected indistinguishable paired read names which will be fixed on-the-fly: ", last_header, " in ", this->fname, "\n");
          if (offset == 0) WARN(this->fname, " is paired but read names are indistinguisable.  example: ", possible_header, "\n");
        } else {
          DIE("Invalid unpaired fastq file that contains identical sequential read names: ", last_header, " in ", this->fname, "\n");
        }
      }
      rtrim(header);
      if (!get_fq_name(header)) continue;

      // now read another three lines, but check that next line is sequence and the third line is a + separator and the fourth line
      // is the same length as the second
      int64_t test_tell = tells[i];
      int seqlen = 0;
      bool record_found = true;
      DBG_VERBOSE("Testing for header: ", header, "\n");
      for (int j = 0; j < 3; j++) {
        if (i + 1 + j >= lines.size()) DIE("Missing record info at pos ", tells[i]);
        string &tmp2 = lines[i + 1 + j];

        if (j == 0) {
          if (tmp2[0] == '@') {
            // previous was quality line, this next line is the header line.
            record_found = false;
            break;
          }
          // sequence line should have only sequence characters
          seqlen = tmp2.size();
          for (int k = 0; k < seqlen; k++) {
            switch (tmp2[k]) {
              case ('a'):
              case ('c'):
              case ('g'):
              case ('t'):
              case ('A'):
              case ('C'):
              case ('G'):
              case ('T'):
              case ('N'): break;  // okay
              case ('\n'):        
                // newline is not expected here, but okay if last char
                record_found = k == seqlen - 1;
                if (record_found) break;
              case ('@'):         // previous was quality line, this next line is the header line.
              default:            // not okay
                DBG_VERBOSE("Found non-seq '", tmp2[k], "' at ", k, " in: ", tmp2, "\n");
                record_found = false;
            }
            if (!record_found) break;
          }
          if (!record_found) break;
        }
        if (j == 1 && tmp2[0] != '+') {
          // this is not a correct boundary for a fastq record - move on
          DBG_VERBOSE("Found non + ", tmp2, "\n");
          record_found = false;
          break;
        }
        if (j == 2 && seqlen != tmp2.size()) {
          // qual should be same length as sequence
          DBG_VERBOSE("Found different len ", seqlen, " vs ", tmp2.size(), " in ", tmp2, "\n");
          record_found = false;
          break;
        }
      }
      if (record_found) {
        // good
        i += 3;
      } else {
        // test next line as a potential header
        continue;
      }

      if (!is_paired()) {
        DBG("Found record for unpaired read file: ", header, "\n");
        last_tell = this_tell;
        break;
      }

      // need this to be the second of the pair
      this_pair = header[header.length() - 1];
      DBG("Found possible pair ", (char)this_pair, ", header: ", header, "\n");
      if (last_pair == this_pair) {
        if (_fix_paired_name) {
          if (last_header.compare(possible_header) == 0) {
            DBG("Second indistinguishable pair, so keep at the first record\n");
            break;
          }
        } else if (header[header.length() - 2] == '/') {
          // one of two split files or a merged file without pairs
          DBG("Same pair, so keep at the first record\n");
          break;
        }
      } else if (last_pair == '1' && this_pair == '2') {
        // proper pair
        DBG("Found proper pair 1&2\n");
        break;
      }

      // did not find valid new start of (possibly paired) record
      last_tell = this_tell;
      last_pair = this_pair;
      last_header = possible_header;
    }
    if (i >= lines.size()) DIE("Could not find a valid line in the fastq file ", fname, ", last line: ", buf);
  }
  io_t.stop();
  DBG("Found record at ", last_tell, " after offset ", offset, "\n");
  if (offset == 0 && last_tell != 0) {
    WARN("First rank to read really should return the first byte for ", fname, " offset=", offset, " last_tell=", last_tell, "\n");
  }
  return last_tell;
}

FastqReader::FastqReader(const string &_fname, upcxx::future<> first_wait)
    : fname(_fname)
    , in(nullptr)
    , max_read_len(0)
    , read_count(0)
    , fqr2(nullptr)
    , first_file(true)
    , _is_paired(true)
    , _fix_paired_name(false)
    , _first_pair(true)
    , io_t("fastq IO for " + fname)
    , dist_prom(world())
    , open_fut(make_future()) {
  Timer construction_timer("FastqReader construct " + get_basename(fname));
  construction_timer.initiate_entrance_reduction();
  bool need_wait = first_wait.ready();
  string fname2;
  size_t pos;
  LOG("Constructed FastqReader with fname=", fname, " first_wait=", (first_wait.ready() ? "Ready" : "NOT READY"), "\n");
  if ((pos = fname.find(':')) != string::npos) {
    if (pos == fname.size() - 1) {
      // unpaired/single file
      _is_paired = false;
      fname = fname.substr(0, pos);
      LOG("New unpaired fname=", fname, "\n");
    } else {
      // colon separating a pair into two files
      fname2 = fname.substr(pos + 1);
      fname = fname.substr(0, pos);
      LOG("New two file paired fname=", fname, " fname2=", fname2, "\n");
    }
  }
  static int rotate_rank = 0;
  int query_rank = rotate_rank++ % rank_n();
  if (rank_me() == query_rank) {
    // only one rank gets the file size, to prevent many hits on metadata
    io_t.start();
    file_size = upcxx_utils::get_file_size(fname);
    in = std::make_unique<ifstream>(fname);
    io_t.stop();
    buf.reserve(BUF_SIZE);
    DBG("Found file_size=", file_size, " for ", fname, "\n");
  }

  future<> file_size_fut = upcxx::broadcast(file_size, query_rank).then([&self = *this](int64_t sz) {
    self.file_size = sz;
    assert(sz>=0);
  });
  file_size_fut = file_size_fut.then([&self = *this]() { self.know_file_size.fulfill_anonymous(1); });

  // continue opening IO operations to find this rank's start record in a separate thread
  open_fut = when_all(open_fut, file_size_fut, first_wait).then([&self = *this]() { return self.continue_open(); });

  if (!fname2.empty()) {
    // this second reader is generally hidden from the user
    LOG("Opening second FastqReader with ", fname2, "\n");
    fqr2 = make_shared<FastqReader>(fname2, open_fut);
    open_fut = when_all(open_fut, fqr2->open_fut);

    // find the positions in the file that correspond to matching pairs
    auto sh_promstartstop1 = make_shared<dist_object<PromStartStop>>(world());
    auto sh_promstartstop2 = make_shared<dist_object<PromStartStop>>(world());
    auto sh_prombarrier = make_shared<upcxx_utils::PromiseBarrier>(world());
    open_fut = when_all(open_fut, fqr2->open_fut)
                   .then([this, sh_promstartstop1, sh_promstartstop2, sh_prombarrier]() {
                     FastqReader &fqr1 = *this;
                     FastqReader &fqr2 = *(this->fqr2);
                     return set_matching_pair(fqr1, fqr2, *sh_promstartstop1, *sh_promstartstop2);
                   })
                   .then([sh_promstartstop1, sh_promstartstop2, sh_prombarrier]() {
                     sh_prombarrier->fulfill();
                     return sh_prombarrier->get_future();
                   });
    fqr2->open_fut = open_fut;
    future<> free_mem =
        open_fut.then([sh_promstartstop1, sh_promstartstop2, sh_prombarrier]() {});  // keep shared_ptrs alive until finished
    upcxx_utils::Timings::set_pending(free_mem);                                     // eventualy will clean up
  }

  DBG("Constructed new FQR on ", fname, " - ", (void *)this, "\n");

  if (need_wait) {
    open_fut.wait();
    if (fqr2) fqr2->open_fut.wait();
  }
}

future<> FastqReader::set_matching_pair(FastqReader &fqr1, FastqReader &fqr2, dist_object<PromStartStop> &dist_start_stop1,
                                        dist_object<PromStartStop> &dist_start_stop2) {

  DBG("Starting matching pair ", fqr1.start_read, " and ", fqr2.start_read, "\n");
  assert(fqr1.in && "FQ 1 is open");
  assert(fqr2.in && "FQ 2 is open");
  assert(fqr1.start_read == fqr1.tellg() && "fqr1 is at the naive start");
  assert(fqr2.start_read == fqr2.tellg() && "fqr2 is at the naive start");
  // disable subsampling temporarily
  auto old_subsample = fqr1.subsample_pct;
  assert(fqr2.subsample_pct == old_subsample);
  fqr1.subsample_pct = 100;
  fqr2.subsample_pct = 100;
  int64_t pos1 = fqr1.start_read, pos2 = fqr2.start_read;
  // allow search to extend past the original block size
  fqr1.end_read = fqr1.file_size;
  fqr2.end_read = fqr2.file_size;
  string id, seq, qual;
  int64_t offset1 = 0, offset2 = 0;
  string read1, read2;
  int64_t bytes1 = 0, bytes2 = 0;
  while (true) {
    offset1 += bytes1;
    bytes1 = fqr1.get_next_fq_record(id, seq, qual, false);
    if (bytes1 == 0) {
      // no match found in entire block. so no reads will be returned
      // this should only happen in the last ranks of a small file in a large job
      offset1 = fqr1.end_read - fqr1.start_read;
      offset2 = fqr2.end_read - fqr2.start_read;
      break;
    }
    rtrim(id);
    DBG("id read1 '", id, "' '", id[id.size() - 2], "' '", id[id.size() - 1], "'\n");
    assert(id.size() >= 3 && id[id.size() - 2] == '/' && id[id.size() - 1] == '1' && "read1 has the expected format");
    id = id.substr(0, id.size() - 2);
    if (read2.compare(id) == 0) {
      offset2 = 0;
      DBG("Found pair read1 ", id, " at ", pos1 + offset1, " matches read2 ", read2, " at ", pos2, " offset1=", offset1, "\n");
      read1 = id;
      break;
    }
    if (read1.empty()) read1 = id;
    offset2 += bytes2;
    bytes2 = fqr1.get_next_fq_record(id, seq, qual, false);  // use fqr1 as it calls fqr2 to mimic interleaving!!!
    if (bytes2 == 0) {
      // no match found in entire block. so no reads will be returned
      // this should only happen in the last ranks of a small file in a large job
      offset1 = fqr1.end_read - fqr1.start_read;
      offset2 = fqr2.end_read - fqr2.start_read;
      break;
    }
    rtrim(id);
    DBG("id read2 ", id, "\n");
    assert(id.size() >= 3 && id[id.size() - 2] == '/' && id[id.size() - 1] == '2' && "read2 has the expected format");
    id = id.substr(0, id.size() - 2);
    if (read1.compare(id) == 0) {
      offset1 = 0;
      DBG("Found pair read1 ", read1, " at ", pos1, " matches read2 ", id, " at ", pos2 + offset2, " offset2=", offset2, "\n");
      read2 = id;
      break;
    }
    if (read2.empty()) read2 = id;
  }
  if (rank_me() == 0) {
    assert(pos1 + offset1 == 0 && pos2 + offset2 == 0 && "Rank0 starts at pos 0");
  }
  LOG("Found matching pair for ", fqr1.fname, " at ", pos1 + offset1, " and ", fqr2.fname, " at ", pos2 + offset2, " - ", read1,
      "\n");
  dist_start_stop1->start_prom.fulfill_result(pos1 + offset1);
  dist_start_stop2->start_prom.fulfill_result(pos2 + offset2);
  if (rank_me() > 0) {
    rpc_ff(
        rank_me() - 1,
        [](dist_object<PromStartStop> &dist_start_stop1, dist_object<PromStartStop> &dist_start_stop2, int64_t end1, int64_t end2) {
          dist_start_stop1->stop_prom.fulfill_result(end1);
          dist_start_stop2->stop_prom.fulfill_result(end2);
        },
        dist_start_stop1, dist_start_stop2, pos1 + offset1, pos2 + offset2);
  }
  if (rank_me() == rank_n() - 1) {
    // special case of eof -- use file_size
    dist_start_stop1->stop_prom.fulfill_result(fqr1.file_size);
    dist_start_stop2->stop_prom.fulfill_result(fqr2.file_size);
  }
  auto fut1 = dist_start_stop1->set(fqr1).then([&fqr1]() {
    fqr1.seek_start();
    fqr1.first_file = true;
  });
  auto fut2 = dist_start_stop2->set(fqr2).then([&fqr2]() { fqr2.seek_start(); });
  return when_all(fut1, fut2).then([&fqr1, &fqr2, old_subsample]() {
    // restore subsampling
    fqr1.subsample_pct = old_subsample;
    fqr2.subsample_pct = old_subsample;
    DBG("Found matching pair ", fqr1.start_read, " and ", fqr2.start_read, "\n");
  });
}

// Find my boundary start and communicate to prev rank for their boundary end (if my start>0)...
upcxx::future<> FastqReader::continue_open() {
  if (block_size == -1) return continue_open_default_per_rank_boundaries();  // the old algorithm
  // use custome block start and block_size
  if (block_size == 0) {
    // this rank does not read this file
    DBG("This rank will not read ", fname, "\n");
    dist_prom->start_prom.fulfill_result(file_size);
    dist_prom->stop_prom.fulfill_result(file_size);
    return dist_prom->set(*this);
  }

  io_t.start();
  if (!in) {
    in.reset(new ifstream(fname));
    LOG("Opened ", fname, " in ", io_t.get_elapsed_since_start(), "s.\n");
  }
  if (!in) {
    SDIE("Could not open file ", fname, ": ", strerror(errno));
  }
  DBG("in.tell=", in->tellg(), "\n");

  io_t.stop();

  using DPSS = dist_object<PromStartStop>;
  assert(block_start <= file_size);
  int64_t my_start = get_fptr_for_next_record(block_start);
  if (rank_me() > 0 && my_start != 0) {
    rpc_ff(rank_me() - 1,
           [](DPSS &dpss, int64_t end_pos) {
             DBG("end_pos=", end_pos, "\n");
             dpss->stop_prom.fulfill_result(end_pos);
           },
           dist_prom, my_start);
  }
  dist_prom->start_prom.fulfill_result(my_start);

  DBG("my_start=", my_start, "\n");

  if (block_start + block_size >= file_size) {
    // next rank will never send their end as it will not even try to open this file
    dist_prom->stop_prom.fulfill_result(file_size);
  } else {
    // expecting an rpc from rank_me + 1
  }
  auto fut_set = dist_prom->set(*this);
  auto fut_seek = fut_set.then([this]() {
    DBG("start_read=", start_read, " end_read=", end_read, "\n");
    if (start_read != end_read) this->seek_start();
  });
  progress();
  return fut_seek;
}

// all ranks open, 1 rank per node finds block boundaries
upcxx::future<> FastqReader::continue_open_default_per_rank_boundaries() {
  assert(upcxx::master_persona().active_with_caller());
  assert(know_file_size.get_future().ready());
  assert(block_size == -1);
  auto sz = INT_CEIL(file_size, rank_n());
  set_block(sz * rank_me(), sz);
  io_t.start();
  if (!in) {
    in.reset(new ifstream(fname));
  }
  if (!in) {
    SDIE("Could not open file ", fname, ": ", strerror(errno));
  }
  DBG("in.tell=", in->tellg(), "\n");
  LOG("Opened ", fname, " in ", io_t.get_elapsed_since_start(), "s.\n");
  io_t.stop();
  if (rank_me() == 0) {
    // special for first rank set start as 0
    DBG_VERBOSE("setting 0 for rank 0\n");
    dist_prom->start_prom.fulfill_result(0);
  }
  if (rank_me() == rank_n() - 1) {
    // special for last rank set end as file size
    DBG_VERBOSE("Setting file_size=", file_size, " for last rank\n");
    dist_prom->stop_prom.fulfill_result(file_size);
  }
  // have all other local ranks delay their seeks and I/O until these seeks are finished.
  promise wait_prom(1);
  if (local_team().rank_me() == local_team().rank_n() - 1) {
    DBG_VERBOSE("Fseeking for the team\n");
    // do all the fseeking for the local team
    using DPSS = dist_object<PromStartStop>;
    int64_t file_bytes = file_size;
    int64_t read_block = INT_CEIL(file_bytes, rank_n());
    auto first_rank = rank_me() + 1 - local_team().rank_n();
    for (auto rank = first_rank; rank < first_rank + local_team().rank_n(); rank++) {
      // just a part of the file is read by this thread
      if (rank == 0) {
        // already done - special already set 0 for first rank
        continue;
      }
      assert(rank > 0);
      auto pos = block_start;
      pos = get_fptr_for_next_record(block_start);
      DBG_VERBOSE("Notifying with rank=", rank, " pos=", pos, " after block_start=", block_start, "\n");
      wait_prom.get_future().then([rank, pos, &dist_prom = this->dist_prom]() {
        DBG_VERBOSE("Sending pos=", pos, " to stop ", rank - 1, " and start ", rank, "\n");
        // send end to prev rank
        rpc_ff(rank - 1, [](DPSS &dpss, int64_t end_pos) { dpss->stop_prom.fulfill_result(end_pos); }, dist_prom, pos);
        // send start to rank
        rpc_ff(rank, [](DPSS &dpss, int64_t start_pos) { dpss->start_prom.fulfill_result(start_pos); }, dist_prom, pos);
      });
    }
  }
  // all the seeks are done, send results to local team
  wait_prom.fulfill_anonymous(1);
  auto fut_set = dist_prom->set(*this);
  auto fut_seek = fut_set.then([this]() { this->seek_start(); });
  progress();
  return fut_seek;
}

void FastqReader::advise(bool will_need) {
  return;  // TODO fix support for ifstream!
#if defined(__APPLE__) && defined(__MACH__)
// TODO
#else
#ifndef MHM2_NO_FADVISE
/* TODO Find a way to support this with ifstream!
  BaseTimer advise_t("Advise " + fname);
  advise_t.start();
  if (will_need) {
    posix_fadvise(fileno(f), start_read, end_read - start_read, POSIX_FADV_SEQUENTIAL);
    LOG("advised ", fname, " POSIX_FADV_SEQUENTIAL in ", advise_t.get_elapsed_since_start(), "s\n");
  } else {
    posix_fadvise(fileno(f), start_read, end_read - start_read, POSIX_FADV_DONTNEED);
    LOG("advised ", fname, " POSIX_FADV_DONTNEED in ", advise_t.get_elapsed_since_start(), "s\n");
  }
  advise_t.end();
  */
#else
  SLOG_VERBOSE("No posix_fadvice to ", fname, "\n");
#endif
#endif
  if (fqr2) fqr2->advise(will_need);  // advise the second file too!
}

void FastqReader::set_block(int64_t start, int64_t size) {
  if (size > 0) LOG("set_block: file=", fname, " start=", start, " size=", size, " file_size=", file_size, "\n");
  block_start = start;
  block_size = size;
}

void FastqReader::seek_start() {
  // seek to first record
  DBG("Seeking to start_read=", start_read, " reading through end_read=", end_read, " my_file_size=", my_file_size(),
      " fname=", fname, "\n");
  io_t.start();
  seekg(start_read);
  if (!in->good()) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
  SLOG_VERBOSE("Reading FASTQ file ", fname, "\n");
  double fseek_t = io_t.get_elapsed_since_start();
  io_t.stop();
  LOG("Reading fastq file ", fname, " at pos ", start_read, "==", tellg(), " to ", end_read, " seek ", fseek_t,
      "s io to open+find+seek ", io_t.get_elapsed(), "s\n");
}

FastqReader::~FastqReader() {
  if (!open_fut.ready()) {
    WARN("Destructor called before opening completed\n");
    open_fut.wait();
  }

  if (in) {
    io_t.start();
    in->close();
    io_t.stop();
  }
  in.reset();

  io_t.done_all_async();  // will print in Timings' order eventually
  FastqReader::overall_io_t += io_t.get_elapsed();
  DBG("Deconstructed FQR on ", fname, " - ", (void *)this, "\n");
}

string FastqReader::get_fname() { return fname; }

size_t FastqReader::my_file_size() {
  size_t size = 0;
  size = end_read - start_read;
  if (fqr2) size += fqr2->my_file_size();
  return size;
}

future<int64_t> FastqReader::get_file_size() const {
  return know_file_size.get_future().then([this]() { return this->file_size; });
}

size_t FastqReader::get_next_fq_record(string &id, string &seq, string &quals, bool wait_open) {
  if (wait_open && !open_fut.ready()) {
    WARN("Attempt to read ", fname, " before it is ready. wait on open_fut first to avoid this warning!\n");
    open_fut.wait();
  }
  if (end_read == start_read) return 0;
  if (subsample_pct != 100) {
    assert(subsample_pct > 0);

    if (fqr2) assert(subsample_pct == fqr2->subsample_pct);
    if (fqr2 && fqr2->read_count != read_count) {
      // return the second mate normally
    } else {
      int modulo = (read_count / (_is_paired ? 2 : 1)) % 100;
      if (modulo >= subsample_pct) {
        // fast forward and discard
        auto tmp_read_count = read_count;
        int skip = (_is_paired ? 2 : 1) * (100 - modulo);
        DBG_VERBOSE("Skipping ", skip, " reads for ", (_is_paired ? "PAIRED " : ""), "subsample ", subsample_pct,
                    "% read_count=", read_count, "\n");
        for (int i = 0; i < skip; i++) {
          read_count = 0;
          auto bytes_read = get_next_fq_record(id, seq, quals, wait_open);
          if (bytes_read == 0) return bytes_read;
        }
        read_count = tmp_read_count + skip;
        assert(read_count % (_is_paired ? 200 : 100) == 0);
      }
      // read normally
      read_count++;
    }
  }
  if (fqr2) {
    // return a single interleaved file
    if (first_file) {
      first_file = false;
    } else {
      first_file = true;
      assert(read_count > fqr2->read_count);
      auto bytes_read = fqr2->get_next_fq_record(id, seq, quals, wait_open);
      assert(read_count == fqr2->read_count);
      if (_fix_paired_name) id[id.length()-1] = '2'; // this is read2
      return bytes_read;
    }
  }
  if (!in || in->eof() || tellg() >= end_read) return 0;
  io_t.start();
  size_t bytes_read = 0;
  id = "";
  char id2 = '\0';
  for (int i = 0; i < 4; i++) {
    std::getline(*in, buf);
    if (!in->good()) DIE("Read record terminated on file ", fname, " before full record at position ", tellg());
    if (i == 0)
      id.assign(buf);
    else if (i == 1)
      seq.assign(buf);
    else if (i == 2)
      id2 = buf[0];
    else if (i == 3)
      quals.assign(buf);
    bytes_read += buf.size();
  }
  rtrim(id);
  rtrim(seq);
  rtrim(quals);
  if (id[0] != '@') DIE("Invalid FASTQ in ", fname, ": expected read name (@), got: id=", id, " at ", tellg(), "\n");
  if (id2 != '+') DIE("Invalid FASTQ in ", fname, ": expected '+', got: '", id2, "' id=", id, " at ", tellg(), "\n");
  // construct universally formatted name (illumina 1 format)
  if (!get_fq_name(id)) DIE("Invalid FASTQ in ", fname, ": incorrect name format '", id, "'");
  // get rid of spaces
  replace_spaces(id);
  if (_is_paired & _fix_paired_name) {
    if (fqr2) {
      id += "/1";
    } else {
      if (_first_pair) id += "/1";
      else id += "/2";
    }
  }
  if (seq.length() != quals.length())
    DIE("Invalid FASTQ in ", fname, ": sequence length ", seq.length(), " != ", quals.length(), " quals length\n", "id:   ", id,
        "\nseq:  ", seq, "\nquals: ", quals);
  if (seq.length() > max_read_len) max_read_len = seq.length();
  io_t.stop();
  DBG_VERBOSE("Read ", id, " bytes=", bytes_read, "\n");
  _first_pair = !_first_pair;
  return bytes_read;
}

int FastqReader::get_max_read_len() { return std::max(max_read_len, fqr2 ? fqr2->get_max_read_len() : 0u); }

double FastqReader::get_io_time() { return overall_io_t; }

void FastqReader::reset() {
  DBG("reset on FQR of ", fname, "\n");
  if (!open_fut.ready()) {
    open_fut.wait();
  }
  if (block_size > 0) {
    if (!in) {
      DIE("Reset called on unopened file\n");
    }
    io_t.start();
    assert(in && "reset called on active file");
    seekg(start_read);
    if (!in->good()) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    io_t.stop();
    read_count = 0;
    first_file = true;
    DBG("reset on ", fname, " tellg=", in->tellg(), "\n");
  }
  if (fqr2) fqr2->reset();
}

//
// FastqReaders
//

FastqReaders::FastqReaders()
    : readers() {}

FastqReaders::~FastqReaders() {
  DBG("Deconstruct FastqReaders!\n");
  close_all();
}
FastqReaders &FastqReaders::getInstance() {
  static FastqReaders _;
  return _;
}

bool FastqReaders::is_open(const string fname) {
  FastqReaders &me = getInstance();
  auto it = me.readers.find(fname);
  return it != me.readers.end() && it->second->is_open();
}

size_t FastqReaders::get_open_file_size(const string fname) {
  assert(is_open(fname));
  FastqReaders &me = getInstance();
  auto it = me.readers.find(fname);
  assert(it != me.readers.end());
  auto fut_sz = it->second->get_file_size();
  assert(fut_sz.ready());
  return fut_sz.wait();
}

FastqReader &FastqReaders::open(const string fname, int subsample_pct, upcxx::future<> first_wait) {
  FastqReaders &me = getInstance();
  auto it = me.readers.find(fname);
  if (it == me.readers.end()) {
    DBG("Opening a new FQR ", fname, "\n");
    upcxx::discharge();  // opening itself may take some time
    auto sh_fqr = make_shared<FastqReader>(fname, first_wait);
    if (subsample_pct < 100) sh_fqr->set_subsample_pct(subsample_pct);
    it = me.readers.insert(it, {fname, sh_fqr});
    upcxx::discharge();  // opening requires a broadcast with to complete
    upcxx::progress();   // opening requires some user progress too
  } else {
    DBG("Re-opened an existing FQR ", fname, "\n");
  }
  assert(it != me.readers.end());
  assert(it->second);
  return *(it->second);
}

FastqReader &FastqReaders::get(const string fname) {
  FastqReader &fqr = open(fname);
  fqr.reset();
  return fqr;
}

void FastqReaders::close(const string fname) {
  FastqReaders &me = getInstance();
  auto it = me.readers.find(fname);
  if (it != me.readers.end()) {
    me.readers.erase(it);
  }
}

void FastqReaders::close_all() {
  DBG("close_all\n");
  FastqReaders &me = getInstance();
  me.readers.clear();
}
