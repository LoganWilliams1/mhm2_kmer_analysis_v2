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
using upcxx::progress;
using upcxx::rank_me;
using upcxx::rpc_ff;

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
  if (_trim_comment) {
    auto space_pos = header.find_first_of(' ');
    if (space_pos != string::npos) header = header.substr(0, space_pos);
  }
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

int64_t FastqReader::tellg() { return in->tellg(); }

void FastqReader::seekg(int64_t pos) {
  assert(pos < file_size);
  if (!(in && in->is_open() && in->good())) {
    DIE("Did not seekg to ", pos, " fname=", get_ifstream_state(), ": ", strerror(errno));
  }
  in->seekg(pos);
  if (!(in && in->is_open() && in->good())) {
    DIE("Could not seekg to ", pos, " fname=", get_ifstream_state(), ": ", strerror(errno));
  }
  auto new_pos = tellg();
  if (new_pos != pos) {
    DIE("seekg failed to change to pos=", pos, " tellg=", new_pos, " fname=", fname, "\n");
  }
}

bool FastqReader::is_open() { return open_fut.ready(); }

bool FastqReader::is_open(bool wait_for_open) {
  if (wait_for_open) {
    assert(!upcxx::in_progress());
    open_fut.wait();
  }
  return is_open();
}

FastqReader &FastqReader::get_fqr2() {
  assert(fqr2);
  return *fqr2;
}

bool FastqReader::is_sep(const char &sep) const {
  return ((sep == '/') | (sep == '.') | (sep == 'R') | (sep == ':'));  // possible paired read separators
}

int64_t FastqReader::get_fptr_for_next_record(int64_t offset) {
  // first record is the first record, include it.  Every other partition will be at least 1 full record after offset.
  // but read the first few lines anyway

  if (!in || !in->is_open()) DIE("Fastq ", fname, " is not open to find the next record after offset=", offset, "\n");
  // eof - do not read anything
  if (offset > 0 && offset + 1024 >= file_size) {
    DBG("offset=", offset, " is less than 1kb from end of file, returning eof file_size=", file_size, "\n");
    return file_size;
  }

  read_io_t.start();
  in->seekg(offset);
  read_io_t.stop();
  if (!in->good() || in->tellg() != offset)
    DIE("Could not seekg to ", offset, " fname=", get_ifstream_state(), ": ", strerror(errno));

  if (offset != 0) {
    // skip first (likely partial) line after this offset to ensure we start at the beginning of a line
    read_io_t.start();
    std::getline(*in, buf);
    read_io_t.stop();
  }
  if (buf.empty() && (in->eof() || in->fail())) {
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
    read_io_t.start();
    std::getline(*in, buf);
    read_io_t.stop();
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
        if (_is_paired && !_fix_paired_name) {
          _fix_paired_name = true;
          LOG("Detected indistinguishable paired read names which will be fixed on-the-fly: ", last_header, " in ", this->fname,
              "\n");
          if (offset == 0) WARN(this->fname, " is paired but read names are indistinguisable.  example: ", possible_header, "\n");
        } else if (!_is_paired) {
          DIE("Invalid unpaired fastq file that contains identical sequential read names: ", last_header, " in ", this->fname,
              " _is_paired=", _is_paired, " _fix_paired_name=", _fix_paired_name, " i=", i, " offset=", offset, "\n");
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
        if (i + 1 + j >= lines.size())
          DIE("Missing record info at ", get_basename(fname), " around pos ", tells[i], " lines: ", lines[0], " header: ", header,
              "\n");
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
              case ('@'):  // previous was quality line, this next line is the header line.
              default:     // not okay
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
      bool has_sep = header.length() > 2 ? is_sep(header[header.length() - 2]) : false;
      this_pair = header[header.length() - 1];
      if (has_sep) DBG("Found possible pair ", (char)this_pair, ", header: ", header, "\n");
      bool is_same_header = last_header.compare(possible_header) == 0;
      auto tpos1 = string::npos, tpos2 = string::npos;
      if (has_sep && last_pair == this_pair) {
        if (_fix_paired_name) {
          if (is_same_header) {
            LOG("Second indistinguishable pair, so keep at the first record\n");
            break;
          }
        } else if (is_interleaved()) {
          WARN("Improper interleaved-paired file (--reads) at ", this_pair, ".  If this actually two fastq files, "
               "use --paired-reads; if this is unpaired, use --unpaired-reads: ", get_basename(fname));
          WARN("Adjacent reads share the same pair id ", last_header, " vs ", header);
          WARN("Changing file ", this->fname, " from interleaved to unpaired");
          _is_paired = false;
          _is_interleaved = false;
        } else if (is_paired() && !is_interleaved()) {
          DBG("Paired files can have any pair id here\n");
          break;
        } else if (!is_paired()) {
          // one of two split files or a merged file without pairs
          DBG("Unpaired so ignore any pair id keep at the first record\n");
          break;
        }
      } else if (has_sep && last_pair == '1' && this_pair == '2') {
        // proper pair
        DBG("Found proper pair 1&2\n");
        break;
      } else if ((tpos1 = possible_header.find_first_of(' ')) != string::npos &&
                 (tpos2 = last_header.find_first_of(' ')) != string::npos &&
                 possible_header.substr(0, tpos1).compare(last_header.substr(0, tpos2)) == 0) {
        if (offset == 0)
          WARN("Ignoring confounding comments in header of ", get_basename(fname), " '", last_header, "' and '", possible_header,
               "' are proper pair\n");
        _trim_comment = true;
        _fix_paired_name = true;
        DBG("Found proper pair - identical ignoring comment\n");
      } else if (!has_sep && is_same_header) {
        LOG("Second indistinguishable pair (with no pair separator), so keep at the first record\n");
        break;
      }
      if (tpos1 != string::npos && tpos2 != string::npos)
        DBG("Comment trimmed: '", last_header.substr(0, tpos2), "' '", possible_header.substr(0, tpos1), "'\n");

      // did not find valid new start of (possibly paired) record
      last_tell = this_tell;
      last_pair = this_pair;
      last_header = possible_header;
    }
    if (i >= lines.size()) DIE("Could not find a valid line in the fastq file ", fname, ", last line: ", buf);
  }
  if (!(in && in->is_open() && in->good())) {
    DIE("get_fptr_for_next_record: file=", get_ifstream_state(), " offset=", offset, " last_tell=", last_tell, ": ",
        strerror(errno), "\n");
  }

  if (offset > 0 && last_tell + 1024 >= file_size) {
    DBG("Found record at less than 1kb from file end. last_tell=", last_tell, " returning eof file_size=", file_size, "\n");
    return file_size;
  }
  DBG("Found record at ", last_tell, " after offset ", offset, "\n");
  if (offset == 0 && last_tell != 0) {
    WARN("First rank to read really should return the first byte for ", fname, " offset=", offset, " last_tell=", last_tell, "\n");
  }
  return last_tell;
}

FastqReader::FastqReader(const string &_fname, future<> first_wait, bool is_second_file)
    : fname(_fname)
    , in(nullptr)
    , max_read_len(0)
    , read_count(0)
    , avg_bytes_per_read(0)
    , num_reads(0)
    , num_pairs(0)
    , num_bases(0)
    , fqr2(nullptr)
    , first_file(true)
    , _is_paired(true)
    , _is_interleaved(false)
    , _fix_paired_name(false)
    , _first_pair(true)
    , _trim_comment(false)
    , _is_bgzf(false)
    , io_t("fastq IO-ops " + fname)
    , read_io_t("fastq Reading " + fname)
    , dist_prom(world())
    , open_fut(make_future()) {
  assert(!upcxx::in_progress());
  Timer construction_timer("FastqReader construct " + get_basename(fname));
  construction_timer.initiate_entrance_reduction();
  bool need_wait = first_wait.ready();
  string fname2;
  size_t pos;
  LOG("Constructed FastqReader with fname=", fname, " first_wait=", (first_wait.ready() ? "Ready" : "NOT READY"), "\n");

  shared_ptr<dist_object<PromStartStop>> sh_promstartstop1, sh_promstartstop2;
  shared_ptr<upcxx_utils::PromiseBarrier> sh_prombarrier;
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
      DBG("New promises will be fulfilled on fname=", fname, "\n");
      sh_promstartstop1 = make_shared<dist_object<PromStartStop>>(world());
      sh_promstartstop2 = make_shared<dist_object<PromStartStop>>(world());
      sh_prombarrier = make_shared<upcxx_utils::PromiseBarrier>(world());
      LOG("New two file paired fname=", fname, " fname2=", fname2, "\n");
    }
  } else if (!is_second_file)
    _is_interleaved = true;

  static int rotate_rank = 0;
  int query_rank = rotate_rank++ % rank_n();
  if (rank_me() == query_rank) {
    // only one rank gets the file size, to prevent many hits on metadata
    io_t.start();
    file_size = upcxx_utils::get_file_size(fname);
    io_t.stop();
    in = std::make_unique<ifstream>(fname);
    buf.reserve(BUF_SIZE);
    DBG("Found file_size=", file_size, " for ", fname, "\n");
  }

  future<> file_size_fut = upcxx::broadcast(file_size, query_rank).then([&self = *this](int64_t sz) {
    self.file_size = sz;
    assert(sz >= 0);
  });
  file_size_fut = file_size_fut.then([&self = *this]() { self.know_file_size.fulfill_anonymous(1); });

  // continue opening IO operations to find this rank's start record in a separate thread
  open_fut = when_all(open_fut, file_size_fut, first_wait).then([&self = *this]() { return self.continue_open(); });

  if (!fname2.empty()) {
    // this second reader is generally hidden from the user
    LOG("Opening second FastqReader with ", fname2, "\n");
    fqr2 = make_shared<FastqReader>(fname2, first_wait, true);
    open_fut = when_all(open_fut, fqr2->open_fut, dist_prom->get_future(), fqr2->dist_prom->get_future());

    // find the positions in the file that correspond to matching pairs
    // use the promises to coordinate the timings
    open_fut = open_fut
                   .then([this, sh_promstartstop1, sh_promstartstop2, sh_prombarrier]() {
                     FastqReader &fqr1 = *this;
                     FastqReader &fqr2 = *(this->fqr2);
                     return set_matching_pair(fqr1, fqr2, *sh_promstartstop1, *sh_promstartstop2);
                   })
                   .then([sh_promstartstop1, sh_promstartstop2, sh_prombarrier, &fqr1 = *this]() {
                     DBG("Done matching on ", fqr1.get_fname(), " waiting on promise barrier\n");
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
  if (fqr1.start_read == fqr1.end_read) {
    assert(fqr2.start_read == fqr2.end_read);
    DBG("Fulfilling No reading of ", fqr1.fname, " or ", fqr2.fname, "\n");
    dist_start_stop1->start_prom.fulfill_result(fqr1.file_size);
    dist_start_stop2->start_prom.fulfill_result(fqr2.file_size);
    dist_start_stop1->stop_prom.fulfill_result(fqr1.file_size);
    dist_start_stop2->stop_prom.fulfill_result(fqr2.file_size);
    return make_future();
  }
  AsyncTimer t_set_matching_pair("set_matching_pair: " + get_basename(fqr2.fname));
  t_set_matching_pair.start();
  DBG("Starting matching pair on ", fqr1.fname, " ,", fqr1.start_read, " and ", fqr2.start_read, "\n");
  assert(fqr1.in && "FQ 1 is open");
  assert(fqr2.in && "FQ 2 is open");
  assert(fqr1.start_read == fqr1.tellg() && "fqr1 is at the naive start");
  assert(fqr2.start_read == fqr2.tellg() && "fqr2 is at the naive start");
  auto target_read_size = fqr1.end_read - fqr1.start_read;
  // disable subsampling temporarily
  auto old_subsample = fqr1.subsample_pct;
  assert(fqr2.subsample_pct == old_subsample);
  fqr1.subsample_pct = 100;
  fqr2.subsample_pct = 100;

  // allow search to extend past the original block size
  fqr1.end_read = fqr1.file_size;
  fqr2.end_read = fqr2.file_size;
  struct PairPositions {
    int64_t pos1, offset1, pos2, offset2;
  };

  int64_t pos1 = fqr1.start_read, pos2 = fqr2.start_read;
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
    DBG("id read1 '", id, "' offset1=", offset1, " '", id[id.size() - 2], "' '", id[id.size() - 1], "'\n");
    assert(id.size() >= 3 && id[id.size() - 2] == '/' && id[id.size() - 1] == '1' && "read1 has the expected format");
    id = id.substr(0, id.size() - 2);
    if (read2.compare(id) == 0) {
      offset2 = 0;
      DBG("Found pair read1 ", id, " at ", pos1 + offset1, " matches read2 ", read2, " at ", pos2, " offset1=", offset1,
          " offset2=", offset2, "\n");
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
    DBG("id read2 ", id, " offset2=", offset2, "\n");
    assert(id.size() >= 3 && id[id.size() - 2] == '/' && id[id.size() - 1] == '2' && "read2 has the expected format");
    id = id.substr(0, id.size() - 2);
    if (read1.compare(id) == 0) {
      offset1 = 0;
      DBG("Found pair read1 ", read1, " at ", pos1, " matches read2 ", id, " at ", pos2 + offset2, " offset2=", offset2,
          " offset1=", offset1, "\n");
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

  DBG("Fulfilling matching self starts fname=", fqr1.fname, " and ", fqr2.fname, "\n");
  dist_start_stop1->start_prom.fulfill_result(pos1 + offset1);
  dist_start_stop2->start_prom.fulfill_result(pos2 + offset2);
  if (pos1 > 0) {
    assert(pos2 > 0);
    assert(rank_me() > 0);
    // tell the previous rank where I am starting
    rpc_ff(
        rank_me() - 1,
        [](dist_object<PromStartStop> &dist_start_stop1, dist_object<PromStartStop> &dist_start_stop2, int64_t end1, int64_t end2,
           string fname) {
          DBG("Fulfilling both matching end from next rank fname=", fname, " and 2, end1=", end1, " end2=", end2, "\n");
          dist_start_stop1->stop_prom.fulfill_result(end1);
          dist_start_stop2->stop_prom.fulfill_result(end2);
        },
        dist_start_stop1, dist_start_stop2, pos1 + offset1, pos2 + offset2, fqr1.fname);
  }
  if (pos1 + offset1 + target_read_size >= fqr1.file_size) {
    // special case of eof -- use file_size
    DBG("Fulfilling matching end of file fname=", fqr1.fname, " and ", fqr2.fname, "\n");
    dist_start_stop1->stop_prom.fulfill_result(fqr1.file_size);
    dist_start_stop2->stop_prom.fulfill_result(fqr2.file_size);
  }
  auto fut1 = dist_start_stop1->set(fqr1).then([&fqr1]() {
    fqr1.seek_start();
    fqr1.first_file = true;
  });
  auto fut2 = dist_start_stop2->set(fqr2).then([&fqr2]() { fqr2.seek_start(); });
  return when_all(fut1, fut2).then([&fqr1, &fqr2, old_subsample, t_set_matching_pair]() {
    // restore subsampling
    fqr1.subsample_pct = old_subsample;
    fqr2.subsample_pct = old_subsample;
    DBG("Found matching pair ", fqr1.start_read, " and ", fqr2.start_read, "\n");
    t_set_matching_pair.stop();
  });
}

// Find my boundary start and communicate to prev rank for their boundary end (if my start>0)...
future<> FastqReader::continue_open() {
  if (block_size == -1) return continue_open_default_per_rank_boundaries();  // the old algorithm
  AsyncTimer t_continue_open("continue_open: " + get_basename(fname));
  t_continue_open.start();
  // use custom block start and block_size
  if (block_size == 0) {
    // this rank does not read this file
    DBG("Fulfilling rank will not read fname=", fname, "\n");
    dist_prom->start_prom.fulfill_result(file_size);
    dist_prom->stop_prom.fulfill_result(file_size);
    return dist_prom->set(*this).then([t_continue_open]() { t_continue_open.stop(); });
  }

  int attempts = 0;
  if (!in) {
    io_t.start();
    while (attempts < 3) {
      in.reset(new ifstream(fname));
      if (in && in->is_open() && in->good()) break;
      attempts++;
      upcxx_utils::ThreadPool::sleep_ns(500000000 + upcxx::local_team().rank_me() * 1000000);
    }
    LOG("Opened ", fname, " in ", io_t.get_elapsed_since_start(), "s.\n");
    io_t.stop();
  }
  if (!(in && in->is_open() && in->good())) {
    DIE("continue_open: file=", get_ifstream_state(), "\n");
  }
  if (attempts > 0) WARN("It took ", attempts, " attempt(s) to open ", fname, "!\n");
  DBG("in.tell=", in->tellg(), "\n");

  assert(block_start <= file_size);

  int64_t my_start = get_fptr_for_next_record(block_start);

  using DPSS = dist_object<PromStartStop>;
  if (rank_me() > 0 && my_start != 0) {
    rpc_ff(
        rank_me() - 1,
        [](DPSS &dpss, int64_t end_pos, const string &fname) {
          DBG("Fulfilling my end from next end_pos=", end_pos, " fname=", fname, "\n");
          dpss->stop_prom.fulfill_result(end_pos);
        },
        dist_prom, my_start, fname);
  }
  DBG("Fulfilling my start=", my_start, " fname=", fname, "\n");
  dist_prom->start_prom.fulfill_result(my_start);

  if (block_start + block_size >= file_size) {
    // next rank will never send their end as it will not even try to open this file
    DBG("Fulfilling my_end is eof fname=", fname, "\n");
    dist_prom->stop_prom.fulfill_result(file_size);
  } else {
    // expecting an rpc from rank_me + 1
  }
  auto fut_set = dist_prom->set(*this);
  auto fut_seek = fut_set.then([&self = *this, t_continue_open]() {
    DBG("start_read=", self.start_read, " end_read=", self.end_read, "\n");
    if (self.start_read != self.end_read) self.seek_start();
    t_continue_open.stop();
  });
  progress();
  return fut_seek;
}

// all ranks open, 1 rank per node finds block boundaries
future<> FastqReader::continue_open_default_per_rank_boundaries() {
  AsyncTimer t_continue_open_default_per_rank_boundaries("continue_open_default_per_rank_boundaries: " + get_basename(fname));
  t_continue_open_default_per_rank_boundaries.start();
  assert(upcxx::master_persona().active_with_caller());
  assert(know_file_size.get_future().ready());
  assert(block_size == -1);
  LOG("Opening ", fname, " over all ranks, not by global blocks - IO performance may suffer\n");
  auto sz = INT_CEIL(file_size, rank_n());
  set_block(sz * rank_me(), sz);

  if (!in) {
    io_t.start();
    in.reset(new ifstream(fname));
    io_t.stop();
  }
  if (!in) {
    SDIE("Could not open file ", fname, ": ", strerror(errno));
  }
  DBG("in.tell=", in->tellg(), "\n");
  LOG("Opened ", fname, " in ", io_t.get_elapsed(), "s.\n");

  if (rank_me() == 0) {
    // special for first rank set start as 0
    DBG("Fulfilling setting 0 for rank 0 fname=", fname, "\n");
    dist_prom->start_prom.fulfill_result(0);
  }
  if (rank_me() == rank_n() - 1) {
    // special for last rank set end as file size
    DBG("Fulfilling file_size=", file_size, " for last rank fname=", fname, "\n");
    dist_prom->stop_prom.fulfill_result(file_size);
  }
  // have all other local ranks delay their seeks and I/O until these seeks are finished.
  promise wait_prom(1);
  if (local_team().rank_me() == local_team().rank_n() - 1) {
    DBG("Fseeking for the team\n");
    // do all the fseeking for the local team
    using DPSS = dist_object<PromStartStop>;
    int64_t file_bytes = file_size;
    int64_t read_block = INT_CEIL(file_bytes, rank_n());
    auto first_rank = rank_me() + 1 - local_team().rank_n();
    auto pos = read_block * first_rank;
    for (auto rank = first_rank; rank < first_rank + local_team().rank_n(); rank++) {
      // just a part of the file is read by this thread
      if (rank == 0) {
        // already done - special already set 0 for first rank
        continue;
      }
      assert(rank > 0);
      auto read_start = read_block * rank;
      pos = get_fptr_for_next_record(read_start);
      DBG("Notifying with rank=", rank, " pos=", pos, " after read_start=", read_start, "\n");
      wait_prom.get_future().then([rank, pos, &dist_prom = this->dist_prom, fname = this->fname]() {
        DBG("Sending pos=", pos, " to stop ", rank - 1, " and start ", rank, "\n");
        // send end to prev rank
        rpc_ff(
            rank - 1,
            [](DPSS &dpss, int64_t end_pos, const string &fname) {
              DBG("Fulfill from leader end=", end_pos, " fname=", fname, "\n");
              dpss->stop_prom.fulfill_result(end_pos);
            },
            dist_prom, pos, fname);
        // send start to rank
        rpc_ff(
            rank,
            [](DPSS &dpss, int64_t start_pos, const string &fname) {
              DBG("Fulfill start leader start=", start_pos, " fname=", fname, "\n");
              dpss->start_prom.fulfill_result(start_pos);
            },
            dist_prom, pos, fname);
      });
    }
  }
  // all the seeks are done, send results to local team
  wait_prom.fulfill_anonymous(1);
  auto fut_set = dist_prom->set(*this);
  auto fut_seek = fut_set.then([this, t_continue_open_default_per_rank_boundaries]() {
    this->seek_start();
    t_continue_open_default_per_rank_boundaries.stop();
  });
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
  DBG("Seeking to start_read=", start_read, " reading through end_read=", end_read, " my_file_size=", my_file_size(false),
      " fname=", fname, "\n");
  io_t.start();
  seekg(start_read);
  double fseek_t = io_t.get_elapsed_since_start();
  io_t.stop();
  if (!in->good()) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
  SLOG_VERBOSE("Reading FASTQ file ", fname, "\n");
  LOG("Reading fastq file ", fname, " at pos ", start_read, "==", tellg(), " to ", end_read, " seek ", fseek_t,
      "s io to open+find+seek ", io_t.get_elapsed(), "s\n");
}

FastqReader::~FastqReader() {
  close();
  DBG("Deconstructed FQR on ", fname, " - ", (void *)this, "\n");
}

void FastqReader::close() {
  DBG("Closing FQR on ", fname, " - ", (void *)this, "\n");
  if (!open_fut.ready()) {
    WARN("Close called before opening completed\n");
    open_fut.wait();
  }

  if (is_first_file()) {
    auto &pr = upcxx_utils::Timings::get_promise_reduce();
    auto fut_num_reads = pr.reduce_one(get_num_reads(), op_fast_add);
    auto fut_num_pairs = pr.reduce_one(get_num_pairs(), op_fast_add);
    auto fut_num_bases = pr.reduce_one(get_num_bases(), op_fast_add);
    auto fut_report = when_all(fut_num_reads, fut_num_pairs, fut_num_bases)
                          .then([fname = this->fname](auto tot_num_reads, auto tot_num_pairs, auto tot_num_bases) {
                            SLOG_VERBOSE("Finished reading ", get_basename(fname), " tot_num_reads=", tot_num_reads,
                                         " tot_num_pairs=", tot_num_pairs, " tot_num_bases=", tot_num_bases, "\n");
                          });
    upcxx_utils::Timings::set_pending(fut_report);
  }

  if (in) {
    io_t.start();
    in->close();
    io_t.stop();
  }

  io_t.done_all_async();  // will print in Timings' order eventually
  read_io_t.done_all_async();

  FastqReader::overall_io_t += io_t.get_elapsed() + read_io_t.get_elapsed();

  in.reset();
}

string FastqReader::get_fname() const { return fname; }

string FastqReader::get_ifstream_state() const {
  std::stringstream ss;
  ss << fname << " is ";
  if (in) {
    if (in->is_open()) {
      if (in->good()) {
        ss << "good";
      } else {
        ss << "NOT GOOD (";
        ss << strerror(errno);
        ss << "): ";
        if (in->eof()) {
          ss << " EOF";
        }
        if (in->fail()) {
          ss << " FAIL";
        }
        if (in->bad()) {
          ss << " BAD";
        }
      }
    } else {
      ss << "NOT OPEN";
    }
  } else {
    ss << "NOT INITIALIZED";
  }
  ss << ".";
  return ss.str();
}

size_t FastqReader::my_file_size(bool include_file_2) const {
  size_t size = 0;
  size = end_read - start_read;
  if (include_file_2 && fqr2) size += fqr2->my_file_size();
  return size;
}

future<int64_t> FastqReader::get_file_size(bool include_file_2) const {
  future<int64_t> fut_f2_size = (include_file_2 && fqr2) ? fqr2->get_file_size() : make_future<int64_t>(0);
  return when_all(fut_f2_size, know_file_size.get_future()).then([this](int64_t f2_size) { return this->file_size + f2_size; });
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
      auto bytes_read = fqr2->get_next_fq_record(id, seq, quals, wait_open);
      assert(read_count == fqr2->read_count);
      if (_fix_paired_name) id[id.length() - 1] = '2';  // this is read2
      return bytes_read;
    }
  }
  if (!in || in->eof()) return 0;
  auto pos = tellg();
  if (pos >= end_read) return 0;
  if (pos < 0 || !in->good()) {
    DIE("get_next_fq_record: file=", get_ifstream_state(), "\n");
  }

  size_t bytes_read = 0;
  id = "";
  char id2 = '\0';
  for (int i = 0; i < 4; i++) {
    read_io_t.start();
    std::getline(*in, buf);
    read_io_t.stop();
    if (!in->good()) {
      // Issue175 try to re-open once!
      static int reopens = 0;
      if (reopens++ < 1) {
        auto newpos = pos + bytes_read;
        WARN("problem with file ", get_ifstream_state(), " block_start=", block_start, " block_size=", block_size,
             " start_read=", start_read, " end_read=", end_read, "... Attempting to reopen it at ", newpos, "\n");
        in.reset(new ifstream(fname));
        seekg(newpos);
        std::getline(*in, buf);
      }
      if (!(in && in->is_open() && in->good()))
        DIE("Read record terminated on file ", get_ifstream_state(), " before full record at position tellg=", tellg(), " i=", i,
            " block_start=", block_start, " block_size=", block_size, " start_read=", start_read, " end_read=", end_read);
    }
    if (i == 0)
      id.assign(buf);
    else if (i == 1)
      seq.assign(buf);
    else if (i == 2)
      id2 = buf[0];
    else if (i == 3)
      quals.assign(buf);
    bytes_read += buf.size() + 1;
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
      if (_first_pair)
        id += "/1";
      else
        id += "/2";
    }
  }
  if (seq.length() != quals.length())
    DIE("Invalid FASTQ in ", fname, ": sequence length ", seq.length(), " != ", quals.length(), " quals length\n", "id:   ", id,
        "\nseq:  ", seq, "\nquals: ", quals);
  set_max_read_len(seq.length());
  num_bases += seq.length();
  num_reads++;
  DBG_VERBOSE("Read ", id, " bytes=", bytes_read, "\n");
  if ((_is_paired && is_first_file() && !_is_interleaved) || (_is_interleaved && _first_pair)) num_pairs++;
  _first_pair = !_first_pair;
  return bytes_read;
}

void FastqReader::set_max_read_len(int len) {
  if (len > max_read_len) max_read_len = len;
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

    assert(in && "reset called on active file");
    io_t.start();
    seekg(start_read);
    io_t.stop();
    if (!in->good()) DIE("Could not fseek on ", fname, " to ", start_read, ": ", strerror(errno));
    read_count = 0;
    first_file = true;
    num_reads = num_pairs = num_bases = 0;
    DBG("reset on ", fname, " tellg=", in->tellg(), "\n");
  }
  _first_pair = true;
  if (fqr2) fqr2->reset();
}

void FastqReader::set_avg_bytes_per_read(unsigned bytes) { avg_bytes_per_read = bytes; }
unsigned FastqReader::get_avg_bytes_per_read() const { return avg_bytes_per_read; }

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

bool FastqReaders::is_open(const string fname, bool wait_for_open) {
  FastqReaders &me = getInstance();
  auto it = me.readers.find(fname);
  return it != me.readers.end() && it->second->is_open(wait_for_open);
}

size_t FastqReaders::get_open_file_size(const string fname, bool include_file_2) {
  assert(is_open(fname));
  FastqReaders &me = getInstance();
  auto it = me.readers.find(fname);
  assert(it != me.readers.end());
  auto fut_sz = it->second->get_file_size(include_file_2);
  assert(fut_sz.ready());
  return fut_sz.wait();
}

FastqReader &FastqReaders::open(const string fname, int subsample_pct, future<> first_wait) {
  assert(!upcxx::in_progress());
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
  DBG("Closing fname=", fname, " ", it != me.readers.end() ? " found and erasing " : " not found ", "\n");
  if (it != me.readers.end()) {
    me.readers.erase(it);
  }
}

void FastqReaders::close_all() {
  DBG("close_all\n");
  FastqReaders &me = getInstance();
  me.readers.clear();
}
