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

#include "upcxx_utils/log.hpp"

#if defined(ENABLE_GASNET_STATS)

#include <iostream>
#include <fstream>
#include <regex>

// We may be compiling with debug-mode GASNet with optimization.
// GASNet has checks to prevent users from blindly doing this,
// because it's a bad idea to run that way in production.
// However in this case we know what we are doing...
#undef NDEBUG
#undef __OPTIMIZE__

#include <gasnetex.h>
#include <gasnet_tools.h>

#include "gasnet_stats.hpp"

#include "upcxx_utils/log.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

static string _current_stats_stage = "";

static vector<string> split_string(const string &content, string in_pattern) {
  vector<string> split_content;
  std::regex pattern(in_pattern);
  copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1), std::sregex_token_iterator(),
       back_inserter(split_content));
  return split_content;
}

struct MsgStats {
  long num;
  long min_size, max_size, tot_size, max_tot_size;

  MsgStats()
      : num(0)
      , min_size(0)
      , max_size(0)
      , tot_size(0)
      , max_tot_size(0) {}

  void set_from_str(const string &num_str, const string &s) {
    long num_here = stol(num_str);
    if (!num_here) return;
    num += num_here;
    auto sizes = split_string(s, "/");
    min_size += stol(sizes[1]);
    max_size += stol(sizes[2]);
    tot_size += stol(sizes[3]);
  }

  void reduce_sizes() {
    auto all_num = reduce_one(num, op_fast_add, 0).wait();
    auto all_min_size = reduce_one(min_size, op_fast_min, 0).wait();
    auto all_max_size = reduce_one(max_size, op_fast_max, 0).wait();
    auto all_max_tot_size = reduce_one(tot_size, op_fast_max, 0).wait();
    auto all_tot_size = reduce_one(tot_size, op_fast_add, 0).wait();
    barrier();
    num = all_num;
    min_size = all_min_size;
    max_size = all_max_size;
    max_tot_size = all_max_tot_size;
    tot_size = all_tot_size;
  }

  static string get_headers() { return "num min_size max_size tot_size avg_size load_balance"; }

  string to_string() {
    ostringstream os;
    double avg_size = (num ? (double)tot_size / num : 0);
    double load_balance = (max_tot_size ? ((double)tot_size / rank_n()) / max_tot_size : 0);
    os << fixed << setprecision(2) << num << " " << min_size << " " << max_size << " " << tot_size << " " << avg_size << " "
       << load_balance << " ";
    return os.str();
  }
};

static void aggregate_stats() {
  string my_fname = "/dev/shm/gasnet_stats." + to_string(rank_me());
  ifstream stats_file(my_fname);
  MsgStats gets_stats, puts_stats, am_med_stats, am_long_stats;
  if (!stats_file) {
    WARN("Could not open ", my_fname, " containing GASNet statistics");
  } else {
    string line;
    // this stage will always be last in the stats file
    while (getline(stats_file, line) && line.find(_current_stats_stage) == string::npos);
    //SLOG(line, "\n");
    while (getline(stats_file, line)) {
      //SLOG(line, "\n");
      auto tokens = split_string(line, "\\s+");
      if (tokens.size() < 3) continue;
      if (tokens[1] == "Total" && tokens[2] == "gets:")
        gets_stats.set_from_str(tokens[3], tokens[7]);
      else if (tokens[1] == "Total" && tokens[2] == "puts:")
        puts_stats.set_from_str(tokens[3], tokens[7]);
      else if (tokens[1] == "AMREQUEST_MEDIUM:")
        am_med_stats.set_from_str(tokens[2], tokens[6]);
      else if (tokens[1] == "PREP_REQUEST_MEDIUM:")
        am_med_stats.set_from_str(tokens[2], tokens[6]);
      else if (tokens[1] == "AMREQUEST_LONG:")
        am_long_stats.set_from_str(tokens[2], tokens[6]);
      else if (tokens[1] == "PREP_REQUEST_LONG:")
        am_long_stats.set_from_str(tokens[2], tokens[6]);
    }
  }
  barrier();
  gets_stats.reduce_sizes();
  puts_stats.reduce_sizes();
  am_med_stats.reduce_sizes();
  am_long_stats.reduce_sizes();
  if (!rank_me()) {
    SLOG(KGREEN, "GASNet communication statistics" KNORM "\n");
    SLOG(KGREEN, "Stage gets puts AM_MEDIUM AM_LONG (", MsgStats::get_headers(), ")" KNORM "\n");
    SLOG(KGREEN, _current_stats_stage, " ", gets_stats.to_string(), puts_stats.to_string(), am_med_stats.to_string(),
         am_long_stats.to_string(), KNORM "\n");
    long overall_bytes = gets_stats.tot_size + puts_stats.tot_size + am_med_stats.tot_size + am_long_stats.tot_size;
    long overall_num = gets_stats.num + puts_stats.num + am_med_stats.num + am_long_stats.num;
    SLOG(KLGREEN, "Overall communication for ", _current_stats_stage, " is ", get_size_str(overall_bytes), " with ", overall_num,
         " messages" KNORM "\n");
  }
}

void begin_gasnet_stats(const string &stage) {
  if (_gasnet_stats) {
    _current_stats_stage = stage;
    GASNETT_STATS_SETMASK("PGAU");
    // SWARN("Collecting communication statistics for ", stage);
    GASNETT_STATS_PRINTF_FORCE("MHM2 stage %s\n", stage.c_str());
  }
}

void end_gasnet_stats() {
  if (_gasnet_stats) {
    barrier();
    GASNETT_STATS_DUMP(1);
    GASNETT_STATS_SETMASK("");
    barrier();
    aggregate_stats();
  }
}

#endif
