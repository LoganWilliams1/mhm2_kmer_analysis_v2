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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CLI11.hpp"
#include "version.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

class Options {
  CLI::App app;

  string config_file = "per_rank/mhm2.config";
  string linked_config_file = "mhm2.config";

  vector<string> splitter(string in_pattern, string &content);

  bool extract_previous_lens(vector<unsigned> &lens, unsigned k);

  double setup_output_dir();

  double setup_log_file();

  static string get_job_id();

 public:
  vector<string> reads_fnames;
  vector<string> paired_fnames;
  vector<string> unpaired_fnames;
  bool adapter_trim = true;
  string adapter_fname;
  unsigned kmer_lens = 21;
  bool default_kmer_lens = true;
  int min_kmer_len = -1;

  int qual_offset = 33;
  bool verbose = false;
  int max_kmer_store_mb = 0;  // per rank - default to use 1% of node memory
  int max_rpcs_in_flight = 100;
  int dmin_thres = 2.0;
  int subsample_fastq_pct = 100;  // percentage of fastq files to read
  bool checkpoint = false;
  bool dump_merged = false;

  bool show_progress = false;
  string pin_by = "numa";
  int max_worker_threads = 3;
  vector<int> insert_size = {0, 0};
  int min_ctg_print_len = 500;
  string output_dir;
  string setup_time;
  bool dump_kmers = false;
  bool use_qf = false;
  // very conservative so as not to drop kmers for low depth datasets
  int sequencing_depth = 4;
  string optimize_for = "default";

  Options();
  ~Options();
  void cleanup();

  bool load(int argc, char **argv);

  void write_config_file();
  void adjust_config_option(const string &opt_name, const string &new_val);

  template <typename T>
  static string vec_to_str(const vector<T> &vec, const string &delimiter = ",") {
    std::ostringstream oss;
    bool is_first = true;
    for (const auto &elem : vec) {
      if (is_first)
        is_first = false;
      else
        oss << delimiter;
      oss << elem;
    }
    return oss.str();
  };
};
