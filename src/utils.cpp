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

#include "utils.hpp"

#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/timers.hpp"

extern char **environ;

using namespace upcxx_utils;

using std::min;
using std::pair;
using std::string;
using std::string_view;
using std::to_string;
using std::unordered_map;
using std::vector;

size_t estimate_hashtable_memory(size_t num_elements, size_t element_size) {
  // get the hashtable load factor
  HASH_TABLE<char, char> tmp;
  double max_load_factor = tmp.max_load_factor();

  // apply the load factor
  size_t expanded_num_elements = num_elements / max_load_factor + 1;

  // get the next power of two
  --expanded_num_elements;
  expanded_num_elements |= expanded_num_elements >> 1;
  expanded_num_elements |= expanded_num_elements >> 2;
  expanded_num_elements |= expanded_num_elements >> 4;
  expanded_num_elements |= expanded_num_elements >> 8;
  expanded_num_elements |= expanded_num_elements >> 16;
  expanded_num_elements |= expanded_num_elements >> 32;
  ++expanded_num_elements;

  size_t num_buckets = expanded_num_elements * max_load_factor;

  return expanded_num_elements * element_size + num_buckets * 8;
}

string revcomp(const string &seq) {
  string seq_rc = "";
  seq_rc.reserve(seq.size());
  for (int i = seq.size() - 1; i >= 0; i--) {
    switch (seq[i]) {
      case 'A':
      case 'a': seq_rc += 'T'; break;
      case 'C':
      case 'c': seq_rc += 'G'; break;
      case 'G':
      case 'g': seq_rc += 'C'; break;
      case 'T':
      case 't': seq_rc += 'A'; break;
      case 'N':
      case 'n': seq_rc += 'N'; break;
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
      case 'V': seq_rc += 'N'; break;
      default:
        DIE("Illegal char at ", i, "'", ((seq[i] >= 32 && seq[i] <= 126) ? seq[i] : ' '), "' (", (int)seq[i], ") in revcomp of '",
            seq, "'");
    }
  }
  return seq_rc;
}

char comp_nucleotide(char ch) {
  switch (ch) {
    case 'A':
    case 'a': return 'T';
    case 'C':
    case 'c': return 'G';
    case 'G':
    case 'g': return 'C';
    case 'T':
    case 't': return 'A';
    case 'N':
    case 'n': return 'N';
    case '0': return '0';
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
    case 'V': return 'N';
    default: DIE("Illegal char '", ((ch >= 32 && ch <= 126) ? ch : ' '), "' (", (int)ch, ") in comp nucleotide");
  }
  return 0;
}

int hamming_dist(string_view s1, string_view s2, bool require_equal_len) {
  if (require_equal_len && s2.size() != s1.size())  // abs((int)(s2.size() - s1.size())) > 1)
    DIE("Hamming distance substring lengths don't match, ", s1.size(), ", ", s2.size(), "\n");
  int d = 0;
  int min_size = min(s1.size(), s2.size());
  for (int i = 0; i < min_size; i++) d += (s1[i] != s2[i]);
  return d;
}

string get_merged_reads_fname(string reads_fname) {
  // always relative to the current working directory
  auto col_pos = reads_fname.find(':');
  if (col_pos != string::npos) {
    if (col_pos == reads_fname.size() - 1) {
      // unpaired read file, remove just the colon
      reads_fname = reads_fname.substr(0, col_pos);
    } else {
      // paired read files
      // remove the first pair, if it exists
      reads_fname = reads_fname.substr(col_pos);
    }
  }
  return upcxx_utils::remove_file_ext(get_basename(reads_fname)) + "-merged.fastq";
}

void switch_orient(int &start, int &stop, int len) {
  int tmp = start;
  start = len - stop;
  stop = len - tmp;
}

void dump_single_file(const string &fname, const string &out_str, bool append) {
  BarrierTimer timer(__FILEFUNC__);
  SLOG_VERBOSE("Writing ", fname, "\n");
  SWARN("This is not the most efficient way to write a file anymore...\n");
  auto fut_tot_bytes_written = upcxx::reduce_one(out_str.size(), upcxx::op_fast_add, 0);
  upcxx_utils::dist_ofstream of(fname, append);
  of << out_str;
  of.close();
  SLOG_VERBOSE("Successfully wrote ", get_size_str(fut_tot_bytes_written.wait()), " bytes to ", fname, "\n");
  assert(rank_me() || of.get_last_known_tellp() == fut_tot_bytes_written.wait());
}

vector<string> get_dir_entries(const string &dname, const string &prefix) {
  vector<string> dir_entries;
  DIR *dir = opendir(dname.c_str());
  if (dir) {
    struct dirent *en;
    while ((en = readdir(dir)) != NULL) {
      string dname(en->d_name);
      if (dname.substr(0, prefix.length()) == prefix) dir_entries.push_back(dname);
    }
    closedir(dir);  // close all directory
  } else {
    SWARN("Could not open ", dname);
  }
  return dir_entries;
}

std::string &left_trim(std::string &str) {
  auto it = std::find_if(str.begin(), str.end(), [](char ch) { return !std::isspace<char>(ch, std::locale()); });
  str.erase(str.begin(), it);
  return str;
}

int pin_clear() {
#if defined(__APPLE__) && defined(__MACH__)
// TODO
#else
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  for (unsigned i = 0; i < sizeof(cpu_set_t) * 8; i++) {
    CPU_SET(i, &cpu_set);
  }
  if (sched_setaffinity(getpid(), sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), getpid());
    return -1;
  }
#endif
  return 0;
}

string get_proc_pin() {
  ifstream f("/proc/self/status");
  string line;
  string prefix = "Cpus_allowed_list:";
  while (getline(f, line)) {
    if (line.substr(0, prefix.length()) == prefix) {
      DBG(line, "\n");
      line = line.substr(prefix.length(), line.length() - prefix.length());
      return left_trim(line);
      break;
    }
  }
  return "";
}

vector<int> get_pinned_cpus() {
  vector<int> cpus;
  stringstream ss(get_proc_pin());
  while (ss.good()) {
    string s;
    getline(ss, s, ',');
    s = left_trim(s);
    auto dash_pos = s.find('-');
    if (dash_pos != string::npos) {
      int start = std::stoi(s.substr(0, dash_pos));
      int stop = std::stoi(s.substr(dash_pos + 1)) + 1;
      for (int i = start; i < stop; i++) cpus.push_back(i);
    } else {
      cpus.push_back(std::stoi(s));
    }
  }
  return cpus;
}

void pin_proc(vector<int> cpus) {
#if defined(__APPLE__) && defined(__MACH__)
// TODO
#else
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  for (auto cpu : cpus) {
    CPU_SET(cpu, &cpu_set);
  }
  if (sched_setaffinity(getpid(), sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), getpid());
  }
#endif
}

void pin_cpu() {
  auto pinned_cpus = get_pinned_cpus();
  pin_proc({pinned_cpus[upcxx::rank_me() % pinned_cpus.size()]});
  SLOG_VERBOSE("Pinning to logical cpus: process 0 on node 0 pinned to cpu ", get_proc_pin(), "\n");
}

void pin_core() {
  string numa_node_dir = "/sys/devices/system/node";
  auto numa_node_entries = get_dir_entries(numa_node_dir, "node");
  if (numa_node_entries.empty()) return;
  vector<int> my_thread_siblings;
  for (auto &entry : numa_node_entries) {
    ifstream f(numa_node_dir + "/" + entry + "/cpulist");
    string buf;
    getline(f, buf);
    f.close();
    int numa_node_i = std::stoi(entry.substr(4));
    auto cpu_entries = get_dir_entries(numa_node_dir + "/" + entry, "cpu");
    for (auto &cpu_entry : cpu_entries) {
      if (cpu_entry != "cpu" + to_string(upcxx::local_team().rank_me())) continue;
      if (cpu_entry == "cpulist" || cpu_entry == "cpumap") continue;
      f.open(numa_node_dir + "/" + entry + "/" + cpu_entry + "/topology/thread_siblings_list");
      getline(f, buf);
      stringstream ss(buf);
      while (ss.good()) {
        string s;
        getline(ss, s, ',');
        my_thread_siblings.push_back(std::stoi(s));
      }
      break;
    }
  }
  if (!my_thread_siblings.empty()) {
    pin_proc(my_thread_siblings);
    SLOG_VERBOSE("Pinning to cores: process 0 on node 0 pinned to cpus ", get_proc_pin(), "\n");
  }
}

void pin_numa(bool round_robin) {
  string numa_node_dir = "/sys/devices/system/node";
  auto numa_node_entries = get_dir_entries(numa_node_dir, "node");
  if (numa_node_entries.empty()) return;
  vector<pair<string, vector<int>>> numa_node_list(numa_node_entries.size(), {"", {}});
  int num_cpus = 0;
  int hdw_threads_per_core = 0;
  for (auto &entry : numa_node_entries) {
    ifstream f(numa_node_dir + "/" + entry + "/cpulist");
    string buf;
    getline(f, buf);
    f.close();
    int numa_node_i = std::stoi(entry.substr(4));
    numa_node_list[numa_node_i].first = buf;
    auto cpu_entries = get_dir_entries(numa_node_dir + "/" + entry, "cpu");
    for (auto &cpu_entry : cpu_entries) {
      if (cpu_entry == "cpulist" || cpu_entry == "cpumap") continue;
      if (!hdw_threads_per_core) {
        f.open(numa_node_dir + "/" + entry + "/" + cpu_entry + "/topology/thread_siblings_list");
        getline(f, buf);
        // assume that the threads are separated by commans - is this always true?
        hdw_threads_per_core = std::count(buf.begin(), buf.end(), ',') + 1;
      }
      numa_node_list[numa_node_i].second.push_back(std::stoi(cpu_entry.substr(3)));
      num_cpus++;
    }
  }
  SLOG_VERBOSE("On node 0, found a total of ", num_cpus, " hardware threads with ", hdw_threads_per_core, " threads per core on ",
               numa_node_list.size(), " NUMA domains\n");
  // pack onto numa nodes
  int hdw_threads_per_numa_node = num_cpus / numa_node_list.size();
  int cores_per_numa_node = hdw_threads_per_numa_node / hdw_threads_per_core;
  int numa_nodes_to_use = (upcxx::local_team().rank_n() + cores_per_numa_node - 1) / cores_per_numa_node;
  if (numa_nodes_to_use > numa_node_list.size()) numa_nodes_to_use = numa_node_list.size();
  if (numa_nodes_to_use == 0) numa_nodes_to_use = 1;
  int my_numa_node = upcxx::local_team().rank_me() / cores_per_numa_node;
  if (round_robin) my_numa_node = upcxx::local_team().rank_me() % numa_node_list.size();
  vector<int> my_cpu_list = numa_node_list[my_numa_node].second;
  sort(my_cpu_list.begin(), my_cpu_list.end());
  pin_proc(my_cpu_list);
  SLOG_VERBOSE("Pinning to ", numa_nodes_to_use, " NUMA domains each with ", cores_per_numa_node, " cores, ",
               hdw_threads_per_numa_node, " cpus: process 0 on node 0 is pinned to cpus ", get_proc_pin(), "\n");
  DBGLOG("Pinned to numa domain ", my_numa_node, ": ", get_proc_pin(), "\n");
}

void log_local(std::string header, std::string msg) {
  // Consolidate a potentially redundant message to the first rank of the node
  upcxx::dist_object<vector<upcxx::promise<string>>> msgs(upcxx::local_team(), upcxx::local_team().rank_n());
  barrier(upcxx::local_team());
  if (upcxx::local_team().rank_me()) {
    rpc(
        upcxx::local_team(), 0,
        [](upcxx::dist_object<vector<upcxx::promise<string>>> &msgs, const string &msg, int from) {
          (*msgs)[from].fulfill_result(msg);
        },
        msgs, msg, upcxx::local_team().rank_me())
        .wait();
  } else {  // local rank 0
    using _MAP = unordered_map<string, vector<int>>;
    _MAP msg_rank_groups;
    (*msgs)[0].fulfill_result(msg);
    for (int i = 0; i < upcxx::local_team().rank_n(); i++) {
      auto cur_msg = (*msgs)[i].get_future().wait();
      msg_rank_groups[cur_msg].push_back(i);
    }
    // order messages by ranks
    vector<string> ordered_msgs;
    for (auto &[ranks_msg, ranks] : msg_rank_groups) {
      ordered_msgs.push_back(ranks_msg);
    }
    sort(ordered_msgs.begin(), ordered_msgs.end(),
         [&msg_rank_groups](const string &a, const string &b) { return msg_rank_groups[a][0] < msg_rank_groups[b][0]; });
    for (auto const &ranks_msg : ordered_msgs) {
      string ranks;  // collapse serial ranks in the message
      bool in_sequence = false, printed = false;
      int last_rank = -2;
      auto it = msg_rank_groups[ranks_msg].begin();
      for (auto &rank : msg_rank_groups[ranks_msg]) {
        if (rank != last_rank + 1) {
          if (in_sequence && last_rank >= 0) ranks += to_string(last_rank);
          in_sequence = false;
          if (!ranks.empty()) ranks += ",";
          ranks += to_string(rank);
          printed = true;
        } else {
          if (!in_sequence) ranks += "-";
          in_sequence = true;
          printed = false;
        }
        last_rank = rank;
      }
      if (!printed && last_rank >= 0) {
        if (!in_sequence && !ranks.empty()) ranks += ",";
        ranks += to_string(last_rank);
      }
      LOG(header, " - local rank(s) ", ranks, ": ", ranks_msg, "\n");
    }
  }
  barrier(upcxx::local_team());
}

void log_pins() { log_local("CPU Pinnings", get_proc_pin()); }

void log_env() {
  string upcxx_msg;
  string job_msg;
  string other_msg;
  for (char **_env = environ; *_env; ++_env) {
    string env(*_env);
    auto pos = env.find("GASNET");
    if (pos == string::npos) pos = env.find("UPCXX");
    if (pos == string::npos) pos = env.find("FI_CXI_");
    if (pos == string::npos) pos = env.find("MHM2");
    if (pos != string::npos) {
      upcxx_msg += env + ",";
      continue;
    }
    if (pos == string::npos) pos = env.find("JOB");
    if (pos == string::npos) pos = env.find("NODE");
    if (pos == string::npos) pos = env.find("STEP");
    if (pos != string::npos) {
      job_msg += env + ", ";
      continue;
    }

#ifndef DEBUG
    // skip cruft unless debugging
    pos = env.find("PATH");  // skip path
    if (pos != string::npos) continue;
    pos = env.find("_DIR");  // skip dir
    if (pos != string::npos) continue;
    pos = env.find("CRAY");  // skip cray
    if (pos != string::npos) continue;
    pos = env.find("PE_");  // skip cray programming envs
    if (pos != string::npos) continue;
    pos = env.find("LANG");  // skip language
    if (pos != string::npos) continue;
    pos = env.find("SSH");  // skip ssh
    if (pos != string::npos) continue;
    pos = env.find("TERM");  // skip terminal
    if (pos != string::npos) continue;
    pos = env.find("TTY");  // skip tty
    if (pos != string::npos) continue;
    pos = env.find("COLOR");  // skip color
    if (pos != string::npos) continue;
#endif
    other_msg += env + ",";
  }
  log_local("GASNET/UPCXX Environment", upcxx_msg);
  log_local("JOB Environment", job_msg);
  LOG("OTHER Environment: ", other_msg, "\n");
}
