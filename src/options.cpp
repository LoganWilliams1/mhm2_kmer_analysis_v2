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

#include "options.hpp"

#include <sys/stat.h>
#include <unistd.h>

#include <cstdlib>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <string_view>
#include <upcxx/upcxx.hpp>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/mem_profile.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"
#include "version.h"

using namespace upcxx_utils;
using namespace std;

#define YES_NO(X) ((X) ? "YES" : "NO")

vector<string> Options::splitter(string in_pattern, string &content) {
  vector<string> split_content;
  std::regex pattern(in_pattern);
  copy(std::sregex_token_iterator(content.begin(), content.end(), pattern, -1), std::sregex_token_iterator(),
       back_inserter(split_content));
  return split_content;
}

bool Options::extract_previous_lens(vector<unsigned> &lens, unsigned k) {
  for (size_t i = 0; i < lens.size(); i++) {
    if (lens[i] == k) {
      lens.erase(lens.begin(), lens.begin() + i + 1);
      return true;
    }
  }
  return false;
}

double Options::setup_output_dir() {
  auto t_start = clock_now();
  if (output_dir.empty()) DIE("Invalid empty ouput_dir");
  if (!upcxx::rank_me()) {
    // create the output directory (and possibly stripe it)

    // always try to make the output_dir
    if (mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */) == -1) {
      // could not create the directory
      if (errno != EEXIST) {
        SDIE("Could not create output directory '", output_dir, "': ", strerror(errno), "\n");
      }
    }

    // always ensure striping is set or reset wide when lustre is available
    int status = std::system("which lfs 2>&1 > /dev/null");
    bool set_lfs_stripe = WIFEXITED(status) & (WEXITSTATUS(status) == 0);
    int num_osts = 0;
    if (set_lfs_stripe) {
      cout << "Trying to set lfs\n";
      // detect the number of OSTs with "lfs df -l $output_dir" and count the entries with OST:
      string cmd("lfs df -l ");
      cmd += output_dir;
      FILE *f_osts = popen(cmd.c_str(), "r");
      if (f_osts) {
        char buf[256];
        while (char *ret = fgets(buf, 256, f_osts)) {
          std::string_view s(buf);
          if (s.find("OST:") != string::npos) num_osts++;
        }
        pclose(f_osts);
      }
      // reduce to the minimum of 90% or rank_n()
      num_osts = std::min(9 * num_osts / 10, std::min((int)144, (int)rank_n()));  // see Issue #70
    }
    set_lfs_stripe &= num_osts > 0;
    if (set_lfs_stripe) {
      // stripe with count -1 to use all the OSTs
      string cmd = "lfs setstripe --stripe-count " + std::to_string(num_osts) + " --stripe-size 16M " + output_dir;
      status = std::system(cmd.c_str());
      if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        cout << "Set Lustre striping on the output directory: count=" << num_osts << ", size=16M\n";
      else
        cout << "Failed to set Lustre striping on output directory: " << WEXITSTATUS(status) << endl;
    }

    // ensure per_rank dir exists (and additionally has stripe 1, if lfs exists)
    string per_rank = output_dir + "/per_rank";
    status = mkdir(per_rank.c_str(), S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
    if (status != 0 && errno != EEXIST) SDIE("Could not create '", per_rank, "'! ", strerror(errno));
    if (set_lfs_stripe) {
      // stripe per_rank directory with count 1
      string cmd = "lfs setstripe --stripe-count 1 " + per_rank;
      status = std::system(cmd.c_str());
      if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        cout << "Set Lustre striping on the per_rank output directory: " << per_rank << "\n";
      else
        cout << "Failed to set Lustre striping on per_rank output directory: " << WEXITSTATUS(status) << endl;
    }
    // this should avoid contention on the filesystem when ranks start racing to creating these top levels
    char basepath[256];
    for (int i = 0; i < rank_n(); i += 1000) {
      sprintf(basepath, "%s/%08d", per_rank.c_str(), i / 1000);
      status = mkdir(basepath, S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
      if (status != 0 && errno != EEXIST) {
        SWARN("Could not create '", basepath, "'! ", strerror(errno));
        break;  // ignore any errors, just stop
      }
    }
    sprintf(basepath, "%s/00000000/%08d", per_rank.c_str(), 0);
    status = mkdir(basepath, S_IRWXU | S_IRWXG | S_IRWXO | S_ISGID /*use default mode/umask */);
    if (status != 0 && errno != EEXIST) SDIE("Could not mkdir rank 0 per thread directory '", basepath, "'! ", strerror(errno));

    cout << "Using output dir: " << output_dir << "\n";  // required for mhm2.py to find the output_dir when not verbose!

    cout.flush();
    cerr.flush();
  }

  upcxx::barrier();
  // after we change to the output directory, relative paths could be incorrect, so make sure we have the correct path of the
  // reads files
  char cwd_str[FILENAME_MAX];
  if (!getcwd(cwd_str, FILENAME_MAX)) DIE("Cannot get current working directory: ", strerror(errno));
  for (auto &fname : reads_fnames) {
    if (fname[0] != '/') {
      string dir = string(cwd_str) + "/";
      auto spos = fname.find_first_of(':');
      if (spos == string::npos || spos == fname.size() - 1) {
        // interleaved (no colon) or single read (colon at the end)
        fname = dir + fname;
      } else {
        // paired read (colon in between)
        fname = dir + fname.substr(0, spos) + ":" + dir + fname.substr(spos + 1);
      }
    }
  }
  if (!adapter_fname.empty() && adapter_fname[0] != '/') {
    // prepend cwd to adapters file
    adapter_fname = string(cwd_str) + "/" + adapter_fname;
  }
  // all change to the output directory
  auto chdir_attempts = 0;
  while (chdir(output_dir.c_str()) != 0) {
    // failed, retry for 5 more seconds - Issue #69
    if (chdir_attempts++ > 10) DIE("Cannot change to output directory ", output_dir, ": ", strerror(errno));
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }
  DBG("Changed dir to ", output_dir, "\n");
  upcxx::barrier();

  chrono::duration<double> t_elapsed = clock_now() - t_start;
  return t_elapsed.count();
}

double Options::setup_log_file() {
  auto t_start = clock_now();
  if (!upcxx::rank_me()) {
    // check to see if mhm2.log exists. If so, and not restarting, rename it
    if (file_exists("mhm2.log")) {
      string new_log_fname = "mhm2-" + setup_time + ".log";
      cerr << KLRED << "WARNING: " << KNORM << output_dir << "/mhm2.log exists. Renaming to " << output_dir << "/" << new_log_fname
           << endl;
      if (rename("mhm2.log", new_log_fname.c_str()) == -1) DIE("Could not rename mhm2.log: ", strerror(errno));
      // also unlink the rank0 per_rank file (a hard link if it exists)
      unlink("per_rank/00000000/00000000/mhm2.log");  // ignore any errors
    }
  }
  upcxx::barrier();
  chrono::duration<double> t_elapsed = clock_now() - t_start;
  return t_elapsed.count();
}

string Options::get_job_id() {
  static const char *env_ids[] = {"SLURM_JOB_ID", "SLURM_JOBID",  "LSB_JOBID", "JOB_ID",
                                  "COBALT_JOBID", "LOAD_STEP_ID", "PBS_JOBID"};
  static string job_id;
  if (job_id.empty()) {
    for (auto env : env_ids) {
      auto env_p = std::getenv(env);
      if (env_p) {
        job_id = env_p;
        break;
      }
    }
    if (job_id.empty()) job_id = std::to_string(getpid());
    auto sz = upcxx::broadcast(job_id.size(), 0, upcxx::world()).wait();
    job_id.resize(sz, ' ');
    upcxx::broadcast(job_id.data(), sz, 0, upcxx::world()).wait();
  }
  return job_id;
}

Options::Options() {
  char buf[32];
  memset(buf, 0, sizeof(buf));
  if (!upcxx::rank_me()) {
    setup_time = get_current_time(true);
    strncpy(buf, setup_time.c_str(), sizeof(buf) - 1);
  }
  upcxx::broadcast(buf, sizeof(buf), 0, world()).wait();
  setup_time = string(buf);
  output_dir = string("mhm2-kmer-analysis-output-n") + to_string(upcxx::rank_n()) + "-N" +
               to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + setup_time;
}

Options::~Options() {
  flush_logger();
  cleanup();
}

void Options::cleanup() {
  // cleanup and close loggers that Options opened in load
  close_logger();
#ifdef DEBUG
  close_dbg();
#endif
}

CLI::Option *add_flag_def(CLI::App &app, string flag, bool &var, string help) {
  return app.add_flag(flag, var, help + " (" + (var ? "true" : "false") + ")");
}

bool Options::load(int argc, char **argv) {
  // MHM2 version v0.1-a0decc6-master (Release) built on 2020-04-08T22:15:40 with g++
  string full_version_str = "MHM2 version " + string(MHM2_VERSION) + "-" + string(MHM2_BRANCH) + " with upcxx-utils " +
                            string(UPCXX_UTILS_VERSION) + " built on " + string(MHM2_BUILD_DATE);
  app.description(full_version_str);
  app.option_defaults()->always_capture_default();
  // basic options - see user guide
  app.add_option("-r, --reads", reads_fnames,
                 "Files containing interleaved paired reads in FASTQ format (comma or space separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile)
      ->group("Basic options");
  app.add_option("-p, --paired-reads", paired_fnames,
                 "Alternating read files containing separate paired reads in FASTQ format (comma or space separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile)
      ->group("Basic options");
  app.add_option("-u, --unpaired-reads", unpaired_fnames, "Unpaired or single reads in FASTQ format (comma or space separated).")
      ->delimiter(',')
      ->check(CLI::ExistingFile)
      ->group("Basic options");
  app.add_option("-k, --kmer-lens", kmer_lens, "kmer lengths (comma separated) for contigging (set to 0 to disable contigging).")
      ->delimiter(',')
      ->group("Basic options");
  auto *output_dir_opt = app.add_option("-o,--output", output_dir, "Output directory.")->group("Basic options");
  add_flag_def(app, "--checkpoint", checkpoint, "Enable checkpointing.")->group("Basic options");
  add_flag_def(app, "--progress", show_progress, "Show progress bars for operations.")->group("Basic options");
  add_flag_def(app, "-v, --verbose", verbose, "Verbose output: lots of detailed information (always available in the log).")
      ->group("Basic options");
  auto *cfg_opt = app.set_config("--config", "", "Load options from a configuration file.")->group("Basic options");

  // advanced options

  // quality tuning
  add_flag_def(app, "--adapter-trim", adapter_trim, "Trim adapters using reference given by --adapter-refs")
      ->group("Quality tuning options");
  app.add_option("--adapter-refs", adapter_fname, "File containing adapter sequences for trimming in FASTA format.")
      ->check(CLI::ExistingFile)
      ->group("Quality tuning options");
  app.add_option("-i, --insert", insert_size, "Insert size (average:stddev) (autodetected by default).")
      ->delimiter(':')
      ->expected(2)
      ->check(CLI::Range(1, 10000))
      ->group("Quality tuning options");
  app.add_option("--min-depth-thres", dmin_thres, "Absolute mininimum depth threshold for DeBruijn graph traversal")
      ->check(CLI::Range(1, 100))
      ->group("Quality tuning options");
  app.add_option("--optimize", optimize_for,
                 "Optimize setting: (contiguity, correctness, default) - improve contiguity at the cost of increased errors; "
                 "reduce errors at the cost of contiguity; default balance between contiguity and correctness")
      ->check(CLI::IsMember({"default", "contiguity", "correctness"}))
      ->group("Quality tuning options");
  // performance trade-offs
  app.add_option("--max-kmer-store", max_kmer_store_mb, "Maximum size for kmer store in MB per rank (set to 0 for auto 1% memory).")
      ->check(CLI::Range(0, 5000))
      ->group("Performance trade-off options");
  app.add_option("--max-rpcs-in-flight", max_rpcs_in_flight,
                 "Maximum number of RPCs in flight, per process (set to 0 for unlimited).")
      ->check(CLI::Range(0, 10000))
      ->group("Performance trade-off options");
  app.add_option("--max-worker-threads", max_worker_threads, "Number of threads in the worker ThreadPool.")
      ->check(CLI::Range(0, 16))
      ->group("Performance trade-off options");
  if (getenv("MHM2_PIN")) pin_by = getenv("MHM2_PIN");  // get default from the environment if it exists
  app.add_option("--pin", pin_by,
                 "Restrict processes according to logical CPUs, cores (groups of hardware threads), "
                 "or NUMA domains (cpu, core, numa, none).")
      ->check(CLI::IsMember({"cpu", "core", "numa", "rr_numa", "none"}))
      ->group("Performance trade-off options");
  app.add_option("--sequencing-depth", sequencing_depth, "Expected average sequencing depth")
      ->check(CLI::Range(1, 100))
      ->group("Performance trade-off options");

  add_flag_def(app, "--use-qf", use_qf,
               "Use quotient filter to reduce memory at the cost of slower processing (only applies to GPUs).")
      ->group("Performance trade-off options");
  app.add_option("--subsample-pct", subsample_fastq_pct, "Percentage of fastq files to read.")
      ->check(CLI::Range(1, 100))
      ->group("Performance trade-off options");
  // miscellaneous
  add_flag_def(app, "--dump-merged", dump_merged, "(debugging option) dumps merged fastq files in the output directory")
      ->multi_option_policy()
      ->group("Other options");
  add_flag_def(app, "--dump-kmers", dump_kmers, "Write kmers out after kmer counting.")->group("Other options");
  app.add_option("-Q, --quality-offset", qual_offset, "Phred encoding offset (auto-detected by default).")
      ->check(CLI::IsMember({0, 33, 64}))
      ->group("Other options");

  try {
    app.parse(argc, argv);
  } catch (const CLI::CallForHelp &e) {
    app.exit(e);
    exit(0);
  } catch (const CLI::CallForAllHelp &e) {
    app.exit(e);
    exit(0);
  } catch (const CLI::ParseError &e) {
    if (upcxx::rank_me() == 0) {
      if (e.get_exit_code() != 0) cerr << "\nError (" << e.get_exit_code() << ") in command line:\n";
      app.exit(e);
    }
    return false;
  }

  {
    std::unordered_set<std::string> check_input_files;
    for (auto name : reads_fnames) {
      if (!check_input_files.insert(name).second) SDIE("Duplicate read file detected: ", name);
    }

    if (!paired_fnames.empty()) {
      // convert pairs to colon ':' separated unpaired/single files for FastqReader to process
      if (paired_fnames.size() % 2 != 0) {
        if (!rank_me()) cerr << "Did not get pairs of files in -p: " << paired_fnames.size() << endl;
        return false;
      }
      while (paired_fnames.size() >= 2) {
        if (!check_input_files.insert(paired_fnames[0]).second)
          SDIE("Duplicate first paired-read file detected: ", paired_fnames[0]);
        if (!check_input_files.insert(paired_fnames[1]).second)
          SDIE("Duplicate second paired-read file detected: ", paired_fnames[1]);

        reads_fnames.push_back(paired_fnames[0] + ":" + paired_fnames[1]);
        paired_fnames.erase(paired_fnames.begin());
        paired_fnames.erase(paired_fnames.begin());
      }
    }

    if (!unpaired_fnames.empty()) {
      // append a ':' to the file name, signaling to FastqReader that this is a unpaired/sinngle file
      // a paired file would have another name after the ':'
      for (auto name : unpaired_fnames) {
        if (!check_input_files.insert(name).second) SDIE("Duplicate unpaired-read file detected: ", name);
        reads_fnames.push_back(name + ":");
      }
      unpaired_fnames.clear();
    }
  }

  if (!adapter_trim) adapter_fname.clear();

  if (reads_fnames.empty()) {
    if (!rank_me()) cerr << "\nError in command line:\nRequire read names\nRun with --help for more information\n";
    return false;
  }

  upcxx::barrier();

  if (!*output_dir_opt) {
    string first_read_fname = reads_fnames[0];
    // strip out the paired or unpaired/single the possible ':' in the name string
    auto colpos = first_read_fname.find_last_of(':');
    if (colpos != string::npos) first_read_fname = first_read_fname.substr(0, colpos);
    first_read_fname = remove_file_ext(get_basename(first_read_fname));
    auto spos = first_read_fname.find_first_of(':');
    if (spos != string::npos) first_read_fname = first_read_fname.substr(0, spos);
    output_dir = string("mhm2-kmer-analysis-output-n") + to_string(upcxx::rank_n()) + "-N" +
                 to_string(upcxx::rank_n() / upcxx::local_team().rank_n()) + "-" + setup_time;
    output_dir_opt->default_val(output_dir);
  }

  double setup_time = 0;
  setup_time += setup_output_dir();
  setup_time += setup_log_file();

  min_kmer_len = kmer_lens;

  if (!app.get_option("--kmer-lens")->empty()) default_kmer_lens = false;

  if (max_kmer_store_mb == 0) {
    // use 1% of the minimum available memory
    max_kmer_store_mb = get_free_mem(true) / 1024 / 1024 / 100;
    max_kmer_store_mb = upcxx::reduce_all(max_kmer_store_mb / upcxx::local_team().rank_n(), upcxx::op_fast_min).wait();
  }

  if (optimize_for == "contiguity") {
  }

  auto logger_t = clock_now();
  if (upcxx::local_team().rank_me() == 0) {
    // open 1 log per node
    // all have logs in per_rank

    init_logger("mhm2.log", verbose, true);
    // if not restarting, hardlink just the rank0 log to the output dir
    // this ensures a stripe count of 1 even when the output dir is striped wide
    if (rank_me() == 0) {
      auto ret = link("per_rank/00000000/00000000/mhm2.log", "mhm2.log");
      if (ret != 0) SWARN("Could not hard link mhm2.log from per_rank/00000000/00000000/mhm2.log: ", strerror(errno), "\n");
    }
  }

  barrier();
  chrono::duration<double> logger_t_elapsed = clock_now() - logger_t;
  SLOG_VERBOSE("init_logger took ", setprecision(2), fixed, logger_t_elapsed.count(), " s at ", get_current_time(), " (",
               setup_time, " s for io)\n");

#ifdef DEBUG
  open_dbg("debug");
#endif

  SLOG(KLBLUE, full_version_str, KNORM, "\n");

// https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/developer-guide-reference/2023-2/use-predefined-macros-to-specify-intel-compilers.html
#if !defined(__INTEL_LLVM_COMPILER)
  if (upcxx::rank_me() == 0) {
    // print out all compiler definitions
    SLOG_VERBOSE(KLBLUE, "_________________________", KNORM, "\n");
    SLOG_VERBOSE(KLBLUE, "Compiler definitions:", KNORM, "\n");
    std::istringstream all_defs_ss(ALL_DEFNS);
    vector<string> all_defs((std::istream_iterator<string>(all_defs_ss)), std::istream_iterator<string>());
    for (auto &def : all_defs) SLOG_VERBOSE("  ", def, "\n");
    SLOG_VERBOSE("_________________________\n");
    SLOG(KLBLUE, "Options:", KNORM, "\n");
    auto all_opts_strings = app.config_to_str_vector(true, false);
    for (auto &opt_str : all_opts_strings) SLOG(KLBLUE, opt_str, KNORM, "\n");
    SLOG(KLBLUE, "_________________________", KNORM, "\n");
  }
#endif
  auto num_nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
  SLOG(KBLUE, "Starting run with ", upcxx::rank_n(), " processes on ", num_nodes, " node", (num_nodes > 1 ? "s" : ""), " at ",
       get_current_time(), KNORM, "\n");
#ifdef DEBUG
  SWARN("Running low-performance debug mode");
#endif
  write_config_file();
  upcxx::barrier();
  return true;
}

void Options::adjust_config_option(const string &opt_name, const string &new_val) {
  app.get_option(opt_name)->clear();
  app.get_option(opt_name)->add_result(new_val);
}

void Options::write_config_file() {
  if (!upcxx::rank_me()) {
    // write out configuration file for restarts
    ofstream ofs(config_file);
    ofs << app.config_to_str(false, true);
    ofs.close();
    unlink(linked_config_file.c_str());  // ignore errors
    auto ret = link(config_file.c_str(), linked_config_file.c_str());
    if (ret != 0) LOG("Could not hard link config file, continuing\n");
  }
}
