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

#include <sys/resource.h>
#include <unistd.h>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <iomanip>

#include "contigging.hpp"
#include "fastq.hpp"
#include "packed_reads.hpp"
#include "stage_timers.hpp"
#include "gasnet_stats.hpp"
#include "upcxx_utils.hpp"
#include "upcxx_utils/thread_pool.hpp"
#include "utils.hpp"

#include "kmer.hpp"

#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#include <limits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#endif

using std::fixed;
using std::setprecision;

using namespace upcxx_utils;

void init_devices();
void done_init_devices();
void teardown_devices();


void print_exec_cmd(int argc, char **argv) {
  string executed = argv[0];
  executed += ".py";  // assume the python wrapper was actually called
  for (int i = 1; i < argc; i++) executed = executed + " " + argv[i];
    SLOG_VERBOSE("Executed as: ", executed, "\n");
}


int merge_reads(vector<string> reads_fname_list, int qual_offset, double &elapsed_write_io_t, PackedReadsList &packed_reads_list,
                bool checkpoint, const string &adapter_fname, int min_kmer_len, int subsample_pct, bool use_blastn_scores);

void set_process_affinity(const string pin_by) {
  SLOG_VERBOSE("Process 0 on node 0 is initially pinned to ", get_proc_pin(), "\n");
  // pin ranks only in production
  if (pin_by == "cpu")
    pin_cpu();
  else if (pin_by == "core")
    pin_core();
  else if (pin_by == "numa")
    pin_numa();
  else if (pin_by == "rr_numa")
    pin_numa(true);
  log_pins();
}

void set_thread_pool(int max_worker_threads) {
  const int num_threads = max_worker_threads;  // reserve up to threads in the singleton thread pool.
  upcxx_utils::ThreadPool::get_single_pool(num_threads);
  // FIXME if (!max_worker_threads) upcxx_utils::FASRPCCounts::use_worker_thread() = false;
  SLOG_VERBOSE("Allowing up to ", num_threads, " extra threads in the thread pool\n");
}

void update_rlimits(int num_input_files) {
  // update rlimits on RLIMIT_NOFILE files if necessary
  if (num_input_files > 1) {
    struct rlimit limits;
    int status = getrlimit(RLIMIT_NOFILE, &limits);
    if (status == 0) {
      limits.rlim_cur = std::min(limits.rlim_cur + num_input_files * 8, limits.rlim_max);
      status = setrlimit(RLIMIT_NOFILE, &limits);
      SLOG_VERBOSE("Set RLIMIT_NOFILE to ", limits.rlim_cur, "\n");
    }
    if (status != 0) SWARN("Could not get/set rlimits for NOFILE\n");
  }
}

void calc_input_files_size(const vector<string> &reads_fnames) {
  auto nodes = upcxx::rank_n() / upcxx::local_team().rank_n();
  auto total_free_mem = get_free_mem(true) * nodes;
  if (!upcxx::rank_me()) {
    // get total file size across all libraries
    double tot_file_size = 0;
    for (auto const &reads_fname : reads_fnames) {
      auto spos = reads_fname.find_first_of(':');  // support paired reads
      if (spos == string::npos) {
        auto sz = get_file_size(reads_fname);
        SLOG(KBLUE, "Reads file ", reads_fname, " is ", get_size_str(sz), KNORM, "\n");
        tot_file_size += sz;
      } else {
        // paired files
        auto r1 = reads_fname.substr(0, spos);
        auto s1 = get_file_size(r1);
        auto r2 = reads_fname.substr(spos + 1);
        auto s2 = get_file_size(r2);
        SLOG("Paired files ", r1, " and ", r2, " are ", get_size_str(s1), " and ", get_size_str(s2), "\n");
        tot_file_size += s1 + s2;
      }
    }
    SOUT(KBLUE, "Total size of ", reads_fnames.size(), " input file", (reads_fnames.size() > 1 ? "s" : ""), " is ",
         get_size_str(tot_file_size), "; ", get_size_str(tot_file_size / rank_n()), " per rank; ",
         get_size_str(local_team().rank_n() * tot_file_size / rank_n()), " per node", KNORM, "\n");

    if (total_free_mem < 3 * tot_file_size)
      SWARN("There may not be enough memory in this job of ", nodes,
            " nodes for this amount of data.\n\tTotal free memory is approx ", get_size_str(total_free_mem),
            " and should be at least 3x the data size of ", get_size_str(tot_file_size), "\n");
  }
}

void run_contigging(Options &options, PackedReadsList &packed_reads_list, int &rlen_limit) {
  BarrierTimer("Start Contigging");

  // contigging loops
  // if (options.kmer_lens.size()) {
  //   for (auto kmer_len : options.kmer_lens) {
  //     if (kmer_len <= 1) continue;  // short circuit to just load reads
  int kmer_len = options.kmer_lens;
  
  auto max_k = (kmer_len / 32 + 1) * 32;
  LOG(upcxx_utils::GasNetVars::getUsedShmMsg(), "\n");



#define CONTIG_K(KMER_LEN) \
  case KMER_LEN: contigging<KMER_LEN>(kmer_len, rlen_limit, packed_reads_list, options); break


      switch (max_k) {
        CONTIG_K(32);
#if MAX_BUILD_KMER >= 64
        CONTIG_K(64);
#endif
#if MAX_BUILD_KMER >= 96
        CONTIG_K(96);
#endif
#if MAX_BUILD_KMER >= 128
        CONTIG_K(128);
#endif
#if MAX_BUILD_KMER >= 160
        CONTIG_K(160);
#endif
        default: DIE("Built for max k = ", MAX_BUILD_KMER, " not k = ", max_k);
      }
#undef CONTIG_K

    // }
  // }

}

void run_pipeline(Options &options, MemoryTrackerThread &memory_tracker, timepoint_t start_t) {
  memory_tracker.start();
  LOG_MEM("Preparing to load reads\n");
  PackedReadsList packed_reads_list;
  for (auto const &reads_fname : options.reads_fnames) {
    packed_reads_list.push_back(new PackedReads(options.qual_offset, get_merged_reads_fname(reads_fname)));
  }
  LOG_MEM("Opened read files\n");
  auto before_merge_mem = get_free_mem(true);

  double elapsed_write_io_t = 0;

  // merge the reads and insert into the packed reads memory cache (always do this even for restarts)
  begin_gasnet_stats("merge_reads");
  stage_timers.merge_reads->start();
  auto avg_read_len =
      merge_reads(options.reads_fnames, options.qual_offset, elapsed_write_io_t, packed_reads_list, options.dump_merged,
                  options.adapter_fname, options.min_kmer_len, options.subsample_fastq_pct, options.optimize_for == "contiguity");
  stage_timers.merge_reads->stop();
  end_gasnet_stats();
  auto after_merge_mem = get_free_mem(true);
  SLOG_VERBOSE(KBLUE, "Cache used ", setprecision(2), fixed, get_size_str(before_merge_mem - after_merge_mem),
               " memory on node 0 for reads", KNORM, "\n");

  // if (avg_read_len < 110 && !options.restart && options.default_kmer_lens) {
  //   options.kmer_lens.pop_back();
  //   if (options.default_scaff_kmer_lens) options.scaff_kmer_lens.front() = options.kmer_lens.back();
  //   SOUT("Average read length is ", avg_read_len, ". Adjusting value of k:\n");
  //   SOUT("  kmer-lens = ", Options::vec_to_str(options.kmer_lens), "\n");
  // }
  options.adjust_config_option("--kmer-lens", to_string(options.kmer_lens));
  // keep track of all the changes in the config file
  options.write_config_file();

  int rlen_limit = 0;
  for (auto packed_reads : packed_reads_list) {
    rlen_limit = max(rlen_limit, (int)packed_reads->get_max_read_len());
    packed_reads->report_size();
  }
  Timings::wait_pending();  // report all I/O stats here

  std::chrono::duration<double> init_t_elapsed = clock_now() - start_t;
  SLOG("\n");
  auto post_init_free_mem = get_free_mem(true);
  SLOG(KBLUE, "Completed initialization in ", setprecision(2), fixed, init_t_elapsed.count(), " s at ", get_current_time(), " (",
       get_size_str(post_init_free_mem), " free memory on node 0)", KNORM, "\n");

#ifndef ENABLE_KOKKOS
  done_init_devices();
#endif 

  run_contigging(options, packed_reads_list, rlen_limit);

  // cleanup
  LOG_MEM("Preparing to close all fastq");
  FastqReaders::close_all();  // needed to cleanup any open files in this singleton
  auto finalization_start_t = clock_now();
  for (auto packed_reads : packed_reads_list) {
    delete packed_reads;
  }
  packed_reads_list.clear();
  LOG_MEM("Closed all fastq");
  // output final assembly
  // SLOG(KBLUE "_________________________", KNORM, "\n");
  // stage_timers.dump_ctgs->start();
  // ctgs.dump_contigs("final_assembly.fasta", options.min_ctg_print_len, "scaffold_");
  // stage_timers.dump_ctgs->stop();

  // SLOG(KBLUE "_________________________", KNORM, "\n");
  // ctgs.print_stats(options.min_ctg_print_len);
  // std::chrono::duration<double> fin_t_elapsed = clock_now() - finalization_start_t;
  // SLOG("\n");
  // auto post_finalize_free_mem = get_free_mem(true);
  // SLOG(KBLUE, "Completed finalization in ", setprecision(2), fixed, fin_t_elapsed.count(), " s at ", get_current_time(), " (",
  //      get_size_str(post_finalize_free_mem), " free memory on node 0)", KNORM, "\n");

  SLOG(KBLUE "_________________________", KNORM, "\n");
  SLOG("Stage timing:\n");
  SLOG("    Initialization: ", init_t_elapsed.count(), "\n");
  // if (!options.restart)
    SLOG("    ", stage_timers.merge_reads->get_final(), "\n");
  // else
  //   SLOG("    ", stage_timers.cache_reads->get_final(), "\n");
  SLOG("    ", stage_timers.analyze_kmers->get_final(), "\n");
  SLOG("    FASTQ total read time: ", FastqReader::get_io_time(), "\n");
  SLOG("    merged FASTQ write time: ", elapsed_write_io_t, "\n");
  SLOG("    Contigs write time: ", stage_timers.dump_ctgs->get_elapsed(), "\n");
  SLOG(KBLUE "_________________________", KNORM, "\n");
  memory_tracker.stop();
  std::chrono::duration<double> t_elapsed = clock_now() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(), " for ", MHM2_VERSION, "\n");

}


string init_upcxx(BaseTimer &total_timer) {
  total_timer.start();
  // capture the free memory and timers before upcxx::init is called
  auto starting_free_mem = get_free_mem();
  char *proc_id = getenv("OMPI_COMM_WORLD_NODE_RANK");
  if (!proc_id) proc_id = getenv("SLURM_PROCID");
  int my_rank = -1;
  if (proc_id) my_rank = atol(proc_id);
  BaseTimer init_timer("upcxx::init");
  BaseTimer first_barrier("FirstBarrier");
  init_timer.start();
  if (!my_rank) {
    char hnbuf[64];
    gethostname(hnbuf, sizeof(hnbuf) - 1);
    std::cout << "Starting Rank0 with " << get_size_str(starting_free_mem) << " on " << hnbuf << " pid=" << getpid() << " at "
              << get_current_time() << std::endl;
  }
  upcxx::init();
  auto init_entry_msm_fut = init_timer.reduce_start();
  init_timer.stop();
  auto init_timings_fut = init_timer.reduce_timings();
  // upcxx::promise<> prom_report_init_timings(1);

  // we wish to have all ranks start at the same time to determine actual timing
  first_barrier.start();
  barrier();
  first_barrier.stop();
  auto post_init_free_mem = get_free_mem();
  barrier(local_team());
  auto msm_starting_free_mem_fut = min_sum_max_reduce_one((float)starting_free_mem / ONE_GB, 0);
  auto msm_post_init_free_mem_fut = min_sum_max_reduce_one((float)post_init_free_mem / ONE_GB, 0);

  auto fut_report_init_timings =
      when_all(/*prom_report_init_timings.get_future(),*/ init_entry_msm_fut, init_timings_fut, first_barrier.reduce_timings(),
               msm_starting_free_mem_fut, msm_post_init_free_mem_fut)
          .then([&total_timer](const upcxx_utils::MinSumMax<double> &entry_msm, upcxx_utils::ShTimings sh_timings,
                               upcxx_utils::ShTimings sh_first_barrier_timings,
                               const upcxx_utils::MinSumMax<float> &starting_mem_msm,
                               const upcxx_utils::MinSumMax<float> &post_init_mem_msm) {
            LOG("Time since start including init, first barrier and reductions: ", total_timer.get_elapsed_since_start(), "\n");
            return string("upcxx::init Before=") + entry_msm.to_string() + string("\nupcxx::init After=") +
                   sh_timings->to_string() + string("\nupcxx::init FirstBarrier=") + sh_first_barrier_timings->to_string() +
                   string("\nupcxx::init Starting RAM=") + starting_mem_msm.to_string() + string(" GB\nupcxx::init Post RAM=") +
                   post_init_mem_msm.to_string() + string(" GB\n");
          });

  // prom_report_init_timings.fulfill_anonymous(1);
  return fut_report_init_timings.wait();
}

void print_log_results(string output_dir, double elapsed_time) {
    // parse proxy results from log file and print
  

    ifstream log_file("mhm2.log");
    if (!log_file) {
      perror("Cannot open mhm2.log");
      return;
    }

    cout << "\n\n\n----------------------------\n\n" <<
    "proxy_results_summary.csv can be found in " << output_dir <<
    "\n\n" << endl;

    string line, token, total_kmers, unique_kmers, mem, read_count;

    std::unordered_map<string, std::pair<int, string*>> results_map = {
      {"tot_num_reads", {6, &read_count}},
      {"Total kmer count sum", {7, &total_kmers}},
      {"Total kmers", {5, &unique_kmers}},
      {"Peak memory", {10, &mem}}
    };

    while (std::getline(log_file, line)) {
      for (const auto& [phrase, pair] : results_map) {
        if (line.find(phrase) != string::npos) {
          std::istringstream iss(line);
          for (int i = 0; i < pair.first; i++) {
            iss >> token;
          }
          *pair.second = token;
          // fix mem token e.g. 4GB -> 4 GB
          // fix reads token
          if (phrase == "Peak memory") { pair.second->erase(pair.second->size() - 2); }
          else if (phrase == "tot_num_reads") { pair.second->erase(0, 14); }
          break;
        }
      }
    }

    //cout << "uniq test: " << unique_kmers << "\ntot test: " << total_kmers << "\nread test: " << read_count << endl;

    double frac = std::stold(unique_kmers) / std::stold(total_kmers);

    ofstream csv("proxy_results_summary.csv");

    // column headers
    csv << "Reads,Unique kmers,Total kmers,Fraction of Unique Kmers,Peak Memory (GB),Timing (seconds)\n";

    // results
    csv << read_count << "," <<
            unique_kmers << "," <<
            total_kmers << "," <<
            std::fixed << std::setprecision(3) << frac << "," <<
            std::setprecision(2) << mem << "," <<
            elapsed_time << "\n";


    csv.close();
 

}








void print_log_results(string output_dir, double elapsed_time) {
    // parse proxy results from log file and print
  

    ifstream log_file("mhm2.log");
    if (!log_file) {
      perror("Cannot open mhm2.log");
      return;
    }

    cout << "\n\n\n----------------------------\n\n" <<
    "proxy_results_summary.csv can be found in " << output_dir <<
    "\n\n" << endl;

    string line, token, total_kmers, unique_kmers, mem, read_count;

    std::unordered_map<string, std::pair<int, string*>> results_map = {
      {"tot_num_reads", {6, &read_count}},
      {"Total kmer count sum", {7, &total_kmers}},
      {"Total kmers", {5, &unique_kmers}},
      {"Peak memory", {10, &mem}}
    };

    while (std::getline(log_file, line)) {
      for (const auto& [phrase, pair] : results_map) {
        if (line.find(phrase) != string::npos) {
          std::istringstream iss(line);
          for (int i = 0; i < pair.first; i++) {
            iss >> token;
          }
          *pair.second = token;
          // fix mem token e.g. 4GB -> 4 GB
          // fix reads token
          if (phrase == "Peak memory") { pair.second->erase(pair.second->size() - 2); }
          else if (phrase == "tot_num_reads") { pair.second->erase(0, 14); }
          break;
        }
      }
    }

    //cout << "uniq test: " << unique_kmers << "\ntot test: " << total_kmers << "\nread test: " << read_count << endl;

    double frac = std::stold(unique_kmers) / std::stold(total_kmers);

    ofstream csv("proxy_results_summary.csv");

    // column headers
    csv << "Reads,Unique kmers,Total kmers,Fraction of Unique Kmers,Peak Memory (GB),Timing (seconds)\n";

    // results
    csv << read_count << "," <<
            unique_kmers << "," <<
            total_kmers << "," <<
            std::fixed << std::setprecision(3) << frac << "," <<
            std::setprecision(2) << mem << "," <<
            elapsed_time << "\n";


    csv.close();
 

}






int main(int argc, char **argv, char **envp) {
  
  BaseTimer total_timer("Total Time", nullptr);  // no PromiseReduce possible
  auto init_timings = init_upcxx(total_timer);
  auto am_root = !rank_me();
  Options options;

#ifdef ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
   double kokkos_elapsed_time;
  Kokkos::Timer kokkos_timer;
  {
#endif

#ifdef CIDS_FROM_HASH
  SWARN("Generating contig IDs with hashing - this could result in duplicate CIDs and should only be used for checking "
        "consistency across small runs");
#endif

  const char *gasnet_statsfile = getenv("GASNET_STATSFILE");

#if defined(ENABLE_GASNET_STATS)
  if (gasnet_statsfile) _gasnet_stats = true;
#else
  if (gasnet_statsfile) SWARN("No GASNet statistics will be collected - use Debug or RelWithDebInfo modes to enable collection.");
#endif

  srand(rank_me() + 10);
  auto start_t = clock_now();
  
  // if we don't load, return "command not found"
  if (!options.load(argc, argv)) return 127;
  print_exec_cmd(argc, argv);
  vector<string> envvars;
  for (char **env = envp; *env != 0; env++) {
    if (!strncmp(*env, "GASNET", strlen("GASNET"))) envvars.push_back(string(*env));
    if (!strncmp(*env, "UPCXX", strlen("UPCXX"))) envvars.push_back(string(*env));
    if (!strncmp(*env, "FI_", strlen("FI_"))) envvars.push_back(string(*env));
  }
  sort(envvars.begin(), envvars.end());
  SLOG_VERBOSE(KBLUE "_________________________", KNORM, "\n");
  SLOG_VERBOSE(KLBLUE "GASNet/UPCXX/OFI environment variables:", KNORM, "\n");
  for (auto envvar : envvars) SLOG_VERBOSE(KBLUE, "  ", envvar, KNORM, "\n");
  SLOG_VERBOSE(KBLUE "_________________________", KNORM, "\n");
  SLOG_VERBOSE(KLCYAN, "Timing reported as min/my/average/max, balance", KNORM, "\n");
  // only write them here to honor the verbose flag in options
  SLOG_VERBOSE(init_timings);
  ProgressBar::SHOW_PROGRESS = options.show_progress;
  set_process_affinity(options.pin_by);
  log_env();
  update_rlimits(options.reads_fnames.size());
  set_thread_pool(options.max_worker_threads);
  calc_input_files_size(options.reads_fnames);

#ifndef ENABLE_KOKKOS
  init_devices();
#endif

  MemoryTrackerThread memory_tracker;  // write only to mhm2.log file(s), not a separate one too

  run_pipeline(options, memory_tracker, start_t);

  LOG("Cleaning up and completing remaining tasks\n");

  upcxx_utils::ThreadPool::join_single_pool();  // cleanup singleton thread pool
  upcxx_utils::Timings::wait_pending();         // ensure all outstanding timing summaries have printed
  LOG("Done waiting for all pending.\n");
  barrier();
  LOG("All ranks done. Flushing logs and finalizing.\n");

  BaseTimer flush_logs_timer("flush_logger", nullptr);  // no PromiseReduce possible
  flush_logs_timer.start();

#ifdef DEBUG
  _dbgstream.flush();
  while (close_dbg());
  LOG("closed DBG.\n");
#endif

  if (am_root)
    upcxx_utils::flush_logger();
  else
    upcxx_utils::close_logger();

  flush_logs_timer.stop();
  auto sh_flush_timings = flush_logs_timer.reduce_timings().wait();
  barrier();
  SLOG_VERBOSE("Total time before close and finalize: ", total_timer.get_elapsed_since_start(), "\n");
  SLOG_VERBOSE("All ranks flushed logs: ", sh_flush_timings->to_string(), "\n");

#ifdef ENABLE_KOKKOS
  //test_kokkos();
  }
#endif

  BaseTimer finalize_timer("upcxx::finalize", nullptr);  // no PromiseReduce possible
  finalize_timer.start();
  upcxx::finalize();

#ifdef ENABLE_KOKKOS
  kokkos_elapsed_time = kokkos_timer.seconds();
  Kokkos::finalize();
#endif

  finalize_timer.stop();
  total_timer.stop();
  if (am_root) {
    cout << "Total time: " << fixed << setprecision(3) << total_timer.get_elapsed() << " s (upcxx::finalize in "
         << finalize_timer.get_elapsed() << " s)" << endl;

#ifdef ENABLE_KOKKOS   
    print_log_results(options.output_dir, kokkos_elapsed_time);
#else
    print_log_results(options.output_dir, total_timer.get_elapsed());
#endif
  }
  

  return 0;
}
