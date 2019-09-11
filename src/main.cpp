//mhm - UPC++ version
//Steven Hofmeyr, LBNL, June 2019

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <chrono>
#include <upcxx/upcxx.hpp>

using namespace std;

#ifndef NDEBUG
#define DEBUG
ofstream _dbgstream;
#endif

#include "utils.hpp"
#include "options.hpp"
#include "merge_reads.hpp"
#include "kcount.hpp"
#include "dbjg_traversal.hpp"
#include "klign.hpp"

using namespace std;
using namespace upcxx;

unsigned int Kmer::k = 0;


int main(int argc, char **argv)
{
  upcxx::init();
  auto start_t = chrono::high_resolution_clock::now();

#ifdef DEBUG
  //time_t curr_t = std::time(nullptr);
  //string dbg_fname = "debug" + to_string(curr_t) + ".log";
  string dbg_fname = "debug.log";
  get_rank_path(dbg_fname, rank_me());
  _dbgstream.open(dbg_fname);
#endif

  SOUT("MHM version ", MHM_VERSION, "\n");
  double start_mem_free = get_free_mem_gb();
  SOUT("Initial free memory on node 0: ", setprecision(3), fixed, start_mem_free, " GB\n");
  SOUT("Running on ", rank_n(), " ranks\n");

  // print out all compiler definitions
  SOUT(KBLUE "_________________________\nCompiler definitions:\n");
  istringstream all_defs_ss(ALL_DEFNS);
  vector<string> all_defs((istream_iterator<string>(all_defs_ss)), istream_iterator<string>());
  for (auto &def : all_defs) SOUT("  ", def, "\n");
  
#ifdef DEBUG
  SOUT(KLRED "WARNING: Running low-performance debug mode\n", KNORM);
#endif
  
  auto options = make_shared<Options>();
  options->load(argc, argv);

  // get total file size across all libraries
  double tot_file_size = 0;
  if (!rank_me()) {
    for (auto const &reads_fname : options->reads_fname_list) tot_file_size += get_file_size(reads_fname);
    SOUT("Total size of ", options->reads_fname_list.size(), " input file", (options->reads_fname_list.size() > 1 ? "s" : ""),
         " is ", (tot_file_size / ONE_GB), " GB\n");
  }

  // first merge reads - the results will go in the per_rank directory
  merge_reads(options->reads_fname_list, options->qual_offset);

  Contigs ctgs;
  for (auto kmer_len : options->kmer_lens) {
    auto loop_start_t = chrono::high_resolution_clock::now();
    SOUT(KBLUE "_________________________\nContig generation k = ", kmer_len, "\n\n", KNORM);
    Kmer::k = kmer_len;
    {
      // scope is to ensure that kmer_dht is freed by destructor
      auto my_num_kmers = estimate_num_kmers(kmer_len, options->reads_fname_list);
      dist_object<KmerDHT> kmer_dht(world(), my_num_kmers, options->max_kmer_store, options->use_bloom);
      barrier();
      analyze_kmers(kmer_len, options->qual_offset, options->reads_fname_list, options->use_bloom, options->min_depth_cutoff,
                    options->dynamic_min_depth, ctgs, kmer_dht);
      barrier();
      traverse_debruijn_graph(kmer_len, kmer_dht, ctgs);
      // FIXME: dump as single file if checkpoint option is specified
      ctgs.dump_contigs("uutigs", kmer_len);
    }
    find_alignments(kmer_len, CONTIG_SEED_SPACE, options->reads_fname_list, options->max_kmer_store, ctgs);
    
    chrono::duration<double> loop_t_elapsed = chrono::high_resolution_clock::now() - loop_start_t;
    SOUT("Completed contig round k = ", kmer_len, " in ", setprecision(2), fixed, loop_t_elapsed.count(), " s at ",
         get_current_time(), ", free memory ", get_free_mem_gb(), " GB\n");
    barrier();
  }
  SOUT(KBLUE "_________________________\nScaffolding\n\n", KNORM);
  SOUT(KBLUE "_________________________\n\n", KNORM);
  double end_mem_free = get_free_mem_gb();
  SOUT("Final free memory on node 0: ", setprecision(3), fixed, end_mem_free,
       " GB (unreclaimed ", (start_mem_free - end_mem_free), " GB)\n");
  
  chrono::duration<double> t_elapsed = chrono::high_resolution_clock::now() - start_t;
  SOUT("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(), "\n"); 
  barrier();

#ifdef DEBUG
  _dbgstream.flush();
  _dbgstream.close();
#endif
  barrier();
  upcxx::finalize();
  return 0;
}

