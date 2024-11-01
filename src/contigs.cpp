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

#include "contigs.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <string>
#include <upcxx/upcxx.hpp>
#include <vector>

#include "upcxx_utils/log.hpp"
#include "upcxx_utils/ofstream.hpp"
#include "upcxx_utils/progress_bar.hpp"
#include "upcxx_utils/timers.hpp"
#include "utils.hpp"

using std::endl;
using std::max;
using std::memory_order_relaxed;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

using upcxx::barrier;
using upcxx::dist_object;
using upcxx::op_fast_add;
using upcxx::op_fast_max;
using upcxx::promise;
using upcxx::rank_me;
using upcxx::rank_n;
using upcxx::reduce_one;
using upcxx::rpc;
using upcxx::rpc_ff;

using namespace upcxx_utils;

void Contigs::clear() {
  contigs.clear();
  vector<Contig>().swap(contigs);
  begin_idx = 0;
  end_idx = 0;
  max_clen = 0;
  tot_length = 0;
}

void Contigs::set_capacity(int64_t sz) { contigs.reserve(sz); }

void Contigs::add_contig(const Contig &contig) {
  contigs.push_back(contig);
  end_idx++;
  tot_length += contig.seq.length();
}

void Contigs::add_contig(Contig &&contig) {
  contigs.emplace_back(std::move(contig));
  end_idx++;
  tot_length += contig.seq.length();
}

size_t Contigs::size() const { return end_idx - begin_idx; }

size_t Contigs::get_length() const { return tot_length; }

void Contigs::print_stats(unsigned min_ctg_len) const {
  BarrierTimer timer(__FILEFUNC__);
  int64_t tot_len = 0, max_len = 0;
  double tot_depth = 0;
  vector<pair<unsigned, uint64_t>> length_sums({{1, 0}, {5, 0}, {10, 0}, {25, 0}, {50, 0}});
  int64_t num_ctgs = 0;
  int64_t num_ns = 0;
  double tot_depth_per_base = 0;
  vector<unsigned> lens;
  lens.reserve(contigs.size());
  for (auto ctg : contigs) {
    auto len = ctg.seq.length();
    if (len < min_ctg_len) continue;
    num_ctgs++;
    tot_len += len;
    tot_depth += ctg.depth;
    tot_depth_per_base += (ctg.depth * len);
    max_len = max(max_len, static_cast<int64_t>(len));
    for (auto &length_sum : length_sums) {
      if (len >= length_sum.first * 1000) length_sum.second += len;
    }
    num_ns += count(ctg.seq.begin(), ctg.seq.end(), 'N');
    lens.push_back(len);
  }
  /*
  // Compute local N50 and then take the median across all of them. This gives a very good approx of the exact N50 and is much
  // cheaper to compute
  // Actually, the approx is only decent if there are not too many ranks
  sort(lens.rbegin(), lens.rend());
  int64_t sum_lens = 0;
  int64_t n50 = 0;
  for (auto l : lens) {
    sum_lens += l;
    if (sum_lens >= tot_len / 2) {
      n50 = l;
      break;
    }
  }
  barrier();
  dist_object<int64_t> n50_val(n50);
  int64_t median_n50 = 0;
  if (!rank_me()) {
    vector<int64_t> n50_vals;
    n50_vals.push_back(n50);
    for (int i = 1; i < rank_n(); i++) {
      n50_vals.push_back(n50_val.fetch(i).wait());
      //SOUT(i, " n50 fetched ", n50_vals.back(), "\n");
    }
    sort(n50_vals.begin(), n50_vals.end());
    median_n50 = n50_vals[rank_n() / 2];
  }
  // barrier to ensure the other ranks dist_objects don't go out of scope before rank 0 is done
  barrier();
  */
  int64_t all_num_ctgs = reduce_one(num_ctgs, op_fast_add, 0).wait();
  int64_t all_tot_len = reduce_one(tot_len, op_fast_add, 0).wait();
  int64_t all_max_len = reduce_one(max_len, op_fast_max, 0).wait();
  double all_tot_depth = reduce_one(tot_depth, op_fast_add, 0).wait();
  int64_t all_num_ns = reduce_one(num_ns, op_fast_add, 0).wait();
  double all_tot_depth_per_base = reduce_one(tot_depth_per_base, op_fast_add, 0).wait() / all_tot_len;
  // int64_t all_n50s = reduce_one(n50, op_fast_add, 0).wait();

  SLOG("Assembly statistics (contig lengths >= ", min_ctg_len, ")\n");
  SLOG("    Number of contigs:       ", all_num_ctgs, "\n");
  SLOG("    Total assembled length:  ", all_tot_len, "\n");
  SLOG("    Average contig depth:    ", all_tot_depth / all_num_ctgs, "\n");
  SLOG("    Average depth per base:  ", all_tot_depth_per_base, "\n");
  SLOG("    Number of Ns/100kbp:     ", (double)all_num_ns * 100000.0 / all_tot_len, " (", all_num_ns, ")", KNORM, "\n");
  // SLOG("    Approx N50 (average):    ", all_n50s / rank_n(), " (rank 0 only ", n50, ")\n");
  // SLOG("    Approx N50:              ", median_n50, "\n");
  SLOG("    Max. contig length:      ", all_max_len, "\n");
  SLOG("    Contig lengths:\n");
  for (auto &length_sum : length_sums) {
    SLOG("        > ", std::left, std::setw(19),
         to_string(length_sum.first) + "kbp:", perc_str(reduce_one(length_sum.second, op_fast_add, 0).wait(), all_tot_len), "\n");
  }
}

void Contigs::dump_contigs(const string &fname, unsigned min_ctg_len, const string &prefix) {
  BarrierTimer timer(__FILEFUNC__);
  dist_ofstream of(fname);
  of << std::setprecision(3);
  size_t my_num_ctgs = 0;
  for (auto it = contigs.begin(); it != contigs.end(); ++it) {
    auto ctg = it;
    if (ctg->seq.length() < min_ctg_len) continue;
    of << ">" << prefix << to_string(ctg->id);
    if (fname != "final_assembly.fasta") of << " " << ctg->depth;
    of << "\n";
    string rc_uutig = revcomp(ctg->seq);
    string seq = (rc_uutig < ctg->seq ? rc_uutig : ctg->seq);
    // for (int64_t i = 0; i < ctg->seq.length(); i += 50) of << ctg->seq.substr(i, 50) << "\n";
    of << seq << "\n";
    my_num_ctgs++;
  }
  of.close();  // sync and output stats
#ifdef DEBUG
  // two important things here.
  // 1: touch and test the load_contigs code when debugging
  // 2: ensure restarts keep identical contigs in the ranks when debugging after load_contigs balances the input
  SLOG_VERBOSE("Reloading contigs from file to rebalance\n");
  auto num_prev_ctgs = reduce_one(my_num_ctgs, op_fast_add, 0).wait();
  load_contigs(fname, prefix);
  // check the correct number of contigs was loaded again
  auto num_ctgs = reduce_one(contigs.size(), op_fast_add, 0).wait();
  if (rank_me() == 0 && num_prev_ctgs != num_ctgs) SDIE("Saved ", num_prev_ctgs, " but only loaded ", num_ctgs, " contigs");
#endif
}

void Contigs::load_contigs(const string &ctgs_fname, const string &prefix) {
  auto get_file_offset_for_rank = [ctgs_fname](ifstream &f, int rank, string &ctg_prefix, size_t file_size) -> size_t {
    if (rank == 0) return 0;
    assert(rank < rank_n());
    size_t offset = file_size / rank_n() * rank;
    f.seekg(offset);
    string line;
    while (getline(f, line)) {
      if (substr_view(line, 0, ctg_prefix.size()) == ctg_prefix) {
        LOG("in file ", ctgs_fname, " found ctg_prefix ", ctg_prefix, " line ", line);
        getline(f, line);
        return f.tellg();
      }
    }
    DIE("Could not find prefix ", ctg_prefix, " in file ", ctgs_fname);
    return 0;
  };

  SLOG_VERBOSE("Loading contigs from fasta file ", ctgs_fname, "\n");
  BarrierTimer timer(__FILEFUNC__);
  clear();
  dist_object<upcxx::promise<size_t>> dist_stop_prom(world());
  string line;
  string ctg_prefix = ">" + prefix;
  string cname, seq, buf;
  size_t bytes_read = 0;
  size_t file_size = 0;

  // broadcast the file size
  if (rank_me() == 0) file_size = upcxx_utils::get_file_size(ctgs_fname);
  file_size = upcxx::broadcast(file_size, 0).wait();
  DBG("Got filesize=", file_size, "\n");
  ifstream ctgs_file(ctgs_fname);
  if (!ctgs_file.is_open()) DIE("Could not open ctgs file '", ctgs_fname, "': ", strerror(errno));

  auto start_offset = get_file_offset_for_rank(ctgs_file, rank_me(), ctg_prefix, file_size);
  if (rank_me() > 0) {
    // notify previous rank of its stop offset
    rpc_ff(
        rank_me() - 1,
        [](dist_object<promise<size_t>> &dist_stop_prom, size_t stop_offset) {
          dist_stop_prom->fulfill_result(stop_offset);
          LOG("received stop_offset=", stop_offset, "\n");
        },
        dist_stop_prom, start_offset);
    LOG("Sent my start_offset to ", rank_me() - 1, " ", start_offset, "\n");
  }
  if (rank_me() == rank_n() - 1) {
    dist_stop_prom->fulfill_result(file_size);
  }
  auto stop_offset = dist_stop_prom->get_future().wait();
  LOG("Got my stop_offset=", stop_offset, "\n");

  size_t tot_len = 0;
  ProgressBar progbar(stop_offset - start_offset, "Parsing contigs");
  // these can be equal if the contigs are very long and there are many ranks so this one doesn't get even a full contig
  ctgs_file.seekg(start_offset);
  int64_t prev_id = -1;
  while (!ctgs_file.eof()) {
    if ((size_t)ctgs_file.tellg() >= stop_offset) break;
    getline(ctgs_file, cname);
    if (cname == "") break;
    getline(ctgs_file, seq);
    if (seq == "") break;
    tot_len += seq.length();
    bytes_read += cname.length() + seq.length();
    progbar.update(bytes_read);
    // extract the id
    char *endptr;
    int64_t id = strtol(cname.c_str() + ctg_prefix.length(), &endptr, 10);
    if (id == prev_id) DIE("Duplicate ids ", prev_id, " in file ", ctgs_fname);
    prev_id = id;
    // depth is the last field in the cname
    double depth = (ctgs_fname == "final_assembly.fasta" ? 0 : strtod(endptr, NULL));
    Contig contig = {.id = id, .seq = seq, .depth = depth};
    add_contig(contig);
    max_clen = max(max_clen, (int)seq.length());
  }
  if (ctgs_file.tellg() < stop_offset)
    DIE("Did not read the entire contigs file from ", start_offset, " to ", stop_offset, " tellg=", ctgs_file.tellg());
  LOG("Got contigs=", contigs.size(), " tot_len=", tot_len, " max_clen=", max_clen, "\n");
  auto &pr = Timings::get_promise_reduce();
  auto fut_tot_contigs = pr.reduce_one(contigs.size(), op_fast_add, 0);
  auto fut_tot_len = pr.reduce_one(tot_len, op_fast_add, 0);
  auto fut_done = progbar.set_done();
  auto fut_report = when_all(fut_tot_contigs, fut_tot_len, fut_done).then([ctgs_fname](uint64_t tot_contigs, uint64_t tot_len) {
    SLOG_VERBOSE("Loaded ", tot_contigs, " contigs (", get_size_str(tot_len), ") from ", ctgs_fname, "\n");
  });
  Timings::set_pending(fut_report);
  // implicit exit barrrier from BarrierTimer
}

size_t Contigs::get_num_ctg_kmers(int kmer_len) const {
  size_t num_ctg_kmers = 0;
  for (auto &ctg : contigs) {
    if (ctg.seq.length() > kmer_len) num_ctg_kmers += (ctg.seq.length() - kmer_len + 1);
  }
  return num_ctg_kmers;
}

int Contigs::get_max_clen() const { return max_clen; }

void Contigs::set_next_slice(int num_slices) {
  if (end_idx == contigs.size() && begin_idx == 0) end_idx = 0;  // first call
  begin_idx = end_idx;
  int expected_slice_size = tot_length / num_slices;
  size_t slice_size = 0;
  for (size_t i = begin_idx; i < contigs.size(); i++) {
    slice_size += contigs[i].seq.length();
    if (slice_size >= expected_slice_size) {
      end_idx = i + 1;
      break;
    }
  }
  if (end_idx == begin_idx) end_idx = contigs.size();
  auto msm_slice_size = min_sum_max_reduce_one(slice_size, 0).wait();
  SLOG_VERBOSE("Rank 0, contig slice: begin idx ", begin_idx, " end idx ", end_idx, " num ctgs ", contigs.size(), "\n");
  SLOG_VERBOSE("Rank 0, expected contig slice size ", expected_slice_size, " actual slice size ", slice_size, " total size ",
               tot_length, "\n");
  SLOG_VERBOSE("All contig slice sizes ", msm_slice_size.to_string(), "\n");
}

void Contigs::clear_slices() {
  begin_idx = 0;
  end_idx = contigs.size();
}

std::vector<Contig>::iterator Contigs::begin() { return contigs.begin() + begin_idx; }

std::vector<Contig>::iterator Contigs::end() { return contigs.begin() + end_idx; }

std::vector<Contig>::const_iterator Contigs::begin() const { return contigs.begin() + begin_idx; }

std::vector<Contig>::const_iterator Contigs::end() const { return contigs.begin() + end_idx; }

void Contigs::sort_by_length() {
  sort(contigs.begin(), contigs.end(), [](const Contig &c1, const Contig &c2) { return c1.seq.length() > c2.seq.length(); });
}
