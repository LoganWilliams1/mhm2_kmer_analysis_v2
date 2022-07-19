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

#include "klign.hpp"
#include "kmer.hpp"
#include "aligner_cpu.hpp"

#include "gpu-utils/gpu_utils.hpp"
#include "adept-sw/driver.hpp"

#ifdef __PPC64__  // FIXME remove after solving Issues #60 #35 #49
#define NO_KLIGN_CPU_WORK_STEAL
#endif

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

static adept_sw::GPUDriver *gpu_driver;

static upcxx::future<> gpu_align_block(shared_ptr<AlignBlockData> aln_block_data, Alns *alns, bool report_cigar,
                                       IntermittentTimer &aln_kernel_timer) {
  future<> fut = upcxx_utils::execute_in_thread_pool([aln_block_data, report_cigar, &aln_kernel_timer] {
    DBG_VERBOSE("Starting _gpu_align_block_kernel of ", aln_block_data->kernel_alns.size(), "\n");
    aln_kernel_timer.start();

    // align query_seqs, ref_seqs, max_query_size, max_ref_size
    gpu_driver->run_kernel_forwards(aln_block_data->read_seqs, aln_block_data->ctg_seqs, aln_block_data->max_rlen,
                                    aln_block_data->max_clen);
    gpu_driver->kernel_block();
    gpu_driver->run_kernel_backwards(aln_block_data->read_seqs, aln_block_data->ctg_seqs, aln_block_data->max_rlen,
                                     aln_block_data->max_clen);
    gpu_driver->kernel_block();
    aln_kernel_timer.stop();

    auto aln_results = gpu_driver->get_aln_results();

    for (int i = 0; i < aln_block_data->kernel_alns.size(); i++) {
      Aln &aln = aln_block_data->kernel_alns[i];
      // FIXME: need to get the second best score
      // FIXME: need to get the number of mismatches
      aln.set(aln_results.ref_begin[i], aln_results.ref_end[i], aln_results.query_begin[i], aln_results.query_end[i],
              aln_results.top_scores[i], 0, 0, aln_block_data->read_group_id);
      if (report_cigar) {
        SWARN("Trying to produce SAM outputs with GPU alignments, which is not supported");
        aln.set_sam_string("*", "*");  // FIXME until there is a valid:ssw_aln.cigar_string);
      }
      aln_block_data->alns->add_aln(aln);
    }
  });
  fut = fut.then([alns = alns, aln_block_data]() {
    DBG_VERBOSE("appending and returning ", aln_block_data->alns->size(), "\n");
    alns->append(*(aln_block_data->alns));
  });

  return fut;
}

void init_aligner(AlnScoring &aln_scoring, int rlen_limit) {
  if (!gpu_utils::gpus_present()) {
    // CPU only
    SWARN("No GPU will be used for alignments");
  } else {
    double init_time;
    gpu_driver = new adept_sw::GPUDriver(local_team().rank_me(), local_team().rank_n(), (short)aln_scoring.match,
                                         (short)-aln_scoring.mismatch, (short)-aln_scoring.gap_opening,
                                         (short)-aln_scoring.gap_extending, rlen_limit, init_time);
    SLOG_VERBOSE("Initialized GPU adept_sw driver in ", init_time, " s\n");
  }
}

void cleanup_aligner() {
  if (gpu_utils::gpus_present()) delete gpu_driver;
}

void kernel_align_block(CPUAligner &cpu_aligner, vector<Aln> &kernel_alns, vector<string> &ctg_seqs, vector<string> &read_seqs,
                        Alns *alns, future<> &active_kernel_fut, int read_group_id, int max_clen, int max_rlen,
                        IntermittentTimer &aln_kernel_timer) {
  BaseTimer steal_t("CPU work steal");
  steal_t.start();
  auto num = kernel_alns.size();
  // steal work from this kernel block if the previous kernel is still active
  // if true, this balances the block size that will be sent to the kernel
  while ((!gpu_utils::gpus_present() || !active_kernel_fut.ready()) && !kernel_alns.empty()) {
    assert(!ctg_seqs.empty());
    assert(!read_seqs.empty());
#ifndef NO_KLIGN_CPU_WORK_STEAL
    // steal one from the block
    cpu_aligner.ssw_align_read(alns, kernel_alns.back(), ctg_seqs.back(), read_seqs.back(), read_group_id);
    kernel_alns.pop_back();
    ctg_seqs.pop_back();
    read_seqs.pop_back();
#endif
    progress();
  }
  steal_t.stop();
  auto steal_secs = steal_t.get_elapsed();
  if (num != kernel_alns.size()) {
    auto num_stole = num - kernel_alns.size();
    LOG("Stole from kernel block ", num_stole, " alignments in ", steal_secs, "s (",
        (steal_secs > 0 ? num_stole / steal_secs : 0.0), " aln/s), while waiting for previous block to complete",
        (kernel_alns.empty() ? " - THE ENTIRE BLOCK" : ""), "\n");
  } else if (steal_secs > 0.01) {
    LOG("Waited ", steal_secs, "s for previous block to complete\n");
  }
  if (!kernel_alns.empty()) {
    assert(active_kernel_fut.ready() && "active_kernel_fut should already be ready");
    active_kernel_fut.wait();  // should be ready already
    shared_ptr<AlignBlockData> aln_block_data =
        make_shared<AlignBlockData>(kernel_alns, ctg_seqs, read_seqs, max_clen, max_rlen, read_group_id, cpu_aligner.aln_scoring);
    assert(kernel_alns.empty());
    // for now, the GPU alignment doesn't support cigars
    if (!cpu_aligner.ssw_filter.report_cigar && gpu_utils::gpus_present()) {
      SLOG_VERBOSE("GPU align block\n");
      active_kernel_fut = gpu_align_block(aln_block_data, alns, cpu_aligner.ssw_filter.report_cigar, aln_kernel_timer);
    } else if (!gpu_utils::gpus_present()) {
      active_kernel_fut = cpu_aligner.ssw_align_block(aln_block_data, alns, aln_kernel_timer);
    } else {
#ifdef __PPC64__
      SWARN("FIXME Issue #49,#60 no cigars for gpu alignments\n");
      active_kernel_fut = gpu_align_block(aln_block_data, alns, cpu_aligner.ssw_filter.report_cigar, aln_kernel_timer);
#else
      active_kernel_fut = cpu_aligner.ssw_align_block(aln_block_data, alns, aln_kernel_timer);
#endif
    }
  }
}
