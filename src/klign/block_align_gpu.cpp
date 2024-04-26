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

#define NO_KLIGN_CPU_WORK_STEAL

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

static adept_sw::GPUDriver *gpu_driver;

static upcxx::future<> gpu_align_block(shared_ptr<AlignBlockData> aln_block_data, Alns *alns, bool report_cigar,
                                       KlignTimers &klign_timers) {
  assert(upcxx::master_persona().active_with_caller());

  future<> fut = upcxx_utils::execute_serially_in_thread_pool([aln_block_data, report_cigar, &klign_timers]() {
    DBG("Starting _gpu_align_block_kernel of ", aln_block_data->kernel_alns.size(), "\n");

    unsigned maxContigSize = aln_block_data->max_clen;
    unsigned maxReadSize = aln_block_data->max_rlen;
    auto &aln_kernel_timer = klign_timers.aln_kernel;
    auto &block_timer = klign_timers.aln_kernel_block;
    double launch_time = 0, mem_time = 0;
    unsigned maxCIGAR = (maxContigSize > maxReadSize) ?
                            3 * maxContigSize :
                            3 * maxReadSize;  // 3* size to eliminate overflow FIXME: Truncate CIGARs that are over maxCIGAR
    // printf("GPU align block: maxContigSize passed in = %d, maxReadSize passed in = %d\n", maxContigSize, maxReadSize);
    if (report_cigar) {
      aln_kernel_timer.start();

      gpu_driver->run_kernel_traceback(aln_block_data->read_seqs, aln_block_data->ctg_seqs, aln_block_data->max_rlen,
                                       aln_block_data->max_clen, launch_time, mem_time);
      block_timer.start();
      gpu_driver->kernel_block_fwd();
      block_timer.stop();

      aln_kernel_timer.stop();
      klign_timers.aln_kernel_launch += launch_time;
      klign_timers.aln_kernel_mem += mem_time;
    } else {
      aln_kernel_timer.start();

      // align query_seqs, ref_seqs, max_query_size, max_ref_size
      gpu_driver->run_kernel_forwards(aln_block_data->read_seqs, aln_block_data->ctg_seqs, aln_block_data->max_rlen,
                                      aln_block_data->max_clen, launch_time, mem_time);
      block_timer.start();
      gpu_driver->kernel_block_fwd();
      block_timer.stop();
      gpu_driver->run_kernel_backwards(aln_block_data->read_seqs, aln_block_data->ctg_seqs, aln_block_data->max_rlen,
                                       aln_block_data->max_clen, launch_time, mem_time);
      block_timer.start();
      gpu_driver->kernel_block_rev();
      block_timer.stop();

      aln_kernel_timer.stop();
      klign_timers.aln_kernel_launch += launch_time;
      klign_timers.aln_kernel_mem += mem_time;
    }

    auto aln_results = gpu_driver->get_aln_results();
    aln_block_data->alns->reserve(aln_block_data->alns->size() + aln_block_data->kernel_alns.size());

    for (int i = 0; i < aln_block_data->kernel_alns.size(); i++) {
      Aln &aln = aln_block_data->kernel_alns[i];
      // FIXME: need to get the second best score
      // FIXME: need to get the number of mismatches
      aln.set(aln_results.ref_begin[i], aln_results.ref_end[i], aln_results.query_begin[i], aln_results.query_end[i],
              aln_results.top_scores[i], 0, 0, aln_block_data->read_group_id);
      if (report_cigar) {
        // std::string cig = "GPU_CIGAR: "; //use this to calculate which percentage of traceback is being done on GPU
        std::string cig = "";
        int k = i * maxCIGAR;
        while (aln_results.cigar[k]) {
          cig += aln_results.cigar[k];
          k++;
        }
        aln.set_sam_string(aln_block_data->read_seqs[i], cig);
        // cig = "GPU_CIGAR: "; //use this to calculate which percentage of traceback is being done on GPU
        cig = "";
      }
      aln_block_data->alns->add_aln(aln);
    }
  });
  fut = fut.then([alns = alns, aln_block_data]() {
    assert(upcxx::master_persona().active_with_caller());
    LOG("appending and returning ", aln_block_data->alns->size(), "\n");
    alns->append(*(aln_block_data->alns));
  });
  return fut;
}

void init_aligner(int match_score, int mismatch_penalty, int gap_opening_penalty, int gap_extending_penalty, int ambiguity_penalty,
                  int rlen_limit, bool compute_cigar) {
  if (!gpu_utils::gpus_present()) {
    // CPU only
    SWARN("No GPU will be used for alignments");
  } else {
    double init_time;
    gpu_driver =
        new adept_sw::GPUDriver(local_team().rank_me(), local_team().rank_n(), (short)match_score, (short)-mismatch_penalty,
                                (short)-gap_opening_penalty, (short)-gap_extending_penalty, rlen_limit, compute_cigar, init_time);
    SLOG_VERBOSE("Initialized GPU adept_sw driver in ", init_time, " s\n");
  }
}

void cleanup_aligner() {
  if (gpu_utils::gpus_present()) delete gpu_driver;
}

void kernel_align_block(CPUAligner &cpu_aligner, vector<Aln> &kernel_alns, vector<string> &ctg_seqs, vector<string> &read_seqs,
                        Alns *alns, future<> &active_kernel_fut, int read_group_id, int max_clen, int max_rlen,
                        KlignTimers &klign_timers) {
  assert(!upcxx::in_progress());
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
        make_shared<AlignBlockData>(kernel_alns, ctg_seqs, read_seqs, max_clen, max_rlen, read_group_id);
    assert(kernel_alns.empty());
    if (gpu_utils::gpus_present()) {
      active_kernel_fut = gpu_align_block(aln_block_data, alns, cpu_aligner.ssw_filter.report_cigar, klign_timers);
    } else {
      active_kernel_fut = cpu_aligner.ssw_align_block(aln_block_data, alns, klign_timers.aln_kernel);
    }
    progress();
  }
}
