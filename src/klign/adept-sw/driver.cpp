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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>

#include "gpu-utils/gpu_compatibility.hpp"
#include "gpu-utils/gpu_common.hpp"
#include "gpu-utils/gpu_utils.hpp"

#include "driver.hpp"
#include "kernel.hpp"

using timepoint_t = std::chrono::time_point<std::chrono::high_resolution_clock>;
#define clock_now std::chrono::high_resolution_clock::now

struct gpu_alignments {
  short* ref_start_gpu;
  short* ref_end_gpu;
  short* query_start_gpu;
  short* query_end_gpu;
  short* scores_gpu;
  unsigned* offset_ref_gpu;
  unsigned* offset_query_gpu;

  gpu_alignments(int max_alignments) {
    ERROR_CHECK(Malloc(&offset_query_gpu, (max_alignments) * sizeof(int)));
    ERROR_CHECK(Malloc(&offset_ref_gpu, (max_alignments) * sizeof(int)));
    ERROR_CHECK(Malloc(&ref_start_gpu, (max_alignments) * sizeof(short)));
    ERROR_CHECK(Malloc(&ref_end_gpu, (max_alignments) * sizeof(short)));
    ERROR_CHECK(Malloc(&query_end_gpu, (max_alignments) * sizeof(short)));
    ERROR_CHECK(Malloc(&query_start_gpu, (max_alignments) * sizeof(short)));
    ERROR_CHECK(Malloc(&scores_gpu, (max_alignments) * sizeof(short)));
  }

  ~gpu_alignments() {
    ERROR_CHECK(Free(offset_ref_gpu));
    ERROR_CHECK(Free(offset_query_gpu));
    ERROR_CHECK(Free(ref_start_gpu));
    ERROR_CHECK(Free(ref_end_gpu));
    ERROR_CHECK(Free(query_start_gpu));
    ERROR_CHECK(Free(query_end_gpu));
    ERROR_CHECK(Free(scores_gpu));
  }
};

struct gpu_alignments_traceback {
  char* longCIGAR_gpu;
  char* CIGAR_gpu;
  char* H_ptr_gpu;

  gpu_alignments_traceback(int max_alignments, int maxCIGAR, unsigned const maxMatrixSize) {
    // printf("new operator called, creating gpu_alignments_traceback data structure\n");
    ERROR_CHECK(Malloc(&CIGAR_gpu, (max_alignments) * sizeof(char) * maxCIGAR));
    // added a buffer because of cuda-error in larger sequences
    ERROR_CHECK(Malloc(&H_ptr_gpu, 1.25 * sizeof(char) * maxMatrixSize * (max_alignments)));
    ERROR_CHECK(Malloc(&longCIGAR_gpu, sizeof(char) * maxCIGAR * (max_alignments)));
  }

  ~gpu_alignments_traceback() {
    // printf("delete operator called, deleting gpu_alignments_traceback data structure\n");
    ERROR_CHECK(Free(CIGAR_gpu));
    ERROR_CHECK(Free(H_ptr_gpu));
    ERROR_CHECK(Free(longCIGAR_gpu));
  }
};

void asynch_mem_copies_htd(gpu_alignments* gpu_data, unsigned* offsetA_h, unsigned* offsetB_h, char* strA, char* strA_d, char* strB,
                           char* strB_d, unsigned half_length_A, unsigned half_length_B, unsigned totalLengthA,
                           unsigned totalLengthB, int sequences_per_stream, int sequences_stream_leftover, Stream_t* streams_cuda,
                           int max_rlen) {
  ERROR_CHECK(
      MemcpyAsync(gpu_data->offset_ref_gpu, offsetA_h, (sequences_per_stream) * sizeof(int), MemcpyHostToDevice, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(gpu_data->offset_ref_gpu + sequences_per_stream, offsetA_h + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(int), MemcpyHostToDevice, streams_cuda[1]));

  ERROR_CHECK(MemcpyAsync(gpu_data->offset_query_gpu, offsetB_h, (sequences_per_stream) * sizeof(int), MemcpyHostToDevice,
                          streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(gpu_data->offset_query_gpu + sequences_per_stream, offsetB_h + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(int), MemcpyHostToDevice, streams_cuda[1]));

  ERROR_CHECK(MemcpyAsync(strA_d, strA, half_length_A * sizeof(char), MemcpyHostToDevice, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(strA_d + half_length_A, strA + half_length_A, (totalLengthA - half_length_A) * sizeof(char),
                          MemcpyHostToDevice, streams_cuda[1]));

  size_t size_strB = sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE;
  if (size_strB < totalLengthB) {
    std::cerr << "<" << __FILE__ << ":" << __LINE__ << "> ERROR: size_strB " << size_strB << " < " << "totalLengthB "
              << totalLengthB << " max rlen " << max_rlen << "\n";
    std::abort();
  }
  if (size_strB < half_length_B) {
    std::cerr << "<" << __FILE__ << ":" << __LINE__ << "> ERROR: size_strB " << size_strB << " < " << "half_length_B "
              << half_length_B << " max rlen " << max_rlen << "\n";
    std::abort();
  }

  ERROR_CHECK(MemcpyAsync(strB_d, strB, half_length_B * sizeof(char), MemcpyHostToDevice, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(strB_d + half_length_B, strB + half_length_B, (totalLengthB - half_length_B) * sizeof(char),
                          MemcpyHostToDevice, streams_cuda[1]));
}

void asynch_mem_copies_dth_mid(gpu_alignments* gpu_data, short* alAend, short* alBend, int sequences_per_stream,
                               int sequences_stream_leftover, Stream_t* streams_cuda) {
  ERROR_CHECK(
      MemcpyAsync(alAend, gpu_data->ref_end_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alAend + sequences_per_stream, gpu_data->ref_end_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(alBend, gpu_data->query_end_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alBend + sequences_per_stream, gpu_data->query_end_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));
}

void asynch_mem_copies_dth(gpu_alignments* gpu_data, short* alAbeg, short* alBbeg, short* top_scores_cpu, int sequences_per_stream,
                           int sequences_stream_leftover, Stream_t* streams_cuda) {
  ERROR_CHECK(
      MemcpyAsync(alAbeg, gpu_data->ref_start_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alAbeg + sequences_per_stream, gpu_data->ref_start_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(alBbeg, gpu_data->query_start_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alBbeg + sequences_per_stream, gpu_data->query_start_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(top_scores_cpu, gpu_data->scores_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(top_scores_cpu + sequences_per_stream, gpu_data->scores_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));
}

void asynch_mem_copies_htd_t(gpu_alignments* gpu_data, gpu_alignments_traceback* gpu_data_traceback, unsigned* offsetA_h,
                             unsigned* offsetB_h, char* strA, char* strA_d, char* strB, char* strB_d, unsigned half_length_A,
                             unsigned half_length_B, unsigned totalLengthA, unsigned totalLengthB, int sequences_per_stream,
                             int sequences_stream_leftover, Stream_t* streams_cuda, int max_rlen) {
  ERROR_CHECK(
      MemcpyAsync(gpu_data->offset_ref_gpu, offsetA_h, (sequences_per_stream) * sizeof(int), MemcpyHostToDevice, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(gpu_data->offset_ref_gpu + sequences_per_stream, offsetA_h + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(int), MemcpyHostToDevice, streams_cuda[1]));

  ERROR_CHECK(MemcpyAsync(gpu_data->offset_query_gpu, offsetB_h, (sequences_per_stream) * sizeof(int), MemcpyHostToDevice,
                          streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(gpu_data->offset_query_gpu + sequences_per_stream, offsetB_h + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(int), MemcpyHostToDevice, streams_cuda[1]));

  ERROR_CHECK(MemcpyAsync(strA_d, strA, half_length_A * sizeof(char), MemcpyHostToDevice, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(strA_d + half_length_A, strA + half_length_A, (totalLengthA - half_length_A) * sizeof(char),
                          MemcpyHostToDevice, streams_cuda[1]));

  size_t size_strB = sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE;
  if (size_strB < totalLengthB) {
    std::cerr << "<" << __FILE__ << ":" << __LINE__ << "> ERROR: size_strB " << size_strB << " < " << "totalLengthB "
              << totalLengthB << " max rlen " << max_rlen << "\n";
    std::abort();
  }
  if (size_strB < half_length_B) {
    std::cerr << "<" << __FILE__ << ":" << __LINE__ << "> ERROR: size_strB " << size_strB << " < " << "half_length_B "
              << half_length_B << " max rlen " << max_rlen << "\n";
    std::abort();
  }
  ERROR_CHECK(MemcpyAsync(strB_d, strB, half_length_B * sizeof(char), MemcpyHostToDevice, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(strB_d + half_length_B, strB + half_length_B, (totalLengthB - half_length_B) * sizeof(char),
                          MemcpyHostToDevice, streams_cuda[1]));
}

void asynch_mem_copies_dth_t(gpu_alignments* gpu_data, gpu_alignments_traceback* gpu_data_traceback, short* alAbeg, short* alBbeg,
                             short* alAend, short* alBend, short* top_scores_cpu, char* cigar_cpu, int maxCIGAR,
                             int sequences_per_stream, int sequences_stream_leftover, Stream_t* streams_cuda) {
  ERROR_CHECK(
      MemcpyAsync(alAbeg, gpu_data->ref_start_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alAbeg + sequences_per_stream, gpu_data->ref_start_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(alBbeg, gpu_data->query_start_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alBbeg + sequences_per_stream, gpu_data->query_start_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(alAend, gpu_data->ref_end_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alAend + sequences_per_stream, gpu_data->ref_end_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(alBend, gpu_data->query_end_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(alBend + sequences_per_stream, gpu_data->query_end_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(
      MemcpyAsync(top_scores_cpu, gpu_data->scores_gpu, sequences_per_stream * sizeof(short), MemcpyDeviceToHost, streams_cuda[0]));
  ERROR_CHECK(MemcpyAsync(top_scores_cpu + sequences_per_stream, gpu_data->scores_gpu + sequences_per_stream,
                          (sequences_per_stream + sequences_stream_leftover) * sizeof(short), MemcpyDeviceToHost, streams_cuda[1]));

  ERROR_CHECK(MemcpyAsync(cigar_cpu, gpu_data_traceback->CIGAR_gpu, sequences_per_stream * sizeof(char) * maxCIGAR,
                          MemcpyDeviceToHost, streams_cuda[0]));

  ERROR_CHECK(MemcpyAsync(
      cigar_cpu + sequences_per_stream * maxCIGAR, gpu_data_traceback->CIGAR_gpu + sequences_per_stream * maxCIGAR,
      (sequences_per_stream + sequences_stream_leftover) * sizeof(char) * maxCIGAR, MemcpyDeviceToHost, streams_cuda[1]));
}

int get_new_min_length(short* alAend, short* alBend, int blocksLaunched) {
  int newMin = 1000;
  int maxA = 0;
  int maxB = 0;
  for (int i = 0; i < blocksLaunched; i++) {
    if (alBend[i] > maxB) maxB = alBend[i];
    if (alAend[i] > maxA) maxA = alAend[i];
  }
  newMin = (maxB > maxA) ? maxA : maxB;
  return newMin;
}

struct adept_sw::DriverState {
  int rank_me;
  Stream_t streams_cuda[NSTREAMS];
  unsigned* offsetA_h;
  unsigned* offsetB_h;
  char *strA_d, *strB_d;
  char* strA;
  char* strB;

  Event_t event_fwd_0, event_fwd_1, event_rev_0, event_rev_1;
  short matchScore, misMatchScore, startGap, extendGap;
  gpu_alignments* gpu_data;
  gpu_alignments_traceback* gpu_data_traceback;
  unsigned half_length_A = 0;
  unsigned half_length_B = 0;
  int max_rlen = 0;
};

adept_sw::GPUDriver::GPUDriver(int upcxx_rank_me, int upcxx_rank_n, short match_score, short mismatch_score,
                               short gap_opening_score, short gap_extending_score, int max_rlen, int max_clen, bool compute_cigar,
                               double& init_time) {
  timepoint_t t = clock_now();

  // calculate if the batch will fit into global memory
  unsigned maxCIGAR = (max_clen > max_rlen) ? 3 * max_clen : 3 * max_rlen;

  int maxMatrixSize = max_rlen * max_clen;
  size_t gpu_mem_avail = gpu_utils::get_gpu_avail_mem();
  size_t tot_mem_req_per_aln =
      max_rlen + (3 * max_rlen) + 2 * sizeof(int) + 6 * sizeof(short) + (1.25 * max_rlen * max_clen) + 2 * (maxCIGAR);
  float factor = 0.75 * 1 / NSTREAMS;
  unsigned max_alns_gpu = ceil(((double)gpu_mem_avail * factor) / tot_mem_req_per_aln);
  unsigned max_alns_sugg = KLIGN_GPU_BLOCK_SIZE;
  // printf("LLCOMMENT: max_rlen = %d, maxCIGAR = %d, batch size = %d, max_alns_gpu = %d\n",max_rlen, maxCIGAR, max_alns_sugg,
  // max_alns_gpu);
  max_alns_gpu = max_alns_gpu > max_alns_sugg ? max_alns_sugg : max_alns_gpu;
  // int its = (max_per_device>max_alns_gpu)?(ceil((double)max_per_device/max_alns_gpu)):1;
  // printf("LLCOMMENT:new max_alns_gpu = %d\n", max_alns_gpu);
  std::chrono::duration<double> elapsed;
  driver_state = new DriverState();
  driver_state->rank_me = upcxx_rank_me;
  // elapsed =  clock_now() - t; os << " new=" << elapsed.count();
  driver_state->matchScore = match_score;
  driver_state->misMatchScore = mismatch_score;
  driver_state->startGap = gap_opening_score;
  driver_state->extendGap = gap_extending_score;
  driver_state->max_rlen = max_rlen;

  // elapsed =  clock_now() - t; os << " get_num_gpus=" << elapsed.count();

  ERROR_CHECK(MallocHost((void**)&(alignments.ref_begin), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&(alignments.ref_end), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&(alignments.query_begin), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&(alignments.query_end), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&(alignments.top_scores), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&(alignments.cigar),
                         sizeof(char) * 3 * max_rlen * KLIGN_GPU_BLOCK_SIZE));  // FIXME LL: estimate here for maxCIGAR
  // elapsed =  clock_now() - t; os << " mallocHost=" << elapsed.count();

  gpu_utils::set_gpu_device(driver_state->rank_me);
  // elapsed =  clock_now() - t; os << " set_device=" << elapsed.count();

  for (int stm = 0; stm < NSTREAMS; stm++) {
    ERROR_CHECK(StreamCreate(&driver_state->streams_cuda[stm]));
  }
  // elapsed =  clock_now() - t; os << " streamcreate=" << elapsed.count();

  ERROR_CHECK(MallocHost((void**)&driver_state->offsetA_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&driver_state->offsetB_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE));
  // elapsed =  clock_now() - t; os << " mallocHost2=" << elapsed.count();

  // FIXME: hack for max contig and read size -> multiplying max_rlen by 2 for contigs
  ERROR_CHECK(Malloc(&driver_state->strA_d, 2 * max_rlen * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));
  ERROR_CHECK(Malloc(&driver_state->strB_d, max_rlen * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));
  // elapsed =  clock_now() - t; os << " mallocs=" << elapsed.count();

  ERROR_CHECK(MallocHost((void**)&driver_state->strA, 2 * sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE));
  ERROR_CHECK(MallocHost((void**)&driver_state->strB, sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE));
  // elapsed =  clock_now() - t; os << " mallocHost3=" << elapsed.count();

  driver_state->gpu_data = new gpu_alignments(KLIGN_GPU_BLOCK_SIZE);  // gpu mallocs
  // elapsed =  clock_now() - t; os << " final=" << elapsed.count();

  if (compute_cigar)
    driver_state->gpu_data_traceback = new gpu_alignments_traceback(KLIGN_GPU_BLOCK_SIZE, maxCIGAR, maxMatrixSize);  // gpu mallocs

  ERROR_CHECK(EventCreateWithFlags(&driver_state->event_fwd_0, EventDisableTiming | EventBlockingSync));
  ERROR_CHECK(EventCreateWithFlags(&driver_state->event_fwd_1, EventDisableTiming | EventBlockingSync));

  ERROR_CHECK(EventCreateWithFlags(&driver_state->event_rev_0, EventDisableTiming | EventBlockingSync));
  ERROR_CHECK(EventCreateWithFlags(&driver_state->event_rev_1, EventDisableTiming | EventBlockingSync));

  elapsed = clock_now() - t;
  init_time = elapsed.count();
}

adept_sw::GPUDriver::~GPUDriver() {
  // won't have been allocated if there was no GPU present
  if (!alignments.ref_begin) return;
  gpu_utils::set_gpu_device(driver_state->rank_me);
  ERROR_CHECK(FreeHost(alignments.ref_begin));
  ERROR_CHECK(FreeHost(alignments.ref_end));
  ERROR_CHECK(FreeHost(alignments.query_begin));
  ERROR_CHECK(FreeHost(alignments.query_end));
  ERROR_CHECK(FreeHost(alignments.top_scores));
  ERROR_CHECK(FreeHost(alignments.cigar));

  ERROR_CHECK(Free(driver_state->strA_d));
  ERROR_CHECK(Free(driver_state->strB_d));
  ERROR_CHECK(FreeHost(driver_state->offsetA_h));
  ERROR_CHECK(FreeHost(driver_state->offsetB_h));
  ERROR_CHECK(FreeHost(driver_state->strA));
  ERROR_CHECK(FreeHost(driver_state->strB));

  ERROR_CHECK(EventDestroy(driver_state->event_rev_0));
  ERROR_CHECK(EventDestroy(driver_state->event_rev_1));
  ERROR_CHECK(EventDestroy(driver_state->event_fwd_0));
  ERROR_CHECK(EventDestroy(driver_state->event_fwd_1));
  for (int i = 0; i < NSTREAMS; i++) ERROR_CHECK(StreamDestroy(driver_state->streams_cuda[i]));
  delete driver_state->gpu_data;
  delete driver_state->gpu_data_traceback;
  delete driver_state;
}

void adept_sw::GPUDriver::kernel_block_fwd() {
  ERROR_CHECK(EventSynchronize(driver_state->event_fwd_0));
  ERROR_CHECK(EventSynchronize(driver_state->event_fwd_1));
}

void adept_sw::GPUDriver::kernel_block_rev() {
  ERROR_CHECK(EventSynchronize(driver_state->event_rev_0));
  ERROR_CHECK(EventSynchronize(driver_state->event_rev_1));
}

void adept_sw::GPUDriver::run_kernel_forwards(std::vector<std::string>& reads, std::vector<std::string>& contigs,
                                              unsigned maxReadSize, unsigned maxContigSize, double& launch_time, double& mem_time) {
  gpu_utils::set_gpu_device(driver_state->rank_me);
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length
  std::chrono::duration<double> elapsed;

  short* alAend = alignments.ref_end;
  short* alBend = alignments.query_end;
  // memory on CPU for copying the results

  int blocksLaunched = totalAlignments;
  std::vector<std::string>::const_iterator beginAVec;
  std::vector<std::string>::const_iterator endAVec;
  std::vector<std::string>::const_iterator beginBVec;
  std::vector<std::string>::const_iterator endBVec;
  beginAVec = contigs.begin();
  endAVec = contigs.begin() + totalAlignments;
  beginBVec = reads.begin();
  endBVec = reads.begin() + totalAlignments;

  std::vector<std::string> sequencesA(beginAVec, endAVec);
  std::vector<std::string> sequencesB(beginBVec, endBVec);
  unsigned running_sum = 0;
  int sequences_per_stream = (blocksLaunched) / NSTREAMS;
  int sequences_stream_leftover = (blocksLaunched) % NSTREAMS;
  driver_state->half_length_A = 0;
  driver_state->half_length_B = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    running_sum += sequencesA[i].size();
    driver_state->offsetA_h[i] = running_sum;  // sequencesA[i].size();
    if (i == sequences_per_stream - 1) {
      driver_state->half_length_A = running_sum;
      running_sum = 0;
    }
  }
  unsigned totalLengthA = driver_state->half_length_A + driver_state->offsetA_h[sequencesA.size() - 1];

  running_sum = 0;
  for (int i = 0; i < (int)sequencesB.size(); i++) {
    running_sum += sequencesB[i].size();
    driver_state->offsetB_h[i] = running_sum;  // sequencesB[i].size();
    if (i == sequences_per_stream - 1) {
      driver_state->half_length_B = running_sum;
      running_sum = 0;
    }
  }
  unsigned totalLengthB = driver_state->half_length_B + driver_state->offsetB_h[sequencesB.size() - 1];

  unsigned offsetSumA = 0;
  unsigned offsetSumB = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    char* seqptrA = driver_state->strA + offsetSumA;
    memcpy(seqptrA, sequencesA[i].c_str(), sequencesA[i].size());
    char* seqptrB = driver_state->strB + offsetSumB;
    memcpy(seqptrB, sequencesB[i].c_str(), sequencesB[i].size());
    offsetSumA += sequencesA[i].size();
    offsetSumB += sequencesB[i].size();
  }

  timepoint_t t = clock_now();
  asynch_mem_copies_htd(driver_state->gpu_data, driver_state->offsetA_h, driver_state->offsetB_h, driver_state->strA,
                        driver_state->strA_d, driver_state->strB, driver_state->strB_d, driver_state->half_length_A,
                        driver_state->half_length_B, totalLengthA, totalLengthB, sequences_per_stream, sequences_stream_leftover,
                        driver_state->streams_cuda, driver_state->max_rlen);
  elapsed = clock_now() - t;
  mem_time += elapsed.count();

  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 3 * (minSize + 1) * sizeof(short);
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad;

  if (ShmemBytes > 48000)
    ERROR_CHECK(FuncSetAttribute((const void*)gpu_bsw::dna_kernel, FuncAttributeMaxDynamicSharedMemorySize, ShmemBytes));

  t = clock_now();
  LaunchKernelGGL(gpu_bsw::dna_kernel, sequences_per_stream, minSize, ShmemBytes, driver_state->streams_cuda[0],
                  driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu,
                  driver_state->gpu_data->offset_query_gpu, driver_state->gpu_data->ref_start_gpu,
                  driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
                  driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu, driver_state->matchScore,
                  driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap, false);

  LaunchKernelGGL(
      gpu_bsw::dna_kernel, sequences_per_stream + sequences_stream_leftover, minSize, ShmemBytes, driver_state->streams_cuda[1],
      driver_state->strA_d + driver_state->half_length_A, driver_state->strB_d + driver_state->half_length_B,
      driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
      driver_state->gpu_data->offset_query_gpu + sequences_per_stream, driver_state->gpu_data->ref_start_gpu + sequences_per_stream,
      driver_state->gpu_data->ref_end_gpu + sequences_per_stream, driver_state->gpu_data->query_start_gpu + sequences_per_stream,
      driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
      driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap, false);

  // does not work without the below stream synchs on AMDGPUs
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[0]));
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[1]));
  elapsed = clock_now() - t;
  launch_time += elapsed.count();

  // copyin back end index so that we can find new min
  t = clock_now();
  asynch_mem_copies_dth_mid(driver_state->gpu_data, alAend, alBend, sequences_per_stream, sequences_stream_leftover,
                            driver_state->streams_cuda);
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[0]));
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[1]));

  ERROR_CHECK(EventRecord(driver_state->event_fwd_0, driver_state->streams_cuda[0]));
  ERROR_CHECK(EventRecord(driver_state->event_fwd_1, driver_state->streams_cuda[1]));
  elapsed = clock_now() - t;
  mem_time += elapsed.count();
}

void adept_sw::GPUDriver::run_kernel_backwards(std::vector<std::string>& reads, std::vector<std::string>& contigs,
                                               unsigned maxReadSize, unsigned maxContigSize, double& launch_time,
                                               double& mem_time) {
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length

  short* alAbeg = alignments.ref_begin;
  short* alBbeg = alignments.query_begin;
  short* alAend = alignments.ref_end;
  short* alBend = alignments.query_end;
  ;  // memory on CPU for copying the results
  short* top_scores_cpu = alignments.top_scores;
  int blocksLaunched = totalAlignments;
  int sequences_per_stream = (blocksLaunched) / NSTREAMS;
  int sequences_stream_leftover = (blocksLaunched) % NSTREAMS;
  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 3 * (minSize + 1) * sizeof(short);
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad;

  int newMin = get_new_min_length(alAend, alBend, blocksLaunched);  // find the new largest of smaller lengths

  if (ShmemBytes > 48000)
    ERROR_CHECK(FuncSetAttribute((const void*)gpu_bsw::dna_kernel, FuncAttributeMaxDynamicSharedMemorySize, ShmemBytes));

  std::chrono::duration<double> elapsed;
  timepoint_t t = clock_now();
  LaunchKernelGGL(gpu_bsw::dna_kernel, sequences_per_stream, newMin, ShmemBytes, driver_state->streams_cuda[0],
                  driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu,
                  driver_state->gpu_data->offset_query_gpu, driver_state->gpu_data->ref_start_gpu,
                  driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
                  driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu, driver_state->matchScore,
                  driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap, true);

  LaunchKernelGGL(
      gpu_bsw::dna_kernel, sequences_per_stream + sequences_stream_leftover, newMin, ShmemBytes, driver_state->streams_cuda[1],
      driver_state->strA_d + driver_state->half_length_A, driver_state->strB_d + driver_state->half_length_B,
      driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
      driver_state->gpu_data->offset_query_gpu + sequences_per_stream, driver_state->gpu_data->ref_start_gpu + sequences_per_stream,
      driver_state->gpu_data->ref_end_gpu + sequences_per_stream, driver_state->gpu_data->query_start_gpu + sequences_per_stream,
      driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
      driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap, true);

  // does not work without the below stream synchs on AMDGPUs
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[0]));
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[1]));
  elapsed = clock_now() - t;
  launch_time += elapsed.count();

  // copyin back end index so that we can find new min
  t = clock_now();
  asynch_mem_copies_dth(driver_state->gpu_data, alAbeg, alBbeg, top_scores_cpu, sequences_per_stream, sequences_stream_leftover,
                        driver_state->streams_cuda);
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[0]));
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[1]));

  ERROR_CHECK(EventRecord(driver_state->event_rev_0, driver_state->streams_cuda[0]));
  ERROR_CHECK(EventRecord(driver_state->event_rev_1, driver_state->streams_cuda[1]));
  elapsed = clock_now() - t;
  mem_time += elapsed.count();
}

void adept_sw::GPUDriver::run_kernel_traceback(std::vector<std::string>& reads, std::vector<std::string>& contigs,
                                               unsigned maxReadSize, unsigned maxContigSize, double& launch_time,
                                               double& mem_time) {
  gpu_utils::set_gpu_device(driver_state->rank_me);

  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length
  unsigned maxCIGAR = (maxContigSize > maxReadSize) ? 3 * maxContigSize : 3 * maxReadSize;  // check if this is now necessary
  unsigned const maxMatrixSize = (maxContigSize + 1) * (maxReadSize + 1);
  // memory on CPU for copying the results
  short* alAbeg = alignments.ref_begin;
  short* alBbeg = alignments.query_begin;
  short* alAend = alignments.ref_end;
  short* alBend = alignments.query_end;
  short* top_scores_cpu = alignments.top_scores;
  char* cigar_cpu = alignments.cigar;

  int blocksLaunched = totalAlignments;
  std::vector<std::string>::const_iterator beginAVec;
  std::vector<std::string>::const_iterator endAVec;
  std::vector<std::string>::const_iterator beginBVec;
  std::vector<std::string>::const_iterator endBVec;
  beginAVec = contigs.begin();
  endAVec = contigs.begin() + totalAlignments;
  beginBVec = reads.begin();
  endBVec = reads.begin() + totalAlignments;

  std::vector<std::string> sequencesA(beginAVec, endAVec);
  std::vector<std::string> sequencesB(beginBVec, endBVec);
  unsigned running_sum = 0;
  int sequences_per_stream = (blocksLaunched) / NSTREAMS;
  int sequences_stream_leftover = (blocksLaunched) % NSTREAMS;
  driver_state->half_length_A = 0;
  driver_state->half_length_B = 0;

  long unsigned largestA = 0;
  long unsigned largestB = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    running_sum += sequencesA[i].size();
    driver_state->offsetA_h[i] = running_sum;
    if (i == sequences_per_stream - 1) {
      driver_state->half_length_A = running_sum;
      running_sum = 0;
    }

    if (sequencesA[i].size() > largestA) {
      largestA = sequencesA[i].size();
    }
  }
  unsigned totalLengthA = driver_state->half_length_A + driver_state->offsetA_h[sequencesA.size() - 1];

  running_sum = 0;
  for (int i = 0; i < (int)sequencesB.size(); i++) {
    running_sum += sequencesB[i].size();
    driver_state->offsetB_h[i] = running_sum;
    if (i == sequences_per_stream - 1) {
      driver_state->half_length_B = running_sum;
      running_sum = 0;
    }

    if (sequencesB[i].size() > largestB) {
      largestB = sequencesB[i].size();
    }
  }

  unsigned totalLengthB = driver_state->half_length_B + driver_state->offsetB_h[sequencesB.size() - 1];

  unsigned offsetSumA = 0;
  unsigned offsetSumB = 0;

  for (int i = 0; i < (int)sequencesA.size(); i++) {
    char* seqptrA = driver_state->strA + offsetSumA;
    memcpy(seqptrA, sequencesA[i].c_str(), sequencesA[i].size());
    char* seqptrB = driver_state->strB + offsetSumB;
    memcpy(seqptrB, sequencesB[i].c_str(), sequencesB[i].size());
    offsetSumA += sequencesA[i].size();
    offsetSumB += sequencesB[i].size();
  }

  std::chrono::duration<double> elapsed;
  timepoint_t t = clock_now();
  asynch_mem_copies_htd_t(driver_state->gpu_data, driver_state->gpu_data_traceback, driver_state->offsetA_h,
                          driver_state->offsetB_h, driver_state->strA, driver_state->strA_d, driver_state->strB,
                          driver_state->strB_d, driver_state->half_length_A, driver_state->half_length_B, totalLengthA,
                          totalLengthB, sequences_per_stream, sequences_stream_leftover, driver_state->streams_cuda,
                          driver_state->max_rlen);
  elapsed = clock_now() - t;
  mem_time += elapsed.count();
  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned maxSize = (maxReadSize > maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 6 * (minSize + 1) * sizeof(short) + 6 * minSize + (minSize + 1) + maxSize;
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad + sizeof(int) * (maxContigSize + maxReadSize + 2);

  if (ShmemBytes > 48000)
    ERROR_CHECK(
        FuncSetAttribute((const void*)gpu_bsw::sequence_dna_kernel_traceback, FuncAttributeMaxDynamicSharedMemorySize, ShmemBytes));

  t = clock_now();
  LaunchKernelGGL(
      gpu_bsw::sequence_dna_kernel_traceback, sequences_per_stream, minSize, ShmemBytes, driver_state->streams_cuda[0],
      driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu, driver_state->gpu_data->offset_query_gpu,
      driver_state->gpu_data->ref_start_gpu, driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
      driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu, driver_state->gpu_data_traceback->longCIGAR_gpu,
      driver_state->gpu_data_traceback->CIGAR_gpu, driver_state->gpu_data_traceback->H_ptr_gpu, maxCIGAR, maxMatrixSize,
      driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  LaunchKernelGGL(
      gpu_bsw::sequence_dna_kernel_traceback, sequences_per_stream + sequences_stream_leftover, minSize, ShmemBytes,
      driver_state->streams_cuda[1], driver_state->strA_d + driver_state->half_length_A,
      driver_state->strB_d + driver_state->half_length_B, driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
      driver_state->gpu_data->offset_query_gpu + sequences_per_stream, driver_state->gpu_data->ref_start_gpu + sequences_per_stream,
      driver_state->gpu_data->ref_end_gpu + sequences_per_stream, driver_state->gpu_data->query_start_gpu + sequences_per_stream,
      driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
      driver_state->gpu_data_traceback->longCIGAR_gpu + sequences_per_stream * maxCIGAR,
      driver_state->gpu_data_traceback->CIGAR_gpu + sequences_per_stream * maxCIGAR,
      driver_state->gpu_data_traceback->H_ptr_gpu + sequences_per_stream * maxMatrixSize, maxCIGAR, maxMatrixSize,
      driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  // does not work without the below stream synchs on AMDGPUs
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[0]));
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[1]));
  elapsed = clock_now() - t;
  launch_time += elapsed.count();

  t = clock_now();
  asynch_mem_copies_dth_t(driver_state->gpu_data, driver_state->gpu_data_traceback, alAbeg, alBbeg, alAend, alBend, top_scores_cpu,
                          cigar_cpu, maxCIGAR, sequences_per_stream, sequences_stream_leftover, driver_state->streams_cuda);

  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[0]));
  ERROR_CHECK(StreamSynchronize(driver_state->streams_cuda[1]));

  ERROR_CHECK(EventRecord(driver_state->event_fwd_0, driver_state->streams_cuda[0]));
  ERROR_CHECK(EventRecord(driver_state->event_fwd_1, driver_state->streams_cuda[1]));
  elapsed = clock_now() - t;
  mem_time += elapsed.count();
}
