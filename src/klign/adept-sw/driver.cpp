
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

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>

#include <cuda_runtime_api.h>
#include <cuda.h>

#include "driver.hpp"
#include "kernel.hpp"
#include "gpu-utils/gpu_common.hpp"
#include "gpu-utils/gpu_utils.hpp"

struct gpu_alignments {
  short* ref_start_gpu;
  short* ref_end_gpu;
  short* query_start_gpu;
  short* query_end_gpu;
  short* scores_gpu;
  unsigned* offset_ref_gpu;
  unsigned* offset_query_gpu;

  gpu_alignments(int max_alignments) {
    cudaErrchk(cudaMalloc(&offset_query_gpu, (max_alignments) * sizeof(int)));
    cudaErrchk(cudaMalloc(&offset_ref_gpu, (max_alignments) * sizeof(int)));
    cudaErrchk(cudaMalloc(&ref_start_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&ref_end_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&query_end_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&query_start_gpu, (max_alignments) * sizeof(short)));
    cudaErrchk(cudaMalloc(&scores_gpu, (max_alignments) * sizeof(short)));
  }

  ~gpu_alignments() {
    cudaErrchk(cudaFree(offset_ref_gpu));
    cudaErrchk(cudaFree(offset_query_gpu));
    cudaErrchk(cudaFree(ref_start_gpu));
    cudaErrchk(cudaFree(ref_end_gpu));
    cudaErrchk(cudaFree(query_start_gpu));
    cudaErrchk(cudaFree(query_end_gpu));
    cudaErrchk(cudaFree(scores_gpu));
  }
};

void asynch_mem_copies_htd(gpu_alignments* gpu_data, unsigned* offsetA_h, unsigned* offsetB_h, char* strA, char* strA_d, char* strB,
                           char* strB_d, unsigned half_length_A, unsigned half_length_B, unsigned totalLengthA,
                           unsigned totalLengthB, int sequences_per_stream, int sequences_stream_leftover,
                           cudaStream_t* streams_cuda, int max_rlen) {
  cudaErrchk(cudaMemcpyAsync(gpu_data->offset_ref_gpu, offsetA_h, (sequences_per_stream) * sizeof(int), cudaMemcpyHostToDevice,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(gpu_data->offset_ref_gpu + sequences_per_stream, offsetA_h + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(int), cudaMemcpyHostToDevice,
                             streams_cuda[1]));

  cudaErrchk(cudaMemcpyAsync(gpu_data->offset_query_gpu, offsetB_h, (sequences_per_stream) * sizeof(int), cudaMemcpyHostToDevice,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(gpu_data->offset_query_gpu + sequences_per_stream, offsetB_h + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(int), cudaMemcpyHostToDevice,
                             streams_cuda[1]));

  cudaErrchk(cudaMemcpyAsync(strA_d, strA, half_length_A * sizeof(char), cudaMemcpyHostToDevice, streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(strA_d + half_length_A, strA + half_length_A, (totalLengthA - half_length_A) * sizeof(char),
                             cudaMemcpyHostToDevice, streams_cuda[1]));

  size_t size_strB = sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE;
  if (size_strB < totalLengthB) {
    std::cerr << "<" << __FILE__ << ":" << __LINE__ << "> ERROR: size_strB " << size_strB << " < "
              << "totalLengthB " << totalLengthB << " max rlen " << max_rlen << "\n";
    std::abort();
  }
  if (size_strB < half_length_B) {
    std::cerr << "<" << __FILE__ << ":" << __LINE__ << "> ERROR: size_strB " << size_strB << " < "
              << "half_length_B " << half_length_B << " max rlen " << max_rlen << "\n";
    std::abort();
  }
  cudaErrchk(cudaMemcpyAsync(strB_d, strB, half_length_B * sizeof(char), cudaMemcpyHostToDevice, streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(strB_d + half_length_B, strB + half_length_B, (totalLengthB - half_length_B) * sizeof(char),
                             cudaMemcpyHostToDevice, streams_cuda[1]));
}

void asynch_mem_copies_dth_mid(gpu_alignments* gpu_data, short* alAend, short* alBend, int sequences_per_stream,
                               int sequences_stream_leftover, cudaStream_t* streams_cuda) {
  cudaErrchk(cudaMemcpyAsync(alAend, gpu_data->ref_end_gpu, sequences_per_stream * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(alAend + sequences_per_stream, gpu_data->ref_end_gpu + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[1]));

  cudaErrchk(cudaMemcpyAsync(alBend, gpu_data->query_end_gpu, sequences_per_stream * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(alBend + sequences_per_stream, gpu_data->query_end_gpu + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[1]));
}

void asynch_mem_copies_dth(gpu_alignments* gpu_data, short* alAbeg, short* alBbeg, short* top_scores_cpu, int sequences_per_stream,
                           int sequences_stream_leftover, cudaStream_t* streams_cuda) {
  cudaErrchk(cudaMemcpyAsync(alAbeg, gpu_data->ref_start_gpu, sequences_per_stream * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(alAbeg + sequences_per_stream, gpu_data->ref_start_gpu + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[1]));

  cudaErrchk(cudaMemcpyAsync(alBbeg, gpu_data->query_start_gpu, sequences_per_stream * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(alBbeg + sequences_per_stream, gpu_data->query_start_gpu + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[1]));

  cudaErrchk(cudaMemcpyAsync(top_scores_cpu, gpu_data->scores_gpu, sequences_per_stream * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[0]));
  cudaErrchk(cudaMemcpyAsync(top_scores_cpu + sequences_per_stream, gpu_data->scores_gpu + sequences_per_stream,
                             (sequences_per_stream + sequences_stream_leftover) * sizeof(short), cudaMemcpyDeviceToHost,
                             streams_cuda[1]));
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
  cudaStream_t streams_cuda[NSTREAMS];
  unsigned* offsetA_h;
  unsigned* offsetB_h;
  char *strA_d, *strB_d;
  char* strA;
  char* strB;
  cudaEvent_t event;
  short matchScore, misMatchScore, startGap, extendGap;
  gpu_alignments* gpu_data;
  unsigned half_length_A = 0;
  unsigned half_length_B = 0;
  int max_rlen = 0;
};

adept_sw::GPUDriver::GPUDriver(int upcxx_rank_me, int upcxx_rank_n, short match_score, short mismatch_score,
                               short gap_opening_score, short gap_extending_score, int max_rlen, double& init_time) {
  using timepoint_t = std::chrono::time_point<std::chrono::high_resolution_clock>;
  timepoint_t t = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed;
  driver_state = new DriverState();
  driver_state->rank_me = upcxx_rank_me;
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " new=" << elapsed.count();
  driver_state->matchScore = match_score;
  driver_state->misMatchScore = mismatch_score;
  driver_state->startGap = gap_opening_score;
  driver_state->extendGap = gap_extending_score;
  driver_state->max_rlen = max_rlen;

  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " get_num_gpus=" << elapsed.count();

  cudaErrchk(cudaMallocHost(&(alignments.ref_begin), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  cudaErrchk(cudaMallocHost(&(alignments.ref_end), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  cudaErrchk(cudaMallocHost(&(alignments.query_begin), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  cudaErrchk(cudaMallocHost(&(alignments.query_end), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  cudaErrchk(cudaMallocHost(&(alignments.top_scores), sizeof(short) * KLIGN_GPU_BLOCK_SIZE));
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " mallocHost=" << elapsed.count();

  gpu_utils::set_gpu_device(driver_state->rank_me);
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " set_device=" << elapsed.count();

  for (int stm = 0; stm < NSTREAMS; stm++) {
    cudaErrchk(cudaStreamCreate(&driver_state->streams_cuda[stm]));
  }
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " streamcreate=" << elapsed.count();

  cudaErrchk(cudaMallocHost(&driver_state->offsetA_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE));
  cudaErrchk(cudaMallocHost(&driver_state->offsetB_h, sizeof(int) * KLIGN_GPU_BLOCK_SIZE));
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " mallocHost2=" << elapsed.count();

  // FIXME: hack for max contig and read size -> multiplying max_rlen by 2 for contigs
  cudaErrchk(cudaMalloc(&driver_state->strA_d, 2 * max_rlen * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));
  cudaErrchk(cudaMalloc(&driver_state->strB_d, max_rlen * KLIGN_GPU_BLOCK_SIZE * sizeof(char)));
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " mallocs=" << elapsed.count();

  cudaErrchk(cudaMallocHost(&driver_state->strA, 2 * sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE));
  cudaErrchk(cudaMallocHost(&driver_state->strB, sizeof(char) * max_rlen * KLIGN_GPU_BLOCK_SIZE));
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " mallocHost3=" << elapsed.count();

  driver_state->gpu_data = new gpu_alignments(KLIGN_GPU_BLOCK_SIZE);  // gpu mallocs
  // elapsed =  std::chrono::high_resolution_clock::now() - t; os << " final=" << elapsed.count();

  elapsed = std::chrono::high_resolution_clock::now() - t;
  init_time = elapsed.count();
}

adept_sw::GPUDriver::~GPUDriver() {
  // won't have been allocated if there was no GPU present
  if (!alignments.ref_begin) return;
  gpu_utils::set_gpu_device(driver_state->rank_me);
  cudaErrchk(cudaFreeHost(alignments.ref_begin));
  cudaErrchk(cudaFreeHost(alignments.ref_end));
  cudaErrchk(cudaFreeHost(alignments.query_begin));
  cudaErrchk(cudaFreeHost(alignments.query_end));
  cudaErrchk(cudaFreeHost(alignments.top_scores));

  cudaErrchk(cudaFree(driver_state->strA_d));
  cudaErrchk(cudaFree(driver_state->strB_d));
  cudaErrchk(cudaFreeHost(driver_state->offsetA_h));
  cudaErrchk(cudaFreeHost(driver_state->offsetB_h));
  cudaErrchk(cudaFreeHost(driver_state->strA));
  cudaErrchk(cudaFreeHost(driver_state->strB));
  for (int i = 0; i < NSTREAMS; i++) cudaErrchk(cudaStreamDestroy(driver_state->streams_cuda[i]));
  delete driver_state->gpu_data;
  delete driver_state;
}

bool adept_sw::GPUDriver::kernel_is_done() {
  if (cudaEventQuery(driver_state->event) != cudaSuccess) return false;
  cudaErrchk(cudaEventDestroy(driver_state->event));
  return true;
}

void adept_sw::GPUDriver::kernel_block() {
  cudaErrchk(cudaEventSynchronize(driver_state->event));
  cudaErrchk(cudaEventDestroy(driver_state->event));
}

void adept_sw::GPUDriver::run_kernel_forwards(std::vector<std::string>& reads, std::vector<std::string>& contigs,
                                              unsigned maxReadSize, unsigned maxContigSize) {
  gpu_utils::set_gpu_device(driver_state->rank_me);
  unsigned totalAlignments = contigs.size();  // assuming that read and contig vectors are same length

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

  cudaErrchk(cudaEventCreateWithFlags(&driver_state->event, cudaEventDisableTiming | cudaEventBlockingSync));

  asynch_mem_copies_htd(driver_state->gpu_data, driver_state->offsetA_h, driver_state->offsetB_h, driver_state->strA,
                        driver_state->strA_d, driver_state->strB, driver_state->strB_d, driver_state->half_length_A,
                        driver_state->half_length_B, totalLengthA, totalLengthB, sequences_per_stream, sequences_stream_leftover,
                        driver_state->streams_cuda, driver_state->max_rlen);
  unsigned minSize = (maxReadSize < maxContigSize) ? maxReadSize : maxContigSize;
  unsigned totShmem = 3 * (minSize + 1) * sizeof(short);
  unsigned alignmentPad = 4 + (4 - totShmem % 4);
  size_t ShmemBytes = totShmem + alignmentPad;
  if (ShmemBytes > 48000)
    cudaErrchk(cudaFuncSetAttribute(gpu_bsw::sequence_dna_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, ShmemBytes));

  gpu_bsw::sequence_dna_kernel<<<sequences_per_stream, minSize, ShmemBytes, driver_state->streams_cuda[0]>>>(
      driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu, driver_state->gpu_data->offset_query_gpu,
      driver_state->gpu_data->ref_start_gpu, driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
      driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu, driver_state->matchScore,
      driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  gpu_bsw::
      sequence_dna_kernel<<<sequences_per_stream + sequences_stream_leftover, minSize, ShmemBytes, driver_state->streams_cuda[1]>>>(
          driver_state->strA_d + driver_state->half_length_A, driver_state->strB_d + driver_state->half_length_B,
          driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
          driver_state->gpu_data->offset_query_gpu + sequences_per_stream,
          driver_state->gpu_data->ref_start_gpu + sequences_per_stream, driver_state->gpu_data->ref_end_gpu + sequences_per_stream,
          driver_state->gpu_data->query_start_gpu + sequences_per_stream,
          driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
          driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  // copyin back end index so that we can find new min
  asynch_mem_copies_dth_mid(driver_state->gpu_data, alAend, alBend, sequences_per_stream, sequences_stream_leftover,
                            driver_state->streams_cuda);

  cudaErrchk(cudaEventRecord(driver_state->event));
}

void adept_sw::GPUDriver::run_kernel_backwards(std::vector<std::string>& reads, std::vector<std::string>& contigs,
                                               unsigned maxReadSize, unsigned maxContigSize) {
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

  cudaErrchk(cudaEventCreateWithFlags(&driver_state->event, cudaEventDisableTiming | cudaEventBlockingSync));
  gpu_bsw::sequence_dna_reverse<<<sequences_per_stream, newMin, ShmemBytes, driver_state->streams_cuda[0]>>>(
      driver_state->strA_d, driver_state->strB_d, driver_state->gpu_data->offset_ref_gpu, driver_state->gpu_data->offset_query_gpu,
      driver_state->gpu_data->ref_start_gpu, driver_state->gpu_data->ref_end_gpu, driver_state->gpu_data->query_start_gpu,
      driver_state->gpu_data->query_end_gpu, driver_state->gpu_data->scores_gpu, driver_state->matchScore,
      driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  gpu_bsw::
      sequence_dna_reverse<<<sequences_per_stream + sequences_stream_leftover, newMin, ShmemBytes, driver_state->streams_cuda[1]>>>(
          driver_state->strA_d + driver_state->half_length_A, driver_state->strB_d + driver_state->half_length_B,
          driver_state->gpu_data->offset_ref_gpu + sequences_per_stream,
          driver_state->gpu_data->offset_query_gpu + sequences_per_stream,
          driver_state->gpu_data->ref_start_gpu + sequences_per_stream, driver_state->gpu_data->ref_end_gpu + sequences_per_stream,
          driver_state->gpu_data->query_start_gpu + sequences_per_stream,
          driver_state->gpu_data->query_end_gpu + sequences_per_stream, driver_state->gpu_data->scores_gpu + sequences_per_stream,
          driver_state->matchScore, driver_state->misMatchScore, driver_state->startGap, driver_state->extendGap);

  asynch_mem_copies_dth(driver_state->gpu_data, alAbeg, alBbeg, top_scores_cpu, sequences_per_stream, sequences_stream_leftover,
                        driver_state->streams_cuda);
  cudaErrchk(cudaEventRecord(driver_state->event));
}
