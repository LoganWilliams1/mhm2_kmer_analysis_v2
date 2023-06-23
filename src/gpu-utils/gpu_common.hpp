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
#include <chrono>

#include "gpu_compatibility.hpp"

namespace gpu_common {

static __constant__ uint64_t GPU_TWINS[256] = {
    0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F, 0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F, 0xFB, 0xBB, 0x7B, 0x3B,
    0xEB, 0xAB, 0x6B, 0x2B, 0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B, 0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
    0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07, 0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23, 0xD3, 0x93, 0x53, 0x13,
    0xC3, 0x83, 0x43, 0x03, 0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E, 0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
    0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A, 0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A, 0xF6, 0xB6, 0x76, 0x36,
    0xE6, 0xA6, 0x66, 0x26, 0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06, 0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
    0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02, 0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D, 0xDD, 0x9D, 0x5D, 0x1D,
    0xCD, 0x8D, 0x4D, 0x0D, 0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29, 0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
    0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25, 0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05, 0xF1, 0xB1, 0x71, 0x31,
    0xE1, 0xA1, 0x61, 0x21, 0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01, 0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
    0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C, 0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28, 0xD8, 0x98, 0x58, 0x18,
    0xC8, 0x88, 0x48, 0x08, 0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24, 0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
    0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20, 0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00};

// Functions that are common to all cuda code; not to be used by upcxx code

#define ERROR_CHECK(ans) gpu_common::gpu_die((ans), __FILE__, __LINE__)

#define GPU_CHECK(ans) gpuAssert((ans), __FILE__, __LINE__);

void gpuAssert(Error_t code, const char *file, int line, bool abort = true);

void gpu_die(Error_t code, const char *file, int line, bool abort = true);

// we are typecasting uint64_t into this, so we need to check them
static_assert(sizeof(unsigned long long) == sizeof(uint64_t));

using timepoint_t = std::chrono::time_point<std::chrono::high_resolution_clock>;

class QuickTimer {
  timepoint_t t;
  double secs = 0;

 public:
  QuickTimer();
  void start();
  void stop();
  void inc(double s);
  double get_elapsed();
};

class GPUTimer {
  Event_t start_event, stop_event;
  float elapsed_t_ms = 0;

 public:
  GPUTimer();
  ~GPUTimer();
  void start();
  void stop();
  double get_elapsed();
};

template <typename T>
inline __device__ T warpReduceSum(T val, int n) {
  unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
#ifdef HIP_GPU
  for (int offset = warpSize / 2; offset > 0; offset /= 2) val += __shfl_down(val, offset);
#endif
#ifdef CUDA_GPU
  unsigned mask = __ballot_sync(0xffffffff, threadid < n);
  for (int offset = warpSize / 2; offset > 0; offset /= 2) val += __shfl_down_sync(mask, val, offset);
#endif
  return val;
};

template <typename T>
inline __device__ T blockReduceSum(T val, int n) {
  static __shared__ T shared[32];  // Shared mem for 32 partial sums
  int lane_id = threadIdx.x % warpSize;
  int warp_id = threadIdx.x / warpSize;

  val = warpReduceSum(val, n);              // Each warp performs partial reduction

  if (lane_id == 0) shared[warp_id] = val;  // Write reduced value to shared memory

  __syncthreads();                          // Wait for all partial reductions

  // read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane_id] : 0;

  if (warp_id == 0) val = warpReduceSum(val, n);  // Final reduce within first warp
  __syncthreads();
  return val;
};

inline __device__ void reduce(uint32_t count, int num, uint32_t *result) {
  uint32_t block_num = blockReduceSum(count, num);
  if (threadIdx.x == 0) atomicAdd(result, block_num);
}

inline __device__ void reduce(uint64_t count, int num, uint64_t *result) {
  unsigned long long int block_num = blockReduceSum(count, num);
  if (threadIdx.x == 0) atomicAdd((unsigned long long int *)result, block_num);
}

template <class T>
inline void get_kernel_config(unsigned max_val, T func, int &gridsize, int &threadblocksize) {
  int mingridsize = 0;
  threadblocksize = 0;  // 1024
  ERROR_CHECK(OccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, func, 0, 0));
  gridsize = (max_val + threadblocksize - 1) / threadblocksize;
};

inline __device__ char comp_nucleotide(char ch) {
  switch (ch) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
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
    default: return 0;
  }
}

inline __device__ void revcomp(uint64_t *longs, uint64_t *rc_longs, int kmer_len, int num_longs) {
  int last_long = (kmer_len + 31) / 32;
  for (size_t i = 0; i < last_long; i++) {
    uint64_t v = longs[i];
    rc_longs[last_long - 1 - i] = (GPU_TWINS[v & 0xFF] << 56) | (GPU_TWINS[(v >> 8) & 0xFF] << 48) |
                                  (GPU_TWINS[(v >> 16) & 0xFF] << 40) | (GPU_TWINS[(v >> 24) & 0xFF] << 32) |
                                  (GPU_TWINS[(v >> 32) & 0xFF] << 24) | (GPU_TWINS[(v >> 40) & 0xFF] << 16) |
                                  (GPU_TWINS[(v >> 48) & 0xFF] << 8) | (GPU_TWINS[(v >> 56)]);
  }
  uint64_t shift = (kmer_len % 32) ? 2 * (32 - (kmer_len % 32)) : 0;
  uint64_t shiftmask = (kmer_len % 32) ? (((((uint64_t)1) << shift) - 1) << (64 - shift)) : ((uint64_t)0);
  rc_longs[0] = rc_longs[0] << shift;
  for (size_t i = 1; i < last_long; i++) {
    rc_longs[i - 1] |= (rc_longs[i] & shiftmask) >> (64 - shift);
    rc_longs[i] = rc_longs[i] << shift;
  }
}

inline __device__ bool pack_seq_to_kmer(char *seqs, int kmer_len, int num_longs, uint64_t *kmer) {
  int l = 0, prev_l = 0;
  uint64_t longs = 0;
  memset(kmer, 0, sizeof(uint64_t) * num_longs);
  // each thread extracts one kmer
  for (int k = 0; k < kmer_len; k++) {
    char s = seqs[k];
    switch (s) {
      case 'a': s = 'A'; break;
      case 'c': s = 'C'; break;
      case 'g': s = 'G'; break;
      case 't': s = 'T'; break;
      case 'A':
      case 'C':
      case 'G':
      case 'T': break;
      default: return false;
    }
    int j = k % 32;
    prev_l = l;
    l = k / 32;
    // we do it this way so we can operate on the variable longs in a register, rather than local memory in the array
    if (l > prev_l) {
      kmer[prev_l] = longs;
      longs = 0;
      prev_l = l;
    }
    uint64_t x = (s & 4) >> 1;
    longs |= ((x + ((x ^ (s & 2)) >> 1)) << (2 * (31 - j)));
  }
  kmer[l] = longs;
  return true;
}

inline __device__ uint16_t atomicAddUint16(uint16_t *address, uint16_t val) {
  unsigned int *base_address = (unsigned int *)((size_t)address & ~2);
  unsigned int long_val = ((size_t)address & 2) ? ((unsigned int)val << 16) : val;
  unsigned int long_old = atomicAdd(base_address, long_val);
  return ((size_t)address & 2) ? (uint16_t)(long_old >> 16) : (uint16_t)(long_old & 0xffff);
}

inline __device__ void atomicAddUint16_thres(uint16_t *address, uint16_t val, uint16_t thres) {
  if (atomicAddUint16(address, 0) < thres - val) atomicAddUint16(address, val);
}

inline __device__ short atomicAndUint16(uint16_t *address, uint16_t val) {
  unsigned int *base_address = (unsigned int *)((char *)address - ((size_t)address & 2));
  int shift = ((size_t)address & 2) << 3;
  unsigned int long_val = ((unsigned int)val << shift) | (0xffff0000u >> shift);
  return (atomicAnd(base_address, long_val) >> shift) & 0xffff;
}

}  // namespace gpu_common
