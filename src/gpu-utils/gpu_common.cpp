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
#include <array>

#include "upcxx_utils/colors.h"
#include "gpu_compatibility.hpp"
#include "gpu_common.hpp"

#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

namespace gpu_common {

void gpuAssert(Error_t code, const char* file, int line, bool abort) {
  if (code != Success) {
    fprintf(stderr, "GPUassert: %s %s %d\n", GetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

void gpu_die(Error_t code, const char* file, int line, bool abort) {
  if (code != Success) {
    std::cerr << KLRED << "<" << file << ":" << line << "> ERROR:" << KNORM << GetErrorString(code) << "\n";
    std::abort();
    // do not throw exceptions -- does not work properly within progress() throw std::runtime_error(outstr);
  }
}

QuickTimer::QuickTimer()
    : secs(0) {}

void QuickTimer::start() { t = clock_now(); }

void QuickTimer::stop() {
  std::chrono::duration<double> t_elapsed = clock_now() - t;
  secs += t_elapsed.count();
}

void QuickTimer::inc(double s) { secs += s; }

double QuickTimer::get_elapsed() { return secs; }

GPUTimer::GPUTimer() {
  ERROR_CHECK(EventCreate(&start_event));
  ERROR_CHECK(EventCreate(&stop_event));
  elapsed_t_ms = 0;
}

GPUTimer::~GPUTimer() {
  ERROR_CHECK(EventDestroy(start_event));
  ERROR_CHECK(EventDestroy(stop_event));
}

void GPUTimer::start() { ERROR_CHECK(EventRecord(start_event, 0)); }

void GPUTimer::stop() {
  ERROR_CHECK(EventRecord(stop_event, 0));
  ERROR_CHECK(EventSynchronize(stop_event));
  float ms;
  ERROR_CHECK(EventElapsedTime(&ms, start_event, stop_event));
  elapsed_t_ms += ms;
}

double GPUTimer::get_elapsed() { return elapsed_t_ms / 1000.0; }



#ifdef ENABLE_KOKKOS
KOKKOS_INLINE_FUNCTION
#endif
#ifndef ENABLE_KOKKOS
inline __device__
#endif
void revcomp(uint64_t *longs, uint64_t *rc_longs, int kmer_len, int num_longs) {
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

#ifdef ENABLE_KOKKOS
KOKKOS_INLINE_FUNCTION
#endif
#ifndef ENABLE_KOKKOS
inline __device__
#endif
bool pack_seq_to_kmer(char *seqs, int kmer_len, int num_longs, uint64_t *kmer) {
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

}  // namespace gpu_common
