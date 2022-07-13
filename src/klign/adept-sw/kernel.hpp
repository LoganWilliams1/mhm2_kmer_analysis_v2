#include "hip/hip_runtime.h"
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

#pragma once
#include "hip/hip_runtime.h"
#include <sys/time.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#define NUM_OF_AA 21
#define ENCOD_MAT_SIZE 91
#define SCORE_MAT_SIZE 576

namespace gpu_bsw {

__device__ short
warpReduceMax_with_index(short val, short& myIndex, short& myIndex2, unsigned lengthSeqB, bool reverse);

__device__ short
warpReduceMax(short val, unsigned lengthSeqB);

__device__ short
blockShuffleReduce_with_index(short myVal, short& myIndex, short& myIndex2, unsigned lengthSeqB, bool reverse);

__device__ short
blockShuffleReduce(short val, unsigned lengthSeqB);

__device__ __host__ short
           findMax(short array[], int length, int* ind);

__device__ __host__ short
            findMaxFour(short first, short second, short third, short fourth);

__device__ void
traceBack(short current_i, short current_j, short* seqA_align_begin,
          short* seqB_align_begin, const char* seqA, const char* seqB, short* I_i,
          short* I_j, unsigned lengthSeqB, unsigned lengthSeqA, unsigned int* diagOffset);

__global__ void
dna_kernel(char* seqA_array, char* seqB_array, unsigned* prefix_lengthA,
                unsigned* prefix_lengthB, short* seqA_align_begin, short* seqA_align_end,
                short* seqB_align_begin, short* seqB_align_end, short* top_scores, short matchScore, short misMatchScore, short startGap, short extendGap, bool reverse);

__global__ void
aa_kernel(char* seqA_array, char* seqB_array, unsigned* prefix_lengthA,unsigned* prefix_lengthB, short* seqA_align_begin, short* seqA_align_end, 
            short* seqB_align_begin, short* seqB_align_end, short* top_scores, short startGap, short extendGap, short* scoring_matrix, short* encoding_matrix, bool reverse);

}  // namespace gpu_bsw
