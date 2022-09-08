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

/* Compatibility between HIP and Cuda Libraries from AMD_HIP_Supported_CUDA_API_Reference_Guide.pdf */

#ifdef HIP_GPU
  #include <hip/hip_runtime_api.h>
  #include <hip/hip_runtime.h>
  #include <hip/hip_profile.h>
  #define MHM2_GPU hip
  #define LaunchKernel(func,blocks,threads_per_block) hipLaunchKernelGGL(func, blocks, threads_per_block, 0,0,
  #define FreeHost hipHostFree
#endif
#ifdef CUDA_GPU
  #include <cuda_runtime_api.h>
  #include <cuda.h>
  #include <cuda_profiler_api.h>
  #define MHM2_GPU cuda
  #define FreeHost cudaFreeHost
  #define LaunchKernel(func,blocks,threads_per_block) func<<<blocks,threads_per_block>>>(
#endif

#define concat(a,b) a##b

#define Success                concat(MHM2_GPU, Success)
#define Error_t                concat(MHM2_GPU, Error_t)
#define GetErrorString         concat(MHM2_GPU, GetErrorString)
#define Stream_t               concat(MHM2_GPU, Stream_t)

#define Alloc                  concat(MHM2_GPU, Alloc)
#define Free                   concat(MHM2_GPU, Free)
#define HostAlloc              concat(MHM2_GPU, HostAlloc)
#define HostFree               concat(MHM2_GPU, HostFree)
#define MallocHost             concat(MHM2_GPU, MallocHost)

#define GetDeviceCount         concat(MHM2_GPU, GetDeviceCount)
#define DeviceProp             concat(MHM2_GPU, DeviceProp)
#define GetDevice              concat(MHM2_GPU, GetDevice)
#define GetDeviceProperties    concat(MHM2_GPU, GetDeviceProperties)
#define SetDevice              concat(MHM2_GPU, SetDevice)
#define DeviceReset            concat(MHM2_GPU, DeviceReset)
#define DeviceSynchronize      concat(MHM2_GPU, DeviceSynchronize)

#define MemGetInfo             concat(MHM2_GPU, MemGetInfo)
#define Memcpy                 concat(MHM2_GPU, Memcpy)
#define MemcpyAsync            concat(MHM2_GPU, MemcpyAsync)
#define MemcpyHostToDevice     concat(MHM2_GPU, MemcpyHostToDevice)
#define MemcpyDeviceToHost     concat(MHM2_GPU, MemcpyDeviceToHost)
#define Memset                 concat(MHM2_GPU, Memset)

#define HostAlloc              concat(MHM2_GPU, HostAlloc)

#define StreamCreate           concat(MHM2_GPU, StreamCreate)
#define StreamDestroy          concat(MHM2_GPU, StreamDestroy)
#define StreamSynchronize      concat(MHM2_GPU, StreamSynchronize)


#define EventCreate            concat(MHM2_GPU, EventCreate)
#define EventCreateWithFlags   concat(MHM2_GPU, EventCreateWithFlags)
#define EventDestroy           concat(MHM2_GPU, EventDestroy)
#define EventRecord            concat(MHM2_GPU, EventRecord)
#define EventQuery             concat(MHM2_GPU, EventQuery)
#define EventSynchronize       concat(MHM2_GPU, EventSynchronize)
#define EventElapsedTime       concat(MHM2_GPU, EventElapsedTime)
#define EventDisableTiming     concat(MHM2_GPU, EventDisableTiming)
#define EventBlockingSync      concat(MHM2_GPU, EventBlockingSync)

#define OccupancyMaxPotentialBlockSize concat(MHM2_GPU, OccupancyMaxPotentialBlockSize)


#undef concat

// Functions that are common to all cuda code; not to be used by upcxx code

#define ERROR_CHECK(ans) gpu_common::gpu_die((ans), __FILE__, __LINE__)

inline void
gpuAssert(Error_t code, const char* file, int line, bool abort = true)
{
    if(code != Success)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", GetErrorString(code), file, line);
        if(abort)
            exit(code);
    }
}

#define GPU_ASSERT(ans) gpuAssert((ans), __FILE__, __LINE__);

// we are typecasting uint64_t into this, so we need to check them
static_assert(sizeof(unsigned long long) == sizeof(uint64_t));
