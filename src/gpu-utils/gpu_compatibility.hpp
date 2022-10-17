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

#define LaunchKernel(func, blocks, threads_per_block, args...) hipLaunchKernelGGL(func, blocks, threads_per_block, 0, 0, args)
#define LaunchKernelGGL(func, blocks, threads_per_block, x, y, args...) \
  hipLaunchKernelGGL(func, blocks, threads_per_block, x, y, args)

#define Success hipSuccess
#define GetErrorString hipGetErrorString
#define Error_t hipError_t
#define Event_t hipEvent_t
#define Stream_t hipStream_t

#define Malloc hipMalloc
#define Free hipFree
#define HostAlloc hipHostAlloc
#define HostFree hipHostFree
#define FreeHost hipHostFree
#define MallocHost hipHostMalloc

#define GetDeviceCount hipGetDeviceCount
#define DeviceProp hipDeviceProp_t
#define GetDevice hipGetDevice
#define GetDeviceProperties hipGetDeviceProperties
#define SetDevice hipSetDevice
#define DeviceReset hipDeviceReset
#define DeviceSynchronize hipDeviceSynchronize

#define MemGetInfo hipMemGetInfo
#define Memcpy hipMemcpy
#define MemcpyAsync hipMemcpyAsync
#define MemcpyHostToDevice hipMemcpyHostToDevice
#define MemcpyDeviceToHost hipMemcpyDeviceToHost
#define Memset hipMemset

#define HostAlloc hipHostAlloc

#define StreamCreate hipStreamCreate
#define StreamDestroy hipStreamDestroy
#define StreamSynchronize hipStreamSynchronize

#define EventCreate hipEventCreate
#define EventCreateWithFlags hipEventCreateWithFlags
#define EventDestroy hipEventDestroy
#define EventRecord hipEventRecord
#define EventQuery hipEventQuery
#define EventSynchronize hipEventSynchronize
#define EventElapsedTime hipEventElapsedTime
#define EventDisableTiming hipEventDisableTiming
#define EventBlockingSync hipEventBlockingSync

#define OccupancyMaxPotentialBlockSize hipOccupancyMaxPotentialBlockSize

#define FuncSetAttribute hipFuncSetAttribute
#define FuncGetAttribute hipFuncGetAttribute
#define FuncAttributeMaxDynamicSharedMemorySize hipFuncAttributeMaxDynamicSharedMemorySize

#endif
#ifdef CUDA_GPU
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cuda_profiler_api.h>
#define MHM2_GPU cuda
#define FreeHost cudaFreeHost
#define LaunchKernel(func, blocks, threads_per_block, args...) func<<<blocks, threads_per_block>>>(args)
#define LaunchKernelGGL(func, blocks, threads_per_block, x, y, args...) func<<<blocks, threads_per_block, x, y>>>(args)
#define Success cudaSuccess
#define GetErrorString cudaGetErrorString
#define Error_t cudaError_t
#define Event_t cudaEvent_t
#define Stream_t cudaStream_t

#define Malloc cudaMalloc
#define Free cudaFree
#define HostAlloc cudaHostAlloc
#define HostFree cudaHostFree
#define MallocHost cudaMallocHost

#define GetDeviceCount cudaGetDeviceCount
#define DeviceProp cudaDeviceProp
#define GetDevice cudaGetDevice
#define GetDeviceProperties cudaGetDeviceProperties
#define SetDevice cudaSetDevice
#define DeviceReset cudaDeviceReset
#define DeviceSynchronize cudaDeviceSynchronize

#define MemGetInfo cudaMemGetInfo
#define Memcpy cudaMemcpy
#define MemcpyAsync cudaMemcpyAsync
#define MemcpyHostToDevice cudaMemcpyHostToDevice
#define MemcpyDeviceToHost cudaMemcpyDeviceToHost
#define Memset cudaMemset

#define HostAlloc cudaHostAlloc

#define StreamCreate cudaStreamCreate
#define StreamDestroy cudaStreamDestroy
#define StreamSynchronize cudaStreamSynchronize

#define EventCreate cudaEventCreate
#define EventCreateWithFlags cudaEventCreateWithFlags
#define EventDestroy cudaEventDestroy
#define EventRecord cudaEventRecord
#define EventQuery cudaEventQuery
#define EventSynchronize cudaEventSynchronize
#define EventElapsedTime cudaEventElapsedTime
#define EventDisableTiming cudaEventDisableTiming
#define EventBlockingSync cudaEventBlockingSync

#define OccupancyMaxPotentialBlockSize cudaOccupancyMaxPotentialBlockSize

#define FuncSetAttribute cudaFuncSetAttribute
#define FuncGetAttribute cudaFuncGetAttribute
#define FuncAttributeMaxDynamicSharedMemorySize cudaFuncAttributeMaxDynamicSharedMemorySize

#endif
