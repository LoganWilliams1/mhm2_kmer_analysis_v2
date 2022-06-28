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
#include <sstream>
#include <chrono>
#include <array>
#include <iomanip>
#include <hip/hip_runtime_api.h>
#include <hip/hip_runtime.h>

#include "gpu_utils.hpp"
#include "upcxx_utils/colors.h"
#include "gpu_common.hpp"

using namespace std;

static int _device_count = 0;
static int _rank_me = -1;

static int get_gpu_device_count() {
  if (!_device_count) {
    auto res = hipGetDeviceCount(&_device_count);
    if (res != hipSuccess) return 0;
  }
  return _device_count;
}

void gpu_utils::set_gpu_device(int rank_me) {
  if (rank_me == -1) {
    std::cerr << KLRED "Cannot set GPU device for rank -1; device is not yet initialized\n" KNORM;
    exit(1);
  }
  int num_devs = get_gpu_device_count();
  ERROR_CHECK(hipSetDevice(rank_me % num_devs));
}

size_t gpu_utils::get_gpu_tot_mem() {
  set_gpu_device(_rank_me);
  hipDeviceProp_t prop;
  ERROR_CHECK(hipGetDeviceProperties(&prop, 0));
  return prop.totalGlobalMem;
}

size_t gpu_utils::get_gpu_avail_mem() {
  set_gpu_device(_rank_me);
  size_t free_mem, tot_mem;
  ERROR_CHECK(hipMemGetInfo(&free_mem, &tot_mem));
  return free_mem;
}

string gpu_utils::get_gpu_device_name() {
  set_gpu_device(_rank_me);
  hipDeviceProp_t prop;
  ERROR_CHECK(hipGetDeviceProperties(&prop, 0));
  return prop.name;
}

static string get_uuid_str(char uuid_bytes[16]) {
  ostringstream os;
  for (int i = 0; i < 16; i++) {
    os << std::setfill('0') << std::setw(2) << std::hex << (0xff & (unsigned int)uuid_bytes[i]);
  }
  return os.str();
}

vector<string> gpu_utils::get_gpu_uuids() {
  vector<string> uuids;
  int num_devs = get_gpu_device_count();
  for (int i = 0; i < num_devs; ++i) {
    hipDeviceProp_t prop;
    ERROR_CHECK(hipGetDeviceProperties(&prop, i));
// #if (CUDA_VERSION >= 10000)
//     uuids.push_back(get_uuid_str(prop.uuid.bytes));
// #else
    ostringstream os;
    os << prop.name << ':' << prop.pciDeviceID << ':' << prop.pciBusID;// << ':' << prop.pciDomainID << prop.multiGpuBoardGroupID;
    uuids.push_back(os.str());
// #endif
  }
  return uuids;
}

string gpu_utils::get_gpu_uuid() {
  auto uuids = get_gpu_uuids();
  return uuids[_rank_me % uuids.size()];
}

bool gpu_utils::gpus_present() { return get_gpu_device_count(); }

void gpu_utils::initialize_gpu(double& time_to_initialize, int rank_me) {
  using timepoint_t = chrono::time_point<chrono::high_resolution_clock>;
  timepoint_t t = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed;

  if (!gpus_present()) return;
  _rank_me = rank_me;
  set_gpu_device(_rank_me);
  ERROR_CHECK(hipDeviceReset());
  elapsed = chrono::high_resolution_clock::now() - t;
  time_to_initialize = elapsed.count();
}

string gpu_utils::get_gpu_device_descriptions() {
  hipDeviceProp_t prop;
  int num_devs = get_gpu_device_count();
  ostringstream os;
  os << "Number of GPU devices visible: " << num_devs << "\n";
  for (int i = 0; i < num_devs; ++i) {
    ERROR_CHECK(hipGetDeviceProperties(&prop, i));

    os << "GPU Device number: " << i << "\n";
    os << "  Device name: " << prop.name << "\n";
    os << "  PCI device ID: " << prop.pciDeviceID << "\n";
    os << "  PCI bus ID: " << prop.pciBusID << "\n";
    // os << "  PCI domainID: " << prop.pciDomainID << "\n";
    // os << "  MultiGPUBoardGroupID: " << prop.multiGpuBoardGroupID << "\n";
// #if (CUDA_VERSION >= 10000)
//     os << "  UUID: " << get_uuid_str(prop.uuid.bytes) << "\n";
// #endif
    os << "  Compute capability: " << prop.major << "." << prop.minor << "\n";
    os << "  Clock Rate: " << prop.clockRate << "kHz\n";
    os << "  Total SMs: " << prop.multiProcessorCount << "\n";
    os << "  Max Shared Memory Per SM: " << prop.maxSharedMemoryPerMultiProcessor << " bytes\n";
    os << "  Registers Per Block: " << prop.regsPerBlock << " 32-bit\n";
    os << "  Max threads per SM: " << prop.maxThreadsPerMultiProcessor << "\n";
    os << "  L2 Cache Size: " << prop.l2CacheSize << " bytes\n";
    os << "  Total Global Memory: " << prop.totalGlobalMem << " bytes\n";
    os << "  Memory Clock Rate: " << prop.memoryClockRate << " kHz\n\n";

    os << "  Max threads per block: " << prop.maxThreadsPerBlock << "\n";
    os << "  Max threads in X-dimension of block: " << prop.maxThreadsDim[0] << "\n";
    os << "  Max threads in Y-dimension of block: " << prop.maxThreadsDim[1] << "\n";
    os << "  Max threads in Z-dimension of block: " << prop.maxThreadsDim[2] << "\n\n";

    os << "  Max blocks in X-dimension of grid: " << prop.maxGridSize[0] << "\n";
    os << "  Max blocks in Y-dimension of grid: " << prop.maxGridSize[1] << "\n";
    os << "  Max blocks in Z-dimension of grid: " << prop.maxGridSize[2] << "\n\n";

    os << "  Shared Memory Per Block: " << prop.sharedMemPerBlock << " bytes\n";
    os << "  Registers Per Block: " << prop.regsPerBlock << " 32-bit\n";
    os << "  Warp size: " << prop.warpSize << "\n\n";
  }
  return os.str();
}
