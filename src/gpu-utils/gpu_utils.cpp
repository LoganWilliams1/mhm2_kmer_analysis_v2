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

#include "upcxx_utils/colors.h"
#include "gpu_compatibility.hpp"
#include "gpu_utils.hpp"
#include "gpu_common.hpp"

using namespace std;

static int _device_count = 0;
static int _rank_me = -1;

static int get_gpu_device_count() {
  static bool tested = false;
  if (!_device_count && !tested) {
    bool success = false;
    auto res = GetDeviceCount(&_device_count);
    if (res == Success) success = true;
    tested = true;
    if (!success) {
      _device_count = 0;
      return 0;
    }
  }
  return _device_count;
}

// singleton to cache and avoid overloading calls to get properties that do not change
// defaults to the current device that is set
static DeviceProp &get_gpu_properties(int device_id = -1) {
  static vector<DeviceProp> _{};
  auto num_devs = get_gpu_device_count();
  if (num_devs <= 0) {
    std::cerr << KLRED "Cannot get GPU properties when there are no GPUs\n" KNORM;
    exit(1);
  }
  if (device_id == -1) {
    ERROR_CHECK(GetDevice(&device_id));
  }
  if (device_id < 0 || device_id >= num_devs) {
    std::cerr << KLRED "Cannot get GPU properties for device " << device_id << " when there are " << num_devs << " GPUs\n" KNORM;
    exit(1);
  }
  if (_.empty()) {
    _.resize(get_gpu_device_count());
    for (int i = 0; i < num_devs; ++i) {
      auto idx = (1 + i + _rank_me) % num_devs;  // stagger access to devices, end on delegated device
      ERROR_CHECK(GetDeviceProperties(&_[idx], idx));
    }
  }
  return _[device_id];
}

void gpu_utils::set_gpu_device(int rank_me) {
#ifdef ENABLE_KOKKOS
  return;
#endif

  if (rank_me == -1) {
    std::cerr << KLRED "Cannot set GPU device for rank -1; device is not yet initialized\n" KNORM;
    exit(1);
  }
  int current_device = -1;
  int num_devs = get_gpu_device_count();
  if (num_devs == 0) return;
  auto mod = rank_me % num_devs;
  ERROR_CHECK(GetDevice(&current_device));
  if (current_device != mod) {
    ERROR_CHECK(SetDevice(mod));
    ERROR_CHECK(GetDevice(&current_device));
    if (current_device != mod) {
      std::cerr << KLRED "Did not set device to " << mod << " rank_me=" << rank_me << "\n" KNORM;
      exit(1);
    }
  }
}

size_t gpu_utils::get_gpu_tot_mem() {
  set_gpu_device(_rank_me);
  return get_gpu_properties().totalGlobalMem;
}

size_t gpu_utils::get_gpu_avail_mem() {
  set_gpu_device(_rank_me);
  size_t free_mem, tot_mem;
  ERROR_CHECK(MemGetInfo(&free_mem, &tot_mem));
  return free_mem;
}

string gpu_utils::get_gpu_device_name() {
  set_gpu_device(_rank_me);
  return get_gpu_properties().name;
}

static string get_uuid_str(char uuid_bytes[16]) {
  ostringstream os;
  for (int i = 0; i < 16; i++) {
    os << std::setfill('0') << std::setw(2) << std::hex << (0xff & (unsigned int)uuid_bytes[i]);
  }
  return os.str();
}

vector<string> gpu_utils::get_gpu_uuids() {
  static vector<string> uuids = []() {
    vector<string> uuids;
    int num_devs = get_gpu_device_count();
    for (int i = 0; i < num_devs; ++i) {
      DeviceProp &prop = get_gpu_properties(i);
      bool set_dev = false;
#ifdef CUDA_GPU
#if (CUDA_VERSION >= 10000)
      uuids.push_back(get_uuid_str(prop.uuid.bytes));
      set_dev = true;
#endif
#endif
      if (!set_dev) {
        ostringstream os;
        os << prop.name << ':' << prop.pciDeviceID << ':' << prop.pciBusID;
        uuids.push_back(os.str());
        set_dev = true;
      }
    }
    return uuids;
  }();
  return uuids;
}

string gpu_utils::get_gpu_uuid() {
  set_gpu_device(_rank_me);
  auto uuids = get_gpu_uuids();
  int current_device = -1;
  ERROR_CHECK(GetDevice(&current_device));
  return uuids[current_device] + " device" + to_string(current_device) + "of" + to_string(get_gpu_device_count());
}

bool gpu_utils::gpus_present() { return get_gpu_device_count(); }

void gpu_utils::initialize_gpu(double &time_to_initialize, int rank_me, int team_rank_n) {
  using timepoint_t = chrono::time_point<chrono::high_resolution_clock>;
  timepoint_t t = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed;

  if (!gpus_present()) return;
  _rank_me = rank_me;
  set_gpu_device(_rank_me);
  ERROR_CHECK(DeviceReset());

  auto gpu_mem_avail_per_rank = get_gpu_avail_mem() / team_rank_n;

  char *buf;
  size_t buf_size = gpu_mem_avail_per_rank * 0.8;
  ERROR_CHECK(Malloc(&buf, buf_size));
  ERROR_CHECK(Memset(buf, 0, buf_size));
  ERROR_CHECK(Free(buf));

  elapsed = chrono::high_resolution_clock::now() - t;
  time_to_initialize = elapsed.count();
}

string gpu_utils::get_gpu_device_descriptions() {
  int num_devs = get_gpu_device_count();
  ostringstream os;
  os << "Number of GPU devices visible: " << num_devs << "\n";
  auto uuids = get_gpu_uuids();
  for (int i = 0; i < num_devs; ++i) {
    DeviceProp &prop = get_gpu_properties(i);

    os << "GPU Device number: " << i << "\n";
    os << "  Device name: " << prop.name << "\n";
    os << "  PCI device ID: " << prop.pciDeviceID << "\n";
    os << "  PCI bus ID: " << prop.pciBusID << "\n";
    os << "  UUID: " << uuids[i] << "\n";
    os << "  PCI domainID: " << prop.pciDomainID << "\n";
    os << "  Compute capability: " << prop.major << "." << prop.minor << "\n";
    os << "  Clock Rate: " << prop.clockRate << "kHz\n";
    os << "  Total SMs: " << prop.multiProcessorCount << "\n";
#ifdef CUDA_GPU
    os << "  MultiGPUBoardGroupID: " << prop.multiGpuBoardGroupID << "\n";
    os << "  Shared Memory Per SM: " << prop.sharedMemPerMultiprocessor << " bytes\n";
    os << "  Registers Per SM: " << prop.regsPerMultiprocessor << " 32-bit\n";
#endif
#ifdef HIP_GPU
    os << "  Max Shared Memory Per SM: " << prop.maxSharedMemoryPerMultiProcessor << " bytes\n";
    os << "  Registers Per Block: " << prop.regsPerBlock << " 32-bit\n";
#endif
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
