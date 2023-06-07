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

#include <upcxx/upcxx.hpp>

#include "upcxx_utils/thread_pool.hpp"
#include "gpu-utils/gpu_utils.hpp"
#include "devices_gpu.hpp"
#include "utils.hpp"

using namespace std;
using namespace upcxx;
using namespace upcxx_utils;

// #define SLOG_GPU(...) SLOG(KLMAGENTA, __VA_ARGS__, KNORM)
#define SLOG_GPU SLOG_VERBOSE

static bool init_gpu_thread = true;
static future<> detect_gpu_fut;
static double gpu_startup_duration = 0;
static int num_gpus_on_node = 0;

size_t get_gpu_avail_mem_per_rank() {
  auto &gpu_team = get_gpu_team();
  barrier(gpu_team);
  auto avail_mem = gpu_utils::get_gpu_avail_mem() / gpu_team.rank_n();
  barrier(gpu_team);
  return avail_mem;
}

// singleton to fetch UUIDs only once as it can overload the gpus if all ranks access the same one at the same time
static vector<string> &get_gpu_uuids() {
  static vector<string> uuids = {};
  if (uuids.empty() && gpu_utils::gpus_present()) {
    uuids = gpu_utils::get_gpu_uuids();
  }
  return uuids;
}

upcxx::team &get_gpu_team() {
  static upcxx::team tm = []() {
    assert(upcxx::master_persona().active_with_caller() && "Called from master persona");
    upcxx::intrank_t color = upcxx::team::color_none;
    if (gpu_utils::gpus_present()) {
      auto my_uuid = gpu_utils::get_gpu_uuid();
      color = std::hash<string>{}(my_uuid) & 0xffffffff;
      if (color < 0) color = -color;
    } else {
      color = 0;  // i.e. just a copy of the local team
    }
    assert(color != upcxx::team::color_none);
    log_local("GPU team color", std::to_string(color));
    return upcxx::local_team().split(color, upcxx::local_team().rank_me());
  }();
  return tm;
}

void init_devices() {
  SLOG("Initializing GPUs\n");
  init_gpu_thread = true;
  // initialize the GPU and first-touch memory and functions in a new thread as this can take many seconds to complete
  detect_gpu_fut = execute_in_thread_pool([]() {
                     DBG("Initializing GPUs\n");
                     gpu_utils::initialize_gpu(gpu_startup_duration, rank_me());
                     stringstream ss;
                     ss << "Done initializing GPU: " << (gpu_utils::gpus_present() ? "Found" : "NOT FOUND");
                     if (gpu_utils::gpus_present()) {
                       auto uuids = get_gpu_uuids();
                       ss << " with " << uuids.size() << " uuids:\t";
                       for (auto &uuid : uuids) ss << uuid << "\t";
                     }
                     ss << "\n";
                     DBG(ss.str());
                   }).then([]() {
    BaseTimer t;
    t.start();
    // also set the device in the master personna thread for any direct calls this should take practically no time
    gpu_utils::set_gpu_device(rank_me());
    t.stop();
    DBG("Set GPU device on master personna in ", t.get_elapsed(), "s \n");
  });
}

void done_init_devices() {
  if (init_gpu_thread) {
    Timer t("Waiting for GPU to be initialized (should be noop)");
    init_gpu_thread = false;
    detect_gpu_fut.wait();
    auto have_gpus = reduce_all(gpu_utils::gpus_present() ? 1 : 0, op_fast_add).wait();
    if (have_gpus && have_gpus != rank_n()) {
      if (!gpu_utils::gpus_present()) WARN("Found no GPUs\n");
      barrier();
      SDIE("Not all ranks found GPUs: ", have_gpus, " out of ", rank_n(), "\n");
    }
    if (have_gpus == rank_n()) {
      barrier(local_team());
      int num_uuids = 0;
      unordered_set<string> unique_ids;
      dist_object<vector<string>> gpu_uuids(get_gpu_uuids(), local_team());
      for (auto uuid : *gpu_uuids) unique_ids.insert(uuid);
      // auto gpu_avail_mem = 0;
      if (!local_team().rank_me()) {
        for (int i = 1; i < local_team().rank_n(); i++) {
          auto gpu_uuids_i = gpu_uuids.fetch(i).wait();
          num_uuids += gpu_uuids_i.size();
          for (auto uuid : gpu_uuids_i) {
            unique_ids.insert(uuid);
          }
        }
        num_gpus_on_node = unique_ids.size();
        // SLOG_GPU("Found GPU UUIDs:\n");
        // for (auto uuid : unique_ids) {
        //  SLOG_GPU("    ", uuid, "\n");
        //}
        // gpu_utils::set_gpu_device(rank_me());
        // gpu_avail_mem = gpu_utils::get_gpu_avail_mem() * num_gpus_on_node;
      }
      barrier(local_team());
      num_gpus_on_node = broadcast(num_gpus_on_node, 0, local_team()).wait();
      auto msm_local = upcxx_utils::min_sum_max_reduce_one(gpu_startup_duration, 0, upcxx::local_team()).wait();
      auto msm_global = upcxx_utils::min_sum_max_reduce_one(gpu_startup_duration, 0).wait();
      // gpu_avail_mem = broadcast(gpu_avail_mem, 0, local_team()).wait();
      SLOG_GPU("Available number of GPUs on this node ", num_gpus_on_node, ". Detected in ", gpu_startup_duration, " s\n");
      if (!local_team().rank_me()) LOG("Initialized on node in ", msm_local.to_string(), "\n");
      if (!rank_me()) LOG("Initialized globally in ", msm_global.to_string(), "\n");
      // SLOG_GPU("Rank 0 is using GPU ", gpu_utils::get_gpu_device_name(), " on node 0, with ", get_size_str(gpu_avail_mem),
      //         " available memory (", get_size_str(gpu_avail_mem / local_team().rank_n()), " per rank). Detected in ",
      //         gpu_startup_duration, " s\n");
      SLOG_GPU(gpu_utils::get_gpu_device_descriptions());
      LOG("Using GPU device: ", gpu_utils::get_gpu_device_name(), " - ", gpu_utils::get_gpu_uuid(), "\n");
      log_gpu_uuid();
      barrier(local_team());
    } else {
      SDIE("No GPUs available - this build requires GPUs");
    }
  }
}

void tear_down_devices() {
  auto &gpu_team = get_gpu_team();
  gpu_team.destroy();
}

void log_gpu_uuid() {
  // Log the UUIDs for the GPUs used on the first node to rank0
  log_local("GPU UUID", gpu_utils::get_gpu_uuid());
}
