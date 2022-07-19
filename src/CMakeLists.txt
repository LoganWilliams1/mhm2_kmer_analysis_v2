# HipMer v 2.0, Copyright (c) 2020, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any
# required approvals from the U.S. Dept. of Energy).  All rights reserved."

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# (1) Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.

# (2) Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# (3) Neither the name of the University of California, Lawrence Berkeley
# National Laboratory, U.S. Dept. of Energy nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# You are under no obligation whatsoever to provide any bug fixes, patches, or
# upgrades to the features, functionality or performance of the source code
# ("Enhancements") to anyone; however, if you choose to make your Enhancements
# available either publicly, or directly to Lawrence Berkeley National
# Laboratory, without imposing a separate written license agreement for such
# Enhancements, then you hereby grant the following license: a  non-exclusive,
# royalty-free perpetual license to install, use, modify, prepare derivative
# works, incorporate into other computer software, distribute, and sublicense
# such enhancements or derivative works thereof, in binary and source code form.

cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
  message(
    STATUS
      "MHM2 is using the UPCXX::upcxx library interface: ${UPCXX_LIBRARIES}")
  find_package(UPCXX)
else()
  find_program(UPCXX_EXEC upcxx)
  set(CMAKE_CXX_COMPILER ${UPCXX_EXEC})
  message(
    STATUS
      "MHM2 is using upcxx directly as the UPCXX::upcxx library interface is not available in this low version of cmake: ${CMAKE_VERSION}"
  )
  unset(UPCXX_LIBRARIES)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -Wno-catch-value"
    # "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare
    # -Wno-unused-function -Wno-deprecated-gpu-targets
    # -Wno-unused-command-line-argument"
)
find_package(Threads REQUIRED)

find_package(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
message(
  "Adding ZLIB for -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)"
)

if(ENABLE_CUDA)
  message(STATUS "Building for GPUs with CUDA")

  set(CMAKE_CUDA_STANDARD
      14
      CACHE STRING "")
  set(CMAKE_CUDA_STANDARD_REQUIRED
      ON
      CACHE BOOL "")
  set(CMAKE_CUDA_EXTENSIONS
      OFF
      CACHE BOOL "")

  find_package(CUDA REQUIRED)

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    cmake_policy(SET CMP0074 NEW) # use the CUDA_ROOT variable
  endif()

  enable_language(CUDA)
  set(MHM2_MIN_CUDA_ARCH 70)

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
      include(FindCUDA/select_compute_arch)
      cuda_detect_installed_gpus(INSTALLED_GPU_CCS_1)
      string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
      string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
      string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
      set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
      message(
        STATUS "Autodetect CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}"
      )
    endif()
    message(STATUS "Using CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}")
  else()
    find_package(CUDA REQUIRED)
    cuda_select_nvcc_arch_flags(ARCH_FLAGS_LIST)
    if(ARCH_FLAGS_LIST)
      string(REPLACE ";" " " ARCH_FLAGS "${ARCH_FLAGS_LIST}")
      message("Found ARCH_FLAGS_LIST: ${ARCH_FLAGS_LIST}")
      message("Found CUDA_GPU_DETECT_OUTPUT: ${CUDA_GPU_DETECT_OUTPUT}")
      if(NOT CMAKE_CUDA_ARCHITECTURES)
        string(REPLACE "." "" CMAKE_CUDA_ARCHITECTURES
                       "${CUDA_GPU_DETECT_OUTPUT}")
      endif()
    else()
      foreach(_CMAKE_CUDA_COMPUTE_CAPABILITY 70 80)
        set(ARCH_FLAGS
            "${ARCH_FLAGS} -gencode arch=compute_${_CMAKE_CUDA_COMPUTE_CAPABILITY},code=sm_${_CMAKE_CUDA_COMPUTE_CAPABILITY}"
        )
      endforeach()
      if(NOT CMAKE_CUDA_COMPUTE_CAPABILITY)
        set(CMAKE_CUDA_COMPUTE_CAPABILITY
            70
            CACHE STRING "")
      endif()
      message(
        WARNING
          "Could not auto-detect the GPU arch flags building for all gpu compute capabilities 3 - 7 and PTX for ${CMAKE_CUDA_COMPUTE_CAPABILITY}"
      )
      set(ARCH_FLAGS
          "${ARCH_FLAGS} -gencode arch=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY},code=sm_${CMAKE_CUDA_COMPUTE_CAPABILITY}"
      )
      set(ARCH_FLAGS
          "${ARCH_FLAGS} -gencode arch=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY},code=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY}"
      )
    endif()
    message(STATUS "ARCH_FLAGS=${ARCH_FLAGS}")
  endif()

  foreach(test_cuda_arch ${CMAKE_CUDA_ARCHITECTURES})
    if(${test_cuda_arch} LESS ${MHM2_MIN_CUDA_ARCH})
      message(
        FATAL_ERROR
          "CUDA architecture ${test_cuda_arch} is incompatible with the minimum ${MHM2_MIN_CUDA_ARCH}.  Try compiling without cuda: -DENABLE_CUDA=OFF"
      )
    endif()
  endforeach()

  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
  set(CMAKE_CUDA_PTX_COMPILATION ON)

  set(CMAKE_CUDA_FLAGS
      "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall ${ARCH_FLAGS} -Wno-deprecated-gpu-targets"
  )

  message(STATUS "CMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}")

else()
  message(STATUS "Building for CPUs")
endif()

# a bit of jiggery-pokery to get a list of all compile time definitions in the
# code
get_directory_property(dir_defs COMPILE_DEFINITIONS)
get_directory_property(dir_compile_opts COMPILE_OPTIONS)
string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" flags_var)
set(all_defs "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ENABLE_CUDA=${ENABLE_CUDA} COMPILE_OPTIONS='${dir_compile_opts}' CMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}' ${flags_var}='${${flags_var}}'")
message("Definitions:")
foreach(d ${dir_defs})
  message("  " ${d})
  set(all_defs "${all_defs} ${d}")
endforeach()

add_definitions(-DALL_DEFNS=\"${all_defs}\")

include_directories("ssw")
add_subdirectory(ssw)
include_directories("kcount")
add_subdirectory(kcount)
include_directories("klign")
add_subdirectory(klign)
include_directories("localassm")
add_subdirectory(localassm)
include_directories("cgraph")
add_subdirectory(cgraph)

set(MHM2_TARGET_OBJECTS)
add_library(hash_funcs OBJECT hash_funcs.c)
list(APPEND MHM2_TARGET_OBJECTS $<TARGET_OBJECTS:hash_funcs>)

# build all the extern template souce code instantiation explosion

foreach(KMER_LENGTH ${MHM2_KMER_LENGTHS})
  foreach(TEMPLATE "contigging" "scaffolding" "kmer")
    set(TEMPLATE_FILE "${TEMPLATE}-extern-template-${KMER_LENGTH}")
    string(TOUPPER ${TEMPLATE} TEMPLATE_UPPER)
    configure_file("extern_template.in.cpp" "${TEMPLATE_FILE}.cpp" @ONLY)
    list(APPEND EXTERN_TEMPLATE_FILES ${TEMPLATE_FILE})
  endforeach()
endforeach()

foreach(
  tgt
  fastq
  alignments
  contigs
  packed_reads
  utils
  options
  merge_reads
  dbjg_traversal
  histogrammer
  aln_depths
  shuffle_reads
  gasnet_stats
  post_assembly
  ${EXTERN_TEMPLATE_FILES})
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND MHM2_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
endforeach()

if(ENABLE_CUDA)
  set(tgt devices_gpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND MHM2_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
else()
  set(tgt devices_cpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND MHM2_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
endif()

set(MHM2_TARGET_OBJECTS
    ${MHM2_TARGET_OBJECTS}
    PARENT_SCOPE)

add_executable(mhm2 main.cpp ${MHM2_TARGET_OBJECTS})

set(MHM2_LINK_LIBRARIES
    Threads::Threads
    SSW_LIBRARY
    KCOUNT_LIBRARY
    KLIGN_LIBRARY
    LOCALASSM_LIBRARY
    CGRAPH_LIBRARY
    ${ZLIB_LIBRARIES}
    MHM2_VERSION_LIB)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_CUDA)
  set_property(TARGET mhm2 PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
  include_directories("gpu-utils")
  add_subdirectory(gpu-utils)
  set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} GPU_UTILS_LIBRARY_static)
endif()

target_link_libraries(mhm2 ${MHM2_LINK_LIBRARIES} ${UPCXX_LIBRARIES}
                      ${UPCXX_UTILS_LIBRARIES})

set(MHM2_LINK_LIBRARIES
    ${MHM2_LINK_LIBRARIES}
    PARENT_SCOPE)

install(TARGETS mhm2 DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(MHM2PY_RUNTIME_ENV "GASNET_BACKTRACE='1',GASNET_BACKTRACE_MT='1', ")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebug")
  set(MHM2PY_RUNTIME_ENV "GASNET_BACKTRACE='1',GASNET_BACKTRACE_MT='1', ")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  set(MHM2PY_RUNTIME_ENV "GASNET_BACKTRACE='1',GASNET_BACKTRACE_MT='1', ")
endif()
configure_file(mhm2.py ${CMAKE_CURRENT_BINARY_DIR}/mhm2.py)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mhm2.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(PROGRAMS mhm2_parse_run_log.pl DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# special spawning script for Perlmutter
install(
  PROGRAMS mps-wrapper-perlmutter.sh
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  RENAME mhm2-mps-wrapper-perlmutter.sh)
