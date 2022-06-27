set(SRCS parse_and_pack.cpp gpu_hash_table.cpp hashutil.cpp partitioned_counter.cpp)

if(MHM2_USE_QF)
  add_definitions(-DUSE_GQF)
  set(SRCS ${SRCS} gqf.cpp)
endif()

add_library(KCOUNT_GPU_LIBRARY_obj OBJECT ${SRCS})
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13)
  target_link_libraries(KCOUNT_GPU_LIBRARY_obj INTERFACE)
  # unsuccessful attempt to suppress warning about recursion
  target_link_options(
    KCOUNT_GPU_LIBRARY_obj PRIVATE $<DEVICE_LINK:--nvlink-options
    -suppress-stack-size-warning>)
endif()

foreach(_SRC ${SRCS})
  set_source_files_properties(${_SRC} PROPERTIES LANGUAGE CUDA LINKER_LANGUAGE CUDA)
endforeach()

option(KCOUNT_GPU_SHARED "kcount-GPU shared library" OFF)
option(KCOUNT_GPU_STATIC "kcount-GPU static library" ON)

if(KCOUNT_GPU_SHARED)
  add_library(KCOUNT_GPU_LIBRARY_shared SHARED
              $<TARGET_OBJECTS:KCOUNT_GPU_LIBRARY_obj>)
  set_property(TARGET KCOUNT_GPU_LIBRARY_shared
               PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  target_link_libraries(KCOUNT_GPU_LIBRARY_shared INTERFACE)
  install(TARGETS KCOUNT_GPU_LIBRARY_shared LIBRARY DESTINATION lib)
endif()
if(KCOUNT_GPU_STATIC)
  add_library(KCOUNT_GPU_LIBRARY_static STATIC
              $<TARGET_OBJECTS:KCOUNT_GPU_LIBRARY_obj>)
  target_link_libraries(KCOUNT_GPU_LIBRARY_static INTERFACE)
  # unsuccessful attempt to suppress warning about recursion
  target_link_options(
    KCOUNT_GPU_LIBRARY_static PRIVATE $<DEVICE_LINK:--nvlink-options
    -suppress-stack-size-warning>)
  set_property(TARGET KCOUNT_GPU_LIBRARY_static
               PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  install(TARGETS KCOUNT_GPU_LIBRARY_static ARCHIVE DESTINATION lib)
endif()
