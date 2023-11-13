
option(LOCALASSM_GPU_LIBRARY_shared "LOCALASSM_GPU shared library" OFF)
option(LOCALASSM_GPU_LIBRARY_static "LOCALASSM_GPU static library" ON)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB HEADER_FILES *.hpp)
add_library(GPU_LOC_obj OBJECT ${SOURCE_FILES})

if(ENABLE_CUDA)

  message(STATUS "building LocalAssembly using CUDA")

  set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA
                                                         LINKER_LANGUAGE CUDA)

  if(LOCALASSM_GPU_LIBRARY_shared)
    add_library(LOCALASSM_GPU_LIBRARY_shared SHARED
                $<TARGET_OBJECTS:GPU_LOC_obj>)
    set_property(TARGET LOCALASSM_GPU_LIBRARY_shared
                 PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    install(TARGETS LOCALASSM_GPU_LIBRARY_shared LIBRARY DESTINATION lib)
  endif()
  if(LOCALASSM_GPU_LIBRARY_static)
    add_library(LOCALASSM_GPU_LIBRARY_static $<TARGET_OBJECTS:GPU_LOC_obj>)
    set_property(TARGET LOCALASSM_GPU_LIBRARY_static
                 PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    install(TARGETS LOCALASSM_GPU_LIBRARY_static ARCHIVE DESTINATION lib)
  endif()

elseif(ENABLE_HIP)

  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  message(STATUS "building LocalAssembly using HIP")

  set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE HIP
                                                         LINKER_LANGUAGE HIP)

  if(LOCALASSM_GPU_LIBRARY_shared)
    add_library(LOCALASSM_GPU_LIBRARY_shared SHARED
                $<TARGET_OBJECTS:GPU_LOC_obj>)
    # set_property(TARGET LOCALASSM_GPU_LIBRARY_shared PROPERTY
    # CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    install(TARGETS LOCALASSM_GPU_LIBRARY_shared LIBRARY DESTINATION lib)
  endif()
  if(LOCALASSM_GPU_LIBRARY_static)
    add_library(LOCALASSM_GPU_LIBRARY_static $<TARGET_OBJECTS:GPU_LOC_obj>)
    # set_property(TARGET LOCALASSM_GPU_LIBRARY_static PROPERTY
    # CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    install(TARGETS LOCALASSM_GPU_LIBRARY_static ARCHIVE DESTINATION lib)
  endif()

else()
  message(
    FATAL_ERROR "Trying to build LocalAssembly-GPU but CUDA/HIP are not enabled"
  )
endif()
