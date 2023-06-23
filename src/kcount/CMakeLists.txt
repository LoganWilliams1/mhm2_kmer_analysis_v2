set(KCOUNT_TARGET_FILES)

foreach(KMER_LENGTH ${MHM2_KMER_LENGTHS})
  foreach(TEMPLATE "kcount")
    set(TEMPLATE_FILE "${TEMPLATE}-extern-template-${KMER_LENGTH}")
    string(TOUPPER ${TEMPLATE} TEMPLATE_UPPER)
    configure_file("${CMAKE_SOURCE_DIR}/src/extern_template.in.cpp"
                   "${TEMPLATE_FILE}.cpp" @ONLY)
    list(APPEND KCOUNT_TARGET_FILES "${TEMPLATE_FILE}.cpp")
  endforeach()
endforeach()

set(KCOUNT_TARGET_OBJECTS)

foreach(tgt kmer_dht)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND KCOUNT_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
endforeach()

if(ENABLE_CUDA)
  set(tgt kcount_gpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND KCOUNT_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
  set(KCOUNT_TARGET_OBJECTS
      ${KCOUNT_TARGET_OBJECTS}
      PARENT_SCOPE)

elseif(ENABLE_HIP)
  set(tgt kcount_gpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  target_compile_options(${tgt} PUBLIC $<$<COMPILE_LANGUAGE:HIP>:-fgpu-rdc>)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND KCOUNT_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
  set(KCOUNT_TARGET_OBJECTS
      ${KCOUNT_TARGET_OBJECTS}
      PARENT_SCOPE)

else()
  set(tgt kcount_cpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND KCOUNT_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
  set(KCOUNT_TARGET_OBJECTS
      ${KCOUNT_TARGET_OBJECTS}
      PARENT_SCOPE)
endif()

add_library(KCOUNT_LIBRARY ${KCOUNT_TARGET_FILES} ${KCOUNT_TARGET_OBJECTS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_CUDA)
  add_definitions(-DUSE_TCF)
  set_property(TARGET KCOUNT_LIBRARY PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
  include_directories("kcount-gpu")
  add_subdirectory(kcount-gpu)
  set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} KCOUNT_GPU_LIBRARY_static)

elseif(ENABLE_HIP)
  # set_property(TARGET KCOUNT_LIBRARY PROPERTY HIP_RESOLVE_DEVICE_SYMBOLS OFF)
  include_directories("kcount-gpu")
  add_subdirectory(kcount-gpu)
  set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} KCOUNT_GPU_LIBRARY_static)

endif()

target_link_libraries(KCOUNT_LIBRARY ${MHM2_LINK_LIBRARIES} ${UPCXX_LIBRARIES}
                      ${UPCXX_UTILS_LIBRARIES})
