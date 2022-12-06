set(LOCALASSM_TARGET_OBJECTS)

foreach(tgt localassm localassm_core)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND LOCALASSM_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
endforeach()

if(ENABLE_CUDA)
  set(tgt localassm_gpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND LOCALASSM_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
  set(LOCALASSM_TARGET_OBJECTS
      ${LOCALASSM_TARGET_OBJECTS}
      PARENT_SCOPE)

elseif(ENABLE_HIP)
  set(tgt localassm_gpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND LOCALASSM_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
  set(LOCALASSM_TARGET_OBJECTS
      ${LOCALASSM_TARGET_OBJECTS}
      PARENT_SCOPE)

else()
  set(tgt localassm_cpu)
  add_library(${tgt} OBJECT ${tgt}.cpp)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    target_link_libraries(${tgt} ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
  endif()
  list(APPEND LOCALASSM_TARGET_OBJECTS $<TARGET_OBJECTS:${tgt}>)
  set(LOCALASSM_TARGET_OBJECTS
      ${LOCALASSM_TARGET_OBJECTS}
      PARENT_SCOPE)
endif()

add_library(LOCALASSM_LIBRARY ${LOCALASSM_TARGET_OBJECTS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_CUDA)
  set_property(TARGET LOCALASSM_LIBRARY PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS
                                                 OFF)
  include_directories("localassm-gpu")
  add_subdirectory(localassm-gpu)
  set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} LOCALASSM_GPU_LIBRARY_static)
endif()

if(ENABLE_HIP)
  include_directories("localassm-gpu")
  add_subdirectory(localassm-gpu)
  set(MHM2_LINK_LIBRARIES ${MHM2_LINK_LIBRARIES} LOCALASSM_GPU_LIBRARY_static)
endif()

target_link_libraries(LOCALASSM_LIBRARY ${MHM2_LINK_LIBRARIES}
                      ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES})
