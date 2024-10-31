set(SRCS kokkos_pnp.cpp)

add_library(KCOUNT_KOKKOS_LIBRARY_obj OBJECT ${SRCS})
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13)
    target_link_libraries(KCOUNT_KOKKOS_LIBRARY_obj INTERFACE)
endif()

add_library(KCOUNT_KOKKOS_LIBRARY_static STATIC
                $<TARGET_OBJECTS:KCOUNT_KOKKOS_LIBRARY_obj>)
target_link_libraries(KCOUNT_KOKKOS_LIBRARY_static INTERFACE Kokkos::kokkos)
install(TARGETS KCOUNT_KOKKOS_LIBRARY_static ARCHIVE DESTINATION lib)

target_link_libraries(KCOUNT_KOKKOS_LIBRARY_obj Kokkos::kokkos)
