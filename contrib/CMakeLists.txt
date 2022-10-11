# MHM2 contrib

#
# The intent here is to include unit tests, like endian correctness for plain
# old C/C++ and not depend on upcxx / upcxx-run to execute them, so make test
# can work on all platforms including multi node without calling upcxx::init()
# and upcxx::finalize()
#

find_package(ZLIB)
add_executable(interleave_fastq interleave_fastq.c)
target_link_libraries(interleave_fastq ZLIB::ZLIB)
install(TARGETS interleave_fastq DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(PROGRAMS trim_fastq.sh
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
