export UPCXX_NETWORK=smp
export MHM2_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CUDA_COMPILER=$(which nvcc)"
