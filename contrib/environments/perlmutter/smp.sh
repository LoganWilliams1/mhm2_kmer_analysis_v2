module load PrgEnv-gnu
module load cmake
module load cpe-cuda
module load cudatoolkit
module swap gcc/11.2.0

module use /global/common/software/m2878/perlmutter/modulefiles
# module load upcxx
module load upcxx/bleeding-edge ; export GASNET_OFI_RECEIVE_BUFF_SIZE=single

module list
which cc
which CC
which g++
which gcc
which nvcc
which upcxx

CC --version
upcxx --version
nvcc --version

export UPCXX_NETWORK=smp
export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DCMAKE_CUDA_COMPILER=$(which nvcc)"
