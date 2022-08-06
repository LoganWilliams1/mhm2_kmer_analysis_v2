module load PrgEnv-cray
module load cmake
#module swap PrgEnv-cray/8.2.0
#module swap cce/13.0.1

module use /global/common/software/m2878/perlmutter/modulefiles
# module load upcxx
module load upcxx/bleeding-edge ; export GASNET_OFI_RECEIVE_BUFF_SIZE=single

module list
which cc
which CC
which g++
which gcc
which upcxx

CC --version
upcxx --version

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=OFF"
