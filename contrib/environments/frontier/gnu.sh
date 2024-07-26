module rm PrgEnv-cray
module load PrgEnv-gnu
module load git
module load cmake

module rm xl
module load gcc-native

module load rocm/6.0.0
#module load rocm
module rm upcxx
module load ums ums014 upcxx
export GASNET_OFI_RECEIVE_BUFF_SIZE=recv
export FI_MR_CACHE_MONITOR=memhooks
export FI_CXI_RX_MATCH_MODE=software

module rm craype-hugepages2M # FIXME when slingshot is fixed

which upcxx
upcxx --version

module load cray-python

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=Off -DENABLE_HIP=On -DUPCXX_UTILS_IO_NO_THREAD=1 -DCMAKE_CXX_FLAGS=-Wno-deprecated-declarations"
export MHM2_BUILD_THREADS=8

