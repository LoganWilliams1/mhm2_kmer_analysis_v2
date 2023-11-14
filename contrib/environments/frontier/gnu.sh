module rm PrgEnv-cray
module load PrgEnv-gnu

module load git
module load cmake

module rm xl
module load gcc
module load rocm/5.2.0
module rm upcxx
module load ums ums014 upcxx
export GASNET_OFI_RECEIVE_BUFF_SIZE=recv
export FI_MR_CACHE_MONITOR=memhooks
export FI_CXI_RX_MATCH_MODE=software
module rm craype-hugepages2M # FIXME when slingshot is fixed

which upcxx
upcxx --version

module load cray-python

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=Off -DENABLE_HIP=On -DUPCXX_UTILS_IO_NO_THREAD=1"
export MHM2_BUILD_THREADS=8

