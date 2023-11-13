module rm PrgEnv-cray
module load PrgEnv-gnu

module load git
module load cmake

module rm xl
module load gcc
module load rocm/5.2.0
# this causes a remap failure when using GPUs
#module load craype-hugepages2M
module use /lustre/orion/csc296/world-shared/frontier/modulefiles
module rm upcxx
module load upcxx/2023.3.0
export GASNET_OFI_RECEIVE_BUFF_SIZE=recv
export FI_OFI_RXM_RX_SIZE=8192
export FI_CXI_DEFAULT_CQ_SIZE=13107200
export FI_MR_CACHE_MONITOR=memhooks
export FI_CXI_RX_MATCH_MODE=software
export FI_CXI_REQ_BUF_MIN_POSTED=10
export FI_CXI_REQ_BUF_SIZE=25165824
module rm craype-hugepages2M # FIXME when slingshot is fixed

#module list
which upcxx
upcxx --version

module load cray-python

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=Off -DENABLE_HIP=On -DUPCXX_UTILS_IO_NO_THREAD=1"
export MHM2_BUILD_THREADS=8

