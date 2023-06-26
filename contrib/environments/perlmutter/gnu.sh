module load PrgEnv-gnu
module load cmake
module load cpe-cuda
module load cudatoolkit
#module swap PrgEnv-gnu/8.2.0
module swap gcc/11.2.0
module remove darshan
#module load craype-hugepages2M
module rm craype-hugepages2M # FIXME when slingshot is fixed

module use /global/common/software/m2878/perlmutter/modulefiles
module rm upcxx
module load upcxx/2023.3.0
export GASNET_OFI_RECEIVE_BUFF_SIZE=recv ; export GASNET_OFI_NUM_RECEIVE_BUFFS=400 # FIXME when slingshot is fixed
export FI_OFI_RXM_RX_SIZE=8192 # FIXME when slingshot is fixed
#export FI_CXI_DEFAULT_CQ_SIZE=13107200 # FIXME when slingshot is fixed
export FI_CXI_DEFAULT_CQ_SIZE=$((13107200 / 5)) # FIXME when slingshot is fixed
export FI_MR_CACHE_MONITOR=memhooks # FIXME when slingshot is fixed
export FI_CXI_RX_MATCH_MODE=software # FIXME when slingshot is fixed
export FI_CXI_REQ_BUF_MIN_POSTED=10 # FIXME when slingshot is fixed
export FI_CXI_REQ_BUF_SIZE=25165824 # FIXME when slingshot is fixed

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

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DCMAKE_CUDA_COMPILER=$(which nvcc)"

echo "Environment for perlmutter GPU"
env | grep '^GASNET\|^UPCXX\|^FI_'
