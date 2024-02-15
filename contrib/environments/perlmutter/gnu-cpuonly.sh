module rm PrgEnv-cray
module rm PrgEnv-intel
module rm PrgEnv-gnu
module load PrgEnv-gnu/8.5.0
module rm gpu
module load cmake
module rm gcc
module load gcc-native
module swap gcc-native/12.3
module remove darshan
module remove craype-hugepages2M # FIXME when slingshot is fixed

module use /global/common/software/m2878/perlmutter/modulefiles
module rm upcxx
module load upcxx/2023.9.0
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
which upcxx

CC --version
upcxx --version

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=OFF"

echo "Environment for perlmutter CPU"
env | grep '^GASNET\|^UPCXX\|^FI_'
