module load PrgEnv-gnu
#module load cpe-cuda
module load gcc/9.3.0
module load cmake
module load cuda
module load cudatoolkit

export PATH=/global/common/software/m2878/shasta2105/upcxx/TESTING-PrgEnv-gnu/bin:$PATH

export FI_PROVIDER='verbs;ofi_rxm'
export UCX_TLS=dc
export UCX_DC_MLX5_NUM_DCI=16

export UPCXX_NETWORK=ofi

module list
which cc
which CC
which g++
which gcc
which nvcc
which upcxx

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DCMAKE_CUDA_COMPILER=$(which nvcc)"
