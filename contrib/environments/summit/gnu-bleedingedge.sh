module use /gpfs/alpine/world-shared/csc296/summit/modulefiles

module load git
module load cmake

module rm xl
module load gcc
module load cuda
module load upcxx
module load python

export GASNET_RBUF_COUNT=32676 # maximum allowed currently support 778 nodes "fully provisioned for simultaneous all-to-all"
export GASNET_AM_CREDITS_SLACK=0 # to disable warnings
export GASNET_AM_CREDITS_PP=1 # for when 1 hca per rank

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) -DCMAKE_CUDA_COMPILER=$(which nvcc)"
export MHM2_BUILD_THREADS=8
