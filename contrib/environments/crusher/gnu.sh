module rm PrgEnv-cray
module load PrgEnv-gnu

module load git
module load cmake

module rm xl
module load gcc
#module load cuda
#module load upcxx
module use /gpfs/alpine/csc296/world-shared/crusher/modulefiles
module load upcxx
module list
which upcxx

#module load python
module load cray-python

#export GASNET_ODP_VERBOSE=0 # disable warnings
#export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which mpicc) -DCMAKE_CXX_COMPILER=$(which mpicxx) -DCMAKE_CUDA_COMPILER=$(which nvcc)"
export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=Off"
export MHM2_BUILD_THREADS=8

# salloc/sbatch with: -A BIF115_crusher --ntasks-per-node=64 
