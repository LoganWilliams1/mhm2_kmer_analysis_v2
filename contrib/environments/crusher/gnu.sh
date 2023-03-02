module rm PrgEnv-cray
module load PrgEnv-gnu

module load git
module load cmake

module rm xl
module load gcc
# need to load > 5.1 to avoid the background HIP thread
module load rocm/5.2.0
# this causes a remap failure when using GPUs
#module load craype-hugepages2M
module use /gpfs/alpine/csc296/world-shared/crusher/modulefiles
module rm upcxx
module load upcxx/nightly; export GASNET_OFI_RECEIVE_BUFF_SIZE=recv
#module list
which upcxx

module load cray-python

export MHM2_CMAKE_EXTRAS="-DCMAKE_C_COMPILER=$(which cc) -DCMAKE_CXX_COMPILER=$(which CC) -DENABLE_CUDA=Off -DENABLE_HIP=On"
export MHM2_BUILD_THREADS=8

# salloc/sbatch with: -A BIF115_crusher --ntasks-per-node=64 
