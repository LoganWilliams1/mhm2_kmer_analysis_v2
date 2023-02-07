module load gcc
module load cmake/3.18.2
module rm openmpi
module load slurm

# UPCXX built with: BUILD_THREADS=16 ./upcxx-utils/contrib/install_upcxx.sh $SCRATCH/install-upcxx-2022.9 --enable-ibv --with-pmi-version=2 --with-default-net=ibv --enable-pmi --with-ibv-spawner=pmi --with-ibv-physmem-max=1/2 --disable-mpi-compat --with-ibv-max-hcas=3

export PMI_HOME=/cm/shared/apps/slurm/current
export PMIRUN_CMD="srun --mpi=pmi2 %V -n %N -- %C"

export GASNET_PHYSMEM_PROBE=0
export GASNET_PHYSMEM_MAX=1/2
export GASNET_MAX_SEGSIZE=0.45/H

export PATH=$SCRATCH/install-upcxx-2022.9/bin:$PATH
export MHM2_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) "
