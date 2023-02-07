module load gcc
module load openmpi
module load cmake/3.18.2
module load slurm

if ! [[ "$PATH" =~ "/home/regan1/install/bin" ]]
then
    PATH="$PATH:$HOME/install/bin"
fi
export PATH

export PMI_HOME=/cm/shared/apps/slurm/current
export PMIRUN_CMD="srun --mpi=pmi2 %V -n %N -- %C"

export GASNET_PHYSMEM_PROBE=0
export GASNET_PHYSMEM_MAX=1/2
export GASNET_MAX_SEGSIZE=0.45/H

export MHM2_CMAKE_EXTRAS="-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which mpicxx)"
