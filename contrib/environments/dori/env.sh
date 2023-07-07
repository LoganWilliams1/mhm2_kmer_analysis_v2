which mpirun 2>/dev/null || export PATH=$PATH:/clusterfs/jgi/scratch/gentech/genome_analysis/rsegan/Software/openmpi-4.1.5/bin/
which cmake 2>/dev/null || export PATH=$PATH:/clusterfs/jgi/scratch/gentech/genome_analysis/rsegan/Software/cmake-3.26.4/bin/
which upcxx 2>/dev/null || export PATH=$PATH:/clusterfs/jgi/scratch/gentech/genome_analysis/rsegan/Software/upcxx-2023.3.0-ibv/bin

export GASNET_PHYSMEM_MAX='100 GB'
export GASNET_PHYSMEM_PROBE=NO
export OMPI_MCA_btl=tcp,vader,self
export GASNET_USE_ODP=0
