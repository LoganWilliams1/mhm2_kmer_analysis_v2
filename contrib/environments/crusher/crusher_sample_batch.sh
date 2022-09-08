#!/bin/bash

#SBATCH -N 1
#SBATCH -A BIF115
#SBATCH --ntasks-per-node=64
#SBATCH --threads-per-core=2
#SBATCH --gpus-per-node=8
#SBATCH --gpu-bind=closest
#SBATCH -t 30


source ../contrib/environments/crusher/gnu.sh

export GASNET_OFI_RECEIVE_BUFF_SIZE=single
# srun -N1 -n 64 --threads-per-core=2 -m block --gpus-per-node=8 --gpu-bind=closest ./install/bin/mhm2 -r /gpfs/alpine/bif115/scratch/mgawan/arcticsynth/arctic_sample_1.fq -o /gpfs/alpine/bif115/scratch/mgawan/arcticsynth/release_crusher_debug_1 --checkpoint=false

# srun -N1 -n 64 -c2 -m block,nopack --gpus-per-node=8 --gpu-bind=closest ./install/bin/mhm2 -r /gpfs/alpine/bif115/scratch/mgawan/arcticsynth/arctic_sample_1.fq,/gpfs/alpine/bif115/scratch/mgawan/arcticsynth/arctic_sample_2.fq,/gpfs/alpine/bif115/scratch/mgawan/arcticsynth/arctic_sample_3.fq,/gpfs/alpine/bif115/scratch/mgawan/arcticsynth/arctic_sample_4.fq -o /gpfs/alpine/bif115/scratch/mgawan/arcticsynth/release_crusher_debug_4 --checkpoint=false

srun -N1 -n 64 -c2 -m block,nopack --gpus-per-node=8 --gpu-bind=closest ./install/bin/mhm2 -r /gpfs/alpine/bif115/scratch/mgawan/arcticsynth/arctic_sample_* -o /gpfs/alpine/bif115/scratch/mgawan/crusher_mhm2_testing_2/GPU_1_NODES/ --checkpoint=false