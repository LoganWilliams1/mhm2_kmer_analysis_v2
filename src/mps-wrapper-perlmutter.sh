#!/bin/bash
# Example mps-wrapper.sh usage:
# > srun [srun args] mps-wrapper.sh [cmd] [cmd args]
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log
# Launch MPS from a single rank per node
if [ $SLURM_LOCALID -eq 0 ]; then
    CUDA_VISIBLE_DEVICES=$SLURM_JOB_GPUS nvidia-cuda-mps-control -d
fi
# Wait for MPS to start
sleep 5
# Run the command
if [ "$1" == "--" ]
then
  shift
fi

# select_cpu_device blocked by rank (not round robin) : # export CUDA_VISIBLE_DEVICES=$(( SLURM_LOCALID % 4 )) 
export CUDA_VISIBLE_DEVICES=$(( 4 * SLURM_LOCALID / ${SLURM_TASKS_PER_NODE%%\(*} )) 

"$@"
# Quit MPS control daemon before exiting
if [ $SLURM_LOCALID -eq 0 ]; then
    echo quit | nvidia-cuda-mps-control 
fi
