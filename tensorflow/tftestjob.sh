#!/bin/bash

#Tensor flow test job batch file

#submit using the command
# sbatch tftestjob.sh
# monitor queue using queue

#SBATCH --account=bdshe01  # Run job under project <project>
#SBATCH --time=1:0:0         # Run for a max of 1 hour

# Node resources:
# (choose between 1-4 gpus per node)

#SBATCH --partition=gpu    # Choose either "gpu" or "infer" node type
#SBATCH --nodes=1          # Resources from a single node
#SBATCH --gres=gpu:1       # One GPU per node (plus 25% of node CPU and RAM per GPU)

#SBATCH --mem=16G
#SBATCH --mail-user=m.griffiths@sheffield.ac.uk

#SBATCH --gpus=1

module load slurm
module load Anaconda3/2020.02
source activate tensorflow

nvidia-smi
# Replace my_program with the name of the program you want to run
python tftest.py

echo test
