#!/bin/bash
#SBATCH -A swat-ml
#SBATCH --partition=swat-ml
#SBATCH --qos=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=180G
#SBATCH --time=140:00:00
#SBATCH --comment=cv1

module load Anaconda3/5.3.0
module load cuDNN/7.6.4.38-gcccuda-2019b
source activate yimin_env

python sunspot_cv1.py

