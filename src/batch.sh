#!/bin/bash
#SBATCH --partition gpu
#SBATCH --output=batch_%j.out         ### File in which to store job output
#SBATCH --error=batch_%j.err          ### File in which to store job error messages
#SBATCH --cpus-per-task 2
#SBATCH --mem 30G
#SBATCH --time 48:00:00
#SBATCH --job-name bigWigAvgOutBatch_
#SBATCH --gres=gpu:p100:1

#/usr/bin/time -v conda install python=3.10
/usr/bin/time -v python /home/groups/CEDAR/franksto/projects/TAD_project/src/bigWigAvgOut_Batch.py