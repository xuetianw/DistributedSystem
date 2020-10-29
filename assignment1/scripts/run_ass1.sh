#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=ass1_res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100
#SBATCH --cpus-per-task=4


srun python /scratch/assignment1/test_scripts/submission_validator.pyc --tarPath=/home/xuetianw/assignment1.tar.gz