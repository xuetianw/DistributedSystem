#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G


python /scratch/assignment1/test_scripts/triangle_counting_tester.pyc --execPath=/home/xuetianw/triangle_counting_parallel


#srun cp /scratch/assignment1/test_scripts/triangle_counting_tester.pyc ./

#srun cp /scratch/input_graphs/*.cs* ./sfuhome/CMPT431