#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100
#SBATCH --cpus-per-task=4

python /scratch/assignment1/test_scripts/pi_calculation_tester.pyc --execPath=/home/xuetianw/pi_calculation_parallel

#ls /scratch/assignment1/test_scripts/*tester.pyc

srun cp /scratch/assignment1/test_scripts/pi_calculation_tester.pyc ./home/xuetianw/
