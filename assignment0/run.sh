#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100


srun /home/$USER/producer_consumer --nItems 100000 --nProducers 2 --nConsumers 4 --bufferSize 100000