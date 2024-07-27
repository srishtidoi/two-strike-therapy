#!/bin/bash

#SBATCH -D <home directory of the user>
#SBATCH --job-name=<job name>
#SBATCH --partition=nodes 
#SBATCH --ntasks 50 <number of tasks to run simultaneously. 1 task = 1 run>
#SBATCH --array=1-100%15 <number of random seeds%number of random seeds to run simultaneously>
#SBATCH --mem=10G <memory allocation>
#SBATCH --time=10:00:00 <time of the simulation, dd-hh:mm:ss>
#SBATCH -e results/%x_%A_%a.e <error file>
#SBATCH -o results/%x_%A_%a.o <output file>

#Edit the lines above to either remove or replace all text in <...>. 
#They are comments to explain the meaning of each variable

#Enable modules
source /opt/flight/etc/setup.sh
flight env activate gridware
module purge
module add gnu
module add python/3.9.7

#Enter number of cores, data directory and parameter to edit
#The number of cores must be the same as --ntasks above
#The name of the param must be the same as the argument name for the main script
#To see all argument names run '>> python3 seed_main.py --help'

dirname=../../simulation_data/test
n_cores=50
param=cost
value=0

#Creating the data directory
mkdir $dirname
echo "created directory for all simulation data"
mkdir $dirname/$SLURM_ARRAY_TASK_ID
echo "created subfolder for the current simulation data"

#Run one simulation. To see all arguments, run '>> python3 main.py --help'
python3 seed_main.py $dirname/$SLURM_ARRAY_TASK_ID --seed $SLURM_ARRAY_TASK_ID --n_cores $n_cores --$param $value


