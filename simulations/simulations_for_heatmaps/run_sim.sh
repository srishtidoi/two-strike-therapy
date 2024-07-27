#!/bin/bash

#SBATCH -D <home directory of the user>
#SBATCH --job-name=<job name>
#SBATCH --partition=nodes
#SBATCH --ntasks 50 <number of tasks to run simultaneously. 1 task = 1 run of a simulation>
#SBATCH --array=1-11 <number of parameter values>
#SBATCH --mem=10G <memory allocation>
#SBATCH --time=01:00:00 <time allocation, dd-hh:mm:ss>
#SBATCH -e results/%x_%A_%a.e <error file>
#SBATCH -o results/%x_%A_%a.o <output file>

#Edit the lines above to either remove or replace all text in <...>

#Enable modules
source /opt/flight/etc/setup.sh
flight env activate gridware
module purge
module add gnu
module add python/3.9.7


#Set job array variables. The config file lists all parameter values to run the simulations for
config=config.txt
value=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
echo param value: $value

#Set number of cores, loop param and data directory nam
#The name of the param must be the same as the argument name for the main script
#To see all param names run '>python3 main.py --help'
#The number of cores must be the same as --ntasks above

param=n_tau
dirname=../../simulation_data/test
n_cores=50

#creating data directory
mkdir $dirname
echo "created directory for all simulation data"
cp $config $dirname/config.txt
echo "copied config file to data directory"
mkdir $dirname/$SLURM_ARRAY_TASK_ID
echo "created subfolder for the current simulation data"

#Change the fixed parameters for the simulations here (if needed)
#The name of the parameter must be the same as the argument name in the main script
#Run: ">python3 main.py --help" to see all parameter names

python3 main.py $dirname/$SLURM_ARRAY_TASK_ID --$param $value --n_runs 500 --n_cores $n_cores


