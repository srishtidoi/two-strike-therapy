#!/bin/bash

#SBATCH -D <home directory of the user>
#SBATCH --job-name=<job name>
#SBATCH --partition=nodes
#SBATCH --ntasks 50 <number of tasks to run simultaneously. 1 task = 1 run of a simulation>
#SBATCH --array=1-121%15 <number of parameter values%number to run simultaneously>
#SBATCH --mem=15G <memory allocation>
#SBATCH --time=10:00:00 <time allocation, dd-hh:mm:ss>
#SBATCH -e results/%x_%A_%a.e <error file>
#SBATCH -o results/%x_%A_%a.o <output file>

#Enable modules 
source /opt/flight/etc/setup.sh
flight env activate gridware
module purge
module add gnu
module add python/3.9.7


#Set job array variables
config=grid_config.txt
value1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
value2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
echo param values: $value1 x $value2

#set number of cores, loop parameter and data directory name
#the name of the param must be the same as the argument name for the main script
#to see all argument names run '>> python3 parallel_main.py --help'
#the number of cores must be the same as --ntasks above

param1=cost
param2=n_tau
dirname=../../simulation_data/test
n_cores=50

#creating data directory
mkdir $dirname
echo "created directory for all simulation data"
cp $config $dirname/config.txt
echo "copied config file to data directory"
mkdir $dirname/$SLURM_ARRAY_TASK_ID
echo "created subfolder for the current simulation data"

#change the parameter and values for simulations here
#the name of the parameter must be the same as the argument name in the main script
#run: ">> python3 parallel_main.py --help" to see all parameter names

python3 main.py $dirname/$SLURM_ARRAY_TASK_ID --$param1 $value1 --$param2 $value2 --n_runs 500 --wt 1000000 --n_max 1500000 --r1 100 --r2 100 --n_cores $n_cores --T 500 --K 1000200


