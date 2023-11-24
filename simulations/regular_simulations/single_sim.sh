#!/bin/bash

#Directory where the data will be saved
dirname=simulation-data
echo dirname: $dirname
mkdir $dirname
echo "created directory for simulation data"

#Change the parameter and values for this run here
#The name of the parameter must be the same as the argument name in the main.py file 
#Run: ">python3 main.py --help" to see all parameter names
 
python3 parallel_main.py $dirname --n_runs 500 --n_cores 8 --n_tau 1500 

