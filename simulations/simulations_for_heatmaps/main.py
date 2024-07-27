import numpy as np
import time
import params
import functions as fn
import csv
from joblib import Parallel, delayed
from numpy.random import default_rng
rng = default_rng()

# simulation data directory given as a command line argument
sim_data_dir = params.data_dir

# execution time
start_time = time.time()
CPU_start_time = time.process_time()

# reading all parameters from params file and command line arguments and saving them to a csv file. See functions.py and params.py for description and default values of parameters.

n_runs = params.n_runs
n = params.n  
w_lambda = params.w_lambda 
mus_vector = params.mus_vector
cost_vector = params.cost_vector 
mutation_probs = params.mutation_probs
alpha_1 = mutation_probs[0]
alpha_2 = mutation_probs[1]
mutation_vector = params.mutation_vector 
therapy_vector = params.therapy_vector
t_res = params.t_res 
T = params.T 
N_tau = params.N_tau 
N_max = params.N_max
n_cores = params.n_cores 
pop_dist_init = np.array(list(params.pop_dist))
K = params.K

# save all parameter values to a csv file
fn.save_params(sim_data_dir+"/params.csv", "n_runs", "n", "w_lambda", "mus_vector", "cost_vector", "mutation_probs", "mutation_vector", "therapy_vector", "t_res", "T", "N_tau", "N_max", "n_cores", "pop_dist_init","K")

# For 2 strikes we have - (w, r1, r2, r12) for wild-type and the 3 resistant variants.

# Running the simulation            
results = Parallel(n_jobs=n_cores)(delayed(fn.one_run)(run) for run in range(n_runs))

# sorting simulation results
pop_history_allruns = sum([results[run][0] for run in range(n_runs)], [])
time_history_allruns = sum([results[run][1] for run in range(n_runs)], [])
outcomes = [results[run][2] for run in range(n_runs)]

# saving simulation results to files
extinction_prob = outcomes.count(0)/n_runs
progression_prob = outcomes.count(1)/n_runs
persistence_prob = outcomes.count(2)/n_runs
with open(sim_data_dir+'/results.txt','w') as f:
    print("n_runs, ", n_runs, file=f)
    print("extinction probability, ", extinction_prob, file=f)
    print("progression probability, ", progression_prob, file=f)
    print("persistence probability, ", persistence_prob, file=f)

if persistence_prob>0.1:
    print("Persistence probability = ", persistence_prob)
    print("Need more time!")

np.savetxt(sim_data_dir+"/lineage_data.csv", pop_history_allruns, delimiter=", ", fmt="% s")    
np.savetxt(sim_data_dir+"/time_data.csv", time_history_allruns, delimiter=", ", fmt="% s")


# calculate and print execution time
end_time = time.time()
CPU_end_time = time.process_time()
execution_time = (end_time-start_time)/60
CPU_time = (CPU_end_time-CPU_start_time)/60
print("Execution time: ", execution_time, " minutes (", execution_time/n_runs, " minutes per run)")
print("CPU time: ", CPU_time, " minutes (", CPU_time/n_runs, " minutes per run)")

