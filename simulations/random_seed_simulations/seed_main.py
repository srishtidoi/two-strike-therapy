import numpy as np
import time
import seed_params as params
import seed_functions as fn
import random as rnd
from joblib import Parallel, delayed
from numpy.random import default_rng
rng = default_rng()

# simulation data directory given as a command line argument
sim_data_dir = params.data_dir

# execution time
start_time = time.time()
CPU_start_time = time.process_time()

# reading all parameters from params file and command line arguments and saving them to a csv file. See the functions and params files for the descriptions and default values of all parameters.

seed = params.seed

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
N_max = params.N_max 
n_cores = params.n_cores 
pop_dist_init = np.array(list(params.pop_dist))
K = params.K 


# the case with no second strike. If the population relapses, we note the N_min and t_min (time at which N_min occurs)
first_run_results = fn.first_run()
N_min = first_run_results[3]
t_min = first_run_results[4]

# if there is no relapse, then we move to the next random seed
if N_min==0:
    print("N_min is 0, move to the next seed")
    quit()

# save parameter values and random seed value to a csv file
fn.save_params(sim_data_dir+"/params.csv", "n", "w_lambda", "mus_vector",
            "cost_vector", "mutation_probs", "mutation_vector", "seed",
            "therapy_vector", "t_res", "T", "N_max", "n_cores", "pop_dist_init","K")


# the multiples of N_min to be tested with the same random seed
N_tau_factors = np.concatenate((np.linspace(20, 1, 30, endpoint=False), np.linspace(1,20, 31)))
N_tau_range = N_min*N_tau_factors
#N_tau_range = np.append(np.linspace(sum(pop_dist_init), N_min, 20, endpoint=False), np.linspace(N_min, sum(pop_dist_init), 21))
threshold_types  = np.append(np.tile("before-min",30), np.tile("after-min",31))
runs = np.arange(1,len(threshold_types)+1)

# Running the simulation            
results = Parallel(n_jobs=n_cores)(delayed(fn.one_run)(run, N_tau, threshold_type, t_min) for run,N_tau,threshold_type in zip(runs, N_tau_range, threshold_types))

results.append(first_run_results[:3])

# sorting simulation results
pop_history_allruns = sum([results[i][0] for i in range(len(N_tau_range)+1)], [])
time_history_allruns = sum([results[i][1] for i in range(len(N_tau_range)+1)], [])
N_tau_factors = np.append(N_tau_factors, 0)
outcomes = [[j ,results[i][2]] for i,j in enumerate(N_tau_factors)]


np.savetxt(sim_data_dir+"/outcome_data.csv", outcomes, delimiter=", ", fmt="% s")    
np.savetxt(sim_data_dir+"/lineage_data.csv", pop_history_allruns, delimiter=", ", fmt="% s")    
np.savetxt(sim_data_dir+"/time_data.csv", time_history_allruns, delimiter=", ", fmt="% s")

with open(sim_data_dir+"/mins.txt",'w') as f:
    print("N_min, ", N_min, file=f)
    print("t_min, ", t_min, file=f)


# calculate and print execution time
end_time = time.time()
CPU_end_time = time.process_time()
execution_time = (end_time-start_time)/60
CPU_time = (CPU_end_time-CPU_start_time)/60
print("Execution time: ", execution_time, " minutes")
print("CPU time: ", CPU_time, " minutes ")

