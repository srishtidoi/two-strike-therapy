import numpy as np
import random as rnd
import csv
import params
import functions as fn
from numpy.random import default_rng
rng = default_rng()

## parameters and initial conditions, imported from the params file
# number of runs of the simulation
n_runs = params.n_runs
# number of types of cells (# strikes + 1)
n = params.n
# base birth rate of wild type
w_lambda = params.w_lambda
# base death rates of all cell types
mus_vector = params.mus_vector
# cost of resistance, must be less than w_lambda
cost_vector = params.cost_vector
# mutation probabilities of acquiring resistance to each treatment
mutation_probs = params.mutation_probs
alpha_1 = mutation_probs[0]
alpha_2 = mutation_probs[1]
# total mutation rates in both the envs, sum of rates of acquiring mutations
mutation_vector = params.mutation_vector
# increase in death rate fur to treatments in each strike
therapy_vector = params.therapy_vector
# resolution of the time series data
t_res = params.t_res
# total time of the simulation
T = params.T
# population size at the switching point to second treatment
N_tau = params.N_tau
# population size at which the simulation stops
N_max = params.N_max
# number of simulations to run simultaneuosly; -1 will use all available cores
n_cores = params.n_cores
# initial sizes of each subpopulation
pop_dist_init = np.array(list(params.pop_dist))
# carrying capacity
K = params.K 

# birth rates after subtractin cost of resistance (sensitive birth rates remain the same)
lambdas_vector = np.tile(w_lambda, n) - cost_vector
# defining all possible events and labelling them
event_types = dict(zip([0,1,2], ["birth", "death", "mutation"]))
# list of all cell type labels; 0: sensitive, 1: resistant to treatment 1, 2: resistant to treatment 2, 3: resistant to both treatments
cell_types = np.arange(n) 

## defining functions

def save_params(filename, *args):
    # Get global dictionary
    glob = globals()
    d = {}
    for v in args:
        # Copy over desired values
        d[v] = glob[v]
    with open(filename, 'w') as f:
        w = csv.DictWriter(f, d.keys())
        w.writeheader()
        w.writerow(d)

def select_event(birth_rates, death_rates, mutation_rates):
    '''select which event happens, to which cell type and when. Gillespie implementation. All the event rates must be in a specific order, as follows- birth rates: S, R1, R2, R12; death rates: S, R1, R2, R12; mutation rates: S->R1, S->R2, R1->R12, R2->R12'''
    
    all_rates = np.concatenate([birth_rates, death_rates, mutation_rates])
    cum_prob = np.cumsum(all_rates)/np.sum(all_rates)
    chosen_event = np.argmax(rnd.random()<=cum_prob)

    chosen_event_type = event_types[chosen_event//4]
    if chosen_event_type=="mutation":
        if chosen_event==8:
            source_cell_type=0
            target_cell_type=1
        elif chosen_event==9:
            source_cell_type=0
            target_cell_type=2
        elif chosen_event==10:
            source_cell_type=1
            target_cell_type=3
        elif chosen_event==11:
            source_cell_type=2
            target_cell_type=3
    else:
        source_cell_type = chosen_event%4
        target_cell_type = source_cell_type
    
    time_int = rng.exponential(1/sum(all_rates),1)
    return [chosen_event_type, source_cell_type, target_cell_type, float(time_int)]

def birth(pop_dist, mother):
    '''input- vector of population distribution, mother cell type- 0:w, 1:r1, 2:r2, 3:r12'''
    pop_dist_new = np.array(list(pop_dist))
    if pop_dist[mother]>0:
        pop_dist_new[mother]+=1

    return pop_dist_new

def death(pop_dist, chosen):
    pop_dist_new = np.array(list(pop_dist))
    pop_dist_new[chosen]=max(0,pop_dist_new[chosen]-1)
    return pop_dist_new

def mutation(pop_dist, source, target):
    pop_dist_new = np.array(list(pop_dist))
    if pop_dist[source]>0:
        pop_dist_new[source]-=1
        pop_dist_new[target]+=1
        
    return pop_dist_new
    
def strike(N, threshold_type, i, pop_dist, frac):
    """changing the treatment to strike 2 after a certain threshold"""
    
    if threshold_type=="pop":
        if N<N_tau: i=1
    
    if threshold_type=="rescue":
        if pop_dist[1]>frac*K: i=1

    return i

def therapy(i):
    D = therapy_vector[i]
    D = np.tile(D,n)
    D[i+1]=0
    D[-1]=0
    return D

def one_run(run):
    print("run ", run, "of ", n_runs-1)
    
    pop_dist = np.array(list(params.pop_dist))
    N_init = sum(pop_dist)

    # initialize data arrays
    pop_history = []
    time_history = []
    pop_history.append(np.insert(pop_dist, 0, run))
    time_history.append([run,0])

    # initialise other variables
    i=0 # current strike number
    t=0 
    t_record = 0
    stop_flag = 0
    N = N_init

    while stop_flag==0:
        
        i = fn.strike(N, 'pop', i, pop_dist, 0.9)
        birth_rates = (lambdas_vector - ((lambdas_vector-mus_vector)*N/K))*pop_dist 
        death_rates = (mus_vector+fn.therapy(i))*pop_dist
        mutation_rates = (mutation_vector[i]*np.array([alpha_1*pop_dist[0], alpha_2*pop_dist[0], alpha_2*pop_dist[1], alpha_1*pop_dist[2]]))
        
        [event_type, source_cell, target_cell, time_int] = fn.select_event(birth_rates, death_rates, mutation_rates)
        
        if event_type=="birth": pop_dist_new = fn.birth(pop_dist, source_cell)
        elif event_type=="death": pop_dist_new = fn.death(pop_dist, source_cell)
        elif event_type=="mutation": pop_dist_new = fn.mutation(pop_dist, source_cell, target_cell)
        
        t+=time_int
        N = sum(pop_dist_new)
        
        if (t-t_record)>=t_res:
            t_record = t
            pop_history.append(np.insert(pop_dist_new, 0, run))
            time_history.append([run, t_record])
        
        pop_dist = pop_dist_new

        # stopping conditions
        if N==0: 
            stop_flag=1
            outcome = 0 # extinction
        if t>100 and N>=0.99*K:
            stop_flag=1
            outcome = 1 # rescue
        if N>=N_max:
            stop_flag=1
            outcome = 1 # rescue
        if t>=T:
            stop_flag=1
            if N>=N_init: outcome = 1 # rescue
            else: outcome = 2 # need more time

    return [pop_history, time_history, outcome]


