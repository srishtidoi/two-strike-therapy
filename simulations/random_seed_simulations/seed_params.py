import numpy as np
import argparse

# define all the parameters and corresponding arguments 

parser = argparse.ArgumentParser()
parser.add_argument("data_dir", help="path of the directory in which simulation data will be saved, relative to ./simulation-data/")
parser.add_argument("--seed", help="start seed for this set of parameters", type=float, default=0)
parser.add_argument("--n_runs", help="number of runs with these parameter values", type=int, default=10)
parser.add_argument("--wt", help="number of wild-type cells at t=0", type=int, default=10000)
parser.add_argument("--r1", help="number of R1 cells at t=0", type=int, default=1)
parser.add_argument("--r2", help="number of R2 cells at t=0", type=int, default=1)
parser.add_argument("--r12", help="number of R12 cells at t=0", type=int, default=0)
parser.add_argument("--n_max", help="threshold population value at which the simulation is stopped", type=int, default=100000)
parser.add_argument("--n_cores", help="number of cores to parallely run the simulation", type=int, default=1)
parser.add_argument("--T", help="max simulation time", type=int, default=100)
parser.add_argument("--K", help="carrying capacity for the total population", type=int, default=10002)

parser.add_argument("--n", help="number of cell types", type=int, default=4)
parser.add_argument("--w_lambda", help="birth rate of wild-type cells", type=float, default=1)
parser.add_argument("--w_mu", help="death rate of wild-type cells", type=float, default=0.1)
parser.add_argument("--cost", help="cost of resistance (given that it is the same for all resistant cell types)", type=float, default=0.5)
parser.add_argument("-m","--mutation_rate", help="mutation rate in all the environments", type=float, default=0.00001)
parser.add_argument("-t","--therapy", help="effect of therapy on death rate (same for all strikes)", type=float, default=1)
parser.add_argument("-t2","--therapy2", help="death rate due to second strike", type=float, default=1)

args=parser.parse_args()

data_dir = args.data_dir
n = args.n
n_wt = args.wt
n_r1 = args.r1
n_r2 = args.r2
n_r12 = args.r12
pop_dist = np.array([n_wt,n_r1,n_r2,n_r12])
w_lambda = args.w_lambda
w_mu = args.w_mu
mus_vector = np.array([w_mu,w_mu,w_mu,w_mu])
cost = args.cost
cost_vector = np.array([0,cost,cost,cost])
m = args.mutation_rate
mutation_vector = np.array([m,m])
t = args.therapy
t2 = args.therapy2
therapy_vector = np.array([t,t2])
N_max = args.n_max
T = args.T
K = args.K
n_cores = args.n_cores
seed = args.seed

# fixed parameters

t_res = 0.01
mutation_probs = np.array([0.25,0.25]) # probability of acquiring each type of mutation


