suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(deSolve)
})

# Calculates extinction probabilities for a range of
# delta_1 and delta_2 values, for multiple switching points

# edit this block 
#################################################################

# path relative to the working directory
main_dir <- "optimal_doses/all_data" # for plots like Fig2A

# treatment level space
D1range <- seq(1,3,0.02) # delta_1, treatment level in E1
D2range <- seq(1,3,0.02) # delta_2, treatment level in E2

# step and max for range of N_tau values
step <- 4000
nmax <- 100000

#################################################################

# default parameters
ord <- 100 # order of initial population relative to 10^4
nr2 <- 1*ord # initial R2 population
nr1 <- 1*ord # initial R1 population
nr12 <- 0 # initial R12 population
ns <- 10000*ord # initial S population
N0 <- ns+nr1+nr2+nr12 # total initial population
K <- N0 # carrying capacity
p2 <- 0.25 # probability of mutating from S -> R2
p1 <- 0.25 # probability of mutating from S -> R1
mu <- 0.00001 # total mutation rate in E1
mu2 <- 0.00001 # total mutation rate in E2
a1 <- mu*p1
a2 <- mu*p2
a1_2 <- mu2*p1
a2_2 <- mu2*p2
b <- 1 # intrinsic birth rate (of S cells)
d <- 0.1 # intrinsic death rate
s <- b-d # growth rate of S cells
cost <- 0.5 # cost of resistance
r <- s-cost # growth rate of resistant cells (R1,R2,R12)
pi2 <- 1-exp(-2*(r)/(b + d - cost)) # establishment probability


# define growth functions
grow_E1 <- function(t, init_state, growth_params){
  with(as.list(c(init_state, growth_params)),{
    dS = S*(gs*(1-(S+R1+R2+R12)/K) - D) - S*(a1+a2)
    dR1 = R1*gr*(1-(S+R1+R2+R12)/K) - R1*a2 + S*a1
    dR2 = R2*(gr*(1-(S+R1+R2+R12)/K) - D) - R2*a1 + S*a2
    dR12 = R12*gr*(1-(S+R1+R2+R12)/K) + R2*a1 + R1*a2
    
    list(c(dS,dR1,dR2,dR12))
  }) # end with(as.list...
}

grow_E2 <- function(t, init_state, growth_params){
  with(as.list(c(init_state, growth_params)),{
    dS = S*(gs*(1-(S+R1)/K) - D) - S*a1
    dR1 = R1*(gr*(1-(S+R1)/K) - D) + S*a1
    
    list(c(dS,dR1))
  }) # end with(as.list...
}

init_state_E1 <- c(S=ns, R1=nr1, R2=nr2, R12=nr12)
times_E1 <- seq(0,500,by=0.01)
grid_data <- expand.grid(D1range, D2range) 
pb = txtProgressBar(min = 0, 
                    max = nrow(grid_data), 
                    initial = 0, style=3) 
stepi <- 0
for(j in seq(nrow(grid_data))){
  D1 <- grid_data[j,1]
  D2 <- grid_data[j,2]
  stepi <- stepi+1
  setTxtProgressBar(pb,stepi)
  #paste("D=",D1,"D2=",D2)
  
  growth_params <- c(K=K, D=D1, gs=s, gr=r, a1=a1, a2=a2)
  out_E1 <- ode(y=init_state_E1, times=times_E1, func=grow_E1, parms=growth_params)
  time <- out_E1[,"time"]
  S <- out_E1[,"S"]
  R1 <- out_E1[,"R1"]
  R2 <- out_E1[,"R2"]
  R12 <- out_E1[,"R12"]
  T_E1 <- S+R1+R2+R12
  
  nmin_idx <- which.min(T_E1)
  T_E1_before <- T_E1[1:nmin_idx]
  T_E1_after <- T_E1[nmin_idx:length(T_E1)]
  nmin_num <- max(0,min(T_E1))
  ntrange <- seq(nmin_num,nmax,step)
  
  PEs_before <- c()
  PEs_after <- c()
  for(N_tau in ntrange){
    idx_before <- which.min(abs(T_E1_before - N_tau))
    idx_after <- which.min(abs(T_E1_after - N_tau))
    
    ## before nmin
    # growth in E2
    growth_params_E2 <- c(K=K, D=D2, gs=s, gr=r, a1=a1_2, a2=a2_2)
    init_state_E2_before <- c(S=S[idx_before], R1=R1[idx_before])
    times_E2 <- seq(0,200,by=0.1)
    
    out_E2_before <- ode(y=init_state_E2_before, times=times_E2, func=grow_E2, parms=growth_params_E2)
    time_E2_before <- out_E2_before[,"time"]
    S_E2_before <- out_E2_before[,"S"]
    R1_E2_before <- out_E2_before[,"R1"]
    
    # calculating the extinction probabilities
    I_S_before <- 0.1*sum(S_E2_before[which(S_E2_before>=1)])
    I_R1_before <- 0.1*sum(R1_E2_before[which(R1_E2_before>=1)])
    R12_estimate_before <- 0.01*(a2*sum(R1[1:idx_before])+a1*sum(R2[1:idx_before]))
    logPE <- -pi2*R2[idx_before] - pi2*a2*I_S_before - pi2*a2*I_R1_before - pi2*R12_estimate_before
    PE <- exp(logPE)
    PEs_before <- c(PEs_before, PE)
    
    ## after nmin
    tau_idx <- length(T_E1_before)+idx_after-1
    # growth in E2
    init_state_E2_after <- c(S=S[tau_idx], R1=R1[tau_idx])
    times_E2 <- seq(0,200,by=0.1)
    
    out_E2_after <- ode(y=init_state_E2_after, times=times_E2, func=grow_E2, parms=growth_params_E2)
    time_E2_after <- out_E2_after[,"time"]
    S_E2_after <- out_E2_after[,"S"]
    R1_E2_after <- out_E2_after[,"R1"]
    
    # calculating the extinction probabilities
    I_S_after <- 0.1*sum(S_E2_after[which(S_E2_after>=1)])
    I_R1_after <- 0.1*sum(R1_E2_after[which(R1_E2_after>=1)])
    R12_estimate_after <- 0.01*(a2*sum(R1[1:tau_idx])+a1*sum(R2[1:tau_idx]))
    logPE <- -pi2*R2[tau_idx] - pi2*a2*I_S_after - pi2*a2*I_R1_after - pi2*R12_estimate_after
    PE <- exp(logPE)
    PEs_after <- c(PEs_after, PE)
  }
  
  PEs_before[which(PE>1)]=1
  PEs_after[which(PE>1)]=1
  PE_data <- data.frame(ntrange, PEs_before, PEs_after)
  colnames(PE_data) <- c("N_tau", "PE_before", "PE_after")
  
  data_filename <- paste("D1",D1,"_D2",D2,".csv", sep = "")
  data_file <- file.path(main_dir,data_filename)
  write.csv(PE_data, file=data_file, row.names = FALSE)
}
close(pb)