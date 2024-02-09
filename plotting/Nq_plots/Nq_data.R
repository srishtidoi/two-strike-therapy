suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
})

# file to save data
filename <- "Nq_plots/data/filename"

# change the relevant parameters here
ord <- 100 # order of initial population relative to 10^4
nr2 <- 1*ord # initial R_2 population
nr1 <- 1*ord # initial R_1 population
nr12 <- 0 # initial R_1,2 population
ns <- 10000*ord # initial S population
N0 <- ns+nr1+nr2+nr12
K <- N0 # carrying capacity
p2 <- 0.25 # probability of mutating from S --> R2 (given mutation occurs)
p1 <- 0.25 # probability of mutating from S --> R1 (given mutation occurs)
mu <- 0.00001 # total mutation probability in E_1
mu2 <- 0.00001 # total mutation probability in E_2
a1 <- mu*p1 
a2 <- mu*p2
a1_2 <- mu2*p1
a2_2 <- mu2*p2
D <- 2.0 # treatment level in E_1
D2 <- D # treatment level in E_2
b <- 1 # intrinsic birth rate
d <- 0.1 # intrinsic death rate
s <- b-d
cost <- 0.5 # cost of resistance
r <- s-cost
pi2 <- 1-exp(-2*(r)/(b + d - cost)) # establishment probability

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

growth_params <- c(K=K, D=D, gs=s, gr=r, a1=a1, a2=a2)
init_state_E1 <- c(S=ns, R1=nr1, R2=nr2, R12=nr12)
times_E1 <- seq(0,500,by=0.01)

out_E1 <- ode(y=init_state_E1, times=times_E1, func=grow_E1, parms=growth_params)
time <- out_E1[,"time"]
S <- out_E1[,"S"]
R1 <- out_E1[,"R1"]
R2 <- out_E1[,"R2"]
R12 <- out_E1[,"R12"]
T_E1 <- S+R1+R2+R12

nmin_idx <- which.min(T_E1)
T_E1_premin <- T_E1[1:nmin_idx] # choosing thresholds before nmin only
nmin_num <- min(T_E1)

step <- 1000
ntrange <- seq(nmin_num,N0,step)

PEs <- c()
for(N_tau in ntrange){
  print(N_tau)
  idx <- which.min(abs(T_E1_premin - N_tau))

  # growth in E2
  growth_params_E2 <- c(K=K, D=D2, gs=s, gr=r, a1=a1_2, a2=a2_2)
  init_state_E2 <- c(S=S[idx], R1=R1[idx])
  times_E2 <- seq(0,200,by=0.1)
  
  out_E2 <- ode(y=init_state_E2, times=times_E2, func=grow_E2, parms=growth_params_E2)
  time_E2 <- out_E2[,"time"]
  S_E2 <- out_E2[,"S"]
  R1_E2 <- out_E2[,"R1"]
  
  I_S <- 0.1*sum(S_E2[which(S_E2>=1)])
  I_R1 <- 0.1*sum(R1_E2[which(R1_E2>=1)])
  
  logPE <- -pi2*R2[idx] - pi2*R12[idx] - pi2*a2_2*I_S - pi2*a2_2*I_R1
  PE <- exp(logPE)
  PEs <- c(PEs, PE)
}

PEs[which(PE>1)]=1
PE_data <- data.frame(ntrange, PEs)
colnames(PE_data) <- c("n_tau", "PE")

nt_thresholds <- c()
for(pe in seq(1,0,-0.01)){
  nt <- max(PE_data$n_tau[which(PE_data$PE>=pe)])
  nt_thresholds <- c(nt_thresholds, nt)
}

write.csv(PE_data, file =paste(filename,".csv", sep=""), row.names = FALSE)
saveRDS(nt_thresholds, file=paste(filename,".RData", sep=""))

