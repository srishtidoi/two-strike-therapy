suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
  library(metR)
  library(pbapply)
})

# In this script, the y-axis parameter is K (carrying capacity)
# The x-axis parameter is always N(tau)

theme_set(theme_linedraw())

# default parameters 
ord <- 100 # order of initial population relative to 10^4
nr2 <- 1*ord # initial R2 population
nr1 <- 1*ord # initial R1 population
nr12 <- 0 # initial R12 population
ns <- 10000*ord # initial S population
N0 <- ns+nr1+nr2+nr12
K <- N0 # carrying capacity
p2 <- 0.25 # probability of mutating from S -> R2 (given mutation occurs)
p1 <- 0.25 # probability of mutating from S -> R1 (given mutation occurs)
mu <- 0.00001 # total mutation probability in E1
mu2 <- 0.00001 # total mutation probability in E2
a1 <- mu*p1 
a2 <- mu*p2
a1_2 <- mu2*p1
a2_2 <- mu2*p2
D <- 2.0 # treatment level in E1
D2 <- 2.0 # treatment level in E2
b <- 1 # intrinsic birth rate (od S cells)
d <- 0.1 # intrinsic death rate
s <- b-d # growth rate of S cells
cost <- 0.5 # cost of resistance
r <- s-cost # growth rate of resistant cells (R1, R2, R12)
pi2 <- 1-exp(-2*(r)/(b + d - cost)) # establishment probability

main_dir <- "heatmaps/analytical_model_data"
plot_filename <- "analytical_model_plot.pdf"
params_filename <- "params.csv"
data_filename <- "out.csv"

# edit this block 
#################################################################

param_to_plot <- "K"
sub_dir <- "heatmap_K"

update_original <- function(param){
  # edit this line with the parameter to plot
  K <<- param
}

prange <- seq(N0,10*N0,200000) # range of y-axis parameter

#################################################################

output_dir <- file.path(main_dir, sub_dir)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Output directory already exists!")
}

# save parameters
var_names <- ls()
variables_to_save <- c("b","d","cost","D","D2","K","nr1","nr2","nr12","ns","N0","a1","a2","a1_2","a2_2","ord")  # add variables to save
variables_to_save <- variables_to_save[variables_to_save!=param_to_plot]
# uncomment if D = D2
#variables_to_save <- variables_to_save[variables_to_save!="D2"]

variables_to_save <- intersect(variables_to_save, var_names)
var_data <- data.frame(t(sapply(variables_to_save, get)))
write.csv(var_data, file = file.path(output_dir,params_filename), row.names = FALSE) 

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

# generate parameter grid
step <- 10000
ntrange <- seq(0,50000,step) # range of N_tau
grid_data <- expand.grid(ntrange, prange)

# initiate progress bar
pb = txtProgressBar(min = 0, max = nrow(grid_data), initial = 0, style=3) 

PEs_before <- c()
PEs_after <- c()
p_current <- -100 # current value of parameter to be varied
stepi <- 0
for(j in seq(nrow(grid_data))){
  param <- grid_data[j,2] 
  update_original(param)
  
  # uncomment if D = D2
  #D2 <- D
  
  # uncomment if the param is birth rate
  #cost <- b/2
  # uncomment if any of the growth parameters are being varied
  #s <- b-d
  #r <- s-cost
  #pi2 <- 1-exp(-2*(r)/(b + d - cost)) # establishment probability
  
  # update progress bar
  setTxtProgressBar(pb,stepi)
  stepi <- stepi+1
  
  N_tau <- grid_data[j,1]
  
  if(param!=p_current){ # 
    init_state_E1 <- c(S=ns, R1=nr1, R2=nr2, R12=nr12)
    times_E1 <- seq(0,500,by=0.01)
    growth_params <- c(K=K, D=D, gs=s, gr=r, a1=a1, a2=a2)
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
    nmin_num <- min(T_E1)
    
    p_current <- param 
  }
  
  if(N_tau < nmin_num){
    PE_before <- 0
    PE_after <- 0
  }else{
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
    logPE <- -pi2*R2[idx_before] - pi2*a2_2*I_S_before - pi2*a2_2*I_R1_before - pi2*(R12_estimate_before) 
    PE_before <- exp(logPE)
    
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
    logPE <- -pi2*R2[tau_idx] - pi2*a2_2*I_S_after - pi2*a2_2*I_R1_after - pi2*(R12_estimate_after)
    PE_after <- exp(logPE)
  }
  PEs_before <- c(PEs_before, PE_before)
  PEs_after <- c(PEs_after, PE_after)
}
close(pb)

PEs_before[which(PEs_before>1)]=1
PEs_after[which(PEs_after>1)]=1
PE_data <- cbind(grid_data,PEs_before, PEs_after)
colnames(PE_data) <- c("N_tau", "param", "before", "after") 
write.csv(PE_data, file = file.path(output_dir,data_filename), row.names = FALSE)

# uncomment to plot analytical heatmap

# heatmap <- ggplot(data=PE_data)+
#   geom_contour_filled(aes(x=ntau, y=param, z=before))+ 
#   xlim(c(0,50000))+
#   #ylim(c(1, 3))+
#   xlab(TeX("$N(tau)$"))+
#   ylab(param_to_plot)+
#   theme(aspect.ratio = 1, 
#         legend.position = "", 
#         text=element_text(size = rel(5)),
#         panel.grid.major = element_blank())
# # +geom_text_contour(aes(x=ntau, y=param,z = pe),
# #                   stroke = 0.15,
# #                   skip = 0,
# #                   nudge_x = 0,
# #                   check_overlap = TRUE)
# 
# pdf(file = file.path(output_dir, plot_filename), width=6, height = 5)
# print(heatmap)
# dev.off()
