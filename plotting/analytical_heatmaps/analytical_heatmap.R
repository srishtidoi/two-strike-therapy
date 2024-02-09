suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
  library(metR)
})

# In this script, the y-axis parameter is K (carrying capacity)
# The x-axis parameter is always N(tau)

theme_set(theme_linedraw())

## change the relevant parameters here
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

# edit this block according to the parameter to plot against N(tau)
#################################################################
#################################################################
#################################################################
param_to_plot <- "carrying capacity"
dirname <- "analytical_heatmaps/data/heatmap_K"
plotname <- "heatmap_K.pdf"

update_original <- function(param){
  # edit this line with the parameter to plot
  K <<- param
}

prange <- seq(N0,10*N0,100000) # range of y-axis parameter

#################################################################
#################################################################
#################################################################

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
step <- 100
ntrange <- seq(0,50000,step) # range of N_tau
grid_data <- expand.grid(ntrange, prange)

values <- c()
p_current <- -100 # current value of parameter to be varied
for(j in seq(nrow(grid_data))){
  param <- grid_data[j,2] #
  update_original(param)
  
  # uncomment if any of the growth parameters are being varied
  #s <- b-d
  #r <- s-cost
  #pi2 <- 1-exp(-2*(r)/(b + d - cost)) # establishment probability

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
    T_E1_premin <- T_E1[1:nmin_idx] # choosing thresholds before nmin only
    nmin_num <- min(T_E1)
  
    p_current <- param 
    print(paste(param_to_plot,"=",param))
  }
  
  if(N_tau < nmin_num){
    PE <- 0
  }else{
    idx <- which.min(abs(T_E1_premin - N_tau))
    # growth in E2
    init_state_E2 <- c(S=S[idx], R1=R1[idx])
    times_E2 <- seq(0,200,by=0.1)
    growth_params_E2 <- c(K=K, D=D2, gs=s, gr=r, a1=a1, a2=a2)
    
    out_E2 <- ode(y=init_state_E2, times=times_E2, func=grow_E2, parms=growth_params_E2)
    time_E2 <- out_E2[,"time"]
    S_E2 <- out_E2[,"S"]
    R1_E2 <- out_E2[,"R1"]

    I_S <- 0.1*sum(S_E2[which(S_E2>=1)])
    I_R1 <- 0.1*sum(R1_E2[which(R1_E2>=1)])
    logPE <- -pi2*R2[idx] - pi2*R12[idx] - pi2*a2*I_S - pi2*a2*I_R1
    PE <- exp(logPE)
  }
  values <- c(values, PE)
}

values[which(values>1)]=1
PE_data <- cbind(grid_data,values)
colnames(PE_data) <- c("ntau", "param", "pe") 

heatmap <- ggplot(data=PE_data)+
  geom_contour_filled(aes(x=ntau, y=param, z=pe))+ 
  xlim(c(0,50000))+
  #ylim(c(1, 3))+
  xlab(TeX("$N(tau)$"))+
  ylab(param_to_plot)+
  theme(aspect.ratio = 1, 
        legend.position = "", 
        text=element_text(size = rel(5)),
        panel.grid.major = element_blank())
  # +geom_text_contour(aes(x=ntau, y=param,z = pe),
  #                   stroke = 0.15,
  #                   skip = 0,
  #                   nudge_x = 0,
  #                   check_overlap = TRUE)



var_names <- ls()
variables_to_save <- c("b","D","D2","cost","d","ord")  # add variables to save
variables_to_save <- intersect(variables_to_save, var_names)
var_data <- data.frame(t(sapply(variables_to_save, get)))

write.csv(var_data, file = paste(dirname,"params.csv",sep="/"), row.names = FALSE) 
write.csv(PE_data, file = paste(dirname,"out.csv",sep="/"), row.names = FALSE)
pdf(file = paste("analytical_heatmaps/plots/", plotname, sep = ""), width=6, height = 5)
print(heatmap)
dev.off()
