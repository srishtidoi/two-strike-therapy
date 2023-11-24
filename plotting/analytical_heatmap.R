suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
  library(metR)
})

#In this script, the y-axis parameter is K (carrying capacity). The x-axis parameter is always N(tau)

theme_set(theme_linedraw())

## change the relevant parameters here
# order of initial population relative to 10^4
ord <- 100
# initial subpopulation sizes
nr2 <- 1*ord 
nr1 <- 1*ord
nr12 <- 0
ns <- 10000*ord
N0 <- ns+nr1+nr2+nr12
# carrying capacity
#K <- N0
# probability of mutating from S --> R2
p2 <- 0.25
# probability of mutating from S --> R1
p1 <- 0.25
# total mutation rate
mu <- 0.00001
# rate of acquiring resistance to treatment 1
a1 <- mu*p1
# rate of acquiring resistance to treatment 2
a2 <- mu*p2
# treatment levels (delta_1 and delta_2)
D <- 1
D2 <- D 
# intrinsic birth rate
b <- 1
# intrinsic death rate
d <- 0.1
# intrinsic growth rate
s <- b-d
# cost of resistance
c <- 0.5
# growth rate of resistance cell types
r <- s-c

# establishment probability
pi2 <- 1-exp(-2*(r)/(b + d - c))


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

prange <- seq(1,3,0.02) # range of y-axis parameter
step <- 100
ntrange <- seq(0,50000,step)

grid_data <- expand.grid(ntrange, prange)

values <- c()
p_current <- -100 # current value of parameter to be varied
for(j in seq(nrow(grid_data))){
  K <- grid_data[j,2]
  N_tau <- grid_data[j,1]
  
  if(K!=p_current){
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
  
    p_current <- K
  }
  
  print(paste("N(tau)=",N_tau,"K=",K))
  
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


PE_data <- cbind(grid_data,values)
colnames(PE_data) <- c("ntau", "K", "pe")
write.csv(PE_data, file = "out.csv", row.names = FALSE)

var_names <- ls()
variables_to_save <- c("K")  # Replace with the names you want to save
variables_to_save <- intersect(variables_to_save, var_names)
var_data <- data.frame(t(sapply(variables_to_save, get)))
write.csv(var_data, file = "params.csv", row.names = FALSE)

heatmap <- ggplot(data=PE_data)+
  geom_contour_filled(aes(x=ntau, y=K, z=pe))+
  xlim(c(0,50000))+
  #ylim(c(1, 3))+
  xlab(TeX("$N_{\\tau}$"))+
  ylab(TeX("$K$"))+
  theme(aspect.ratio = 1, 
        legend.position = "", 
        text=element_text(size = rel(5)),
        panel.grid.major = element_blank())
  # geom_text_contour(aes(x=ntau, y=K,z = pe),
  #                   stroke = 0.15,
  #                   skip = 0,
  #                   nudge_x = 0,
  #                   check_overlap = TRUE)

pdf(file = "heatmap.pdf", width=6, height = 5)
print(heatmap)
dev.off()
