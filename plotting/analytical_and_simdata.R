suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
})
# setting the theme for plots
theme_set(theme_linedraw())

# directory with simulation data
dirname <- "/path/to/simulation/data"

# number of runs in the simulations, can be checked in the params.csv files
n_runs <- 500

# extracting config data (list of all parameters with ArrayTaskIDs)
data <- read.table(file = paste(dirname,"config.txt",sep = "/"), sep = '\t', header = TRUE)

# extracting extinction probabilities from results file for each parameter set
results <- data.frame()
for(i in data$ArrayTaskID){
  probs <- read.csv(paste(dirname, i, "results.txt", sep = "/"), header = FALSE)
  probs <- probs %>% remove_rownames %>% column_to_rownames(var="V1") 
  results <- rbind(results, c(probs[2,], probs[3,], probs[4,]))
}
colnames(results) <- c("extinction", "progression", "persistence")

# organising data and calculating standard error
data <- cbind(data, results)
data$SE <- sqrt(data$extinction*(1 - data$extinction)/n_runs) 

# point of max extinction probability
maxn <- data[which.max(data$extinction),]$n_tau

# analytical model params
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
c <- 0.5 # cost of resistance
r <- s-c
pi2 <- 1-exp(-2*(r)/(b + d - c)) # establishment probability

# growth equations
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

# growth in E1
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
T_E1_premin <- T_E1[1:nmin_idx]
nmin_num <- min(T_E1)

results <- c()
step <- 10000
ntrange <- seq(nmin_num,N0,step)

# calculating probabilities of extinction (P_E)
PEs <- c()
for(N_tau in ntrange){
  #print(N_tau)
  idx <- which.min(abs(T_E1_premin - N_tau))

  # growth in E2
  init_state_E2 <- c(S=S[idx], R1=R1[idx])
  times_E2 <- seq(0,200,by=0.1)
  
  out_E2 <- ode(y=init_state_E2, times=times_E2, func=grow_E2, parms=growth_params)
  time_E2 <- out_E2[,"time"]
  S_E2 <- out_E2[,"S"]
  R1_E2 <- out_E2[,"R1"]
  
  I_S <- 0.1*sum(S_E2[which(S_E2>=1)])
  I_R1 <- 0.1*sum(R1_E2[which(R1_E2>=1)])
  
  logPE <- -pi2*R2[idx] - pi2*R12[idx] - pi2*a2*I_S - pi2*a2*I_R1
  PE <- exp(logPE)
  PEs <- c(PEs, PE)
}

PEs[which(PE>1)]=1
PE_data <- data.frame(ntrange, PEs)
colnames(PE_data) <- c("n_tau", "PE")

nmin <- min(T_E1)
print(paste("N_min = ",nmin))

# plotting simulation data
nmax <- 100000
PEplot <- ggplot(data=data, aes(y=extinction))+
  geom_point(aes(x=n_tau), size=3)+
  #geom_vline(xintercept = maxn, lty=2,lwd=2,alpha=0.5)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        text=element_text(size=rel(6)),
        aspect.ratio = 1)+
  labs(x=TeX("N(tau)"), y=TeX("$P_{E}$"))+
  ylim(c(0,1))+
  geom_errorbar(aes(x=n_tau, ymin = extinction - 1.96 * SE, ymax = extinction + 1.96 * SE)) # binomial proportion 95% confidence interval

# plotting the analytical prediction
PEplot <- PEplot+
  geom_line(data=PE_data[which(PE_data$n_tau>-1),],
            aes(x=n_tau, y=PE), lwd=1, alpha=0.7)+
  geom_vline(xintercept = nmin, lty=2, lwd=1.5, color='red', alpha=0.7)+
  geom_hline(yintercept = 0, col='grey50')+
  scale_x_continuous(limits = c(0, nmax), breaks = c(0, 0.5*nmax, nmax))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1, legend.position = "",
        text = element_text(size=rel(6)))

pdf(paste(dirname,"PEplot.pdf",sep = "/"), height = 5, width = 6)
print(epplot)
dev.off()
