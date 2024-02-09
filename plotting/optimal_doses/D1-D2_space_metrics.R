suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
  library(viridis)
  library(metR)
})

D1range <- seq(1,3,0.02)
D2range <- seq(1,3,0.02)
grid_data <- expand.grid(D1range, D2range) 

norm_aucs <- c()
for(j in seq(nrow(grid_data))){
  D1 <- grid_data[j,1]
  D2 <- grid_data[j,2]
  
  print(paste("D1=",D1,"D2=",D2))
  filename <- paste("optimal_doses/all_data/","D1",D1,"_D2",D2, sep = "")
  PE_data <- read.csv(file =paste(filename,".csv", sep=""))
  
  total_auc <- 500000
  norm_auc <- sum(-diff(PE_data$PE) * (PE_data$n_tau[-1] + PE_data$n_tau[-length(PE_data$n_tau)]) / 2)/total_auc
  norm_aucs <- c(norm_aucs, norm_auc)
  
  ### other metrics
  ## mean of N_0.1 and N_0.9
  #mid <- mean(c(max(PE_data$n_tau[which(PE_data$PE>=0.1)]),max(PE_data$n_tau[which(PE_data$PE>=0.9)])))
  #mids <- c(mids, mid)
  ## value of N_0.5
  #N_0.5 <- max(PE_data$n_tau[which(PE_data$PE>=0.5)])
  #N_0.5s <- c(N_0.5s, N_0.5)
  
}

data <- cbind(grid_data, norm_aucs)
colnames(data)<- c("D1","D2","norm_auc")

write.csv(data, file = "optimal_doses/normalised_aucs.csv", row.names = FALSE)

# plot the metric in the D1-D2 space
my_palette <- viridis_pal(option = "rocket")(12)
theme_set(theme_linedraw())

# uncomment to load data to plot other metrics
#data <- read.csv(file="optimal_doses/normalised_aucs.csv")

# uncomment to normalise the N_0.5 metric
#data <- data %>% filter(N_0.5>0)
#data$N_0.5 <- data$N_0.5/N0

doseplot <- ggplot(data=data, aes(x=D1,y=D2))+
  geom_contour_filled(aes(z=norm_auc))+
  scale_fill_manual(values = my_palette)+
  # geom_text_contour(aes(x=D1, y=D2, z=norm_auc, fontface="bold"),
  #                   colour="black",
  #                   nudge_y=-0.04,
  #                   nudge_x=0.05,
  #                   skip=1)+
  geom_line(aes(x=D1, y=D1))+
  #ylim(c(0,1.5))+
  labs(x=TeX("$delta_{1}$"), y=TeX("$delta_{2}$"))+
  theme(legend.position = "none",
        text = element_text(size=rel(3.7)))

pdf(file="optimal_doses/plots/normalised_aucs.pdf", height=3, width = 3)
print(doseplot)
dev.off()
