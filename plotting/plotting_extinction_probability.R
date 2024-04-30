suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
})
# setting theme for plots
theme_set(theme_linedraw())

# edit this block 
#################################################################

# simulation data directory path relative to the working directory
main_simdir <- "../simulation_data"
sub_simdir <- "default"
simdir <- file.path(main_simdir, sub_simdir)

# analytical model data directory path relative to the working directory
main_amdir <- "extinction_probability_plots"
sub_amdir <- "default"
amdir <- file.path(main_amdir, sub_amdir)

#################################################################

# simulation data
all_seeds <- list.dirs(simdir, recursive = FALSE)
all_data <- read.csv(file.path(all_seeds[1],"outcome_data.csv"), header=FALSE)[1]
sim_nmins <- c()
for(seed_dir in all_seeds){
  data <- read.csv(file.path(seed_dir,"outcome_data.csv"), header=FALSE)
  nmin_data <- read.csv(file.path(seed_dir,"mins.txt"), header=FALSE)
  all_data <- cbind(all_data, data$V2)
  sim_nmins <- nmin_data[1,2]
}

avg_nmin <- mean(sim_nmins)
sums <- apply(all_data[,2:ncol(all_data)], 1, function(x) length(which(x==0)))
sums_p <- apply(all_data[,2:ncol(all_data)], 1, function(x) length(which(x==2)))
PP <- sums_p/(ncol(all_data)-1) # persistence probability
EP <- sums/(ncol(all_data)-1) # extinction probability
plot_data <- cbind(all_data[1], EP, PP)
# grouping before and after nadir switching points
types <- c(rep("before min", (nrow(all_data)-2)/2), "min", rep("after min",(nrow(all_data)-2)/2), "no 2nd strike")
plot_data <- cbind(plot_data, types)
# error bars
plot_data$SE <- sqrt(plot_data$EP * (1 - plot_data$EP) / 100)
# omitting before N_min switching points
plot_data[which(plot_data$V1<avg_nmin & (plot_data$types=="after min" | plot_data$types=="min")),4]="omit"

# analytical model data
PE_data <- read.csv(file = file.path(amdir,"out.csv"))
colnames(PE_data) <- c("N_tau", "PE_before", "PE_after")
nmax <- 100000 
nmin <- PE_data[1,1]

epplot <- ggplot(data = plot_data,aes(x=V1, color=types))+
  geom_point(aes(y=EP), size=3, alpha=1)+
  labs(x=TeX("$N(tau)$"), y=TeX("$P_{E}$"))+
  scale_colour_manual(values=c("white", "#56B4E9","#56B4E9","grey20","#CC79A7"),name="",
                      breaks=c("omit", "min" ,"after min", "before min",  "no 2nd strike"))+
  geom_errorbar(aes(ymin = EP - 1.96 * SE, ymax = EP + 1.96 * SE), width = 0.2)+
  geom_line(data=PE_data,
            aes(x=N_tau, y=PE_before), lwd=1, color='black')+
  geom_line(data=PE_data,
            aes(x=N_tau, y=PE_after), lwd=1, color='blue')+
  geom_vline(xintercept = nmin, lty=2, lwd=1.5, color='red', alpha=0.7)+
  geom_hline(yintercept = 0, col='grey50')+
  scale_x_continuous(limits = c(0, nmax), breaks = c(0, 0.5*nmax, nmax))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1, 
        legend.position = "",
        text = element_text(size=rel(5)))


pdf(file.path(amdir,"plot_test.pdf"), height = 5, width = 6)
print(epplot)
dev.off()
