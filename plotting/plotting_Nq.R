suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
})
# setting theme for plots
theme_set(theme_linedraw())

main_dir <- "Nq_plots/data"
plots_dir <- "Nq_plots/plots"

# edit this block 
#################################################################

dir_default <- file.path(main_dir,"default") # solid line
dir1 <- file.path(main_dir,"cost_0") # dashed line
dir2 <- file.path(main_dir,"cost_0.2") # dotted line
dir3 <- file.path(main_dir,"cost_0.7") # mixed line
dir4 <- file.path(main_dir,"cost_0.9") # long dashed line

plot_filename <- "cost.pdf"

#################################################################

data_default <- read.csv(file = file.path(dir_default,"out.csv"))
data_first <- read.csv(file = file.path(dir1,"out.csv"))
data_second <- read.csv(file = file.path(dir2,"out.csv"))
data_third <- read.csv(file = file.path(dir3,"out.csv"))
data_fourth <- read.csv(file = file.path(dir4,"out.csv"))

params_default <- read.csv(file = file.path(dir_default,"params.csv"))
params_first <- read.csv(file = file.path(dir1,"params.csv"))
params_second <- read.csv(file = file.path(dir2,"params.csv"))
params_third <- read.csv(file = file.path(dir3,"params.csv"))
params_fourth <- read.csv(file = file.path(dir4,"params.csv"))

compute_thresholds <- function(PE_data, N0){
  befores <- c()
  afters <- c()
  for(pe in seq(1,0,-0.01)){
    nt_before <- max(PE_data$N_tau[which(PE_data$PE_before>=pe)])
    nt_after <- max(PE_data$N_tau[which(PE_data$PE_after>=pe)])
    befores <- c(befores, nt_before)
    afters <- c(afters, nt_after)
  }
  data.frame(befores, afters)/N0
}

default <- compute_thresholds(data_default, params_default$N0)
first <- compute_thresholds(data_first, params_first$N0)
second <- compute_thresholds(data_second, params_second$N0)
third <- compute_thresholds(data_third, params_third$N0)
fourth <- compute_thresholds(data_fourth, params_fourth$N0)

y_breaks <- c(1, 0.1, 0.01) #c(1e+04, 1e+06,1e+08)
y_lims <- c(0.001,1) #c(1e+03,1e+08)

nqplot <- ggplot()+
  geom_line(aes(x=seq(1,0,-0.01), y=default$befores))+
  geom_line(aes(x=seq(1,0,-0.01), y=first$befores), linetype=2)+
  geom_line(aes(x=seq(1,0,-0.01), y=second$befores), linetype=3)+
  geom_line(aes(x=seq(1,0,-0.01), y=third$befores), linetype=4)+
  geom_line(aes(x=seq(1,0,-0.01), y=fourth$befores), linetype=5)+
  geom_line(aes(x=seq(1,0,-0.01), y=default$afters),colour="#56B4E9", alpha=0.8)+
  geom_line(aes(x=seq(1,0,-0.01), y=first$afters),colour="#56B4E9", linetype=2, alpha=0.8)+
  geom_line(aes(x=seq(1,0,-0.01), y=second$afters),colour="#56B4E9", linetype=3, alpha=0.8)+
  geom_line(aes(x=seq(1,0,-0.01), y=third$afters),colour="#56B4E9", linetype=4, alpha=0.8)+
  geom_line(aes(x=seq(1,0,-0.01), y=fourth$afters),colour="#56B4E9", linetype=5, alpha=0.5)+
  scale_y_continuous(trans="log", breaks=y_breaks, limits = y_lims)+
  ylab(TeX("$N_{q}/N(0)$"))+xlab(TeX("$q$"))+
  #ylab("")+xlab("")+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = rel(4.4)))

pdf(file.path(plots_dir,plot_filename), height = 3, width = 3)
print(nqplot)
dev.off()
