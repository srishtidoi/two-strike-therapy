suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(deSolve)
})

# This script plots the normalised N_q vs q plot for varying K
# using the analytical model results

theme_set(theme_linedraw())

# initial population size for normalisation
N0 <- 1000200

default <- readRDS(file = "Nq_plots/data/default.RData")/N0
first <- readRDS(file = "Nq_plots/data/K_2N0.RData")/N0
second <- readRDS(file = "Nq_plots/data/K_100N0.RData")/N0
#third <- readRDS(file = "Nq_plots/data/_")/N0
#fourth <- readRDS(file = "Nq_plots/data/_")/N0

data <- data.frame(default, first, second)

Nq_plot <- ggplot(data=data)+
  geom_line(aes(x=seq(1,0,-0.01), y=default))+
  geom_line(aes(x=seq(1,0,-0.01), y=first), linetype=2)+
  geom_line(aes(x=seq(1,0,-0.01), y=second), linetype=3)+
  #geom_line(aes(x=seq(1,0,-0.01), y=third), linetype=4)+
  #geom_line(aes(x=seq(1,0,-0.01), y=fourth), linetype=5)+
  scale_y_continuous(trans="log", breaks = c(1, 0.1, 0.01, 0.001))+#, limits = c(0.001,1))+
  ylab(TeX("$N_{q}$"))+xlab(TeX("$q$"))+
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = rel(5.8)))

pdf("Nq_plots/plots/K.pdf", height = 5, width = 6)
print(Nq_plot)
dev.off()

