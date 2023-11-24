suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)
dirname <- args[1]
#dirname <- "simulation-data/test/n_tau1000"

theme_set(theme_linedraw())
dir <- paste("~/padhai/thesis/code/ET-new/", dirname, sep='')

# read lineage and time data for all runs
lineage_data <- read.csv(paste(dir,"lineage_data.csv",sep = "/"))
time_data <- read.csv(paste(dir,"time_data.csv",sep = "/"))
params <- read.csv(paste(dir,"params.csv",sep = "/"))

# defining number of cell types and number of runs
n <- params$n
n_runs <- params$n_runs

# compiling all data by adding population size and run number as columns
size_data <- rowSums(lineage_data[,1:n+1])
colnames(lineage_data) <- c('run',0, seq(n-1))
all_data <- cbind(time_data, size_data, lineage_data[,1:n+1])

# initializing plot
lineage_plot <- ggplot()

# looping over all runs and plotting run data
for(run in seq(0,n_runs-1)){
  run_data <- all_data[all_data$run==run,]
  plot_data <- melt(run_data[,2:ncol(run_data)], id.vars=c('t'))
  lineage_plot <- lineage_plot +
    geom_line(data=plot_data, aes(x=t, y=value, colour=variable),size=0.7)
}

# customizing plot (colours, scale, etc.)  
lineage_plot <- lineage_plot+
  scale_color_manual(values=c('grey60','limegreen','tomato', 'red2'), labels=c("Total population","sensitive","resistant (strike 1)", "resistant (strike 2)"))+
  scale_y_continuous(trans="log10")+
  geom_hline(yintercept=params$N_tau, linetype=2, color="grey40")+
  ylab("Number of cells")+
  xlab("t(days)")+
  labs(color="")+
  theme(legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major=element_line(colour="grey70"),panel.grid.minor=element_line(colour="grey70"))

pdf(paste(dir,"Rplots.pdf",sep = "/"))
suppressWarnings(print(lineage_plot))
dev.off()
