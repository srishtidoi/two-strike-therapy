suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  require(scales)
  library(latex2exp)
  library(metR)
})

theme_set(theme_linedraw())

# load all simulation data with 100 random seeds
dirname <- "path/to/data/directory"
all_seeds <- list.dirs(dirname, recursive = FALSE)
all_data <- read.csv(paste(all_seeds[1],"outcome_data.csv", sep = '/'), header=FALSE)[1]

for(seed_dir in all_seeds){
  data <- read.csv(paste(seed_dir,"outcome_data.csv", sep = '/'), header=FALSE)
  all_data <- cbind(all_data, data$V2)
}

sums <- apply(all_data[,2:ncol(all_data)], 1, function(x) length(which(x==0)))
sums_p <- apply(all_data[,2:ncol(all_data)], 1, function(x) length(which(x==2)))
PP <- sums_p/(ncol(all_data)-1)
EP <- sums/(ncol(all_data)-1)
plot_data <- cbind(all_data[1], EP, PP)
types <- c(rep("before min", (nrow(all_data)-2)/2), "min", rep("after min",(nrow(all_data)-2)/2), "no 2nd strike")
plot_data <- cbind(plot_data, types)
plot_data$SE <- sqrt(plot_data$EP * (1 - plot_data$EP) / 100)

seedsplot <- ggplot(data = plot_data, aes(x=V1, color=types))+
  geom_vline(xintercept = 1, lty=2, lwd=2, col="grey35")+
  geom_point(aes(y=EP), size=3, alpha=1)+
  labs(x=TeX(""), y="")+
  scale_colour_manual(values=c( "grey20",'orange', "#56B4E9","#CC79A7"),name="",
                      breaks=c("before min", "min" ,"after min", "no 2nd strike"),
                      labels=c(TeX("before $N_{min}$"), TeX("$N_{min}$"), TeX("after $N_{min}$"), TeX("No strike 2")))+
  geom_errorbar(aes(ymin = EP - 1.96 * SE, ymax = EP + 1.96 * SE), width = 0.2)+
  theme(legend.text.align = 0)+
  ylim(c(-0.01,1.01))+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        legend.text = element_text(size=rel(3)),
        legend.key.size = unit(0.8,"cm"),
        text=element_text(size=rel(6)))

pdf(paste(dirname, "random_seeds_plot.pdf",sep = '/'), height = 5, width = 7)
print(seedsplot)
dev.off()