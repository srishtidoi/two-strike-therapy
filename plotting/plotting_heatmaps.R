suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  require(scales)
  library(latex2exp)
  library(metR)
  library(viridis)
})

theme_set(theme_linedraw())

# path to analytical and simulation data relative to the working directory
main_simdir <- "../simulation_data"
main_amdir <- "heatmaps/analytical_model_data"

# edit this block 
#################################################################

param_to_plot <- "K"
sub_simdir <- "heatmap_K"
sub_amdir <- "heatmap_K"
plot_filename <- "heatmap_K.pdf"

#################################################################

simdir <- file.path(main_simdir,sub_simdir)
amdir <- file.path(main_amdir,sub_amdir)

simdata <- read.table(file = file.path(simdir,"config.txt"), sep = '\t', header = TRUE)
amdata <- read.csv(file = file.path(amdir,"out.csv"))

results <- data.frame()
for(i in simdata$ArrayTaskID){
  probs <- read.csv(file.path(simdir, i, "results.txt"), header = FALSE)
  probs <- probs %>% remove_rownames %>% column_to_rownames(var="V1") 
  results <- rbind(results, c(probs[2,], probs[3,], probs[4,]))
}

colnames(results) <- c("extinction", "progression", "persistence")
colnames(simdata) <- c("sim no.","N_tau","param")
colnames(amdata) <- c("N_tau","param","PE_before","PE_after")
simdata <- cbind(simdata, results)

heatmap <- ggplot(data=simdata)+
  geom_contour_filled(aes(x=N_tau, y=param, z=extinction))+
  geom_contour(data=amdata, aes(x=N_tau, y=param, z=PE_before), colour='white')+
  labs(y=param_to_plot, x=TeX("$N(tau)$"))+
  #xlim(c(0,50000))+
  #ylim(c(10^6,10^7+2000))+
  theme(aspect.ratio = 1, 
        text=element_text(size = rel(4)),
        panel.grid.major = element_blank(),
        legend.title=element_blank(), 
        legend.text = element_text(size = rel(2)),
        legend.key.size = unit(0.8,"cm"))+
  geom_text_contour(data=amdata, aes(x=N_tau, y=param, z = PE_before),
                  colour='white',
                  skip = 1,
                  #nudge_x = 1500,
                  check_overlap = TRUE)

pdf(file.path("heatmaps/plots",plot_filename), height = 5, width=6)
print(heatmap)
dev.off()
