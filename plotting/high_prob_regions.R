suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(tidyverse)
  library(latex2exp)
  library(viridis)
  library(metR)
})
# Extracts data and created heatmaps showing 
# high extinction probability regions in the
# delta_1-delta_2 space (uncommented) and the
# b-d space (commented)


# setting theme for plots
theme_set(theme_linedraw())

# edit this block 
#################################################################

main_dir <- "optimal_doses/data_R2_0/"
# main_dir <- "b-d_space/data_default"
plots_dir <- "optimal_doses/plots"
# plots_dir <- "b-d_space/plots"

output_file <- "heatmap_data_cost_0.csv"
plot_file <- "heatmap_cost_0_after.pdf"

D1range <- seq(1,3,0.02)
D2range <- seq(1,3,0.02)
grid_data <- expand.grid(D1range, D2range)
# brange <- seq(0.05,3,0.05)
# drange <- seq(0,3,0.05)
# rel_cost <- read.csv(file=file.path(main_dir,"params.csv"))$rel_cost[1]
# grid_data <- expand.grid(brange, drange) %>% filter(Var2<=Var1*(1-rel_cost))

hp <- 0.8 # high probability threshold

#################################################################

before <- c() 
after <- c()
optimal_ntaus <- c()
nmins <- c()
optimal_pes <- c()
maxtau_before <- c()
maxtau_after <- c()
ba_flags <- c()

# initiate text progress bar
pb = txtProgressBar(min = 0, max = nrow(grid_data), initial = 0, style=3) 
stepi <- 0
for(j in seq(nrow(grid_data))){
  stepi <- stepi+1
  setTxtProgressBar(pb,stepi)
  
  ######### edit this block for b-d space plots #############
  D1 <- grid_data[j,1]
  D2 <- grid_data[j,2]
  # b <- grid_data[j,1]
  # d <- grid_data[j,2]

  data_filename <- paste("D1",D1,"_D2",D2,".csv",sep = "")
  # data_filename <- paste("b",b,"_d",d,".csv", sep = "")
  ###########################################################  
  
  PE_data <- read.csv(file=file.path(main_dir,data_filename))
  
  # high extinction probability range before N(min)
  highs_before <- PE_data[which(PE_data$PE_before>hp),1]
  min_value_before <- ifelse(length(highs_before) == 0, 0, min(highs_before))
  max_value_before <- ifelse(length(highs_before) == 0, 0, max(highs_before))
  high_range_before <- max_value_before - min_value_before
  before <- c(before, high_range_before)
  maxtau_before <- c(maxtau_before,max_value_before)
  
  # high extinction probability range after N(min)
  highs_after <- PE_data[which(PE_data$PE_after>hp),1]
  min_value_after <- ifelse(length(highs_after) == 0, 0, min(highs_after))
  max_value_after <- ifelse(length(highs_after) == 0, 0, max(highs_after))
  high_range_after <- max_value_after - min_value_after
  after <- c(after, high_range_after)
  maxtau_after <- c(maxtau_after, max_value_after)
  
  # max extinction probability; ba_flag tells whether it lies
  # before or after the nadir (0:at nadir, 1: after, -1:before)
  max_PE <- max(c(PE_data$PE_before, PE_data$PE_after))
  max_row <- which(PE_data==max_PE, arr.ind = TRUE)[1,1]
  max_col <- which(PE_data==max_PE, arr.ind = TRUE)[1,2]
  optimal_pes <- c(optimal_pes, max_PE)
  ba_flag <- ifelse(max_row==1,0, ifelse(max_col==2, -1, 1))
  ba_flags <- c(ba_flags, ba_flag)
  
  # optimal switching point and N(min)
  optimal_ntaus <- c(optimal_ntaus, PE_data$n_tau[max_row])
  nmins <- c(nmins,PE_data[1,1])
}
close(pb)

data <- cbind(grid_data, before, after, optimal_ntaus, nmins, ba_flags, optimal_pes, maxtau_before, maxtau_after)
colnames(data)<- c("X","Y","PE_before","PE_after", "optimal_N_tau","N_min","ba_flag","PE_optimal","max_N_tau_before","max_N_tau_after")
write.csv(data, file = file.path(plots_dir,output_file), row.names = FALSE)


# plotting high extinction probability region heatmap
data <- read.csv(file=file.path(plots_dir,output_file))
my_palette <- viridis_pal(option = "rocket")(13)
hpr_plot <- ggplot(data=data, aes(x=X,y=Y))+
  # edit this line to plot other quantities
  geom_contour_filled(aes(z=PE_after), breaks = c(0,20000,30000,35000,65000))+#,breaks=c(0,3000,7000,12000,16000))+
  scale_fill_manual(values = my_palette[c(1,6,8,10,12)])+
  geom_contour(aes(z=PE_optimal), col="white", breaks = c(0.9,0.94,0.96,0.98,0.99))+
  #geom_contour(aes(z=optimal_N_tau), col="yellow", breaks = c(1000,5000,10000,30000,50000))+
  #### edit the next four lines for b-d space plots
  geom_line(aes(x=X,y=X))+
  labs(x=TeX("$delta_{1}$"), y=TeX("$delta_{2}$"))+
  # geom_line(aes(x=X, y=X*(1-rel_cost)))+ # where d=b-cost
  # geom_line(aes(x=X, y=X-0.9), linetype=2)+ # where b-d=0.9
  # ylim(c(0,1.5))+    
  ####
  guides(fill = guide_legend(reverse = TRUE))+
  theme(legend.title=element_blank(), 
        aspect.ratio = 1,
        legend.text = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8,"cm"),
        text = element_text(size=rel(4.4)))

#pdf(file=file.path(plots_dir,plot_file), height=3, width = 3)
print(hpr_plot)
#dev.off()
