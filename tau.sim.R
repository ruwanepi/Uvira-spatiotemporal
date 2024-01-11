# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script performs a simulation of tau for 100 points that are either geotagged
# to a centroid of a polygon/avenue or randomly-jittered to simuate household locations
# Author: R Ratnayake, 2023, based on A Azman, 2018 (code: https://osf.io/4fsnc/)

# Simulation of household locations

# load packages
x1 <- c("here", "tidyverse", "lubridate", "skimr", "ggplot2", "sf",
        "tmap", "tsibble", "EpiEstim", "rsatscan", "xlsx", "rgdal",
        "cowplot", "ggmap", "lubridate", "incidence2", "raster",
        "exactextractr", "R0", "SciViews", "imputeTS", "feasts",
        "mgcv", "IDSpatialStats", "slider")
# Install any packages not yet installed
x2 <- x1 %in% row.names(installed.packages())
if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
lapply(x1, library, character.only = TRUE)

#set working directory
setwd("C:/R-projects/Uvira-spatial-epid/")

# load pc dataset
source(here('tau.utils.R'))
pcall <- load_RDTpos_cases20162020()
#pc <- read_csv("cleaned_simplified_points_uvira3.csv") #2020 only
str(pcall)
unique(pcall$x, pcall$y)

plot(pcall$x, pcall$y, 
     xlab = "X", ylab = "Y",
     pch = 20, col='red', frame = FALSE)

uvira.ave.shape <- 
  here("C:/R-projects/Uvira-spatial-epid/Uvira_shp/Uvira_avenues.shp") %>% 
  st_read()

plot(uvira.ave.shape)
plot(uvira.ave.shape.2.cent)

ggplot() + 
  geom_sf(data = uvira.ave.shape, fill = "grey", group='UGA') +
  geom_sf(data = uvira.ave.shape.2.cent) +
  coord_sf()

# dataset 1: this dataset is the centroid-tagged dataset
pc_dt1 <- as.data.frame(#pc[1:200, ]
  pcall)
unique(pc_dt1$x, pc_dt1$y)


# dataset 2: retain these points and add a random jitter with a large factor
# to the point to simulate individual household locations
set.seed(1981)

pc_dt2 <- as.data.frame(cbind(jitter(pc_dt1$x, factor=500), 
                              jitter(pc_dt1$y, factor=500), 
                              pc_dt1$time))
colnames(pc_dt2) <- c('x','y','time')
unique(pc_dt2$x, pc_dt1$y)

# dataset 3: use another jittering methods with normal distribution 
# and arbritrary SD of 100

pc_dt3 <- as.data.frame(cbind((pc_dt1$x+rnorm(500*pc_dt1$x, sd=100)), 
                              (pc_dt1$y+rnorm(500*pc_dt1$y, sd=100)),
                              pc_dt1$time))
colnames(pc_dt3) <- c('x','y','time')
unique(pc_dt3$x, pc_dt3$y)

# plot each dataset to see differences
tiff("rdt.sim.tiff", units="in", width=11.69, height=8.27, res=300)
sim.plot(pc_dt1)
dev.off

sim.plot(pc_dt2)
sim.plot(pc_dt3)

tiff("compare.sim.tiff", units="in", height=11.69, width=8.27, res=300)
par(mfrow=c(2,2))
plot(pc_dt1$x, pc_dt1$y, 
     main = "Using centroids of avenues",
     xlab = "X", ylab = "Y",
     pch = 20, col='black', cex = 1, frame = FALSE)
#plot(pc_dt2$x, pc_dt2$y, 
#     main = "Using locations 1 (jitter)",
#     xlab = "X", ylab = "Y",
#     pch = 1, col='red', frame = FALSE)
plot(pc_dt3$x, pc_dt3$y, 
     main = "Using simulated household locations",
     xlab = "X", ylab = "Y",
     pch = 20, col='blue', cex = 0.4, frame = FALSE)
dev.off()

# carry out 2 tau analyses of each simulated dataset

source(here('tau.utils2.R'))

rerun_analyses <- TRUE

## these are the min and max (and midpoints) for the tau-distance windows
r.maxs = seq(100,5000,10)
r.mins = (r.maxs - 50)
r.mids <- (r.maxs + r.mins)/2

## time windows
d.mins<- c(0,5,10,15,20,25,15,1,0)
d.maxs<- c(5,10,15,20,25,30,30,5,2)
d.mids<- (d.maxs+d.mins)/2

## these functions run the primary analyses and save the outputs in the <filename>
(start_time <- Sys.time()) #log start time  
if(rerun_analyses){
  run_main_analysis(point_mat=pc_dt3 %>% as.matrix,
                    filename="pc_dt3.test.rds")
}
end_time <- Sys.time() #log end time  
end_time - start_time #calculate computation time  
#  })

#load data and create moving averages
#pc_dt1_test <- readRDS("C:/R-projects/Uvira-spatial-epid/pc_dt1.rds")
pc_dt1_test <- readRDS("C:/R-projects/Uvira-spatial-epid/pcall.rds")

pc_dt1_test_0_4 <- as.data.frame(cbind(dist = r.mids,
                                       PE = pc_dt1_test[[1]][[2]][["pt.est"]],
                                       LCI = pc_dt1_test[[1]][[2]][["ci.low"]], 
                                       UCI = pc_dt1_test[[1]][[2]][["ci.high"]]))

#pc_dt2_test <- readRDS("C:/R-projects/Uvira-spatial-epid/pc_dt2.rds")

#pc_dt2_test_0_4 <- as.data.frame(cbind(dist = r.mids,
#                                       PE = pc_dt2_test[[1]][[2]][["pt.est"]],
#                                       LCI = pc_dt2_test[[1]][[2]][["ci.low"]], 
#                                       UCI = pc_dt2_test[[1]][[2]][["ci.high"]]))

#pc_dt3_test <- readRDS("C:/R-projects/Uvira-spatial-epid/pc_dt3.rds")
pc_dt3_test <- readRDS("C:/R-projects/Uvira-spatial-epid/pc_dt3.test.rds")

pc_dt3_test_0_4 <- as.data.frame(cbind(dist = r.mids,
                                       PE = pc_dt3_test[[1]][[2]][["pt.est"]],
                                       LCI = pc_dt3_test[[1]][[2]][["ci.low"]], 
                                       UCI = pc_dt3_test[[1]][[2]][["ci.high"]]))

#pruning
pc_dt1_test_0_4_d <- pc_dt1_test_0_4 %>% 
  #filter(pc_dt1_test_0_4$dist >= 95) %>% 
  slice(-c(492:496))
rm(pc_dt1_test_0_4)

#pc_dt2_test_0_4_d <- pc_dt2_test_0_4 %>% 
  #  filter(pc_dt2_test_0_4$dist >= 95) %>% 
#  slice(-c(492:496))
#rm(pc_dt2_test_0_4)

pc_dt3_test_0_4_d <- pc_dt3_test_0_4 %>% 
  #  filter(pc_dt3_test_0_4$dist >= 95) %>% 
  slice(-c(492:496))
rm(pc_dt3_test_0_4)

#create moving averages
pc_dt1_test_0_4_d <- pc_dt1_test_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

#pc_dt2_test_0_4_d <- pc_dt2_test_0_4_d %>% 
#  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
#  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
#  mutate(ma.UCI = slide_mean(UCI, before = 10))

pc_dt3_test_0_4_d <- pc_dt3_test_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

#determine when ma.PE for RR=1 is crossed for each dataset
pc_dt1_test_0_4_d %>% filter(ma.PE < 1 & dist>0) %>% slice(1:3) 
#1645m, 1655m, 1665m*
#pc_dt2_test_0_4_d %>% filter(ma.PE < 1 & dist>0) %>% slice(1:3) 
pc_dt3_test_0_4_d %>% filter(ma.PE < 1 & dist>0) %>% slice(1:3) 
#1795m, 1805m, 1815m*

#determine when ma.LCI for RR=1 is crossed for each dataset
pc_dt1_test_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:3) 
#1085m, 1095m, 1105m
#pc_dt2_test_0_4_d %>% filter(ma.LCI < 1 & dist>0) %>% slice(1:3) 
pc_dt3_test_0_4_d %>% filter(ma.LCI < 1 & dist>0) %>% slice(1:3) 
#1395m, 1405m, 1415m


# Create the 3 plots based on moving averages

## dataset 1: centroid-tagged cases
# RDT-positive cases: plot days 0-4 with moving averages
(dt1.g <- ggplot(pc_dt1_test_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500), limits = c(0,2500)) + ylim(0, 3) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("  Days 0 to 4, Centroid-tagged cases") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "black", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="black", alpha=0.2) +
    geom_vline(xintercept=585, linetype=3, colour="darkgrey", size=0.5) +
    #geom_line(aes(dist, PE), color = "grey", linetype=3, size=0.4) +
    geom_vline(xintercept=1915, linetype=3, colour="black", size=0.5)
)

## dataset 2: jitter function

(dt2.g <- ggplot(pc_dt2_test_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500), limits = c(0,2500)) + ylim(0, 3) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("  Days 0 to 4, Household locations 1") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "red", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="red", alpha=0.2) +
    geom_vline(xintercept=425, linetype=3, colour="darkgrey", size=0.5) +
    #geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
    geom_vline(xintercept=2625, linetype=3, colour="red", size=0.5)
)

## dataset 3: rnorm function

(dt3.g <- ggplot(pc_dt3_test_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500), limits = c(0,2500)) + ylim(0, 3) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("  Days 0 to 4, Simulated household locations") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "cornflowerblue", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="cornflowerblue", alpha=0.2) +
    geom_vline(xintercept=445, linetype=3, colour="darkgrey", size=0.5) +
    #geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
    geom_vline(xintercept=2045, linetype=3, colour="cornflowerblue", size=0.5)
)

# Plot combined graph (moving averages only)

pdf('sidexside.pdf')
plot_grid(dt1.g, #dt2.g, 
          dt3.g, nrow=2)
dev.off() 

# combined graphs

# cutoff points <420m for centroid graph dataset
pc_dt1_test_0_4_d_420 <- pc_dt1_test_0_4_d %>% 
  filter(dist>415)

## point estimates overlap (dt1 and dt3)
tiff("overlap.sim.tiff", units="in", width=11.69, height=8.27, res=300)
(dtc13.g <- ggplot(pc_dt3_test_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500), limits = c(75,2500)) + ylim(0, 4) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    #ggtitle("  Days 0 to 4, Centroid and simulated household location datasets") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 1.5) +
    geom_line(aes(dist, ma.PE), color = "cornflowerblue", linetype=1, size = 1.5) +
    #geom_vline(xintercept=1915, linetype=3, colour="black", size=0.5) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), #color="cornflowerblue", 
                alpha=0.2, fill="cornflowerblue") +
    geom_vline(xintercept=1415, linetype=3, colour="cornflowerblue", size=1) +
    # dt1
    geom_line(data=pc_dt1_test_0_4_d_420, aes(dist, ma.PE), color = "black", linetype=1, size = 1) +
    #geom_vline(xintercept=2045, linetype=3, colour="cornflowerblue", size=0.5) +
    geom_ribbon(data=pc_dt1_test_0_4_d_420, aes(ymin=ma.LCI, ymax=ma.UCI), #color="black", 
                alpha=0.08, fill="black") +
    geom_vline(xintercept=1105, linetype=3, colour="black", size=1)
)
dev.off()

(dtc12.g <- ggplot(pc_dt2_test_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(500,1000,1500,2000,2500,3000), limits = c(400,3000)) + ylim(0, 4) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("  Days 0 to 4, Centroid and household location 2 datasets") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "red", linetype=1) +
    geom_vline(xintercept=2615, linetype=3, colour="red", size=0.5) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), #color="red", 
                alpha=0.2, fill="red"
    ) +
    # dt1
    geom_line(data=pc_dt1_test_0_4_d, aes(dist, ma.PE), color = "black", linetype=2, size=0.5) +
    geom_vline(xintercept=2035, linetype=3, colour="black", size=0.5) +
    geom_ribbon(data=pc_dt1_test_0_4_d, aes(ymin=ma.LCI, ymax=ma.UCI), #color="black", 
                alpha=0.05, fill="black")
  
  
)

# Plot combined graph (moving averages only)
pdf("compare2.pdf")
plot_grid(dtc12.g, dtc13.g, ncol=2)
dev.off()


############## Evaluate similarity between distance-series ###################

#min and max tau moving average estimates
range(pc_dt1_test_0_4_d$ma.PE)
range(pc_dt2_test_0_4_d$ma.PE)
range(pc_dt3_test_0_4_d$ma.PE)
mean(pc_dt1_test_0_4_d$ma.PE)
mean(pc_dt2_test_0_4_d$ma.PE)
mean(pc_dt3_test_0_4_d$ma.PE)


# % differences with dataset1

# difference between MA point estimates that cross 1
((1795-1645)/1645)*100
# difference between MA LCI that cross 1
((1395-1085)/1085)*100

((2615-1905)/1905)*100
((2035-1905)/1905)*100
((415-575)/575)*100
((435-575)/575)*100

# evaluate correlation coefficients

(cor1 <- cor.test(pc_dt1_test_0_4_d$dist, pc_dt1_test_0_4_d$ma.PE,
                  method = "pearson")) #corr coef, -0.752, 95% CI -0.788,-0.711, p-value < 2.2e-16
#(cor2 <- cor.test(pc_dt2_test_0_4_d$dist, pc_dt2_test_0_4_d$ma.PE,
#                  method = "pearson")) #corr coef, -0.912, 95% CI -0.926,-0.90, p-value < 2.2e-16
(cor3 <- cor.test(pc_dt3_test_0_4_d$dist, pc_dt3_test_0_4_d$ma.PE,
                  method = "pearson")) #corr coeft, -0.818, 95% CI -0.845,-0.786, p-value < 2.2e-16
# very similar correlation coefficients (range 75.2% to 91.2%, and significant)

# calculate Euclidean distances between dt2-dt1 and dt3-dt1

## dt1-dt1
sqrt(sum((pc_dt1_test_0_4_d - pc_dt1_test_0_4_d)^2))
## 0

## dt3-dt2
sqrt(sum((pc_dt3_test_0_4_d - pc_dt2_test_0_4_d)^2))
## 14.56478

## dt2-dt1
sqrt(sum((pc_dt2_test_0_4_d - pc_dt1_test_0_4_d)^2))
## 17.16677

## dt3-dt1
sqrt(sum((pc_dt3_test_0_4_d - pc_dt1_test_0_4_d)^2))
## 15.30335



#################### OLD GRAPHS ################################################

(  pc_dt1_test.g <- ggplot(pc_dt1_test_0_4_d, aes(dist, PE)) + 
     geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.15, fill="cornflowerblue") +
     xlim(0, 3000) + 
     ylim(0, 5) +
     geom_hline(yintercept = 1, linetype = "dashed", colour = 'orange2', size = 0.5) +
     labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
     ggtitle("  Days 0 to 4, centroid dataset") +
     theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
           legend.position = "none", text = element_text(size = 10),
           panel.border = element_rect(colour = "black", fill=NA, size=0.5),
           axis.title.y = element_blank(),
           axis.title.x = element_blank()) +
     geom_ma(ma_fun = SMA, n = 10, color = "red")
)

#plot days 0-4, centroid dataset1
(  pc_dt2_test.g <- ggplot(pc_dt2_test_0_4_d, aes(dist, PE)) + 
    geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.15, fill="cornflowerblue") +
    xlim(0, 3000)  
  + ylim(0, 5) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = 'orange2', size = 0.5) +
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("  Days 0 to 4, random location dataset") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    geom_ma(ma_fun = SMA, n = 10, color = "red") 
)

#plot days 0-4, centroid dataset1
(  pc_dt3_test.g <- ggplot(pc_dt3_test_0_4_d, aes(dist, PE)) + 
    geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.15, fill="cornflowerblue") +
    xlim(0, 3000)  
  + ylim(0, 5) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = 'orange2', size = 0.5) +
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("  Days 0 to 4, random location dataset") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    geom_ma(ma_fun = SMA, n = 10, color = "red") 
)

#plot graphs
pc_dt1_test.g
pc_dt2_test.g
pc_dt3_test.g

## combined graph using different moving averages
m = c(10,25,50,100)

for (i in m) {
  
  # Create a variable name to which the plot will be assigned  
  plot_var_name <- str_c(c("ggplot", i), collapse = "_")
  print(plot_var_name)
  
  # Construct the plots
  temp_plot = ggplot(data=pc_dt1_test_0_4_d, aes(dist, PE)) +
    scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 200)) +
    ylim(0, 4) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = 'orange2', size = 0.5) +
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle(i,"metres, Days 0 to 4, centroid with 95% CI (black), HH loc1 with 95% CI (red), HH loc2 with 95% CI (blue), RR=1 (orange)") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "right", text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    geom_ma(ma_fun = SMA, n = i, color = "black", linetype = "dashed") +
    geom_ma(data=pc_dt2_test_0_4_d, inherit.aes = TRUE, 
            ma_fun = SMA, n = i, color = "red", linetype = "solid") +
    geom_ma(data=pc_dt3_test_0_4_d, inherit.aes = TRUE, 
            ma_fun = SMA, n = i, color = "blue", linetype = "solid") +
    geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.15, fill="lightgrey") +
    geom_ribbon(data=pc_dt2_test_0_4_d, aes(ymin=LCI, ymax=UCI), alpha=0.05, fill="red") +
    geom_ribbon(data=pc_dt3_test_0_4_d, aes(ymin=LCI, ymax=UCI), alpha=0.05, fill="blue")
  
  # Assign the plot to plot name
  assign(plot_var_name, temp_plot)
  #ggsave(temp_plot, file=paste0("tau_mvg_avg_", i,".png"), width = 14, height = 10, units = "cm")
}  

# Plot combined graph (moving averages only)
plot_grid(ggplot_10, ggplot_25, ggplot_50, ncol=1)


## two graphs using point estimates and 95% CIs

raw.c.g <- ggplot(data=pc_dt1_test_0_4_d, aes(dist, PE)) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 200)) + ylim(0, 6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = 'orange2', size = 0.5) +
  labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
  ggtitle("Days 0 to 4, centroid with 95% CI (black)") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "right", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  geom_line(size=.5) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.3, fill="grey") #+
#geom_line(data=pc_dt2_test_0_4_d, colour="red", linetype="dashed", size=.5) +
#geom_ribbon(data=pc_dt2_test_0_4_d, aes(ymin=LCI, ymax=UCI), alpha=0.05, fill="red") +
#geom_line(data=pc_dt3_test_0_4_d, colour="blue", linetype="dashed", size=.5) #+
#geom_ribbon(data=pc_dt3_test_0_4_d, aes(ymin=LCI, ymax=UCI), alpha=0.05, fill="blue")

raw.hhloc.g <- ggplot(data=pc_dt1_test_0_4_d, aes(dist, PE)) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 200)) + ylim(0, 6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = 'orange2', size = 0.5) +
  labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
  ggtitle("Days 0 to 4, HH loc1 with 95% CI (red), HH loc2 with 95% CI (blue), RR=1 (orange)") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "right", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  #geom_line(size=.5) +
  #geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.15, fill="lightgrey") +
  geom_line(data=pc_dt2_test_0_4_d, colour="red", size=.5) +
  geom_ribbon(data=pc_dt2_test_0_4_d, aes(ymin=LCI, ymax=UCI), alpha=0.05, fill="red") +
  geom_line(data=pc_dt3_test_0_4_d, colour="blue", size=.5) +
  geom_ribbon(data=pc_dt3_test_0_4_d, aes(ymin=LCI, ymax=UCI), alpha=0.05, fill="blue")

#combined raw graphs
plot_grid(raw.c.g, raw.hhloc.g, nrow=2)

#determine when PE for RR=1 is crossed for each dataset
pc_dt1_test_0_4_d %>% filter(PE < 1) %>% slice(1:5) #145m, 155m, 165m, 175m, 185m
pc_dt2_test_0_4_d %>% filter(PE < 1) %>% slice(1:5) #1335m,1345m, 1355m, 1585m, 1595m
pc_dt3_test_0_4_d %>% filter(PE < 1) %>% slice(1:5) #725m, 865m, 1085m, 1445m, 1455m

#determine when LCI for RR=1 is crossed for each dataset
pc_dt1_test_0_4_d %>% filter(LCI < 1) %>% slice(1:5) #145m, 155m, 165m, 175m, 185m
pc_dt2_test_0_4_d %>% filter(LCI < 1) %>% slice(1:5) #95m, 105m, #115m, 125m, 135m
pc_dt3_test_0_4_d %>% filter(LCI < 1) %>% slice(1:5) #95m, 105m, #115m, 125m, 135m