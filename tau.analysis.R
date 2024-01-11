# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script carries out the analyses for tau.
# Author: R Ratnayake, 2023, based on A Azman, 2018 (code: https://osf.io/4fsnc/)

# Analyses for tau

source(here('tau.utils.R'))

rerun_analyses <- TRUE

#set working directory
setwd("C:/R-projects/Uvira-spatial-epid/")

#load data
sc2020 <- load_suspected_cases2020()
pc2020 <- load_RDTpos_cases2020()
sc2019 <- load_suspected_cases2019()
pc2019 <- load_RDTpos_cases2019()
pc2018 <- load_RDTpos_cases2018()
pc2017 <- load_RDTpos_cases2017()
pc2016 <- load_RDTpos_cases2016()
pcall <- load_RDTpos_cases20162020()

## these are the min and max (and midpoints) for the tau-distance windows
r.maxs = seq(100,5000,10)
r.mins = (r.maxs - 50)
r.mids <- (r.maxs + r.mins)/2

## time windows
d.mins<- c(0,5,10,15,20,25,15,1,0)
d.maxs<- c(5,10,15,20,25,30,30,5,2)
d.mids<- (d.maxs+d.mins)/2

## these functions run the primary analyses and save the outputs in the <filename>
#profvis({
(start_time <- Sys.time()) #log start time  
if(rerun_analyses){
  run_main_analysis(point_mat=pcall %>% as.matrix,
                    filename="pcall.rds")
}
end_time <- Sys.time() #log end time  
end_time - start_time #calculate computation time  
#  })


# Prepare data and moving averages for graphs 

## RDT-positive cases (whole time series, 2016-2020)

pc_all1000 <- readRDS("C:/R-projects/Uvira-spatial-epid/pcall.rds")

pc_all1000_0_4 <- as.data.frame(cbind(dist = r.mids,
                                      PE = pc_all1000[[1]][[2]][["pt.est"]],
                                      LCI = pc_all1000[[1]][[2]][["ci.low"]], 
                                      UCI = pc_all1000[[1]][[2]][["ci.high"]]))

pc_all1000_1_4 <- as.data.frame(cbind(dist = r.mids,
                                      PE = pc_all1000[[8]][[2]][["pt.est"]],
                                      LCI = pc_all1000[[8]][[2]][["ci.low"]], 
                                      UCI = pc_all1000[[8]][[2]][["ci.high"]]))

#pruning
pc_all1000_0_4_d <- pc_all1000_0_4 %>% 
  filter(pc_all1000_0_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc_all1000_0_4)

pc_all1000_1_4_d <- pc_all1000_1_4 %>% 
  filter(pc_all1000_1_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc_all1000_1_4)

# Create moving averages for PE, LCI, UCI
pc_all1000_0_4_d <- pc_all1000_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

pc_all1000_1_4_d <- pc_all1000_1_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

## RDT-positive cases (2020)

pc_main1000 <- readRDS("C:/R-projects/Uvira-spatial-epid/pc1000.rds")

pc_main1000_0_4 <- as.data.frame(cbind(dist = r.mids,
                                       PE = pc_main1000[[1]][[2]][["pt.est"]],
                                       LCI = pc_main1000[[1]][[2]][["ci.low"]], 
                                       UCI = pc_main1000[[1]][[2]][["ci.high"]]))

pc_main1000_1_4 <- as.data.frame(cbind(dist = r.mids,
                                       PE = pc_main1000[[8]][[2]][["pt.est"]],
                                       LCI = pc_main1000[[8]][[2]][["ci.low"]], 
                                       UCI = pc_main1000[[8]][[2]][["ci.high"]]))

#pruning
pc_main1000_0_4_d <- pc_main1000_0_4 %>% 
  filter(pc_main1000_0_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc_main1000_0_4)

pc_main1000_1_4_d <- pc_main1000_1_4 %>% 
  filter(pc_main1000_1_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc_main1000_1_4)

# Create moving averages for PE, LCI, UCI
pc_main1000_0_4_d <- pc_main1000_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

pc_main1000_1_4_d <- pc_main1000_1_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

## Suspected cases

sc_main1000 <- readRDS("C:/R-projects/Uvira-spatial-epid/sc1000.rds")

sc_main1000_0_4 <- as.data.frame(cbind(dist = r.mids,
                                       PE = sc_main1000[[1]][[2]][["pt.est"]],
                                       LCI = sc_main1000[[1]][[2]][["ci.low"]], 
                                       UCI = sc_main1000[[1]][[2]][["ci.high"]]))

sc_main1000_1_4 <- as.data.frame(cbind(dist = r.mids,
                                       PE = sc_main1000[[8]][[2]][["pt.est"]],
                                       LCI = sc_main1000[[8]][[2]][["ci.low"]], 
                                       UCI = sc_main1000[[8]][[2]][["ci.high"]]))

#pruning
sc_main1000_0_4_d <- sc_main1000_0_4 %>% 
  filter(dist >= 95) %>% 
  slice(-c(490:492))
rm(sc_main1000_0_4)

sc_main1000_1_4_d <- sc_main1000_1_4 %>% 
  filter(dist >= 95) %>% 
  slice(-c(490:492))
rm(sc_main1000_1_4)

# create moving averages for PE, LCI, UCI
sc_main1000_0_4_d <- sc_main1000_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

sc_main1000_1_4_d <- sc_main1000_1_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

## RDT-positive cases (whole time series, 2016-2020)
#determine when ma.LCI for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_all1000_0_4_d %>% filter(ma.LCI <= 1 & dist>400) %>% slice(1:3) 
# 1085, 1095, **1105
#determine when ma.PE for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_all1000_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# 1645, 1655, 1665**
#determine when ma.LCI for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_all1000_1_4_d %>% filter(ma.LCI <= 1 & dist>400) %>% slice(1:3) 
# 1085, 1095, **1105
#determine when ma.PE for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_all1000_1_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:4) 
# 1635, 1645, 1655**

## RDT-positive cases (2020 alone)
#determine when ma.LCI for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_main1000_0_4_d %>% filter(ma.LCI <= 1 & dist>400) %>% slice(1:3) 
# 565, 575, **585
#determine when ma.PE for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_main1000_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# 1895, 1905, 1915*
#determine when ma.LCI for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_main1000_1_4_d %>% filter(ma.LCI <= 1 & dist>400) %>% slice(1:3) 
# 405, 415, 425**
#determine when ma.PE for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
pc_main1000_1_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# 1895, 1905, 1915**

## Suspected cases, whole period
#determine when ma.LCI for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
sc_main1000_0_4_d %>% filter(ma.LCI <= 1 & dist>400) %>% slice(1:4) 
# 625, 635, **645 
#determine when ma.PE for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
sc_main1000_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# 2055, 2065, 2075**
#determine when ma.LCI for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
sc_main1000_1_4_d %>% filter(ma.LCI <= 1 & dist>400) %>% slice(1:5) 
#determine when ma.PE for RR=1 is crossed 3x consecutively (i.e. >=30m consecutive)
sc_main1000_1_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:4) 
# 2055, 2065, 2075**

# Plot two graphs for entire time series (2016-2020)

# RDT-positive cases: plot days 0-4 with moving averages
(pc.all.04 <- ggplot(pc_all1000_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(500,1000,1500,2000,2500), limits = c(400,2500)) + ylim(0, 3) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("RDT-positive cases (2016-2020)", subtitle="  Days 0 to 4") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face='bold')) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
    geom_vline(xintercept=1105, linetype=3, colour="orange2", size=0.5) +
    geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
    geom_vline(xintercept=1665, linetype=3, colour="orange2", size=0.5)
)

tiff("tau_pc04_all.tiff", units="in", width=11.69, height=8.27, res=300)
pc.all.04
dev.off()

# RDT-positive cases: plot days 1-4 with moving averages
(pc.all.14 <- ggplot(pc_all1000_1_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(500,1000,1500,2000,2500), limits = c(400,2500)) + ylim(0, 3) + 
    #labs(x = "Distance (metres)", y = "Relative risk of cholera (tau)") +
    ggtitle(" ", subtitle="Days 1 to 4") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
    geom_vline(xintercept=1105, linetype=3, colour="orange2", size=0.5) +
    geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
    geom_vline(xintercept=1655, linetype=3, colour="orange2", size=0.5)
)

tiff("tau_pc14_all.tiff", units="in", width=11.69, height=8.27, res=300)
pc.all.14
dev.off()

# 2020 only

# RDT-positive cases: plot days 0-4 with moving averages
(pc04 <- ggplot(pc_main1000_0_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(500,1000,1500,2000,2500), limits = c(400,2500)) + ylim(0, 3) + 
    labs(x = "distance (metres)", y = "Relative risk of cholera (tau)") + 
    ggtitle("RDT-positive cases (2020)", subtitle="  Days 0 to 4") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face='bold')) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
    geom_vline(xintercept=585, linetype=3, colour="orange2", size=0.5) +
    geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
    geom_vline(xintercept=1915, linetype=3, colour="orange2", size=0.5)
)

tiff("tau_pc04.tiff", units="in", width=11.69, height=8.27, res=300)
pc04
dev.off()

# RDT-positive cases: plot days 1-4 with moving averages
(pc14 <- ggplot(pc_main1000_1_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(500,1000,1500,2000,2500), limits = c(400,2500)) + ylim(0, 3) + 
    #labs(x = "Distance (metres)", y = "Relative risk of cholera (tau)") +
    ggtitle(" ", subtitle="Days 1 to 4") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
    geom_vline(xintercept=425, linetype=3, colour="orange2", size=0.5) +
    geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
    geom_vline(xintercept=1915, linetype=3, colour="orange2", size=0.5)
)

########################################################################

# Graphs of suspected cases

# Suspected cases: plot days 0-4 with moving averages

(sc04 <- ggplot(sc_main1000_0_4_d, aes(dist, ma.PE)) + 
   scale_x_continuous(breaks = c(500,1000,1500,2000,2500), limits = c(400,2500)) + ylim(0, 3) + 
   #labs(#x = "distance (metres)", y = "Relative risk of cholera (tau)") +
   ggtitle("Suspected cases (2020)", subtitle="  Days 0 to 4") +
   theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
         legend.position = "none", text = element_text(size = 11),
         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         plot.title = element_text(hjust = 0.5, face='bold')) +
   geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
   geom_line(aes(dist, ma.PE), color = "cornflowerblue", linetype=1) +
   geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="cornflowerblue", alpha=0.2) +
   geom_vline(xintercept=645, linetype=3, colour="cornflowerblue", size=0.5) +
   geom_line(aes(dist, PE), color = "cornflowerblue", linetype=3, size=0.4) +
   geom_vline(xintercept=2075, linetype=3, colour="cornflowerblue", size=0.5)
)

#plot days 1-4

(sc14 <- ggplot(sc_main1000_1_4_d, aes(dist, ma.PE)) + 
    scale_x_continuous(breaks = c(500,1000,1500,2000,2500), limits = c(400,2500)) + ylim(0, 3) + 
    #labs(x = "Distance (metres)"#, y = "Relative risk of cholera (tau)") +
    ggtitle(" ", subtitle="Days 1 to 4") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none", text = element_text(size = 11),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()
    ) +
    geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
    geom_line(aes(dist, ma.PE), color = "cornflowerblue", linetype=1) +
    geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="cornflowerblue", alpha=0.2) +
    geom_vline(xintercept=635, linetype=3, colour="darkgrey") +
    geom_vline(xintercept=1155, linetype=3, colour="cornflowerblue") +
    geom_line(aes(dist, PE), color = "cornflowerblue", linetype=3, size=0.4) +
    geom_vline(xintercept=2075, linetype=3, colour="cornflowerblue")
)


#create common x and y labels

y.grob <- textGrob("Relative risk of cholera (tau)", gp=gpar(fontsize=11), rot=90)
x.grob <- textGrob("Distance (metres)", gp=gpar(fontsize=11))

## Combined graphs for 2020 only
tiff("tau.tiff", units="in", width=11.69, height=8.27, res=300)
plot <- plot_grid(pc.all.04, pc04, sc04,   
                  pc.all.14, pc14, sc14,   
                  labels = c('A','C','E',
                             'B', 'D', 'F'), ncol = 3)
#plot <- plot_grid(pc04, sc04, pc.all.04,  
#                  pc14, sc14, pc.all.14,  
#                  labels = c('A','C','E','B', 'D', 'F'), ncol = 3)
grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
dev.off()

tiff("tau_all.tiff", units="in", width=11.69, height=8.27, res=300)
plot <- plot_grid(pc.all.04, pc.all.14, 
                  labels = c('A','B'), ncol = 2)
grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
dev.off()



#### 5-panel vertical plot with RDT+ results for 2016-2020

#for loop to load and prune annual data
for (i in 2016:2019) {
  filename <- paste0("pc",i)
  wd <- paste0("pc", i, ".rds")
  assign(filename, readRDS(wd))
}

pc2016_0_4 <- as.data.frame(cbind(dist = r.mids,
                                  PE = pc2016[[1]][[2]][["pt.est"]],
                                  LCI = pc2016[[1]][[2]][["ci.low"]], 
                                  UCI = pc2016[[1]][[2]][["ci.high"]]))
pc2017_0_4 <- as.data.frame(cbind(dist = r.mids,
                                  PE = pc2017[[1]][[2]][["pt.est"]],
                                  LCI = pc2017[[1]][[2]][["ci.low"]], 
                                  UCI = pc2017[[1]][[2]][["ci.high"]]))
pc2018_0_4 <- as.data.frame(cbind(dist = r.mids,
                                  PE = pc2018[[1]][[2]][["pt.est"]],
                                  LCI = pc2018[[1]][[2]][["ci.low"]], 
                                  UCI = pc2018[[1]][[2]][["ci.high"]]))
pc2019_0_4 <- as.data.frame(cbind(dist = r.mids,
                                  PE = pc2019[[1]][[2]][["pt.est"]],
                                  LCI = pc2019[[1]][[2]][["ci.low"]], 
                                  UCI = pc2019[[1]][[2]][["ci.high"]]))
#pruning

pc2016_0_4_d <- pc2016_0_4 %>% 
  filter(pc2016_0_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc2016_0_4)

pc2017_0_4_d <- pc2017_0_4 %>% 
  filter(pc2017_0_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc2017_0_4)

pc2018_0_4_d <- pc2018_0_4 %>% 
  filter(pc2018_0_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc2018_0_4)

pc2019_0_4_d <- pc2019_0_4 %>% 
  filter(pc2019_0_4$dist >= 95) %>% 
  slice(-c(490,491,492))
rm(pc2019_0_4)

# Create moving averages for PE, LCI, UCI
pc2016_0_4_d <- pc2016_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))
pc2017_0_4_d <- pc2017_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))
pc2018_0_4_d <- pc2018_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))
pc2019_0_4_d <- pc2019_0_4_d %>% 
  mutate(ma.PE = slide_mean(PE, before = 10)) %>% 
  mutate(ma.LCI = slide_mean(LCI, before = 10)) %>% 
  mutate(ma.UCI = slide_mean(UCI, before = 10))

#thresholds

#determine when ma.LCI for RR=1 is crossed for each dataset
pc2016_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:3) 
# ***405, 415, ***425
pc2017_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:5) 
# 405, ***855, 865, ***875
pc2018_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:3) 
# ***405, 415, ***425
pc2019_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:3) 
# ***405, 415, ***425
pc_main1000_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:3) 
# ***565, 575, ***585
pc_all1000_0_4_d %>% filter(ma.LCI < 1 & dist>400) %>% slice(1:3) 
# 1085, 1095, 1105

#determine when ma.PE for RR=1 is crossed for each dataset
pc2016_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# ***1355, 1365, ***1375
pc2017_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# ***1465, 1475, ***1485
pc2018_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# ***1135, 1145, ***1155
pc2019_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:4) 
# 915, ***1105, 1115, 1125***
pc_main1000_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# ***1895, 1905, ***1915
pc_all1000_0_4_d %>% filter(ma.PE < 1 & dist>400) %>% slice(1:3) 
# 1645, 1655, 1665

#graphs

g2016 <- ggplot(pc2016_0_4_d, aes(dist, ma.PE)) + 
  scale_x_continuous(breaks = c(500,1000,1500,2000), limits = c(400,2000)) + ylim(0, 4) + 
  labs(#x = "distance (metres)", 
    y = "Relative risk (tau)") + 
  ggtitle("  Days 0 to 4, RDT-positive cases, 2016") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
  geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
  geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
  #geom_vline(xintercept=415, linetype=3, colour="orange2", size=0.5) +
  geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
  geom_vline(xintercept=1365, linetype=3, colour="red3", size=0.5)

g2017 <- ggplot(pc2017_0_4_d, aes(dist, ma.PE)) + 
  scale_x_continuous(breaks = c(500,1000,1500,2000), limits = c(400,2000)) + ylim(0, 4) + 
  labs(#x = "distance (metres)", 
    y = "Relative risk (tau)") + 
  ggtitle("  Days 0 to 4, RDT-positive cases, 2017") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
  geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
  geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
  #geom_vline(xintercept=865, linetype=3, colour="orange2", size=0.5) +
  geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
  geom_vline(xintercept=1475, linetype=3, colour="red3", size=0.5)

g2018 <- ggplot(pc2018_0_4_d, aes(dist, ma.PE)) + 
  scale_x_continuous(breaks = c(500,1000,1500,2000), limits = c(400,2000)) + ylim(0, 4) + 
  labs(#x = "distance (metres)", 
    y = "Relative risk (tau)") + 
  ggtitle("  Days 0 to 4, RDT-positive cases, 2018") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
  geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
  geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
  #geom_vline(xintercept=415, linetype=3, colour="orange2", size=0.5) +
  geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
  geom_vline(xintercept=1145, linetype=3, colour="red3", size=0.5)

g2019 <- ggplot(pc2019_0_4_d, aes(dist, ma.PE)) + 
  scale_x_continuous(breaks = c(500,1000,1500,2000), limits = c(400,2000)) + ylim(0, 4) + 
  labs(#x = "distance (metres)", 
    y = "Relative risk (tau)") + 
  ggtitle("  Days 0 to 4, RDT-positive cases, 2019") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
  geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
  geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
  #geom_vline(xintercept=415, linetype=3, colour="orange2", size=0.5) +
  geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
  geom_vline(xintercept=1115, linetype=3, colour="red3", size=0.5)

g2020 <- ggplot(pc_main1000_0_4_d, aes(dist, ma.PE)) + 
  scale_x_continuous(breaks = c(500,1000,1500,2000), limits = c(400,2000)) + ylim(0, 4) + 
  labs(x = "distance (metres)", y = "Relative risk (tau)") + 
  ggtitle("  Days 0 to 4, RDT-positive cases, 2020") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 0.5) +
  geom_line(aes(dist, ma.PE), color = "orange2", linetype=1) +
  geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange2", alpha=0.2) +
  #geom_vline(xintercept=575, linetype=3, colour="orange2", size=0.5) +
  geom_line(aes(dist, PE), color = "orange2", linetype=3, size=0.4) +
  geom_vline(xintercept=1905, linetype=3, colour="red3", size=0.5)

#produce an overlapping graph

tiff("tau.annual.tiff", units="in", width=11.69, height=8.27, res=300)
ggplot(pc_main1000_0_4_d, aes(dist, ma.PE)) +
  scale_x_continuous(breaks = c(500,1000,1500,2000), limits = c(400,2000)) + ylim(0,5) + 
  labs(x = "Distance (metres)", y = "Relative risk of cholera (tau)") +
  #ggtitle("  Days 0 to 4, RDT-positive cases, annual estimates for 2016-2020") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none", text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.7)) +
  geom_hline(yintercept = 1, linetype = 2, colour = 'grey', size = 1) +
  #2020
  geom_line(aes(dist, ma.PE), color = "red4", linetype=1, size=1) +
  geom_ribbon(aes(ymin=ma.LCI, ymax=ma.UCI), fill="red4", alpha=0.07) +
  #geom_vline(xintercept=585, linetype=3, colour="red4", size=1) +
  #geom_vline(xintercept=1915, linetype=3, colour="red4", size=1) +
  #2019
  geom_line(data=pc2019_0_4_d, aes(dist, ma.PE), color = "cornflowerblue", linetype=1, size=1) +
  geom_ribbon(data=pc2019_0_4_d, aes(ymin=ma.LCI, ymax=ma.UCI), fill="cornflowerblue", alpha=0.07) +
  #geom_vline(xintercept=425, linetype=3, colour="cornflowerblue", size=1) +
  #geom_vline(xintercept=1125, linetype=3, colour="cornflowerblue", size=1) +
  #2018
  geom_line(data=pc2018_0_4_d, aes(dist, ma.PE), color = "green3", linetype=1, size=1) +
  geom_ribbon(data=pc2018_0_4_d, aes(ymin=ma.LCI, ymax=ma.UCI), fill="orange3", alpha=0.07) +
  #geom_vline(xintercept=425, linetype=3, colour="green3", size=1) +
  #geom_vline(xintercept=1155, linetype=3, colour="green3", size=1) +
  #2017
  geom_line(data=pc2017_0_4_d, aes(dist, ma.PE), color = "orange3", linetype=1, size=1) +
  geom_ribbon(data=pc2017_0_4_d, aes(ymin=ma.LCI, ymax=ma.UCI), fill="yellow3", alpha=0.07) +
  #geom_vline(xintercept=875, linetype=3, colour="orange3", size=0.5) +
  #geom_vline(xintercept=1485, linetype=3, colour="orange3", size=1) +
  #2016
  geom_line(data=pc2016_0_4_d, aes(dist, ma.PE), color = "purple3", linetype=1, size=1) +
  geom_ribbon(data=pc2016_0_4_d, aes(ymin=ma.LCI, ymax=ma.UCI), fill="purple3", alpha=0.07) +
  #geom_vline(xintercept=425, linetype=3, colour="purple", size=1) +
  #geom_vline(xintercept=1375, linetype=3, colour="purple3", size=1) 
  geom_line(data=pc_all1000_0_4_d, aes(dist, ma.PE), color = "black", linetype=1, size=1.5) +
  geom_ribbon(data=pc_all1000_0_4_d, aes(ymin=ma.LCI, ymax=ma.UCI), fill="black", alpha=0.1)
dev.off()








