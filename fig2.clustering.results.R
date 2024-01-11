# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script produces Figure 2 (clustering results)
# Author: R Ratnayake, 2023

# Outcomes and graphs script for local clustering

# RDT-positive cases (main analysis)
result_uvira.rdt.pos.2016 <- readRDS('result_uvira.rdt.pos.2016.rds')
result_uvira.rdt.pos.2017 <- readRDS('result_uvira.rdt.pos.2017.rds')
result_uvira.rdt.pos.2018 <- readRDS('result_uvira.rdt.pos.2018.rds')
result_uvira.rdt.pos.2019 <- readRDS('result_uvira.rdt.pos.2019.rds')
result_uvira.rdt.pos.2020 <- readRDS('result_uvira.rdt.pos.2020.rds')

cluster.rdt.pos.results <- bind_rows(result_uvira.rdt.pos.2016$col,
                                     result_uvira.rdt.pos.2017$col,
                                     result_uvira.rdt.pos.2018$col,
                                     result_uvira.rdt.pos.2019$col,
                                     result_uvira.rdt.pos.2020$col)

cluster.rdt.pos.results <- cluster.rdt.pos.results %>% 
  mutate(START_DATE = lubridate::ymd(START_DATE)) %>% 
  mutate(END_DATE = lubridate::ymd(END_DATE)) %>% 
  #mutate(YEAR=year(ymd(START_DATE))) %>% 
  #relocate(YEAR, .before=CLUSTER) %>% 
  relocate(RADIUS, .after=CLUSTER) %>% 
  relocate(POPULATION, .after=RADIUS) %>%
  relocate(OBSERVED, .after=POPULATION) %>%
  relocate(EXPECTED, .after=OBSERVED) %>% 
  relocate(REL_RISK, .after=EXPECTED) %>%
  relocate(P_VALUE, .after=REL_RISK) %>% 
  relocate(START_DATE, .after=P_VALUE) %>% 
  relocate(END_DATE, .after=START_DATE) %>% 
  dplyr::select(CLUSTER,RADIUS,POPULATION,OBSERVED,EXPECTED,
                REL_RISK,P_VALUE,LLR,START_DATE,END_DATE) %>% 
  mutate(RADIUS2 = RADIUS*1000)

cluster.rdt.pos.results <- cluster.rdt.pos.results %>% 
  mutate(duration = END_DATE - START_DATE)

write.csv(cluster.rdt.pos.results,
          "C:/R-projects/Uvira-spatial-epid/clusterresults_rdtpos.csv", 
          row.names = FALSE)

#mean and range of radius (for P<0.05)
mean(cluster.rdt.pos.results$RADIUS2[cluster.rdt.pos.results$P_VALUE<0.05])
range(cluster.rdt.pos.results$RADIUS2[cluster.rdt.pos.results$P_VALUE<0.05])
mean(cluster.rdt.pos.results$OBSERVED[cluster.rdt.pos.results$P_VALUE<0.05])
range(cluster.rdt.pos.results$OBSERVED[cluster.rdt.pos.results$P_VALUE<0.05])
mean(cluster.rdt.pos.results$duration[cluster.rdt.pos.results$P_VALUE<0.05])
range(cluster.rdt.pos.results$duration[cluster.rdt.pos.results$P_VALUE<0.05])
# 652.0179
# 307.5843 1581.6291
# 19.65385
# 4 48
# 24.76923 days
# 0 58 **corrected to 1 day as it cannot be zero.

#cowplot of all graphs (significant clusters only)
tiff("satscan.rdt.tiff", units="in", width=20, height=10, res=300)
plot_grid(clusters.rdt.pos.2016, clusters.rdt.pos.2017, 
          clusters.rdt.pos.2018, clusters.rdt.pos.2019,
          clusters.rdt.pos.2020, clusters.rdt.pos.2016.2020, 
          ncol=3)
dev.off()

ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\satscan_clusters.rdt.pos.report.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off() 

# Sensitivity analysis of suspected cases

result_uvira.2016 <- readRDS('result_uvira.2016.rds')
result_uvira.2017 <- readRDS('result_uvira.2017.rds')
result_uvira.2018 <- readRDS('result_uvira.2018.rds')
result_uvira.2019 <- readRDS('result_uvira.2019.rds')
result_uvira.2020 <- readRDS('result_uvira.2020.rds')

cluster.results <- bind_rows(result_uvira.2016$col,
                             result_uvira.2017$col,
                             result_uvira.2018$col,
                             result_uvira.2019$col,
                             result_uvira.2020$col)

cluster.results <- cluster.results %>% 
  mutate(START_DATE = lubridate::ymd(START_DATE)) %>% 
  mutate(END_DATE = lubridate::ymd(END_DATE)) %>% 
  #mutate(YEAR=year(ymd(START_DATE))) %>% 
  #relocate(YEAR, .before=CLUSTER) %>% 
  relocate(RADIUS, .after=CLUSTER) %>% 
  relocate(POPULATION, .after=RADIUS) %>%
  relocate(OBSERVED, .after=POPULATION) %>%
  relocate(EXPECTED, .after=OBSERVED) %>% 
  relocate(REL_RISK, .after=EXPECTED) %>%
  relocate(P_VALUE, .after=REL_RISK) %>% 
  relocate(START_DATE, .after=P_VALUE) %>% 
  relocate(END_DATE, .after=START_DATE) %>% 
  dplyr::select(CLUSTER,RADIUS,POPULATION,OBSERVED,EXPECTED,
           REL_RISK,P_VALUE,LLR,START_DATE,END_DATE) %>% 
  mutate(RADIUS2 = RADIUS*1000)

write.csv(cluster.results,"C:/R-projects/Uvira-spatial-epid/clusterresults.csv", row.names = FALSE)

#mean and range of radius
mean(cluster.results$RADIUS2)
range(cluster.results$RADIUS2)
mean(cluster.results$OBSERVED)
range(cluster.results$OBSERVED)

#read in saved clustering files
clusters.2016 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.2016.shp")
clusters.2017 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.2017.shp") 
clusters.2018 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.2018.shp") 
clusters.2019 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.2019.shp") 
clusters.2020 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.2020.shp") 
clusters.2016.2020 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.2016.2020.shp") 
clusters.rdt.pos.2016 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2016.shp") 
clusters.rdt.pos.2017 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2017.shp") 
clusters.rdt.pos.2018 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2018.shp") 
clusters.rdt.pos.2019 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2019.shp") 
clusters.rdt.pos.2020 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2020.shp") 
clusters.rdt.pos.2016.2020 <- read_sf("C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2016.2020.shp") 

library(cowplot)

#cowplot of all graphs (significant clusters only)

tiff("satscan.susp.tiff", units="in", width=20, height=10, res=300)
plot_grid(clusters.2016, clusters.2017, 
          clusters.2018, clusters.2019,
          clusters.2020, clusters.2016.2020, ncol=3)
dev.off()

ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\satscan_clusters.report.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off() 