# This script interacts with the SaTScan program to evaluate 
# local clustering of cholera for 2016
#
# A spatiotemporal retrospective analysis is used, with space-time
# permutation of discrete scan statistics

#reference case data 
#in format <Serial ID>, <Ave name>, <Date/time>, <# of cases>
cases <- read.csv(file = "C:/R-projects/Uvira-spatial-epid/uvira.case-2016.csv",
                  header=TRUE, fileEncoding="UTF-8-BOM")
cases <- cases %>% dplyr::select(id,cases,date)

class(cases)
str(cases) 
#1208  cases total

#reference coordinates data
#in format <Serial ID>, <Ave name>, <longitude>, <latitude>
coordinates <- read.csv(file = "C:/R-projects/Uvira-spatial-epid/uvira.coord-2016.csv",
                        header = TRUE, fileEncoding="UTF-8-BOM")
coordinates <- coordinates %>% dplyr::select(id,x,y)

class(coordinates)
str(coordinates)
#1208   matched coordinates total

#reference population data
#in format <Ave name>, <dummy year>, <population>
population <- read.csv(file = "C:/R-projects/Uvira-spatial-epid/uvira.ave.pop-2016.csv",
                       header = TRUE, fileEncoding="UTF-8-BOM")
class(population)
str(population)
#confirmed 1208 pop entries

# Load centroids of the neighbourhoods (our spatial unit analysis)
#ss.avenues <- read_sf('C:/R-projects/Uvira-spatial-epid/Uvira_shp/Uvira_avenues.shp')
ss.avenues <- uvira.ave.shape.2.cent
summary(ss.avenues)

# Load SaTScan
ss.local <- "/Program Files/SaTScan" # change location accordingly 

# Clear and load the parameters for SaTScan
invisible(ss.options(reset=TRUE))
ss.options(c("StartDate=2016/01/01","EndDate=2016/12/31"))
ss.options(list(CaseFile="cholera.cas", 
                PrecisionCaseTimes=3, #day
                PopulationFile="cholera.pop",
                CoordinatesFile="cholera.geo", 
                CoordinatesType=1, #latitude/longitude
                AnalysisType=3, #retrospective space-time 
                ModelType=0, #discrete Poisson
                TimeAggregationUnits=3, #day
                TimeAggregationLength=1, 
                OutputShapefiles="y",
                OutputGoogleEarthKML="n",
                ResultsFile='C:/R-projects/Uvira-spatial-epid/',
                StudyPeriodCheckType=0,
                MaxSpatialSizeInPopulationAtRisk=10, #added
                MinimumTemporalClusterSize=7, #at least 7 days
                MaxTemporalSize=60, #at most 60 days
                MaxTemporalSizeInterpretation=1, #interpret as days
                #MinimumCasesInHighRateClusters=2, #unable to do this for some reason
                #RiskLimitHighClusters="n", #unable to do this for some reason
                MonteCarloReps=999,
                ProspectiveStartDate="2016/01/01",
                ReportHierarchicalClusters="y",
                OutputTemporalGraphHTML="y"))

#check parameter file
head(ss.options(),3)

#write SaTScan-readable formatting files 
td = ("C:/R-projects/Uvira-spatial-epid/")
write.ss.prm(td, "cholera")
write.cas(cases, td, "cholera")
write.geo(coordinates, td, "cholera")
write.pop(population, td, "cholera")

#find the location of Satscan program
setwd("/Program Files/SaTScan")

ss.local <- "C:/Program Files/SaTScan"

#Run analysis
result_uvira.2016 <- satscan(td, "cholera", sslocation="/Program Files/SaTScan")

saveRDS(result_uvira.2016, file = "C:/R-projects/Uvira-spatial-epid/result_uvira.2016.rds") 

#Summarize cluster results succinctly
summary(result_uvira.2016)
summary.default(result_uvira.2016)
result_uvira.2016$main

write.table(result_uvira.2016$main, file = "C:/R-projects/Uvira-spatial-epid/retrores2016.txt")
write.csv(result_uvira.2016$col, 'C:/R-projects/Uvira-spatial-epid/retro2016.csv', row.names=FALSE)

#Summarize detailed cluster information and sig clusters
result_uvira.2016$col
result_uvira.2016$col$OBSERVED[result_uvira.2016$col$P_VALUE<0.05]
result_uvira.2016$col$EXPECTED[result_uvira.2016$col$P_VALUE<0.05]
result_uvira.2016$col$POPULATION[result_uvira.2016$col$P_VALUE<0.05]
result_uvira.2016$col$RADIUS[result_uvira.2016$col$P_VALUE<0.05]*1000
result_uvira.2016$col$REL_RISK[result_uvira.2016$col$P_VALUE<0.05]
result_uvira.2016$col$P_VALUE[result_uvira.2016$col$P_VALUE<0.05]

#Create an output of locations belonging to a cluster
cluster.case.loc <- result_uvira.2016$gis 
#remove clusters with non signifcant p-value
cluster.case.loc <- cluster.case.loc[cluster.case.loc$P_VALUE <=0.049, ]  
#assess #s of cases per cluster
cluster.case.loc %>% group_by(CLUSTER) %>% summarise(n = n())
##note that where mutiple cases are flagged from one location they
##have been combined into one location. So, there are more cases
##than indicated in the output file as follows:

#match dates and lat/long to cluster.cases file
##bind cases and coordinates file
case.coord <- merge(cases, coordinates, by="id")
##bind case.coord with cluster file
cluster.case.loc <- merge(case.coord, cluster.case.loc, 
                          by.x ='id', by.y ='LOC_ID')
##clean up file (remove extraneous variables) and save as RDS
cluster.cases.2016 <- cluster.case.loc %>% 
  dplyr::select(id, cases, date, x, y, CLUSTER, P_VALUE, LOC_LAT, LOC_LONG) %>% 
  relocate(date, .after = id) %>% relocate(CLUSTER, .after = date) %>% 
  rename(cluster=CLUSTER, pval=P_VALUE, cluster_x=LOC_LAT, cluster_y=LOC_LONG)

saveRDS(cluster.cases.2016, file = "C:/R-projects/Uvira-spatial-epid/cluster.cases.2016.rds")         

#map the clusters
#pull in cluster shapefile from SaTScan run in Windows and plot
#clusters on Uvira map

#significant clusters only
shapefile(result_uvira.2016$shapeclust[result_uvira.2016$shapeclust$P_VALUE <=0.049, ], 
          filename='C:/R-projects/Uvira-spatial-epid/clusters.2016.shp', overwrite=TRUE)

clusters.2016.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.2016.shp')

(clusters.2016 <- ggplot() + 
    geom_sf(data=uvira.ave.shape) +
    geom_sf(data=clusters.2016.sig, colour = "darkorange1", 
            size = 0.9, fill = alpha("#2C77BF", .5)) +
    geom_sf_text(data=clusters.2016.sig, 
                 aes(label = OBSERVED), size=3, colour="white") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size = 10)) +
    xlab("Latitude") + ylab("Longitude") +
    labs(subtitle = "A: 2016")
)

ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\satscan_clusters.2016.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off()  

saveRDS(clusters.2016, file = "C:/R-projects/Uvira-spatial-epid/clusters.2016.rds") 

#To save file
sink(file='C:/R-projects/Uvira-spatial-epid/stss_result_2016.txt')
satscan(td, "cholera", sslocation="/Program Files/SaTScan")
sink() 
