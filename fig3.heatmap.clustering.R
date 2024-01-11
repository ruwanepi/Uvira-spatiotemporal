# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script produces Figure 3 (clustering heatmap),
# Author: R Ratnayake, 2023

# Epicurve and cluster alarms and combined heatmap of Uvira

#Combine raster with shape files from 2017 to 2020

############################ RDT-positive cases only ########################

#Combine raster with shape files from 2017 to 2020

##Clusters as shp files
clusters.rdt.pos.2016.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2016.shp')
clusters.rdt.pos.2017.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2017.shp')
clusters.rdt.pos.2018.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2018.shp')
clusters.rdt.pos.2019.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2019.shp')
clusters.rdt.pos.2020.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.rdt.pos.2020.shp')

#In QGIS, calculate the proportion of overlap of clusters each year
# using the Overlap Analysis algorithm to calculate the area and % cover
# by which features from an input layer (shpfile of avenues) are overlapped 
# by features from a selection of overlay layer (clusters in 2016-2020)

#st_write(clusters.rdt.pos.2016.sig, "clusters2016.shp", append=FALSE)
#st_write(clusters.rdt.pos.2017.sig, "clusters2017.shp", append=FALSE)
#st_write(clusters.rdt.pos.2018.sig, "clusters2018.shp", append=FALSE)
#st_write(clusters.rdt.pos.2019.sig, "clusters2019.shp", append=FALSE)
#st_write(clusters.rdt.pos.2020.sig, "clusters2020.shp", append=FALSE)
#st_write(uvira.ave, "uvira_ave.shp", append=FALSE)

overlap2 <- read_sf('C:/R-projects/Uvira-spatial-epid/overlap_rdt.shp')

overlap2 <- overlap2 %>% 
  rename(clusters2016_area=clusters.r,
         clusters2016_pc=clusters_1,
         clusters2017_area=clusters_2,
         clusters2017_pc=clusters_3,
         clusters2018_area=clusters_4,
         clusters2018_pc=clusters_5,
         clusters2019_area=clusters_6,
         clusters2019_pc=clusters_7,
         clusters2020_area=clusters_8,
         clusters2020_pc=clusters_9) %>%
  dplyr::select(UGA, clusters2016_pc, clusters2017_pc,
                clusters2018_pc, clusters2019_pc, clusters2020_pc)

#v1: mean proportion of avenue affected (out of 2016-2020)
overlap2$prop <- ((overlap2$clusters2016_pc + 
                    overlap2$clusters2017_pc +
                    overlap2$clusters2018_pc + 
                    overlap2$clusters2019_pc + 
                    overlap2$clusters2020_pc)/5)

#v2: replace proportions affected with binary values ('ever-affected' out of 2016-2020)
overlap2$clusters2016_pc[overlap2$clusters2016_pc>0] <- 1
overlap2$clusters2016_pc[overlap2$clusters2016_pc==0] <- 0
overlap2$clusters2017_pc[overlap2$clusters2017_pc>0] <- 1
overlap2$clusters2017_pc[overlap2$clusters2017_pc==0] <- 0
overlap2$clusters2018_pc[overlap2$clusters2018_pc>0] <- 1
overlap2$clusters2018_pc[overlap2$clusters2018_pc==0] <- 0
overlap2$clusters2019_pc[overlap2$clusters2019_pc>0] <- 1
overlap2$clusters2019_pc[overlap2$clusters2019_pc==0] <- 0
overlap2$clusters2020_pc[overlap2$clusters2020_pc>0] <- 1
overlap2$clusters2020_pc[overlap2$clusters2020_pc==0] <- 0

#combine 'ever-affected' into a total score
overlap2$index <- (overlap2$clusters2016_pc + 
                    overlap2$clusters2017_pc +
                    overlap2$clusters2018_pc + 
                    overlap2$clusters2019_pc + 
                    overlap2$clusters2020_pc)

#rearrange variables
overlap2 <- overlap2 %>% 
  relocate(prop, index, .before = clusters2016_pc) %>% 
  dplyr::select(UGA,prop,index,geometry)

#combine overlap index with map file
uvira.index2 <- cbind(uvira.ave.shape.2, overlap2, by='UGA') %>% 
  dplyr::select(UGA,prop,index,geometry)

#specify the location of the Uvira general hospital
ctc <- data.frame(x=29.1383, y=-3.39206)
ctu <- data.frame(x=29.129915, y=-3.420542)

# This heatmap shows the mean proportion of each avenue as affected by 
# clustering between 2016-2020
(heatmap.rdt <-
ggplot(uvira.index2) +
  geom_sf(aes(fill=index)) +
  scale_fill_gradientn(colours=rev(brewer.pal(11,"YlOrRd")),trans='reverse') +
  labs(fill="No. of clusters by avenue (0-5)") +
  theme(panel.background=element_rect(fill="grey95"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.background = element_rect(fill="grey95", size=.5, linetype="dotted"),
        legend.position=c(.75,.1),
        legend.box = "horizontal",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1, label.position = "bottom"),
         size = guide_legend(title.position="top", title.hjust = 0.5)) +
  geom_point(data=ctc, aes(x=x,y=y), pch = 24, fill='blue', size=2) +
  geom_point(data=ctu, aes(x=x,y=y), pch = 24, fill='blue', size=2) +
  annotation_scale(bar_cols = c("grey60", "white"), location = "bl")
)

ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\CRI_rdt.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off()  

## Combine map with cluster alarm visual

# Epicurve on top
tiff("alarms_recureence.tiff", 
     units="in", width=8.27, height=11.69, res=300)

plot_grid(as_grob(alarms), as_grob(heatmap.rdt), 
          ncol=1, labels="AUTO", label_size=14, 
          label_x = 0, label_y = 0, hjust = -0.5, vjust = -0.5,
          rel_heights=c(0.5,1)) + 
  theme(plot.margin = margin(1,1,1,1,"cm"))

dev.off()

######### map of population density

## match shpfile with ave population
shp.pop <- merge(uvira.ave.shape.2, uvira.ave.pop.2017,
                 by.x=c("UGA"), by.y=c("ave"))

(population.density <-
    ggplot(shp.pop) +
    geom_sf(aes(fill=pop),color = NA)  +
    geom_sf_label(aes(label=pop), size=2) +
    scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlBu"))) +
    theme(panel.background=element_rect(fill="grey95"),
          legend.title=element_text(size=16),
          legend.text=element_text(size=14),
          legend.spacing.x = unit(1.0, 'cm'),
          legend.position=c(.75,.12),
          plot.caption=element_text(size=14)) +
    geom_point(data=hopital, aes(x=x,y=y), fill="red", size=6, shape=25)
)

# For suspected cases

## First, for suspected cases, then for suspected cases

##Clusters as shp files
clusters.2016.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.2016.shp')
clusters.2017.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.2017.shp')
clusters.2018.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.2018.shp')
clusters.2019.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.2019.shp')
clusters.2020.sig <- read_sf('C:/R-projects/Uvira-spatial-epid/clusters.2020.shp')

#In QGIS, calculate the proportion of overlap of clusters each year
# using the Overlap Analysis algorithm to calculate the area and % cover
# by which features from an input layer (shpfile of avenues) are overlapped 
# by features from a selection of overlay layer (clusters in 2016-2020)

#st_write(clusters.2016.sig, "clusters2016.shp", append=FALSE)
#st_write(clusters.2017.sig, "clusters2017.shp", append=FALSE)
#st_write(clusters.2018.sig, "clusters2018.shp", append=FALSE)
#st_write(clusters.2019.sig, "clusters2019.shp", append=FALSE)
#st_write(clusters.2020.sig, "clusters2020.shp", append=FALSE)
#st_write(uvira.ave, "uvira_ave.shp", append=FALSE)

overlap <- read_sf('C:/R-projects/Uvira-spatial-epid/overlap2.shp')

overlap <- overlap %>% 
  rename(clusters2016_area=clusters20,
         clusters2016_pc=clusters_1,
         clusters2017_area=clusters_2,
         clusters2017_pc=clusters_3,
         clusters2018_area=clusters_4,
         clusters2018_pc=clusters_5,
         clusters2019_area=clusters_6,
         clusters2019_pc=clusters_7,
         clusters2020_area=clusters_8,
         clusters2020_pc=clusters_9) %>%
  dplyr::select(UGA, clusters2016_pc, clusters2017_pc,
                clusters2018_pc, clusters2019_pc, clusters2020_pc)

#v1: mean proportion affected of avenue
overlap$prop <- ((overlap$clusters2016_pc + 
                    overlap$clusters2017_pc +
                    overlap$clusters2018_pc + 
                    overlap$clusters2019_pc + 
                    overlap$clusters2020_pc)/5)

#v2: replace proportions affected with binary values
overlap$clusters2016_pc[overlap$clusters2016_pc>0] <- 1
overlap$clusters2016_pc[overlap$clusters2016_pc==0] <- 0
overlap$clusters2017_pc[overlap$clusters2017_pc>0] <- 1
overlap$clusters2017_pc[overlap$clusters2017_pc==0] <- 0
overlap$clusters2018_pc[overlap$clusters2018_pc>0] <- 1
overlap$clusters2018_pc[overlap$clusters2018_pc==0] <- 0
overlap$clusters2019_pc[overlap$clusters2019_pc>0] <- 1
overlap$clusters2019_pc[overlap$clusters2019_pc==0] <- 0
overlap$clusters2020_pc[overlap$clusters2020_pc>0] <- 1
overlap$clusters2020_pc[overlap$clusters2020_pc==0] <- 0
#combine annual values into a total score
overlap$index <- (overlap$clusters2016_pc + 
                    overlap$clusters2017_pc +
                    overlap$clusters2018_pc + 
                    overlap$clusters2019_pc + 
                    overlap$clusters2020_pc)

#rearrange variables
overlap <- overlap %>% 
  relocate(prop, index, .before = clusters2016_pc) %>% 
  dplyr::select(UGA,prop,index,geometry)

#combine overlap index with map file
uvira.index <- cbind(uvira.ave.shape.2, overlap, by='UGA') %>% 
  dplyr::select(UGA,prop,index,geometry)

#produce heatmap
library(RColorBrewer)
pal <- brewer.pal(6, "OrRd")  
class(pal)

plot(uvira.index["prop"],  
     breaks="quantile", nbreaks=6, pal=pal,
     border=NA, bgc='gray95', main=NULL)

plot(uvira.index["index"],  
     breaks="quantile", nbreaks=6, pal=pal,
     border=NA, bgc='gray95', main=NULL)

(heatmap.suspected <-
    ggplot(uvira.index) +
    geom_sf(aes(fill=prop),color = NA) +
    #scale_fill_viridis_c(option = "inferno") +
    scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlBu"))) +
    labs(fill="Cluster recurrence index [%]") +
    theme(panel.background=element_rect(fill="grey95"),
          legend.title=element_text(size=16),
          legend.text=element_text(size=14),
          legend.spacing.x = unit(1.0, 'cm'),
          legend.position=c(.75,.12))
)

ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\CRI.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off()  

#take clustering down to the 1km x 1km pixel
library(viridis)

uvira.index.2 <- cbind(uvira.index, uvira.ave.shape.2.cent, by='UGA') 

uvira.index.2$X <- round(uvira.index.2$X, 3)
uvira.index.2$Y <- round(uvira.index.2$Y, 3)

ggplot(uvira.index.2) + 
  geom_raster(aes(x=X, y=Y, fill=prop)) + 
  coord_fixed(ratio = 1) +
  scale_fill_viridis(direction = -1) +
  theme_bw() 


