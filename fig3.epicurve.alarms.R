# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script produces Figure 3 (clustering heatmap),
# specifically the timing of the alarms.
# Author: R Ratnayake, 2023

# Cluster alarms

## Define cluster dates
cluster.points <- read.csv(file = "C:/R-projects/Uvira-spatial-epid/clusterstartrdtpos.csv")

##convert character to dates
cluster.points$START_DATE <- as.Date(cluster.points$START_DATE, 
                                     format="%m/%d/%Y")
cluster.points$END_DATE <- as.Date(cluster.points$END_DATE, 
                                   format="%m/%d/%Y")

cluster.points <- cluster.points %>% 
  mutate(START_DATE=lubridate::ymd(START_DATE)) %>% 
  mutate(END_DATE=lubridate::ymd(END_DATE))  

# Define sequence of weekly breaks
weekly_breaks_central <- seq.Date(from=floor_date(min(uvira.cholera.new$date_adm, na.rm=T),   "week", week_start = 1), 
                                  to=ceiling_date(max(uvira.cholera.new$date_adm, na.rm=T), "week", week_start = 1), 
                                  by="week") 

#prune uvira.cholera.new to RDT-positive only
uvira.cholera.new.rdt <- uvira.cholera.new  %>% 
  filter(result=="Positive")

tiff("alarms.tiff", units="in", width=20, height=10, res=300)
(alarms<- 
    ggplot(data = uvira.cholera.new.rdt) + 
    geom_histogram(mapping = aes(x = date_adm),
                   breaks = weekly_breaks_central,  
                   closed = "left", fill = "gray80", size=0.3)+ 
    scale_x_date(expand = c(0,0), date_breaks="6 months", 
                 date_labels = "%b %y")+ # date labels format
    scale_y_continuous(expand = c(0,0))+
    theme_minimal()+ theme(axis.title = element_text(size=14),
                           axis.text = element_text(size=11),            
                           axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank()) +    # axis titles in bold
    #labs(x="Week of case presentation", y="Weekly cases") +
    xlab("\nWeek of case presentation") + ylab("Weekly RDT-positive cases\n") +
    #annotate("rect",ymin=0,ymax=50,alpha=0.3,fill="orange",
    #         xmin  = as.Date(cluster.points$START_DATE[1]), 
    #         xmax  = as.Date(cluster.points$END_DATE[1])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[1])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[2])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[3])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[4])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[5])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[6])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[7])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[8])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[9])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[10])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[11])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[12])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[13])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[14])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[15])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[16])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[17])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[18])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[19])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[20])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[21])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[22])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[23])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[24])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[25])) +
    annotate("pointrange",y=0,ymin=0,ymax=10,colour="orangered3",alpha=3,size=.5,
             x  = as.Date(cluster.points$START_DATE[26])))  
dev.off()  