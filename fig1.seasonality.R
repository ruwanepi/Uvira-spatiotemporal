# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script produces Figure 1 (epidemic curve),
# specifically the seasonality assessment graph
# Author: R Ratnayake, 2023

# Seasonality and decomposition analysis

# Decomposition analysis 

## Trim dataset to 2016 to 2020
uvira.ts <- uvira.cholera.new %>%
  rename(date = date_adm) %>% 
  filter(date >= as.Date("2016-01-01") & date <= as.Date("2020-12-31"))  
## 5447 cases total

### now for RDT-positive cases only
uvira.ts.rdtpos <- uvira.rdt.pos %>%
  rename(date = date_adm) %>% 
  filter(date >= as.Date("2016-01-01") & date <= as.Date("2020-12-31")) 
## 1493 cases

## Group by epiweek
uvira.ts <- uvira.ts %>% 
  group_by(epiweek, year) %>% 
  summarize(n=n()) %>% 
  rename(totalcas=n)
# 264 weeks total 

### now for RDT-positive cases only
uvira.ts.rdtpos <- uvira.ts.rdtpos %>% 
  group_by(epiweek, year) %>% 
  summarize(n=n()) %>% 
  rename(totalcas=n)
# 172 weeks total 
  
## Setup the time series analysis
are_duplicated(uvira.ts, index = epiweek) 
dup<-duplicates(uvira.ts, index = epiweek)
rm(dup)
uvira.ts <- uvira.ts %>% distinct(epiweek, .keep_all = TRUE) 
uvira.ts <- tsibble(uvira.ts, index = epiweek)

# fill gaps in time
uvira.ts <- fill_gaps(uvira.ts, .full = TRUE)
# NA found for 2019-W14 (fill in with year=2019, totalcas=0)
uvira.ts[171, 2] = 2019
uvira.ts[171, 3] = 0

## now for RDT-positive cases
are_duplicated(uvira.ts.rdtpos, index = epiweek) 
dup<-duplicates(uvira.ts.rdtpos, index = epiweek)
rm(dup)
uvira.ts.rdtpos <- uvira.ts.rdtpos %>% distinct(epiweek, .keep_all = TRUE) 
uvira.ts.rdtpos <- tsibble(uvira.ts.rdtpos, index = epiweek)

# fill gaps in time
uvira.ts.rdtpos <- fill_gaps(uvira.ts.rdtpos, .full = TRUE)
# NAs found 
# Fill in missing year
uvira.ts.rdtpos$year <- substr(uvira.ts.rdtpos$epiweek, 1, 4)
# There are two options for dealing with NAs: (1) replace with 0,
# or (2) use multiple imputation (MICE)
# Replace with zero method
#uvira.ts.rdtpos["totalcas"][is.na(uvira.ts.rdtpos["totalcas"])] <- 0

#Use multiple imputation
library(mice)
uvira.ts.rdtpos.imp <- mice(uvira.ts.rdtpos, 
                            m=1, 
                            maxit = 50, 
                            method = 'pmm', 
                            seed = 500)

# check that imputed values are there  
uvira.ts.rdtpos.imp$imp$totalcas

# save imputed dataset
uvira.ts.rdtpos.imputed <- complete(uvira.ts.rdtpos.imp)
# change to TS format
uvira.ts.rdtpos.imputed <- tsibble(uvira.ts.rdtpos.imputed, index = epiweek)

#grah rdt-pos time series (non-imputed)
ggplot(uvira.ts.rdtpos, aes(x = epiweek, y = totalcas)) +
  geom_line(linetype = "twodash", colour="steelblue") +
  geom_point(colour="steelblue") +
  labs(y="Suspected cholera cases", x="Epiweek, 2016-2020") +
  #theme(legend.position = "bottom", text = element_text(size = 20)) + 
  ggtitle("Uvira - seasonality, 2016 to 2020") +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.border = element_blank(),
                     text = element_text(size = 16),
                     axis.text = element_text(size = 16))

#grah rdt-pos time series (imputed)
ggplot(uvira.ts.rdtpos.imputed, aes(x = epiweek, y = totalcas)) +
  geom_line(linetype = "twodash", colour="steelblue") +
  geom_point(colour="steelblue") +
  labs(y="Suspected cholera cases", x="Epiweek, 2016-2020") +
  #theme(legend.position = "bottom", text = element_text(size = 20)) + 
  ggtitle("Uvira - seasonality, 2016 to 2020") +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.border = element_blank(),
                     text = element_text(size = 16),
                     axis.text = element_text(size = 16))

#Seasonal decomposition for all suspected cases

# using Seasonal and Trend decomposition using LOESS with a seasonal 
# window of all obs in series and trend lags of 5 (trend not so smooth)
# and robust LOESS regression
season <- uvira.ts %>% 
  model(STL(totalcas ~ trend(window=14) + season(window='periodic'), 
            robust=TRUE)) %>% 
  components() 

autoplot(season)
ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\decomposition.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off()  

# extract components
comp <- uvira.ts %>% 
  model(STL(totalcas ~ trend(window=14) + season(window='periodic'), 
            robust=TRUE)) %>% 
  components()

# check residuals of remainder component for normality (indicating that they are just noise)
par(mfrow = c(1,2))
plot(comp$remainder, col = "azure4", type = "l", lwd = 2, main = 'Remainder', ylab = "")
qqnorm(comp$remainder)
qqline(comp$remainder)

shapiro.test(comp$remainder)
# the Q-Q plot looks heavy-tailed and is not representative of a normal
# distribution of residuals

# insert serial number into comp
comp$serial <- 1:nrow(comp)

#now for RDT-positive cases only
season2 <- uvira.ts.rdtpos %>% 
  model(STL(totalcas ~ trend(window=14) + season(window='periodic'), 
            robust=TRUE)) %>% 
  components() 

autoplot(season2)
ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\decomposition2.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off()  

# extract components - RDT-pos
comp2 <- uvira.ts.rdtpos %>% 
  model(STL(totalcas ~ trend(window=14) + season(window='periodic'), 
            robust=TRUE)) %>% 
  components()

# check residuals of remainder component for normality (indicating that they are just noise)
par(mfrow = c(1,2))
plot(comp2$remainder, col = "azure4", type = "l", lwd = 2, main = 'Remainder', ylab = "")
qqnorm(comp2$remainder)
qqline(comp2$remainder)

shapiro.test(comp2$remainder)
# the Q-Q plot looks heavy-tailed and is not representative of a normal
# distribution of residuals

# insert serial number into comp
comp2$serial <- 1:nrow(comp2)


# for imputed RDT data

season3 <- uvira.ts.rdtpos.imputed %>% 
  model(STL(totalcas ~ trend(window=14) + season(window='periodic'), 
            robust=TRUE)) %>% 
  components() 

autoplot(season3)
ggsave("C:\\R-projects\\Uvira-spatial-epid\\Uvira_graphs\\decomposition3.jpeg", width = 11.69, height = 8.27,
       dpi="print")
dev.off()  

# extract components - RDT-pos
comp3 <- uvira.ts.rdtpos.imputed %>% 
  model(STL(totalcas ~ trend(window=14) + season(window='periodic'), 
            robust=TRUE)) %>% 
  components()

# check residuals of remainder component for normality (indicating that they are just noise)
par(mfrow = c(1,2))
plot(comp3$remainder, col = "azure4", type = "l", lwd = 2, main = 'Remainder', ylab = "")
qqnorm(comp3$remainder)
qqline(comp3$remainder)

shapiro.test(comp3$remainder)
# the Q-Q plot looks heavy-tailed and is not representative of a normal
# distribution of residuals

# insert serial number into comp
comp3$serial <- 1:nrow(comp3)

# plot seasonal component only

library("grid")
library("ggplotify")

# Suspected cases (all years)
(p <- ~plot.ts(x=comp$serial, y=comp$season_year,
        type = "l", lwd = 2, col = "azure4",
        ylab = "Weekly suspected cases", xlab = "Week of case presentation (2016-2020)",
        cex.main=0.8, cex.lab=0.8, cex.axis=0.8))

# RDT-positive cases (all years)
(p2 <- ~plot.ts(x=comp2$serial, y=comp2$season_year,
               type = "l", lwd = 2, col = "azure4",
               ylab = "Weekly RDT-positive cases", 
               xlab = "Week of case presentation (2016-2020)",
               cex.main=0.8, cex.lab=0.8, cex.axis=0.8))

# Suspected cases (single year)
comp <- comp %>% slice(54:105)
  # add sequential (generic) dates
  # defining start date
  date <- as.Date("2020/01/01")
  # defining length of range 
  len <- 52
  # generating range of dates
  date.g <- seq(date, by = "week", length.out = len)
  # add to comp2 dataframe
  comp$date <- date.g

# RDT-positive cases (single year)
comp2 <- comp2 %>% slice(40:91)
# add sequential (generic) dates
  # defining start date
  date <- as.Date("2020/01/01")
  # defining length of range 
  len <- 52
  # generating range of dates
  date.g <- seq(date, by = "week", length.out = len)
  # add to comp2 dataframe
  comp2$date <- date.g

# Imputed RDT-positive cases (single year)
comp3 <- comp3 %>% slice(40:91)
  # add sequential (generic) dates
  # defining start date
  date <- as.Date("2020/01/01")
  # defining length of range 
  len <- 52
  # generating range of dates
  date.g <- seq(date, by = "week", length.out = len)
  # add to comp2 dataframe
  comp3$date <- date.g

  
## ggplot version of imputed RDT-positive graph
(p3 <- ggplot(comp3, aes(date, season_year)) +
  geom_line(colour="#C93312", size=1) +
  xlab("\nMonth of case presentation") +
  ylab("\nSeasonality") +
  ggtitle("Seasonal trend of RDT-positive cases") +
  ylim(-6, 15) +
  scale_x_date(date_labels = "%m",
               date_breaks = "1 month",
               limits = as.Date(c("2020-01-01", "2020-12-23"))) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.9, 0.9),
        text = element_text(size = 8))
)
  
# non-imputed RDT graph
(p2 <- ggplot(comp2, aes(date, season_year)) +
  geom_line(colour="#C93312", size=1) +
  xlab("\nMonth of case presentation") +
  ylab("\nRDT-positive cases") +
  ggtitle("Seasonal trend") +
  ylim(-6, 15) +
  scale_x_date(date_labels = "%m",
               date_breaks = "1 month") +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.9, 0.9),
        text = element_text(size = 9)))
  
  # imputed RDT graph with March start
  
  comp4 <- comp3 %>%
    arrange(factor(date, levels = LETTERS[c(3, 1, 2)]), desc(Res), desc(Pop))
  
  
  (p2 <- ggplot(comp2, aes(date, season_year)) +
      geom_line(colour="#C93312", size=1) +
      xlab("\nMonth of case presentation") +
      ylab("\nRDT-positive cases") +
      ggtitle("Seasonal trend (cholera year)") +
      ylim(-6, 15) +
      scale_x_date(date_labels = "%m",
                   date_breaks = "1 month") +
      theme(axis.line = element_line(colour = "black"),
            #panel.grid.major = element_blank(), 
            #panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = c(0.9, 0.9),
            text = element_text(size = 9)))