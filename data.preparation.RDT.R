# Cleaning script prepares RDT datasets for 2016, 2017, 2018, 2019, 2020
#
# Main products include
# R objects: uvira.cholera.new (cases), uvira.ave.shape.2.cent (shp),
# uvira.ss (SATSCAN version of cases)
# Files on disk: uvira.cases.2016.2020.rds (cases), uvira.ave.shp (shp)
#
# Ruwan Ratnayake, LSHTM, August 2023

# Create version which filters out cases that were not tested
uvira.rdt.pos <- uvira %>% 
  filter(tested=="Tested") %>% 
  filter(result=="Positive")
# 1508 RDT-tested cases retained

# Create a version that has only the data needs for the 3 SATSCAN
# files

uvira.rdt.pos <- as.data.frame(uvira.rdt.pos)
class(uvira.rdt.pos)

uvira.rdt.pos.ss <- uvira.rdt.pos %>% 
  group_by(ave, date_adm, X, Y) %>% 
  summarize(count_by_avedate_adm=n()) %>% 
  rename(cases=count_by_avedate_adm,
         y=X,x=Y,date=date_adm) %>% 
  arrange(date,ave)
# 1436 cases retained

#create serial id
uvira.rdt.pos.ss$id <- seq.int(nrow(uvira.rdt.pos.ss))
uvira.rdt.pos.ss <- uvira.rdt.pos.ss %>% 
  relocate(id, .before = date) %>% 
  relocate(ave, .after = date)

######################################################################
##
## Create files for each year
##
######################################################################

#### 2016 to 2020

#create case file for entire period (2016-2020)
uvira.case.rdt.pos.2016.2020 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2020-12-31"))
write.csv(uvira.case.rdt.pos.2016.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case.rdt.pos-2016-20.csv', row.names=FALSE)
# 1424 rows

#create coordinates file for entire period (2016-2020)
uvira.coord.rdt.pos.2016.2020 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2020-12-31")) 
write.csv(uvira.coord.rdt.pos.2016.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord.rdt.pos-2016-20.csv', row.names=FALSE)
# 1424 rows

#Create pop file for Satscan with locationID
uvira.pop.rdt.pos.2016.2020 <- merge(uvira.ave.pop.2017, uvira.case.rdt.pos.2016.2020, by='ave')
uvira.pop.rdt.pos.2016.2020 <- uvira.pop.rdt.pos.2016.2020 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.rdt.pos.2016.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop.rdt.pos-2016-20.csv', row.names=FALSE)
# 1424 rows

############################ SINGLE YEARS #########################################

#### 2016

#case file for 2016
uvira.case.rdt.pos.2016 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2016-12-31"))
write.csv(uvira.case.rdt.pos.2016,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case.rdt.pos-2016.csv', row.names=FALSE)
# 213 rows

#coordinates file for 2016
uvira.coord.rdt.pos.2016 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2016-12-31")) 
write.csv(uvira.coord.rdt.pos.2016,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord.rdt.pos-2016.csv', row.names=FALSE)
# 213 rows

#pop file for 2016
uvira.pop.rdt.pos.2016 <- merge(uvira.ave.pop.2017, uvira.case.rdt.pos.2016, by='ave')
uvira.pop.rdt.pos.2016 <- uvira.pop.rdt.pos.2016 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.rdt.pos.2016,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop.rdt.pos-2016.csv', row.names=FALSE)
# 213 rows


#### 2017

#case file for 2017
uvira.case.rdt.pos.2017 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2017-01-01") &
           date <= as.Date("2017-12-31"))
write.csv(uvira.case.rdt.pos.2017,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case.rdt.pos-2017.csv', row.names=FALSE)
# 354 rows

#coordinates file for 2017
uvira.coord.rdt.pos.2017 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2017-01-01") &
           date <= as.Date("2017-12-31")) 
write.csv(uvira.coord.rdt.pos.2017,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord.rdt.pos-2017.csv', row.names=FALSE)
# 354 rows

#pop file for 2017
uvira.pop.rdt.pos.2017 <- merge(uvira.ave.pop.2017, uvira.case.rdt.pos.2017, by='ave')
uvira.pop.rdt.pos.2017 <- uvira.pop.rdt.pos.2017 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.rdt.pos.2017,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop.rdt.pos-2017.csv', row.names=FALSE)
# 354 rows

#### 2018

#case file for 2018
uvira.case.rdt.pos.2018 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2018-01-01") &
           date <= as.Date("2018-12-31"))
write.csv(uvira.case.rdt.pos.2018,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case.rdt.pos-2018.csv', row.names=FALSE)
# 224 rows

#coordinates file for 2018
uvira.coord.rdt.pos.2018 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2018-01-01") &
           date <= as.Date("2018-12-31")) 
write.csv(uvira.coord.rdt.pos.2018,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord.rdt.pos-2018.csv', row.names=FALSE)
# 224 rows

#pop file for 2018
uvira.pop.rdt.pos.2018 <- merge(uvira.ave.pop.2017, uvira.case.rdt.pos.2018, by='ave')
uvira.pop.rdt.pos.2018 <- uvira.pop.rdt.pos.2018 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.rdt.pos.2018,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop.rdt.pos-2018.csv', row.names=FALSE)
# 224 rows


#### 2019

#case file for 2019
uvira.case.rdt.pos.2019 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2019-01-01") &
           date <= as.Date("2019-12-31"))
write.csv(uvira.case.rdt.pos.2019,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case.rdt.pos-2019.csv', row.names=FALSE)
# 253 rows

#coordinates file for 2019
uvira.coord.rdt.pos.2019 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2019-01-01") &
           date <= as.Date("2019-12-31")) 
write.csv(uvira.coord.rdt.pos.2019,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord.rdt.pos-2019.csv', row.names=FALSE)
# 253 rows

#pop file for 2019
uvira.pop.rdt.pos.2019 <- merge(uvira.ave.pop.2017, uvira.case.rdt.pos.2019, by='ave')
uvira.pop.rdt.pos.2019 <- uvira.pop.rdt.pos.2019 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.rdt.pos.2019,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop.rdt.pos-2019.csv', row.names=FALSE)
# 253 rows


#### 2020

#case file for 2020
uvira.case.rdt.pos.2020 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2020-01-01") &
           date <= as.Date("2020-12-31"))
write.csv(uvira.case.rdt.pos.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case.rdt.pos-2020.csv', row.names=FALSE)
# 380 rows

#coordinates file for 2020
uvira.coord.rdt.pos.2020 <- uvira.rdt.pos.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2020-01-01") &
           date <= as.Date("2020-12-31")) 
write.csv(uvira.coord.rdt.pos.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord.rdt.pos-2020.csv', row.names=FALSE)
# 380 rows

#pop file for 2020
uvira.pop.rdt.pos.2020 <- merge(uvira.ave.pop.2017, uvira.case.rdt.pos.2020, by='ave')
uvira.pop.rdt.pos.2020 <- uvira.pop.rdt.pos.2020 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.rdt.pos.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop.rdt.pos-2020.csv', row.names=FALSE)
# 380 rows

# Create weekly dataset for 2016-2020

## Add epiweeks, year, and filter for 2016-2020
uvira.rdt.pos.w <- uvira.rdt.pos %>% 
  mutate(epiweek = yearweek(date_adm, week_start = 1)) %>% 
  mutate(year = lubridate::year(date_adm)) %>% 
  filter(year<=2020)

#check 
uvira.rdt.pos.w %>% count(epiweek) %>% mutate(percent = n / sum(n)*100)