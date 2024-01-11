# Cleaning script prepares datasets for 2016, 2017, 2018, 2019, 2020
#
# Ruwan Ratnayake, LSHTM, August 2023

# Load packages (also look at source.R for description of key analysis scripts)
source('source.R')

# Import cholera surveillance data into R 
uvira.cholera <- read.csv("Uvira_data_2009-2021.csv")

# Put dates in correct format, add epiweeks, and add year
uvira.cholera <- uvira.cholera %>% 
  mutate(dt_case_admission = lubridate::dmy(dt_case_admission)) %>% 
  mutate(dt_case_end = lubridate::dmy(dt_case_end)) %>% 
  mutate(epiweek = yearweek(dt_case_admission, week_start = 1))  

# Extract year 
uvira.cholera <- uvira.cholera %>% 
  mutate(year = lubridate::year(uvira.cholera$dt_case_admission))

# Check whether RDT was done and % by result
prop.test <- uvira.cholera %>% 
  group_by(cat_case_result, year) %>% summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))

## Reorder, rename, create new variables for death and severity, and 
## (otherwise) drop extraneous variables
uvira.cholera <- uvira.cholera %>% 
  dplyr::select(id_case, dt_case_admission, amt_case_age, cat_case_age,
         ind_case_sex, cat_case_dehydration, cat_case_end_status,
         cat_case_quarter, cat_case_airedesante, cat_case_avenue, 
         cat_case_test, cat_case_result, epiweek, year) %>% 
  rename(id = id_case, date_adm = dt_case_admission,
         age = amt_case_age, agecat = cat_case_age, 
         sex = ind_case_sex, dehyd = cat_case_dehydration,
         outcome = cat_case_end_status, quartier = cat_case_quarter,  
         ads = cat_case_airedesante, ave = cat_case_avenue, tested = 
           cat_case_test, result = cat_case_result) %>% 
  mutate(death = ifelse(outcome=="Deceased", "1", "0")) %>% 
  mutate(severe = ifelse(dehyd=="Severe", "1",
                         ifelse(dehyd=="Moderate", "2", 
                                ifelse(dehyd=="Simple", "3", "0")))) %>% 
  relocate(year, .after = id, ) %>% 
  relocate(epiweek, .after = year) %>% 
  relocate(death, .after = outcome) %>% 
  relocate(severe, .after = dehyd)

## Explore missing in data and change character to factor class
skim(uvira.cholera) #29.3%, 98%, 75.4% missing age, dehyd, result
uvira.cholera$agecat <- as.factor(uvira.cholera$agecat)  
uvira.cholera$sex <- as.factor(uvira.cholera$sex)  
uvira.cholera$dehyd <- as.factor(uvira.cholera$dehyd) 

## Sort the data
uvira.cholera <- uvira.cholera %>% arrange(uvira.cholera, date_adm, ave)

## Remove any cases with unidentifiable ave
uvira.cholera <- uvira.cholera %>% 
  filter(ave != 'autre') %>% 
  filter(ave != 'Hors_Uvira') %>% 
  filter(ave != 'hors_zone') %>% 
  filter(ave != 'Prison Centrale') %>% 
  filter(ave != 'Uvira_UG_a_incertaine')

# Import shapefile with GIS coordinates
uvira.ave.shape <- 
  here("C:/R-projects/Uvira-spatial-epid/Uvira_shp/Uvira_avenues.shp") %>% 
  st_read()

# Investigate linking variable
uga <- uvira.ave.shape %>% 
  group_by(UGA) %>%
  dplyr::select(UGA) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))

uvira.cholera.new <- uvira.cholera %>% 
  mutate(ave = replace(ave, ave=='av.alpha', 'av.alpha.sg')) %>% 
  mutate(ave = replace(ave, ave=='av.apollo_1.mul', 'av.apollo')) %>% 
  mutate(ave = replace(ave, ave=='av.apollo_2.mul', 'av.apollo')) %>% 
  mutate(ave = replace(ave, ave=='av.azuhuri', 'av.azuhuri.kas')) %>% 
  mutate(ave = replace(ave, ave=='av.budota', 'av.budota.kas')) %>% 
  mutate(ave = replace(ave, ave=='av.bukavu', 'av.bukavu.kas')) %>% 
  mutate(ave = replace(ave, ave=='av.centre_commercial', 'av.centre_commercial.kal')) %>%
  mutate(ave = replace(ave, ave=='av.cinq_chantiers', 'av.cinq_chantiers.kak')) %>%
  mutate(ave = replace(ave, ave=='av.conforti', 'av.conforti.kas')) %>%
  mutate(ave = replace(ave, ave=='av.de_la_cite', 'av.de_la_cite.mul')) %>%
  mutate(ave = replace(ave, ave=='av.de_l_authenticite', 'av.de_l_authenticite.rb1')) %>%
  mutate(ave = replace(ave, ave=='av.de_la_mission', 'av.de_la_mission.sg')) %>%
  mutate(ave = replace(ave, ave=='av.de_la_paix', 'av.de_la_paix.kim')) %>%
  mutate(ave = replace(ave, ave=='av.du_24_novembre', 'av.du_24_novembre.rb1')) %>%
  mutate(ave = replace(ave, ave=='av.du_30_juin', 'av.du_27_octobre.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.du_lac', 'av.du_lac.kav')) %>%
  mutate(ave = replace(ave, ave=='av.du_marche', 'av.du_marche.kav')) %>% 
  mutate(ave = replace(ave, ave=='av.du_projet', 'av.du_projet.kav')) %>%
  mutate(ave = replace(ave, ave=='av.du_stade', 'av.du_stade.kim')) %>%
  mutate(ave = replace(ave, ave=='av.du_stade', 'av.du_stade.sg')) %>%
  mutate(ave = replace(ave, ave=='av.fizi', 'av.fizi.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.haut_congo', 'av.haut_congo.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.idjwi_1.rb2', 'av.idjwi.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.idjwi_2.rb2', 'av.idjwi.rb2')) %>% 
  mutate(ave = replace(ave, ave=='av.idjwi_3.rb2', 'av.idjwi.rb2')) %>% 
  mutate(ave = replace(ave, ave=='av.kabomboza.kas', 'av.kabomboza')) %>%
  mutate(ave = replace(ave, ave=='av.kabungulu_1', 'av.kabungulu.kim')) %>%
  mutate(ave = replace(ave, ave=='av.kagenge', 'av.kagenge.kal')) %>%
  mutate(ave = replace(ave, ave=='av.kakamba', 'av.kakamba.kal')) %>%
  mutate(ave = replace(ave, ave=='av.kakungwe_1.rb1', 'av.kakungwe.rb1')) %>%
  mutate(ave = replace(ave, ave=='av.kakungwe_1.rb2', 'av.kakungwe.rb1')) %>% 
  mutate(ave = replace(ave, ave=='av.kalehe', 'av.kakungwe.rb1')) %>%
  mutate(ave = replace(ave, ave=='av.kalembelembe', 'av.kalembelembe.kav')) %>%
  mutate(ave = replace(ave, ave=='av.kalimabenge', 'av.kalimabenge.kab')) %>%
  mutate(ave = replace(ave, ave=='av.kamongola', 'av.kamongola.kal')) %>%
  mutate(ave = replace(ave, ave=='av.karigo', 'av.karigo.kal')) %>%
  mutate(ave = replace(ave, ave=='av.kashekebwe', 'av.kashekebwe.kab')) %>%
  mutate(ave = replace(ave, ave=='av.kasia', 'av.kasia.kal')) %>%
  mutate(ave = replace(ave, ave=='av.kigongo', 'av.kigongo.kak')) %>%
  mutate(ave = replace(ave, ave=='av.kijaga', 'av.kijaga.kab')) %>%
  mutate(ave = replace(ave, ave=='av.kimbangu' & quartier=='kabindula', 
                       'av.kimbangu.kab')) %>%
  mutate(ave = replace(ave, ave=='av.kinaga', 'av.kinaga.rug')) %>%
  mutate(ave = replace(ave, ave=='av.kinogono', 'av.kinogono.kal')) %>%
  mutate(ave = replace(ave, ave=='av.kinogono_1.kal', 'av.kinogono.kal')) %>%
  mutate(ave = replace(ave, ave=='av.kitundu', 'av.kitundu.kak')) %>%
  mutate(ave = replace(ave, ave=='av.kitunya', 'av.kitunya.kib')) %>%
  mutate(ave = replace(ave, ave=='av.kyonga', 'av.kyonga.kil')) %>%
  mutate(ave = replace(ave, ave=='av.lenghe_III' & quartier=='kakombe', 
                       'av.lenghe_III.kak')) %>%
  mutate(ave = replace(ave, ave=='av.lenghe_III' & quartier=='kilibula', 
                       'av.lenghe_III.kil')) %>%
  mutate(ave = replace(ave, ave=='av.lumbulumbu', 'av.lumbulumbu.nya')) %>%
  mutate(ave = replace(ave, ave=='av.lumumba', 'av.lumumba.kav')) %>%
  mutate(ave = replace(ave, ave=='av.maendeleo' & quartier=='kasenga', 
                       'av.maendeleo.kas')) %>%
  mutate(ave = replace(ave, ave=='av.maendeleo' & quartier=='kilibula', 
                       'av.maendeleo.kil')) %>%
  mutate(ave = replace(ave, ave=='av.maendeleo' & quartier=='rugenge', 
                       'av.maendeleo.rug')) %>%
  mutate(ave = replace(ave, ave=='av.majengo', 'av.majengo.kas')) %>%
  mutate(ave = replace(ave, ave=='av.major_vangu', 'av.major_vangu.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.makobola', 'av.makobola.mul')) %>%
  mutate(ave = replace(ave, ave=='av.mangondo', 'av.mangondo.kas')) %>%
  mutate(ave = replace(ave, ave=='av.maombi', 'av.maombi.kab')) %>%
  mutate(ave = replace(ave, ave=='av.mapendo' & quartier=='kasenga', 
                       'av.maombi.kas')) %>%
  mutate(ave = replace(ave, ave=='av.mapendo' & quartier=='kavimvira', 
                       'av.maombi.kav')) %>%
  mutate(ave = replace(ave, ave=='av.mapinduzi', 'av.mapinduzi.sg')) %>%
  mutate(ave = replace(ave, ave=='av.matadi_2', 'av.matadi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.mombasa', 'av.mombasa.kil')) %>%
  mutate(ave = replace(ave, ave=='av.munanira', 'av.munanira.nya')) %>%
  mutate(ave = replace(ave, ave=='av.muranvya', 'av.muranvya.kas')) %>%
  mutate(ave = replace(ave, ave=='av.musabwa', 'av.musabwa.kab')) %>%
  mutate(ave = replace(ave, ave=='av.musheru', 'av.musheru.kas')) %>%
  mutate(ave = replace(ave, ave=='av.musulmane', 'av.musulmane.kak')) %>%
  mutate(ave = replace(ave, ave=='av.musumba', 'av.musumba.nya')) %>%
  mutate(ave = replace(ave, ave=='av.mutarure', 'av.mutarure.kal')) %>%
  mutate(ave = replace(ave, ave=='av.ngovi_mgja', 'av.ngovi_mgja.kal')) %>%
  mutate(ave = replace(ave, ave=='av.nyakyoya', 'av.nyakyoya.kak')) %>%
  mutate(ave = replace(ave, ave=='av.nyamianda_1.kim', 'av.nyamianda.kim')) %>%
  mutate(ave = replace(ave, ave=='av.nyangara', 'av.nyangara.kav')) %>%
  mutate(ave = replace(ave, ave=='av.nyatutwa', 'av.nyatutwa.kak')) %>%
  mutate(ave = replace(ave, ave=='av.nyoroka', 'av.nyoroka.kal')) %>%
  mutate(ave = replace(ave, ave=='av.nyoroka_1', 'av.nyoroka.kal')) %>%
  mutate(ave = replace(ave, ave=='av.nyoroka_1.kal', 'av.nyoroka.kal')) %>%
  mutate(ave = replace(ave, ave=='av.nyoroka_2', 'av.nyoroka.kal')) %>%
  mutate(ave = replace(ave, ave=='av.nyoroka_2.kal', 'av.nyoroka.kal')) %>%
  mutate(ave = replace(ave, ave=='av.petrocongo', 'av.petrocongo.rug')) %>%
  mutate(ave = replace(ave, ave=='av.reboisement' & quartier=='kakombe', 
                       'av.reboisement.kak')) %>%
  mutate(ave = replace(ave, ave=='av.reboisement' & quartier=='kasenga', 
                       'av.reboisement.kas')) %>%
  mutate(ave = replace(ave, ave=='av.reboisement.kakkib','av.reboisement.kak')) %>%
  mutate(ave = replace(ave, ave=='av.rond_point', 'av.rond_point.kav')) %>%
  mutate(ave = replace(ave, ave=='av.rugembe', 'av.rugembe.kal')) %>%
  mutate(ave = replace(ave, ave=='av.rugenge_nord', 'av.rugembe.rug')) %>%
  mutate(ave = replace(ave, ave=='av.rugenge_nord.rug', 'av.rugembe.rug')) %>%
  mutate(ave = replace(ave, ave=='av.rugenge_sud', 'av.rugembe.rug')) %>%
  mutate(ave = replace(ave, ave=='av.rugenge_sud.rug', 'av.rugembe.rug')) %>%
  mutate(ave = replace(ave, ave=='av.shaba', 'av.shaba.kak')) %>%
  mutate(ave = replace(ave, ave=='av.shishi', 'av.shishi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.shishi_1', 'av.shishi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.shishi_1.mul', 'av.shishi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.shishi_2.mul', 'av.shishi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.shishi_4', 'av.shishi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.shishi_4.mul', 'av.shishi.mul')) %>%
  mutate(ave = replace(ave, ave=='av.tanganyika', 'av.tanganyika.kav')) %>%
  mutate(ave = replace(ave, ave=='av.tupendane', 'av.tupendane.kav')) %>%
  mutate(ave = replace(ave, ave=='av.ushirika', 'av.ushirika.rug')) %>%
  mutate(ave = replace(ave, ave=='av.uvira', 'av.uvira.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.virunga', 'av.virunga.kak')) %>%
  mutate(ave = replace(ave, ave=='av.yohana', 'av.yohana.mul')) %>%
  mutate(ave = replace(ave, ave=='av.du_port' & quartier=="kalundu", 
                       'av.du_port.kal')) %>% 
  mutate(ave = replace(ave, ave=='av.bongisa', 'av.bongisa.kal')) %>%
  mutate(ave = replace(ave, ave=='av.apollo', 'av.apollo.mul')) %>%
  mutate(ave = replace(ave, ave=='av.du_15_decembre', 'av.du_15_decembre.rb2')) %>%
  mutate(ave = replace(ave, ave=='av.kakungwe_2.rb1', 'av.kakungwe.rb1')) %>%
  mutate(ave = replace(ave, ave=='av.kakungwe_2.rb1', 'av.kakungwe.rb1')) %>%
  mutate(ave = replace(ave, ave=='av.kimbangu' & quartier=='kavimvira', 
                       'av.kimbangu.kav')) %>%
    mutate(ave = replace(ave, ave=='av.kimbangu' & quartier=='kabindula', 
                       'av.kimbangu.kab')) %>%
  mutate(ave = replace(ave, ave=='av.mitumba_1_2_3', 'av.mitumba_1_2_3.mul')) %>%
  mutate(ave = replace(ave, ave=='av.mitumba_1_2_3', 'av.mitumba_1_2_3.mul')) %>%
  mutate(ave = replace(ave, ave=='av.rubenga', 'av.rubenga.kav')) %>%
  mutate(ave = replace(ave, ave=='av.solange', 'av.solange.kal')) %>% 
  mutate(ave = replace(ave, ave=='av.cpgl.kav', 'av.CPGL.kav')) %>%
  mutate(ave = replace(ave, ave=='av.kabomboza', 'av.kabomboza.kas')) 

uvira.cholera.new <- uvira.cholera.new %>% 
  mutate(ave = replace(ave, ave=='av.rugembe.rug', 'av.rugembe.kal')) %>% 
  mutate(ave = replace(ave, ave=='av.kivu_nord', 'av.kivu.kas')) %>% 
  mutate(ave = replace(ave, ave=='av.maombi.kas', 'av.maombi.kab')) %>% 
  mutate(ave = replace(ave, ave=='av.maombi.kav', 'av.maombi.kab')) %>% 
  mutate(ave = replace(ave, ave=='av.namyanda', 'av.nyamianda.kim')) %>% 
  mutate(ave = replace(ave, ave=='av.reboisement.kas', 'av.reboisement.kak'))
  
# check # of cases for different time periods
uvira.cholera.new %>% 
  group_by(year) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1)) 
# RESULT: returns 5447 for 2016-2020 period (1341+1134+1000+922+1050)

  aves <- uvira.cholera.new %>% 
    group_by(ave) %>%
    summarise (n=n()) %>%
    mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1)) 

#ensure minimal missing aves in each dataset
setdiff(aves$ave, uga$UGA)

#clean up
rm(uvira.cholera, prop.test, aves, uga)

#sort again and put key variables together for easy viewing
uvira.cholera.new <- uvira.cholera.new %>% 
  arrange(uvira.cholera.new, date_adm, ave) %>% 
  relocate(ave, .after = date_adm)

#merge polygons with same UGA into a single polygon
uvira.ave.shape.2 <- uvira.ave.shape %>%
  group_by(UGA) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

#check average size of avenue
plot(uvira.ave.shape.2)
uvira.ave.shape.2$area.m2 <- st_area(uvira.ave.shape.2) #in m2
uvira.ave.shape.2$area.km2 <- st_area(uvira.ave.shape.2)/1000000 #in km2
mean(uvira.ave.shape.2$area.m2) ##77213.24 m^2
mean(uvira.ave.shape.2$area.km2) ##0.07721324 km^2

# save Uvira shape file to disk
st_write(uvira.ave.shape.2, "uvira.ave.shp")

# Add centroid X and Y coordinates
uvira.ave.shape.2.cent <- st_centroid(uvira.ave.shape.2)

# plot to check for centroids
ggplot() + geom_sf(data = uvira.ave.shape.2, fill = 'white') + 
  geom_sf(data = uvira.ave.shape.2.cent, color = 'red') 

# separate out the coordinates
coods <- uvira.ave.shape.2.cent %>% st_coordinates()
uvira.ave.shape.2.cent <- cbind(uvira.ave.shape.2.cent, coods)  
uvira.ave.shape.2.cent <- as.data.frame(uvira.ave.shape.2.cent)
class(uvira.ave.shape.2.cent)

# Merge case file and shapefile
uvira <- merge(uvira.cholera.new, uvira.ave.shape.2.cent, 
  by.x = 'ave', by.y = 'UGA', all.x = TRUE
)

# check # of cases for different time periods
uvira.cholera.new %>% 
  filter(year>=2016 & year<=2020) %>% 
  group_by(year) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))
# RESULT: returns 5447 cases (1341 + 1134 + 1000 + 922 + 1050)

# check # of cases tested (or not)
uvira.cholera.new %>% 
  filter(year>=2016 & year<=2020) %>% 
  group_by(tested) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))
# RESULT: 3456 + 1 were tested

# check result, if tested
uvira.cholera.new %>% 
  filter(year>=2016 & year<=2020) %>% 
  filter(tested == 'Tested') %>% 
  group_by(result) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))
# Result, 1493 43.2% were positive 

# check # of cases for different time periods
uvira %>% 
  filter(year>=2016 & year<=2020) %>% 
  group_by(year) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))
# returns 5447 cases
# RESULT (table 1)
#1  2016  1341 24.6%   
#2  2017  1134 20.8%   
#3  2018  1000 18.4%   
#4  2019   922 16.9%   
#5  2020  1050 19.3%  

# check # of cases tested for different time periods
uvira %>% 
  filter(year>=2016 & year<=2020) %>% 
  filter(tested=="Tested") %>% 
  group_by(year) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))
# RESULT: 
#1  2016   617 17.9%   
#2  2017   857 24.8%   
#3  2018   533 15.4%   
#4  2019   597 17.3%   
#5  2020   852 24.7% 

# RDT result by year
uvira %>% 
  filter(year>=2016 & year<=2020) %>% 
  filter(tested=="Tested") %>% 
  group_by(year, result) %>% 
  summarise (n=n()) %>%
  mutate(rel.freq =  scales::percent(n/sum(n), accuracy = 0.1))
# RESULT
#year result       n rel.freq
#1  2016 Negative   391 63.4%   
#2  2016 Positive   226 36.6%   
#3  2017 Negative   483 56.4%   
#4  2017 Positive   374 43.6%   
#5  2018 Negative   300 56.3%   
#6  2018 Positive   233 43.7%   
#7  2019 Negative   337 56.4%   
#8  2019 Positive   260 43.6%   
#9  2020 Negative   452 53.1%   
#10  2020 Positive   400 46.9%  


# save Uvira case file to disk
saveRDS(uvira, file="uvira.cases.2016.2020.rds")

# Create a version that has only the data needs for the 3 SATSCAN
# files

uvira <- as.data.frame(uvira)
class(uvira)

### ** This is the departure point for the RDT-positive only analyses too

#group cases by ave and date
# results in 13 165 rows of 
uvira.ss <- uvira %>% 
  group_by(ave, date_adm, X, Y) %>% 
  summarize(count_by_avedate_adm=n()) %>% 
  rename(cases=count_by_avedate_adm,
         y=X,x=Y,date=date_adm) %>% 
  arrange(date,ave) 

sum(uvira.ss$cases)
#14 400 cases total 

#create serial id
uvira.ss$id <- seq.int(nrow(uvira.ss))
uvira.ss <- uvira.ss %>% 
  relocate(id, .before = date) %>% 
  relocate(ave, .after = date)

# check # of cases for different time periods
uvira.ss %>% 
  filter(date >= as.Date("2016-01-01") & date <= as.Date("2020-12-31"))
sum(uvira.ss$cases)
#13 163 rows total with 13 163 rows (to be processed with Satscan)


##########################################################################################################################################################################################
##
## Use Worldpop representation of Uvira to setup the population of
## Uvira in raster format
## shp file for Sud Kivu: https://data.humdata.org/dataset/a4cc05e3-57ce-42ac-a4e1-5ac891845762/resource/8e0e1982-e12d-4275-bb91-17b713848a51/download/south_kivu_drc_adm1-3.zip
## raster file for population: https://www.worldpop.org/geodata/summary?id=6348
##
## Guidance from: https://www.dante-project.org/vignettes/exactextractr-pop
##
#########################################################################################################################################################################################

# import pop file with 2017 census estimates (from Aurelie's thesis)
uvira.ave.pop.2017 <- read.csv('C:\\R-projects\\Uvira-spatial-epid\\uvpop.csv')
uvira.ave.pop.2017 <- uvira.ave.pop.2017 %>% 
  dplyr::select(ave, pop2, date2) %>% 
  rename(pop=pop2, date=date2)
# 216 avenues total

# calculate mean and range of population sizes of avenues
round(mean(uvira.ave.pop.2017$pop)) # 1177.949
range(uvira.ave.pop.2017$pop) # 180 - 5711
sum(uvira.ave.pop.2017$pop)  # 254 437

# make minor correction to uvira.ss file as there is a 2021 case without an avenue
uvira.ss <- uvira.ss[!(uvira.ss$ave=="av.village_kabindula"),]

# creates distance matrix between avenue centroids (in metres) to check min/max/mean 
pixdistmat <- distm(cbind(uvira.ave.shape.2.cent$X, 
                          uvira.ave.shape.2.cent$Y))
colnames(pixdistmat) <- seq_along(1:216)

# min, max, mean distance between patches
round(min(pixdistmat[pixdistmat > 0]), 1) #43.7 m 
round(max(pixdistmat, na.rm = TRUE), 1) #12435.9 m / 12.4 km
round(mean(pixdistmat, na.rm = TRUE), 1) #3111.2 / 3.1 km
# other measures of dispersal
hist(pixdistmat, breaks=32, col="lightblue1")
lines(density(pixdistmat), col="dodgerblue3", lwd=2)

round(quantile(pixdistmat, probs = 
                 c(0.05, 0.06, 0.1, 0.15, 0.2, 
                   0.25, 0.5, 0.75, 0.95, 1)), digits=2)
# 5%      10%      15%      20%      25%      50%      75%      95% 
# 419.51   648.85   854.90  1068.15  1302.95  2708.98  4466.04  7307.41 
# 100% 
# 12435.92
# Therefore, 5% of distance measures between avenues are <420m
# i.e., =0.05*15876 measures (794 measures)

######################################################################
##
## Create files for each year
##
######################################################################

#### 2016 to 2020

#create case file for entire period (2016-2020)
uvira.case.2016.2020 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2020-12-31"))
write.csv(uvira.case.2016.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case-2016-20.csv', row.names=FALSE)
# 5447 cases total
# 4997 rows total
sum(uvira.case.2016.2020$cases)

#create coordinates file for entire period (2016-2020)
uvira.coord.2016.2020 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2020-12-31")) 
write.csv(uvira.coord.2016.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord-2016-20.csv', row.names=FALSE)
#4997 rows total

#Create pop file for Satscan with locationID
uvira.pop.2016.2020 <- merge(uvira.ave.pop.2017, uvira.case.2016.2020, by='ave')
uvira.pop.2016.2020 <- uvira.pop.2016.2020 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.2016.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop-2016-20.csv', row.names=FALSE)
#4997 rows total
 
############################ SINGLE YEARS #########################################

#### 2016

#case file for 2016
# results in 1 208 cases
uvira.case.2016 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2016-12-31"))
write.csv(uvira.case.2016,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case-2016.csv', row.names=FALSE)
# 1341 cases total
# 1208 rows total
sum(uvira.case.2016$cases)

#coordinates file for 2016
uvira.coord.2016 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2016-01-01") &
           date <= as.Date("2016-12-31")) 
write.csv(uvira.coord.2016,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord-2016.csv', row.names=FALSE)
#1208 rows total

#pop file for 2016
uvira.pop.2016 <- merge(uvira.ave.pop.2017, uvira.case.2016, by='ave')
uvira.pop.2016 <- uvira.pop.2016 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.2016,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop-2016.csv', row.names=FALSE)
#1208 rows total

#### 2017

#case file for 2017
uvira.case.2017 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2017-01-01") &
           date <= as.Date("2017-12-31"))
write.csv(uvira.case.2017,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case-2017.csv', row.names=FALSE)
# 1134 cases total
# 1031 rows total
sum(uvira.case.2017$cases)


#coordinates file for 2017
uvira.coord.2017 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2017-01-01") &
           date <= as.Date("2017-12-31")) 
write.csv(uvira.coord.2017,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord-2017.csv', row.names=FALSE)
#1031 rows total

#pop file for 2017
uvira.pop.2017 <- merge(uvira.ave.pop.2017, uvira.case.2017, by='ave')
uvira.pop.2017 <- uvira.pop.2017 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.2017,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop-2017.csv', row.names=FALSE)
#1031 rows total


#### 2018

#case file for 2018
uvira.case.2018 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2018-01-01") &
           date <= as.Date("2018-12-31"))
write.csv(uvira.case.2018,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case-2018.csv', row.names=FALSE)
# 1000 cases total
# 951 rows total
sum(uvira.case.2018$cases)

#coordinates file for 2018
uvira.coord.2018 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2018-01-01") &
           date <= as.Date("2018-12-31")) 
write.csv(uvira.coord.2018,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord-2018.csv', row.names=FALSE)
#951 rows total

#pop file for 2018
uvira.pop.2018 <- merge(uvira.ave.pop.2017, uvira.case.2018, by='ave')
uvira.pop.2018 <- uvira.pop.2018 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.2018,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop-2018.csv', row.names=FALSE)
#951 rows total

#### 2019

#case file for 2019
uvira.case.2019 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2019-01-01") &
           date <= as.Date("2019-12-31"))
write.csv(uvira.case.2019,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case-2019.csv', row.names=FALSE)
# 922 cases total
# 855 rows total
sum(uvira.case.2019$cases)

#coordinates file for 2019
uvira.coord.2019 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2019-01-01") &
           date <= as.Date("2019-12-31")) 
write.csv(uvira.coord.2019,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord-2019.csv', row.names=FALSE)
# 855 rows total


#pop file for 2019
uvira.pop.2019 <- merge(uvira.ave.pop.2017, uvira.case.2019, by='ave'#,
                        #all.x = TRUE
                        )
uvira.pop.2019 <- uvira.pop.2019 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x)
write.csv(uvira.pop.2019,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop-2019.csv', row.names=FALSE)
# 855 rows total

#case file for 2020
uvira.case.2020 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,cases,date,ave) %>% 
  filter(date >= as.Date("2020-01-01") &
           date <= as.Date("2020-12-31"))
write.csv(uvira.case.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.case-2020.csv', row.names=FALSE)
# 1050 cases total
# 952 rows total
sum(uvira.case.2020$cases)

#coordinates file for 2020
uvira.coord.2020 <- uvira.ss %>% 
  data.frame() %>%
  dplyr::select(id,x,y,date,ave) %>% 
  filter(date >= as.Date("2020-01-01") &
           date <= as.Date("2020-12-31")) 
write.csv(uvira.coord.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.coord-2020.csv', row.names=FALSE)
# 952 rows total

#pop file for 2020
uvira.pop.2020 <- merge(uvira.ave.pop.2017, uvira.case.2020, by='ave')
uvira.pop.2020 <- uvira.pop.2020 %>% 
  dplyr::select(id,date.x,pop) %>% rename(date=date.x) 
write.csv(uvira.pop.2020,
          'C:\\R-projects\\Uvira-spatial-epid\\uvira.ave.pop-2020.csv', row.names=FALSE)
# 952 rows total