############################################################################
## This code runs the spatial meta-population model of case-area targeted ##
## intervention against cholera.                                          ##
##                                                                        ##    
## Author: Ruwan Ratnayake, May 2023                                      ##
############################################################################

# This script:
# 1. Sets up case data, population data, and the meta-populations on a map
# 2. Collates the fixed parameters
# 3. Runs the stochastic simulation
#
# The outputs include:
# 1. Daily counts of each model state for every patch (and summed across system)
# 2. Coverage table showing which patches were covered with which intervention 
#    on which days.
# 3. Reduction in outbreak size at time, t.
# 4. Reduction of the proportion of patches affected by time, t.


######################################################################
# 1. Data preparation and meta-population map and matrix construction
######################################################################

# Set working directory and load packages
setwd("C:/R-projects/Cholera-metapop-model")
source('packages.R') 

# Load case data (and trim to the 2020 cholera season, roughly 2019-2020)
#source('uvdat.R')
uvdat <- readRDS('uvdat.rds')
uvdat2020 <- uvdat[2323:3116, ]

# Load population data and raster layer
source('uvpop.R')
sum(total.pop)
plot(uvpop)

# Process population raster into patches (50m, 100m, 500m, 200m, 500m, 1000m) 
source('popprocess.R')

# Plot 50x50m and 100x100m patch maps
#spplot(population_raster50mshp, "layer", col.regions = bpy.colors(100, cutoff.tails = 0.1, alpha = 1.0), 
#       lwd=0, main = "Uvira population by patch (8012 50m x 50m patches)", cex.main=0.8) 
tiff('C:/Users/User/Documents/LSHTM/- PhD Projects/4 Mathematical modelling (MOD)/Work/Model documentation/popdensity.tiff', 
     units="in", width=32, height=40, res=300)

uvira_proj <- spTransform(population_raster100mshp, CRS("+init=EPSG:4326"))
uvira_proj@bbox  # Check dimensions to help guide offset choices

#     min       max
#x   4493182   4503578
#y -19551405 -19536822

scale = list("SpatialPolygonsRescale", 
             layout.scale.bar(), scale = 25000, 
             fill = c("transparent", "black"), 
             offset = c(5003182, -18001405))

# The scale argument sets length of bar in map units
text1 = list("sp.text", c(41000, 4104000), "0")
text2 = list("sp.text", c(66000, 4104000), "25 km")

arrow = list("SpatialPolygonsRescale", 
             layout.north.arrow(), 
             offset = c(10, 10), scale = 10000)

spplot(uvira_proj, #"layer", 
       sp.layout = list(arrow))

spplot(population_raster100mshp, "layer", col.regions = bpy.colors(100, cutoff.tails = 0.1, alpha = 1.0), 
       lwd=0, par.settings = list(axis.line = list(col = "transparent")), cex.main=0.8) 
dev.off()

# From here on, we use 100m pixel system (that equates to 2003 occupied patches), 
# which means that we can only test CATI radii of sizes of multiples of 50m 
# (50m, 150m, 250m)

# Convert patches into a dataframe of 2003 occupied patches with patchID (1-2003), 
# pixel IDs, latitude, longitude, population:
source('makeunipix.R')

# To save the patch locations for later reference:
# Update raster file to include both patchID and population
population_raster100m[population_raster100m@data@values>0]=unipix100$patchID
# plot out with values into super large TIFF
tiff('C:/Users/User/Documents/LSHTM/- PhD Projects/4 Mathematical modelling (MOD)/Work/Model documentation/patchmap.tiff', 
     units="in", width=32, height=40, res=300)
plot(population_raster100m)
text(population_raster100m, cex=0.5)
dev.off()

# Creates distance matrix between patches 
pixdistmat <- distm(cbind(unipix100$x, unipix100$y))
colnames(pixdistmat) <- seq_along(1:nrow(unipix100))
# Calculate min, max, mean distance between patches
min(pixdistmat[pixdistmat > 0]); max(pixdistmat, na.rm = TRUE); mean(pixdistmat, na.rm = TRUE) 
# min = 92.1 m, max = 14,208.0 m, mean = 4,282.6 m

# Bind patchIDs to case file by applying the pix.id.find function to x and y 
# columns of uvdat; binds to uvdat case dataframe
source('pixidfind.R') 
uvdat <- data.frame(uvdat, patchID = apply(cbind(uvdat[, 2:3]), 1, pix.id.find, unipix100))

# Set options for start date of model, as epiweeks representing different start
# dates of the seasonal outbreaks in 2020 
weekdates = c(175, 176, 177, 178, 179, 180, 202)

# Show epicurve of the outbreak in 2020
#cholera_epi_week <- incidence(x = uvdat2020, date_index = Date, interval = "week")
#plot(cholera_epi_week)

# Load function to identify patches to be targeted with CATI (by size of CATI radius)
source('cholspots.R')  
source('massspots.R')

# Create a "not in" operator for intervention coding and blank "(NUll/NA/na/0)" 
# function for special case of ACP zeroing out
# see https://stackoverflow.com/questions/19655579/a-function-that-returns-true-on-na-null-nan-in-r
`%notin%` <- Negate(`%in%`)
is.blank <- function(x, false.triggers=FALSE){if(is.function(x)) return(FALSE)
  return(is.null(x) || length(x) == 0 || all(is.na(x)) || all(x=="") || (false.triggers && all(!x)))}

# Set seed for random number generation
set.seed(1981)

# Setup parameters for spatial isotropic kernel / force of infection
min(pixdistmat[pixdistmat>0]) # 92.15 m
max(pixdistmat[pixdistmat>0]) # 14207.95 m
median(pixdistmat[pixdistmat>0]) # 3824.11 m
quantile(pixdistmat[pixdistmat>0], prob = 0.75) # 6299.15
D1 = 50 # half of the grid size used (100m x 100m), minimum distance, d_0
D2 = 1000 # max distance of transmission risk 
#D2 = 6299 # 75th percentile of distances, maximum distance, d_max
a = 2.1 # normalization constant for isotropic kernel
# Load function to calculate spatial kernels between patches
source('kernel.R')

################################
# 2. Functions to run the model
################################

# ! Run the model function below first

# model run with CATI
system.time(
cholmod.sim.optim <- CATI.model(weekdates, uvdat, uvpop, unipix100, 
                                pixdistmat, steps, paramsList = list(optim = optim)))

# model run with CATI (no OCV)
system.time(
cholmod.sim.nOCV <- CATI.model(weekdates, uvdat, uvpop, unipix100, 
                               pixdistmat, steps, paramsList = list(nOCV = nOCV)))

# model run with mass intervention
system.time(
cholmod.sim.mass <- CATI.model(weekdates, uvdat, uvpop, unipix100, 
                               pixdistmat, steps, paramsList = list(mass = mass)))

# model run without CATI
rm(WASH_Eff, ACP_Eff, ACP_Dur, OCV_Eff, OCV_Lag)
system.time(
cholmod.sim.nocati <- CATI.model(weekdates, uvdat, uvpop, unipix100, 
                                 pixdistmat, steps, paramsList = list(none = none)))

 ###########################
# 3. Load fixed parameters  
###########################

# There are two free parameters to calibrate to later:
# beta (individual exposure parameter) and alpha (shape of kernel)

IP = 1.4 # median incubation period is 1.4d (Azman, 2013)
IP.g = sort(sample(rgamma(1000, shape=IP, rate=0.5), 8), decreasing=FALSE) # distribution
infper = 5 # duration of shedding/infectiousness is 5d (Weil, 2014)
infper.acp = infper-2.74 # infectious period is reduced by 2.74 days (Lewnard, 2016)
pdetect = 0.8  # assumption about probability of case detection by surveillance system
alpha = 150 # calculated spatial zone of infection risk around a case (Ratnayake, 2023)
asymp.prop = 0.5 # Nelson et al, 2009

# calculation of beta
# assume Re = 1.7 (from Uvira cholera data in 2020)
Re = 1.5 # Re= beta * duration of infectiousness
(beta = Re/infper) 
# when Re = 1.5, beta = 0.3
# when Re = 2.0, beta = 0.4
# when Re = 3.0, beta = 0.6


# 4. Main model function

CATI.model <- function(
    simN,       # added for simulation number (**turn off if running single instance)
    weekdates,  # start and end weeks for simulation
    uvdat,      # case data to seed and fit the model
    uvpop,      # raster of population per pixel
    unipix,     # universal pixel lookup table for given raster
    pixdistmat, # patch distance matrix for given raster
    steps,      # days simulation should run
    paramsList = NULL){ # allocate all CATI interventions

# 4.1. Apply parameters 
  
  if(("optim" %in% names(paramsList)) | ("nOCV" %in% names(paramsList)) |
     ("mass" %in% names(paramsList))){paramsList = c(params = list(data.frame(
      CATI_Rad = 150, ACP_Rad = 50, CATI_Del = 2, OCV_Eff = (1-0.87), OCV_Lag = 7, OCV_Cov = 0.8,
      WASH_Eff = 0.66, ACP_Eff = (1-0.955), ACP_Dur = 2, MASS_Del = 14, OCV_Del = 30,
      MASS_Rad = 1000)), paramsList)
     # Load dynamic and intervention parameters 
     CATI_Rad = as.numeric(paramsList$params[1])
     ACP_Rad = as.numeric(paramsList$params[2])
     CATI_Del = as.numeric(paramsList$params[3])
     OCV_Eff = as.numeric(paramsList$params[4])
     OCV_Lag = as.numeric(paramsList$params[5])
     OCV_Cov = as.numeric(paramsList$params[6])
     WASH_Eff = as.numeric(paramsList$params[7])
     ACP_Eff = as.numeric(paramsList$params[8])
     ACP_Dur = as.numeric(paramsList$params[9]) 
     MASS_Del = as.numeric(paramsList$params[10])
     OCV_Del = as.numeric(paramsList$params[11])
     MASS_Rad = as.numeric(paramsList$params[12])}else 
        if("none" %in% names(paramsList)){paramsList = c(params = list(data.frame(
          CATI_Rad = 150, ACP_Rad = 50)), paramsList)
        # if no intervention
        CATI_Rad = as.numeric(paramsList$params[1])
        ACP_Rad = as.numeric(paramsList$params[2])}
  
# 4.2. Set up SEIR compartments

# (a) S, E, I, R, D compartments (D = patches with newly-detected cases)
# (b) Assume that 10-25% of Susceptibles are immune due to prior OCV vaccination 
# or infection during the past year. They are moved directly to Recovered compartment.
# (c) Exposed individuals remain for 2 days before moving to Infected compartment. 
# (d) Infected individuals remain infectious for 5 days before moving to Recovered compartment.

sing_S = round(unipix100$pop, 0)
sing_E = matrix(0, ncol = nrow(unipix100), nrow = 2) # 2 exposed days
sing_I = matrix(0, ncol = nrow(unipix100), nrow = 5) # 5 infectious days
sing_R = rep(0, nrow(unipix100))
sing_D <- matrix(0, nrow = nrow(sing_I), ncol = ncol(sing_I))

# Assume 25-50% are vaccinated or have infection-acquired immunity 
(prior.imm <- runif(1, 0.25, 0.50))*100
immune = round(sing_S*prior.imm, 0)
# Move immune to Recovered compartment
sing_S = sing_S - immune
sing_R = sing_R + immune 

#check totals
sum(sing_S)
sum(sing_E)
sum(sing_I)
sum(sing_R)
sum(sing_S) + sum(sing_I) + sum(sing_E) + sum(sing_R)


# 4.2. Set up initial conditions

# Calculate the total initial numbers in each state, during the first week that 
# cholera is seeded (defined as first vector in weekdates,and the first day of 
# that week). Here, we want to inform the model with the first clusters from the 
# spatiotemporal analysis.

# Initial dataframe of seeding cases by number of cases and location, 
# for the first week with cases 
seed_inf = data.frame(number=uvdat$Number_of_cases[uvdat$Week==weekdates[3]], #Num cases for the stated week
                      location=uvdat$patchID[uvdat$Week==weekdates[3]]) #Patch location for the stated week

seed_inf = data.frame(number = seed_inf$number,
  number_true = ceiling(ceiling((sum(seed_inf$number)/pdetect))/nrow(seed_inf)),
  location = seed_inf$location,
  timing = sample(1:5, nrow(seed_inf), replace = T),
  pop = unipix100$pop[match(seed_inf$location, unipix100$patchID)])

# Allocate cases to the infectious compartment, by location and sampled timing
for(i in 1:nrow(seed_inf)){
  sing_I[cbind(seed_inf$timing[i], seed_inf$location[i])] = 
    sing_I[cbind(seed_inf$timing[i], seed_inf$location[i])] + seed_inf$number_true[i]
  # sample timings of cases within 5-day infectious period
  timsam <- sample(1:5, sum(seed_inf$number[i]), replace = T)
  # allocate undetected infectious cases
  for(k in 1:length(timsam)){sing_I[cbind(timsam[k], seed_inf$location[k])] = 
    sing_I[cbind(timsam[k], seed_inf$location[k])] + 1}
}


#check totals
sum(sing_S) + sum(sing_I) + sum(sing_E) + sum(sing_R)

# Check total # of current infections 
sum(rowSums(sing_I))
start.inf = cbind(unipix100[colSums(sing_I) > 0,], 
                  numinf = colSums(sing_I)[colSums(sing_I) > 0])
start.inf$numinf = start.inf$numinf 

# 4.3. Calculate the initial spatial force of infection 

# Make list of 2003 patches by infection status
patches_status <- as.data.frame(cbind(1:nrow(unipix100), 0))
colnames(patches_status) <- c('patchID','inf')
# Match patches with infections to patch list
patches_infected <- as.data.frame(cbind(start.inf$patchID, start.inf$numinf))
colnames(patches_infected) <- c('patchID','inf')
patches_status$inf <- patches_infected[match(patches_status$patchID, patches_infected$patchID), 2] %>% replace(is.na(.), 0)

# The force of infection on each patch with be proportional to the values 
# calculated, based on the range, alpha and an exponential functional form. 
# This sums up the foi on each patch, from all the other patches, at a given timestep.
# ! Note, I divided the foi by 100 to render it a proportion (to scale the impact of
# infections from neighboring patches on the current patch)

# Dispersal takes an exponential (the probability of dispersing a distance 
# d ~ exp(−d/a), where a is the range). The exponential model arises if we 
# assume dispersal happens in a constant direction with a constant stopping rate.

#foi.exp = as.data.frame(cbind(unipix100$patchID,
#        round((apply(exp(-(pixdistmat/alpha))*patches_status$inf, 2, sum)), digits=4)-patches_status$inf))
#colnames(foi.exp) <- c('patchID', 'foi')
#foi.exp$foi = foi.exp$foi/100 # convert to a proportion (of beta)
#foi.exp$foi[foi.exp$foi>1]=1

# calculate new matrix of spatial FoI as beta*sum(patches_j with infecteds)*kernel.plk

# 4.3. Setup main model run

# Duration of model (30-60 days for short-term impact)
steps = 60

# For 30-day period, total numbers of persons per state
totals = data.frame(S = rep(NA, steps), 
                    E = rep(NA, steps),
                    I = rep(NA, steps),
                    R = rep(NA, steps),
                    D = rep(NA, steps), #temp compartment for detection of new cases  
                    A = rep(NA, steps),
                    exposed = rep(NA, steps),
                    infected = rep(NA, steps),
                    recovered = rep(NA, steps),
                    CATI0 = rep(NA, steps),
                    CATI1 = rep(NA, steps),
                    CATI2 = rep(NA, steps),
                    CATI3 = rep(NA, steps),
                    CATI4 = rep(NA, steps), 
                    CATI5 = rep(NA, steps),
                    CATI.patches = rep(NA, steps),
                    CATI.num = rep(NA, steps)) 

# For a single timestep, total numbers of persons per state, per patch
hold = data.frame(S = rep(NA, nrow(unipix100)),
                  E = rep(NA, nrow(unipix100)),
                  I = rep(NA, nrow(unipix100)),
                  R = rep(NA, nrow(unipix100)),
                  D = rep(NA, nrow(unipix100)),
                  A = rep(NA, nrow(unipix100)),
                  infected = rep(NA, nrow(unipix100)), # added to document new infections per patch
                  kernel = rep(NA, nrow(unipix100))) 

# Master list of each timestep and each patch 
detailed_totals <- list()
for(i in 1:steps){detailed_totals[[i]] = hold}

# Tracks cholspots and which patches have been treated, by timestep
activespots_log <- vector("list", steps)
cholspots_log <- vector("list", steps)
cholspots_ACP_log <- vector("list", steps)
prevcholspots_log <- vector("list", steps)
treatlog_ACP <- vector("list", steps)
treatlog_WASH <- vector("list", steps)
treatlog_OCV <- vector("list", steps) # for ongoing vaccination effect against symptomatic infection
treatlog_OCV2 <- vector("list", steps) # for one-time vaccination effect to reduce susceptibility

# Tracks the total population size = 280,032
(tpop = sum(sing_S, sing_E, sing_I, sing_R))


# 5. Run the model

# Step 1: Calculate the stage transition probabilities based on dynamics and interventions
# Step 2: Update the states
# Step 3: Document newly-detected cases for CATI-targeting in the next timestep

# Initializes the progress bar and start time

pb <- txtProgressBar(min=0, max=steps, style = 3, width = 50, char = "=")

#foreach(i = 1:steps, 
#        .export = c('%notin%', 'find.cholspots', 'is.blank', 'kernel.plk.fun', 
#                    'make.unipix', 'pix.id.find', 'pop.process', 'pop.process2'),
#        .packages = c("tidyverse", "doParallel", "raster", "sf", "rgdal", "tmap",
#                      "lubridate", "tsibble", "incidence2", "truncnorm",
#                      "maptools", "rgeos")) %do% {

for(i in 1:steps){  # This is the main timestep loop

  # Sets the progress bar to the current state
  setTxtProgressBar(pb, i)
  
                       
###################################################################
# Step 1: calculate stage transition probabilities and CATI targets
###################################################################
  
  #################################################
  # 5.1. Calculate time-varying force of infection
  #################################################
  
#  if(i == 1){foi.exp = foi.exp} # end loop for foi, if i==1
#  else{
  patches_infected <- as.data.frame(cbind(patches_status$patchID[colSums(sing_I)>=1], 1))
  colnames(patches_infected) <- c('patchID','inf')
  # Update patches_status for patches with new infections
  patches_status$inf <- patches_infected[match(patches_status$patchID,
                                               patches_infected$patchID), 2] %>% replace(is.na(.), 0)
  # Update the spatial foi
#  foi.exp$foi = (round((apply(exp(-(pixdistmat/alpha))*patches_status$inf, 
#                             2, sum)), digits=4)-patches_status$inf)}
#  foi.exp$foi[foi.exp$foi>1]=1
#  foi.exp$foi = foi.exp$foi/100
  
  ######################################################
  # 5.2. Find cholspots (patches to be targeted by CATI)
  ######################################################
  
  # Cholspots finds patches within the radius specified by CATI_Rad (50/150/250m). 
  # Note that antibiotics are given to a smaller radius of 50m, so require their 
  # own cholspots calculation. Cholspots are cumulative, in that once a patch 
  # receives CATI treatment, this lasts for the study duration of 30-60 days.

  if("optim" %in% names(paramsList) | "nOCV" %in% names(paramsList) #| "mass" %in% names(paramsList)
     ){
    
  # General cholspots (for WASH, OCV)
  if(i==1){totalspots=0}
  spotlist <- find.cholspots(i, sing_I, sing_D, pixdistmat, radius = CATI_Rad,
                             cumulspots = totalspots)
  
  totalspots <- unlist(spotlist[1]) # cumulative spots
  cholspots <- unlist(spotlist[2]) # new spots
  prevcholspots <- unlist(spotlist[3]) # previous spots

  # Smaller radius cholspots (for ACP)
  if(i==1){totalspots.ACP=0}
  spotlist.ACP <- find.cholspots(i, sing_I, sing_D, pixdistmat, radius = ACP_Rad, 
                                 cumulspots = totalspots.ACP)
  totalspots.ACP <- unlist(spotlist.ACP[1]) # cumulative spots
  cholspots.ACP <- unlist(spotlist.ACP[2]) # new spots

  # Log the cholspots
  prevcholspots_log[[i]] = prevcholspots
  cholspots_log[[i]] = totalspots
  cholspots_ACP_log[[i]] = totalspots.ACP

  # (For monitoring the progress as the model ticks through the calculations...)
  #print(paste("timestep",i,sep=': '))
  #print(paste("total infections", sum(sing_I), sep = ': ')) 
  #print("Cholera hotspots detected"); print(cholspots)
  

  }
  
  # If mass intervention, we are targeting a very large radius around the
  # affected area, and WASH and OCV are implemented there after delays of
  # 5 days (WASH) and 30 days (OCV)
  if("mass" %in% names(paramsList)){
    if(i==1){massspots=0}
    massspots <- find.massspots(i, sing_I, sing_D, pixdistmat, radius.mass = MASS_Rad)
  }
  
  #########################################
  # 5.3. Setup the implementation calendar
  #########################################
  
  # Assignment of interventions (only for CATI or mass intervention)
  if("optim" %in% names(paramsList) | "nOCV" %in% names(paramsList) |
     "mass" %in% names(paramsList)){ 
  
  # For CATI, delay to team getting to field (CATI_Del) ~ 
  if("optim" %in% names(paramsList) | "nOCV" %in% names(paramsList)){
    implem.day = i + CATI_Del # at least a 2 day delay
    # Record which cholspots are targeted and on what day
    treatlog_WASH[[implem.day]] = totalspots
    } 
    
  # For mass intervention, delay to initially getting mass campaign mobilized
  if(("mass" %in% names(paramsList) & i > MASS_Del)){
    implem.day = i # Median delay ~ 6 days (Ratnayake, 2020, time review)
    # Record which cholspots are targeted and on what day
    treatlog_WASH[[implem.day]] = massspots}
  
  # Only include OCV with optimal CATI implementation and mass intervention
  # OCV peak vibriocidal response occurs 7-11 days after administration
  if("optim" %in% names(paramsList)){
  # For reduction of symptomatic infection (applied at first and subsequent timesteps)
    treatlog_OCV[[implem.day + OCV_Lag]] = totalspots}
  # For the one-time susceptibility-reducing effect of vaccination, trim treatlog_OCV
  # to include only the first instance of vaccination in a patch
  else if("mass" %in% names(paramsList) & i > OCV_Del){
    implem.day = i
    treatlog_OCV[[implem.day + OCV_Lag]] = massspots}
  
  # Adjust for one-time vaccination effect
  if("optim" %in% names(paramsList)){
    if(i >= (1 + CATI_Del + OCV_Lag)){ # therefore, eligible for a previous OCV
    prev.OCV = as.data.frame((1 + CATI_Del + OCV_Lag):(i-1))
    colnames(prev.OCV)<-c('timestep')
    for(e in prev.OCV$timestep){
      if(!is.blank(treatlog_OCV[[e]])){treatlog_OCV2[[i]] =
        treatlog_OCV[[i]][(treatlog_OCV[[i]] %notin% treatlog_OCV[[e]])]}}
    # Make sure first instance of OCV is included
    treatlog_OCV2[[1 + CATI_Del + OCV_Lag]] = treatlog_OCV[[1 + CATI_Del + OCV_Lag]]}}
  
  # Adjust for one-time vaccination effect
  else if("mass" %in% names(paramsList)){
    if(i >= (1 + OCV_Del + OCV_Lag)){ # therefore, eligible for a previous OCV
    prev.OCV = as.data.frame((1 + OCV_Del + OCV_Lag):(i-1))
    colnames(prev.OCV)<-c('timestep')
    for(e in prev.OCV$timestep){
      if(!is.blank(treatlog_OCV[[e]])){treatlog_OCV2[[i]] = 
        treatlog_OCV[[i]][treatlog_OCV[[i]] %notin% treatlog_OCV[[e]]]}}
    # Make sure first instance of OCV is included
    treatlog_OCV2[[1 + OCV_Del + OCV_Lag]] = treatlog_OCV[[1 + OCV_Del + OCV_Lag]]}}

  # Only include ACP for CATI
  if("optim" %in% names(paramsList) | "nOCV" %in% names(paramsList)){
    if("optim" %in% names(paramsList)){treatlog_ACP[[implem.day]] = totalspots.ACP}
    if("nOCV" %in% names(paramsList)){treatlog_ACP[[implem.day]] = totalspots}
    
  # Clear ACP effect from patches after 2 days when concentration of 
  # antibiotic is no longer sufficient to have an effect. Also, patches that
  # were treated with antibiotics once, cannot be treated again, so clear patches
  # from treatlog_ACP[i] that appeared more than 2 days ago in the list.
  #
  # ! Note that some patches may be treated with WASH and OCV and only later 
  # have a case that requires ACP. These patches can therefore cycle through
  # WASH+OCV, then ACP+WASH+OCV for 2 days and back to WASH+OCV
  
  if(i > CATI_Del+2){ # +2 here because the first 2 timesteps is in the first time a full ACP course can run
    # Vector of days to check to which patches antibiotics were given earlier
    prev.ACP = as.data.frame((CATI_Del+1):(i-2)) 
    colnames(prev.ACP)<-c('timestep')
    # Check if antibiotics were given earlier, and clear patch from current list
    for(d in prev.ACP$timestep){
      if(!is.blank(treatlog_ACP[[d]])){treatlog_ACP[[i]] = 
          treatlog_ACP[[i]][(treatlog_ACP[[i]] %notin% treatlog_ACP[[d]])]}}}
  } # end ACP assignment loop
  }#else{print("no CATIs were implemented")}
  
###################################################################
# Step 2: Update the states (disease dynamics and interventions)
###################################################################

# We use a Poisson tau-leap algorithm to sample the process at times. The number 
# of events of each time within an interval of duration is  approximately Poisson 
# distributed with mean equal to the rate of the transition at the start of the 
# interval divided by the duration of the interval. 
#
# The number of events is minimized (as the method can give rise to more events 
# than is possible i.e., > # of susceptible persons) using "min"
#
# To couple the meta-populations, we consider that infective individuals in 
# one patch can infect susceptible individuals in another patch, but there is no 
# explicit movement of individuals. This involves calculating the sum of infectious 
# persons in patches_j within the range of alpha of patch_i and scaling their 
# force of infection by the pre-calculated spatial force of infection. 
# As per Lloyd and May (1996), the force of infection on a given patch_i equals
# foi_i = sum(B_i_j * I_j), and is expressed here as intra-patch transmission +
# inter-patch transmission

# Setup empty dataframes to collect new exposed, infected, recovered individuals
# per timestep
#kernels <- data.frame(1)[,`:=`(c("p", "kernel.tot"),NA)][,V1:=NULL][.0]
kernels <- data.frame()
vaccinated <- data.frame()
exposed <- data.frame()
infected <- data.frame()
recovered <- data.frame()
catistatus <- data.frame()
sing_E.temp = matrix(0, ncol = nrow(unipix100), nrow = 2)
sing_I.temp = matrix(0, ncol = nrow(unipix100), nrow = 5)

# To sum up infectious pressure, identify patches_j which have infections
patches_infected <- data.frame(cbind(
  patchID = patches_status$patchID[colSums(sing_I)>=1],
  infecteds = colSums(sing_I)[colSums(sing_I)>=1]))
  
# Patch-transmission loop goes through the 2003 patches to figure out state 
# transitions given stochastic model and impact of some or all of the interventions
#
# Transmission is modeled as follows:
# ! Note that foi_spatial scales the impact of incident cases in other patches
# equivalent to dS/dt = -(beta*S_i/N + foi_spatial$S_i*sum(I_j)) * RR_POUWT
# equivalent to dE/dt =  (beta*S_i/N + foi_spatial$S_i*sum(I_j)) * inc.per*E * RR_OCV * RR_ACP)
# equivalent to dI/dt = inc.per*E - (1/inf.per or inf.per.acp)*I   
# equivalent to dR/dt = (1/inf.per or inf.per.acp)*I

# Run the patch loop across 2003 patches x 60 timesteps in parallel
#foreach(p = 1:nrow(patches_status),
#        .packages = c("tidyverse")) %dopar% {

#system.time(
for(p in 1:nrow(patches_status)){
  
  # To sum up infectious pressure, identify patches_j within the range, alpha (150m),
  # and exclude patch_i itself.
  #patches_j = (1:ncol(pixdistmat))[pixdistmat[ ,p] < alpha] 
  #patches_j = patches_j[patches_j %notin% p]  
  # sum up case counts in patches_j  
  #sing_I_j = sum(sing_I[ ,patches_j])
  
  # Calculate the power-law-distribution shaped isotropic kernel
  # To sum up infectious pressure, identify patches_j which have infections
  # and sum up the infectious pressure from those j patches in patch p
  
  # Patch-specific kernels to be summed up (reset for every patch)
  pkernel = data.frame()
  
  # Calculate the kernels from all patches_j with infections on patch_p (into pkernel)
  
  #pkernel <- foreach(j = 1:nrow(patches_infected), 
  #                   .packages = "dplyr",
  #                   .export = c('kernel.plk', 'patches_infected'),
  #                   .combine = 'bind_rows')%dopar%{
  
  for(j in 1:nrow(patches_infected)){
    pkernel.n <- data.frame(cbind(
      patchID = patches_infected$patchID[j], 
      infecteds = patches_infected$infected[j],
      kernel = kernel.plk[p, patches_infected$patchID[j]]))
    pkernel = bind_rows(pkernel, pkernel.n)
    }
  
  # Calculate master list of summed kernels from infected patches on a given patch, p
  kernels.n = data.frame(cbind(p, kernel.tot = sum(pkernel[,'kernel'])))
  # Per timestep, collect summed kernels for all patches_p
  kernels = bind_rows(kernels, kernels.n)
  
  # Note: Infection includes intra-patch and inter-patch transmission, i.e.:
  # 
  # exposed = foi*sing_S_i/N_i
  # (beta * foi.exp$foi[p] * sing_I_j)/(sing_S[p]+sum(sing_E[ ,p])+sum(sing_I[ ,p]) + sing_R[p])
  # *sing_S[p]/(sing_S[p]+sum(sing_E[ ,p])+sum(sing_I[ ,p]) + sing_R[p])
  
  
  # Cycle through options for each patch (ignore numbering):
  #
  # CATI
  # 0. No CATI
  # 1. Days 1-2: ACP + WASH (50m radius only)
  # 2. Days 3-6: WASH only (50m radius only)
  #    Days 1-6: WASH only
  # 3. Day  7 only: WASH + OCV (one-time susceptibility-reducing effect)
  # 4. Days 7-END: WASH + OCV (extended severe infection-reducing effect)
  # 5. Days 7-END: ACP + WASH + OCV (patch  previously in 150m radius has new case)
  #
  # Mass intervention
  # 0. No intervention
  # (1. No ACP / not applicable)
  # 2. Days 1-6: WASH only
  # 3. Day  7 only: WASH + OCV (one-time susceptibility-reducing effect)
  # 4. Days 7-END: WASH + OCV (extended severe infection-reducing effect)
  # (5. No ACP / not applicable)
  
  # Option 1: 50m, days 1-2 - treat with ACP and WASH  
  if(p%in%treatlog_ACP[[i]] & p%in%treatlog_WASH[[i]] & p%notin%treatlog_OCV[[i]]
     & p%notin%treatlog_OCV2[[i]]){
    
    #placeholder for no new vaccination
    vaccinated.n = data.frame(cbind(patchID=p, vaccinated=0))
    
    exposed.n = data.frame(cbind(patchID=p, exposed=min(sing_S[p], 
        rpois(1, (sing_S[p]*(beta*sum(sing_I[,p]) + kernels$kernel.tot[p]*sum(sing_I))/(sing_S[p] + sum(sing_E[ ,p]) + sum(sing_I[ ,p]) + sing_R[p])
        )*WASH_Eff))))
    infected.n = data.frame(cbind(patchID=p, infected=min(sing_E[2,p], 
        rpois(1, # symptomatic infection risk reduced by effects of ACP 
              1/IP*sing_E[2,p]*ACP_Eff)))) 
    recovered.n = data.frame(cbind(patchID=p, recovered=min(sing_I[5,p], 
        rpois(1, ((1/infper.acp)*sing_I[5,p])))))
    catistatus.n = data.frame(cbind(patchID=p, catistatus=1))
    
    vaccinated = bind_rows(vaccinated, vaccinated.n)
    exposed = bind_rows(exposed, exposed.n)
    infected = bind_rows(infected, infected.n)
    recovered = bind_rows(recovered, recovered.n)
    catistatus = bind_rows(catistatus, catistatus.n)}else

  # Option 2: WASH only; 50m, days 2-6 (antibiotic effect removed) OR 
  #           150m, days 1-6, treat with WASH only
  if(p%in%treatlog_WASH[[i]] & p%notin%treatlog_OCV[[i]] & p%notin%treatlog_ACP[[i]]
     & p%notin%treatlog_OCV2[[i]]){
        
    #placeholder for no new vaccination
    vaccinated.n = data.frame(cbind(patchID=p, vaccinated=0))
    
    exposed.n = data.frame(cbind(patchID=p, exposed=min(sing_S[p], 
        rpois(1, (sing_S[p]*(beta*sum(sing_I[,p]) + kernels$kernel.tot[p]*sum(sing_I))/(sing_S[p] + sum(sing_E[ ,p]) + sum(sing_I[ ,p]) + sing_R[p])
        )*WASH_Eff))))
    infected.n = data.frame(cbind(patchID=p, infected=min(sing_E[2,p], 
      rpois(1, 1/IP*sing_E[2,p])))) 
    recovered.n = data.frame(cbind(patchID=p, recovered=min(sing_I[5,p], 
      rpois(1, ((1/infper)*sing_I[5,p])))))
      catistatus.n = data.frame(cbind(patchID=p, catistatus=2))
        
    vaccinated = bind_rows(vaccinated, vaccinated.n)
    exposed = bind_rows(exposed, exposed.n)
    infected = bind_rows(infected, infected.n)
    recovered = bind_rows(recovered, recovered.n)
    catistatus = bind_rows(catistatus, catistatus.n)}else
  
  # Option 3: 50/150m, days 7+, treat with WASH and OCV (and one-time susceptibility-reducing effect)
  if(p%notin%treatlog_ACP[[i]] & p%in%treatlog_WASH[[i]] & p%in%treatlog_OCV[[i]]
    & p%in%treatlog_OCV2[[i]]){
        
        # Log the new vaccination
    vaccinated.n = data.frame(cbind(patchID=p, 
                   vaccinated=round(sing_S[p]*(1-OCV_Eff))))
    
    exposed.n = data.frame(cbind(patchID=p, exposed=min(sing_S[p], 
      rpois(1, (sing_S[p]*(beta*sum(sing_I[,p]) + kernels$kernel.tot[p]*sum(sing_I))/(sing_S[p] + sum(sing_E[ ,p]) + sum(sing_I[ ,p]) + sing_R[p])
      )*WASH_Eff))))
    infected.n = data.frame(cbind(patchID=p, infected=min(sing_E[2,p], 
      rpois(1, # symptomatic infection process reduced by effect of OCV
      1/IP*sing_E[2,p]*OCV_Eff)))) 
    recovered.n = data.frame(cbind(patchID=p, recovered=min(sing_I[5,p], 
      rpois(1, ((1/infper)*sing_I[5,p])))))
    catistatus.n = data.frame(cbind(patchID=p, catistatus=3))
        
    vaccinated = bind_rows(vaccinated, vaccinated.n)
    exposed = bind_rows(exposed, exposed.n)
    infected = bind_rows(infected, infected.n)
    recovered = bind_rows(recovered, recovered.n)
    catistatus = bind_rows(catistatus, catistatus.n)}else 
  
  # Option 4: 50/150m, days 7+, treat with WASH and OCV 
  if(p%notin%treatlog_ACP[[i]] & p%in%treatlog_WASH[[i]] & p%in%treatlog_OCV[[i]]
     & p%notin%treatlog_OCV2[[i]]){
    
    #placeholder for no new vaccination
    vaccinated.n = data.frame(cbind(patchID=p, vaccinated=0))
    
    exposed.n = data.frame(cbind(patchID=p, exposed=min(sing_S[p], 
        rpois(1, (sing_S[p]*(beta*sum(sing_I[,p]) + kernels$kernel.tot[p]*sum(sing_I))/(sing_S[p] + sum(sing_E[ ,p]) + sum(sing_I[ ,p]) + sing_R[p])
        )*WASH_Eff))))
    infected.n = data.frame(cbind(patchID=p, infected=min(sing_E[2,p], 
        rpois(1, # symptomatic infection process reduced by effect of OCV
              1/IP*sing_E[2,p]*OCV_Eff)))) 
    recovered.n = data.frame(cbind(patchID=p, recovered=min(sing_I[5,p], 
        rpois(1, ((1/infper)*sing_I[5,p])))))
    catistatus.n = data.frame(cbind(patchID=p, catistatus=4))
  
    vaccinated = bind_rows(vaccinated, vaccinated.n)
    exposed = bind_rows(exposed, exposed.n)
    infected = bind_rows(infected, infected.n)
    recovered = bind_rows(recovered, recovered.n)
    catistatus = bind_rows(catistatus, catistatus.n)}else 
  
    
  # Option 5: 50m, new antibiotic treatment in already treated patch (gets everything)
  if(p%in%treatlog_ACP[[i]] & p%in%treatlog_WASH[[i]] & p%in%treatlog_OCV[[i]]
     & p%notin%treatlog_OCV2[[i]]){
    
    #placeholder for no new vaccination
    vaccinated.n = data.frame(cbind(patchID=p, vaccinated=0))
    
    exposed.n = data.frame(cbind(patchID=p, exposed=min(sing_S[p], 
         rpois(1, (sing_S[p]*(beta*sum(sing_I[,p]) + kernels$kernel.tot[p]*sum(sing_I))/(sing_S[p] + sum(sing_E[ ,p]) + sum(sing_I[ ,p]) + sing_R[p])
                 )*WASH_Eff))))
    infected.n = data.frame(cbind(patchID=p, infected=min(sing_E[2,p], 
         rpois(1, # symptomatic infection process reduced by effect of OCV and ACP
               1/IP*sing_E[2,p]*ACP_Eff*OCV_Eff)))) 
    recovered.n = data.frame(cbind(patchID=p, recovered=min(sing_I[5,p], 
         rpois(1, ((1/infper.acp)*sing_I[5,p])))))
    
    catistatus.n = data.frame(cbind(patchID=p, catistatus=5))
  
    vaccinated = bind_rows(vaccinated, vaccinated.n)
    exposed = bind_rows(exposed, exposed.n)
    infected = bind_rows(infected, infected.n)
    recovered = bind_rows(recovered, recovered.n)
    catistatus = bind_rows(catistatus, catistatus.n)}else{
 
      
  # Option 0 -- no targeting/CATI
  
    #placeholder for no new vaccination
    vaccinated.n = data.frame(cbind(patchID=p, vaccinated=0))
      
    exposed.n = data.frame(cbind(patchID=p, exposed=min(sing_S[p], 
        rpois(1, sing_S[p]*(beta*sum(sing_I[,p]) + kernels$kernel.tot[p]*sum(sing_I))/
                (sing_S[p] + sum(sing_E[ ,p]) + sum(sing_I[ ,p]) + sing_R[p])))))
    infected.n = data.frame(cbind(patchID=p, infected=min(sing_E[2,p], 
         rpois(1, 1/IP*sing_E[2,p])))) 
    recovered.n = data.frame(cbind(patchID=p, recovered=min(sing_I[5,p], 
         rpois(1, ((1/infper)*sing_I[5,p])))))   
    catistatus.n = data.frame(cbind(patchID=p, catistatus=0))
    
    vaccinated = bind_rows(vaccinated, vaccinated.n)
    exposed = bind_rows(exposed, exposed.n)
    infected = bind_rows(infected, infected.n)
    recovered = bind_rows(recovered, recovered.n)
    catistatus = bind_rows(catistatus, catistatus.n)
    }
 
  } # Patch-transmission loop ends 
#)
# Update the states, per patch, per timestep

# Temp vectors
# Update susceptible compartments
sing_S.temp = sing_S - exposed$exposed - vaccinated$vaccinated
# Once exposed is removed from Susc, calculation and transfer asymptomatics 
# from exposed to recovered
infected.asymp = data.frame(infected.asymp = round(exposed$exposed * asymp.prop))
exposed$exposed = exposed$exposed - infected.asymp$infected.asymp
# Update exposed compartments
sing_E.temp[1,] = sing_E[1,] + exposed$exposed # add newly-exposed to d1
sing_E.temp[1,] = sing_E.temp[1,] - sing_E[1,] # minus d1 exposed from d1 (shuffle)
sing_E.temp[2,] = sing_E[2,] + sing_E[1,]      # add d1 exposed to d2 (shuffle)
sing_E.temp[2,] = sing_E.temp[2,] - infected$infected # minus newly-infected from d2
# Update infected compartments
sing_I.temp[1,] = sing_I[1,] + infected$infected # add newly-infected to d1 infected
sing_I.temp[1,] = sing_I.temp[1,] - sing_I[1,]   # minus d1 infected  
#
sing_I.temp[2,] = sing_I[2,] + sing_I[1,]        # add d1 infected to d2
sing_I.temp[2,] = sing_I.temp[2,] - sing_I[2,]   # minus d2 infected
#
sing_I.temp[3,] = sing_I[3,] + sing_I[2,]        # add d2 infected to d3
sing_I.temp[3,] = sing_I.temp[3,] - sing_I[3,]   # minus d3 infected 
#
sing_I.temp[4,] = sing_I[4,] + sing_I[3,]        # add d3 infected to d4
sing_I.temp[4,] = sing_I.temp[4,] - sing_I[4,]   # minus d4 infected
#
sing_I.temp[5,] = sing_I[5,] + sing_I[4,]        # add d4 infected to d5
sing_I.temp[5,] = sing_I.temp[5,] - recovered$recovered # minus recovered
# Update recovered compartment with asymptomatically infected and newly vaccinated
sing_R.temp = sing_R + recovered$recovered + infected.asymp$infected.asymp +
  vaccinated$vaccinated  # **checked and works

# Update the original compartments for this timestep with the temp files   
sing_S     = sing_S.temp
sing_E[1,] = sing_E.temp[1,]
sing_E[2,] = sing_E.temp[2,]
sing_I[1,] = sing_I.temp[1,]
sing_I[2,] = sing_I.temp[2,]
sing_I[3,] = sing_I.temp[3,]
sing_I[4,] = sing_I.temp[4,]
sing_I[5,] = sing_I.temp[5,]
sing_R     = sing_R.temp

# Update CATI status
catistatus.cum = catistatus

############################################################
## Step 3. Update on newly-infected individuals and patches
############################################################

# Identify newly-infected persons (in Infectious state, day 1) and 
# patches with new infections

new_detections <- as.data.frame(cbind((1:ncol(sing_I))[colSums(sing_I) > 0],
                                      (colSums(sing_I))[colSums(sing_I) > 0]))
colnames(new_detections) <- c("patchID", "cases")

# For t=1, substitute sing_I for sing_D
if(i==1){sing_D = sing_I}
# For t>1, indicate newly-affected patches with sing_D
else{
  sing_newI <- as.data.frame(cbind(location = seq(1,ncol(sing_I)), newcases = 0))
  sing_newI$newcases = sing_I[1,] 
  # Trim to only patches with cases
  sing_newI <- as.data.frame(cbind(location=sing_newI$location[sing_newI$newcases > 0],
                                   newcases=as.integer(sing_newI$newcases[sing_newI$newcases > 0])))
  # Disaggregate rows into as many rows as there are cases 
  sing_newI <- sing_newI %>% uncount(newcases, .remove = FALSE) %>%
    group_by(location) %>% mutate(newcases = 1) %>% ungroup
  sum(sing_newI$newcases)
  
  # Leave option open to randomly sample pdetect (%) of these infections 
  
  # This is non-sampled version where surveillance detects ALL new cases 
  #sing_newD <- as.data.frame(
  #  sing_newI[sample(nrow(sing_newI), round(sum(sing_newI$newcases)*1.0, digits=0)), ])
  
  # This is sampled version where surveillance detects 75% of new cases only
  sing_newD <- as.data.frame(
    sing_newI[sample(nrow(sing_newI), round(sum(sing_newI$newcases)*0.75, digits=0)), ])
  
  colnames(sing_newD) <- c('location', 'detected')
  # Match to sing_D matrix to indicate patches with cases detected; may be smaller than sum(sing_newD) if a patch has >1 infection
  sing_newD2 <- as.data.frame(cbind(location = seq(1,ncol(sing_I)), detected=0))
  sing_newD2$detected = sing_newD$detected[match(sing_newD2$location, sing_newD$location)]
  sing_newD2[is.na(sing_newD2)] = 0
  sing_D[1,] = sing_newD2$detected 
  } 

# For monitoring progress during execution.    
print(paste("timestep",i,"of", steps, sep = ' '))
#if("optim" %in% names(paramsList) | "nOCV" %in% names(paramsList)){
#  print(paste("new CATIs", round(sum(sing_D)/9),sep = ': '))}

########################
# STEP 4. Store results
########################

# System totals
totals[i, ] = c(sum(sing_S),
                sum(sing_E),
                sum(sing_I), 
                sum(sing_R), 
                sum(sing_D),
                sum(infected.asymp),
                sum(exposed$exposed),
                sum(infected$infected),
                sum(recovered$recovered),
                sum(with(catistatus.cum, catistatus==0)),
                sum(with(catistatus.cum, catistatus==1)),
                sum(with(catistatus.cum, catistatus==2)),
                sum(with(catistatus.cum, catistatus==3)),
                sum(with(catistatus.cum, catistatus==4)),
                sum(with(catistatus.cum, catistatus==5)),
                sum(with(catistatus.cum, catistatus==1 | catistatus==2 | catistatus==3 | catistatus==4 | catistatus==5)),
                round(sum(with(catistatus.cum, catistatus==1 | catistatus==2 | catistatus==3 | catistatus==4 | catistatus==5))/9))

# Patch totals
detailed_totals[[i]][ ,1:8] = cbind(sing_S,
                                    colSums(sing_E),
                                    colSums(sing_I), 
                                    sing_R, 
                                    colSums(sing_D),
                                    colSums(infected.asymp),
                                    colSums(infected), # added to document new infections per patch
                                    kernels$kernel.tot)

# Break timestep loop if sing_I drops to zero
if((sum(sing_I) == 0)){
  print(paste("failed at timestep", i, "due to sing_I dropping to zero"))
  break}

# Break timestep loop if new cases drops to zero
#if((sing_newI$newcases == 0)){
#  print(paste("failed at timestep", i, "due to no new cases"))
#  break}

} # End timestep loop at t = steps   

close(pb)


# Compile treatlog
treatlog = list(cholspots = cholspots_log,
                cholspots.ACP = cholspots_ACP_log,
                prevcholspots = prevcholspots_log,
                ACP = treatlog_ACP,
                WASH = treatlog_WASH,
                OCV = treatlog_OCV,
                OCV2 = treatlog_OCV2)

# Compile final results
templist <- list(sing_S, sing_E, sing_I, sing_R, sing_D, infected.asymp,
                 exposed, infected, recovered,
                 catistatus, totals, detailed_totals, treatlog)
names(templist) = c("sing_S", "sing_E", "sing_I", "sing_R", "sing_D", "asymp",
                    "exposed", "infected", "recovered",
                    "catistatus", "totals", "detailed_totals", "treatlog")
mainmodrun <- templist

if(length(mainmodrun) > 1){
  model_outlist <- list()
  model_outlist[[1]] <- mainmodrun$totals
  model_outlist[[2]] <- mainmodrun$detailed_totals
  model_outlist[[3]] = mainmodrun$treatlog
  names(model_outlist) = c("Poptotals", "Patchtotals", "treatlog")
  return(model_outlist)
  }else{return(NA)}
} # end:: model function



saveRDS(mainmodrun, file="mainmodrun.RData")
check <- readRDS("mainmodrun.RData")