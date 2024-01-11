# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script sets up the utility functios for tau.
# Author: R Ratnayake, 2023, based on A Azman, 2018 (code: https://osf.io/4fsnc/)

# Utility functions for tau
# Note: Running tau.analysis.R will automatically run tau.utils.R.

## loads point data from all suspected cases (2020)
load_suspected_cases2020 <- function(){
  sc2020 <- read_csv("cleaned_simplified_points_uvira2a.csv")
  return(sc)
}

## loads point data from RDT-positive cases alone (2020)
load_RDTpos_cases2020 <- function(){
  pc2020 <- read_csv("cleaned_simplified_points_uvira3.csv")
  return(pc2020)
}

## loads point data from all suspected cases (2019)
load_suspected_cases2019 <- function(){
  sc2019 <- read_csv("cleaned_simplified_points_uvira2019.csv")
  return(sc2019)
}

## loads point data from RDT-positive cases alone (2019)
load_RDTpos_cases2019 <- function(){
  pc2019 <- read_csv("cleaned_simplified_points_uvira2019rdt.csv")
  return(pc2019)
}

## loads point data from RDT-positive cases alone (2018)
load_RDTpos_cases2018 <- function(){
  pc2018 <- read_csv("cleaned_simplified_points_uvira2018rdt.csv")
  return(pc2018)
}

## loads point data from RDT-positive cases alone (2017)
load_RDTpos_cases2017 <- function(){
  pc2017 <- read_csv("cleaned_simplified_points_uvira2017rdt.csv")
  return(pc2017)
}

## loads point data from RDT-positive cases alone (2016)
load_RDTpos_cases2016 <- function(){
  pc2016 <- read_csv("cleaned_simplified_points_uvira2016rdt.csv")
  return(pc2016)
}

## loads point data from ALL RDT-positive cases alone (2016-2020)
load_RDTpos_cases20162020 <- function(){
  pcall <- read_csv("cleaned_simplified_points_uvira20162020rdt.csv")
  return(pcall)
}

##' function for main analysis
##' @param r.mins vector of minimum distances for tau windows
##' @param r.maxs vector of maximum distances for tau windows
##' @param d.mins vector of min time windows for tau 
##' @param d.maxs vector of max time windows for tau 
##' @param point_mat matrix of x,y and time
##' @param r.mids vector of mid-point of distance windows for tau
##' @param d.mids vector of mid-point time windows for tau 

run_main_analysis <- function(#r.mins=c(0, seq(5,950,10)),
                              #r.maxs=c(40, seq(55,1000,10)),
                              r.mins = (r.maxs - 50),
                              r.maxs = seq(50,5000,10),
                              d.mins=c(0,5,10,15,20,25,15,1,0),
                              d.maxs=c(5,10,15,20,25,30,30,5,2),
                              point_mat,
                              filename){
  
  #convert df into matrix if needed
  if(!is.matrix(pcall)) {pcall <- as.matrix(pcall)} 
  #if(!is.matrix(pc2016)) {pc2016 <- as.matrix(pc2016)} 
  
  set.seed(1981)
  
  r.mids <- (r.maxs + r.mins)/2 #calculate midpoints for distances
  d.mids <- (d.maxs + d.mins)/2 #calculate midpoints for times
  rc <- list() 
  
  for(i in seq_along(d.mins)){
    #combined function to produce tau and CIs with bootstrapping
    rc[[i]] <- get_tau_and_cis(day_min=d.mins[i],     
                               day_max=d.maxs[i],
                               my.mat=point_mat,
                               r.max=r.maxs,
                               r.min=r.mins,
                               n_boots=1000,
                               infectious.process=TRUE #considers only infectious cases as those after the index case
    )
    print(i)
  }
  
  saveRDS(rc,file = filename)  
}

##' gets tau estimates and confidence intervals using 
##' IDDSpatial Stats implementation of tau
##' @param n_boots number of bootstrap iterations
##' @param infectious.process flag if we only want to compute tau considering related cases as those that occur **after** an *index* case
##' @param full_boots flag to return full bootstrap replicate matrix (TRUE) or relevant quantiles for confidence intervals
get_tau_and_cis <- function(day_min,
                            day_max,
                            my.mat,
                            r.max,
                            r.min,
                            n_boots,
                            infectious.process=TRUE,
                            full_boots=FALSE){
  
  
  ## make new relation function 
  if(infectious.process){
    new.gen.func <- function(r1,r2,
                             within.range=c(day_min,day_max)){
      is.in.generation(r1,r2,
                       within.range,
                       infectious=TRUE)
    }
  } else {
    
    new.gen.func <- function(r1,r2,
                             within.range=c(day_min,day_max)){
      is.in.generation(r1,r2,
                       within.range,
                       infectious=FALSE)
    }
  }
  
  
  ## get tau estimate
  tau <- get.tau(posmat=my.mat,
                 fun=new.gen.func,
                 r = r.max,
                 r.low=r.min,
                 comparison.type="independent")
  
  if(!full_boots){
    ## get cis rather than full boot
    boot <- get.tau.ci(
      posmat=my.mat,
      fun=new.gen.func,
      r = r.max,
      r.low=r.min,
      boot.iter=n_boots,
      comparison.type="independent",
      ci.low=0.025,
      ci.high=0.975
    )
  } else {
    ## get full bootstrap estimates matrix
    boot <- get.tau.bootstrap(
      posmat=my.mat,
      fun=new.gen.func,
      r = r.max,
      r.low=r.min,
      boot.iter=n_boots,
      comparison.type="independent")
    
  }
  
  return(list(tau,boot))
}

## the tau function needs a 'relation function' to
## determine how pairs are related
##' @param row1 - case line 1 (expecting 3rd element to be time in vector)
##' @param row2 - case line 2
##' @param within.range - days that case should occur within (lower and upper)
##' @param infectious - if infectious, then we only look at time in one direction
##' @return
##' @author asa
is.in.generation <- function(row1,
                             row2,
                             within.range=c(10,20), # 2 x serial interval
                             infectious=TRUE){
  if(infectious){
    
    if((row1[3]-row2[3]) >= within.range[1] & (row1[3]-row2[3]) < within.range[2]){rc=1}
    else{rc=2}
    
  } else {
    if(abs(row1[3]-row2[3]) >= within.range[1] & abs(row1[3]-row2[3]) < within.range[2]){rc=1}
    else {rc=2}
  }
  
  return(rc)
}