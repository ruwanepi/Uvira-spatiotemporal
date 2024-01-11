# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script sets up the packages for the analysis, so
# run it first.
# Author: R Ratnayake, 2023

# Load packages

x1 <- c("here", "tidyverse", "lubridate", "skimr", "ggplot2", "sf",
        "tmap", "tsibble", "EpiEstim", "rsatscan", "xlsx", "rgdal",
        "cowplot", "ggmap", "lubridate", "incidence2", "raster",
        "exactextractr", "R0", "SciViews", "imputeTS", "feasts",
        "mgcv", "IDSpatialStats", "tidyquant", "slider", "surveillance",
        "ggspatial", "grid", "gridExtra", "scales", "geosphere", "readr")
# Install any packages not yet installed
x2 <- x1 %in% row.names(installed.packages())
if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
lapply(x1, library, character.only = TRUE)

# Set working directory
setwd("C:/R-projects/Uvira-spatial-epid/")