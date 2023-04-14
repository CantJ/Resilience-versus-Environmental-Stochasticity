# This script is for extracting environmental data (Temperature and precipitation) from CHELSA for terrestrial populations
# Date last modified: December 2020
# Author: James Cant

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(data.table)
library(raster)
library(dplyr)
library(tiff)
library(parrallel)
library(rgdal)

# set working directory
setwd("FILE_PATH")
store <- "FILE_PATH" #this is for the chelsa extraction

# load in data with species names automatically set as row names, and for which marine abiotic data has already been downloaded.
data <- read.csv("CSVFILE.csv", stringsAsFactors = FALSE)

# Separate out coordinates of interest 
lat_lon <- subset(data, data$Realm %in% c("Terrestrial", "Aquatic"))
# shrink the data file
lat_lon <- lat_lon[, c("Lat", "Lon", "StartYear", "EndYear", "SpeciesAccepted")]
# remove entries missing Year information
lat_lon <- na.omit(lat_lon)

# extract only site coordinates 
site_coords <- dplyr::select(lat_lon, Lon, Lat) %>%  #this creates a matrix of just the date and location info for each population of interest - this is the only information the CHELSA estraction will require.
  as.matrix( dimnames = list(rep('value',nrow(lat_lon)),
                             c('Long', 'Lat')) )

# CHELSA extraction set up 
# set up reading directory
read_dir <- 'https://envidatrepo.wsl.ch/uploads/chelsa/chelsa_V1/timeseries/' #this searches the more accurate CHELSA time series for between the years 1979 and 2015
read_dir2 <- 'https://envidatrepo.wsl.ch/uploads/chelsa/chelsa_V1/chelsa_cruts/' #this searches the CHELSAcrusts time series which has data (though at a lower accuracy for 1901 - 1978)

# create functions for downloading spatial files from the internet 

# 1. produce file name based on index 'ii'
produce_file_name <- function(ii){
  
  if(chelsa_df$year[ii] >= 1979) {
    file_root  <- paste0(chelsa_df$variable[ii],'/CHELSA_',chelsa_df$variable[ii],'_',
                         chelsa_df$year[ii],'_',
                         chelsa_df$month[ii],'_V1.2.1.tif')
  } else {
    file_root  <- paste0(chelsa_df$variable[ii],'/CHELSAcruts_',chelsa_df$variable[ii],'_',
                         as.numeric(chelsa_df$month[ii]),'_',
                         chelsa_df$year[ii],'_V.1.0.tif')  
  }
  
  return(file_root)
  
}

# 2. function for extracting year and monthly data
extract_year_month <- function(ii){
  
  # extract with archive 
  file_path <- file_links[ii]
  download.file(file_path, destfile = file_dest[ii], mode = "wb")
  
  ##### now get climate information from the downloaded .tif file
  
  # This section is for if the package rgdal is accessible ----------------------
  # read raster
  raster_file <- grep('.tif|.sdat$',list.files(store), value=T)
  repP        <- raster(paste0(store,'/',raster_file))
  # -----------------------------------------------------------------------------
  
  # This section is for if rgdal isnt working -----------------------------------
  # Identify relevant file name
  raster_file <- grep('.tif|.sdat$',list.files(store), value=T)
  # readTIFF opens the file as a raster but in a way that makes it unreadable by spatial R packages (hence why it is being repacked as a matrix for re-rastering)
  repP        <- readTIFF(paste0(store,'/',raster_file), as.is = TRUE) #this will open the tiff file using unsigned pixels
  repP[repP >= 32768] <- -65536 + repP[repP >= 32768] # this re signs them.
  repP[repP < -30798] <- NA  # in chelsa extremely low negative values denote NAs
  # convert the native raster to a matrix for re-rastering
  repP        <- as.matrix(repP)
  # and re-raster
  repP        <- raster(repP, xmn = -180.0001, xmx = 179.9999, ymn = -90.00014, ymx = 83.99986, crs = "+proj=longdat +datum=WGS84 +no_defs")
  # -----------------------------------------------------------------------------
  
  # extract info 
  values_clim <- raster::extract(repP, site_coords)
  clim_df     <- lat_lon %>% 
    mutate( variable = chelsa_df$variable[ii],
            year     = chelsa_df$year[ii],
            month    = chelsa_df$month[ii],
            value    = values_clim)
  
  file.remove( file_dest[ii] )
  
  print(ii)
  
  return(clim_df)
  
}

# Extract CHELSA DATA ------------------------------------------------------------------------

# what do I need from CHELSA for this population?
chelsa_df <- expand.grid( variable = c('tmax','tmin','prec'), # when submitting a job to the ARC cluster there is too much data for this to all be done at once so these variables need to be submitted seperatly.
                          year     = c((min(lat_lon$StartYear)-50):1957), # equally even within each variable there is too much data needed so the time frame for extraction needs seperating accordingly.
                         # year     = c(1958:max(lat_lon$EndYear)),
                          month    = c(paste0('0',1:9),paste0(10:12)),
                          stringsAsFactors = F) %>% 
  arrange(variable,year,month)
chelsa_df <- chelsa_df[which(chelsa_df$year <= 2013),] # there is no data after 2013

# get all CHELSA file links (from file name)
file_names <- lapply(1:nrow(chelsa_df), produce_file_name) %>% unlist #determine file name
file_links <- rep(NA, length(file_names)) #create blank links vector

# calculate the links based on the file name and the year needed.
for (z in 1:length(file_names)) {
  if (chelsa_df$year[z] < 1979) {
    file_links[z] <- paste0(read_dir2,file_names[z])
  } else {
    file_links[z] <- paste0(read_dir,file_names[z])
  }
}

#assign file destination names
file_dest  <- gsub("tmin/|tmax/|prec/","",file_names)

# Now run the loop to extract all climate data for the areas identified in the chelsa_df object.
start <- Sys.time() #just a little loop timer
range_clim  <- 1:dim(chelsa_df)[1]

# run the function on all identified month year combinations for the GPS sites needed.
try(climate_all <- lapply(range_clim, extract_year_month))
Sys.time() - start #how long was taken?

# and export file - the line of code used will depend on the data being extracted
write.csv(climate_all = "Chelsa max Temperature data 1901_1957.csv")
write.csv(climate_all = "Chelsa max Temperature data 1957_2013.csv")
write.csv(climate_all = "Chelsa min Temperature data 1901_1957.csv")
write.csv(climate_all = "Chelsa min Temperature data 1957_2013.csv")
write.csv(climate_all = "Chelsa Precipitation data 1901_1957.csv")
write.csv(climate_all = "Chelsa Precipitation data 1957_2013.csv")

#**************************** End of Code ****************************** 
