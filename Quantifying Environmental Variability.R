# This script is for estimating values of environmental variability using the precipitation and temperature records extracted from CHELSA.
# this extraction will be done in the form of a spectrum analysis to determine the spectrum of noise in the annual variation
# Date last modified: December 2020
# Author: James Cant

#Clear the workspace
rm(list=ls(all=TRUE))

# load required packages
library(raster)
library(rgdal)
library(ncdf4)
library(stringr) 
library(dplyr)
library(colorednoise)
library(stats)
library(data.table)

# set working directory
setwd("FILE_PATH")

# and set pathway for reading in extracted CHELSA data.
ARC.store <- "PATH_TO_EXTRACTED_CHELSA_DATA"

# load in data with species names automatically set as row names.
data <- read.csv("RawData.csv", stringsAsFactors = FALSE)

###### --------------------------------------------------------------------
# Because of the way I am merging the climate data into the main dataframe below I need to here remove the marine environmental data already in the raw data file
# and reattach it along with the terrestrial data - it streamlines the code doing it this way.
# add columns for additional precipitation variables
data$Prec.Frequency.full <- NA
data$Prec.Autocorrelation.full <- NA
data$Prec.Frequency.inter <- NA
data$Prec.Autocorrelation.inter <- NA
data$Prec.Frequency.intra <- NA
data$Prec.Autocorrelation.intra <- NA
# Extract the environmental data already attached for marine populations
climate_marine <- data[which(data$Realm == "Marine"), c("SpeciesAccepted","StartYear","EndYear","Lat",                       
                                                        "Lon","Frequency.full","Autocorrelation.full","range.full",           
                                                        "Frequency.inter","Autocorrelation.inter","range.inter","Frequency.intra",      
                                                        "Autocorrelation.intra","range.intra","Prec.Frequency.full","Prec.Autocorrelation.full", 
                                                        "Prec.Frequency.inter","Prec.Autocorrelation.inter","Prec.Frequency.intra","Prec.Autocorrelation.intra")]
# match the colnames in climate marine with what they will be called in the overall climate output
colnames(climate_marine) <-  c("SpeciesAccepted", "StartYear", "EndYear", "Lat", "Lon", 
                               "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", 
                               "Temp.Frequency.inter", "Temp.Autocorrelation.inter", "Temp.range.inter", 
                               "Temp.Frequency.intra", "Temp.Autocorrelation.intra", "Temp.range.intra", 
                               "Prec.Frequency.full",  "Prec.Autocorrelation.full", "Prec.Frequency.inter", 
                               "Prec.Autocorrelation.inter", "Prec.Frequency.intra",  "Prec.Autocorrelation.intra")
# And remove the environmental columns from the raw data
data[, c("X.1", "X", "Frequency.full","Autocorrelation.full","range.full", "Frequency.inter","Autocorrelation.inter",
         "range.inter","Frequency.intra", "Autocorrelation.intra","range.intra","Prec.Frequency.full", 
         "Prec.Autocorrelation.full", "Prec.Frequency.inter","Prec.Autocorrelation.inter","Prec.Frequency.intra",
         "Prec.Autocorrelation.intra")] <- NULL
###### --------------------------------------------------------------------

# load in required climate data
prec1 <- read.csv(paste0(ARC.store,"/Chelsa Precipitation Data 1901_1957.csv", sep = ""), row.names = 1)
prec2 <- read.csv(paste0(ARC.store,"/Chelsa Precipitation Data 1958_2013.csv", sep = ""), row.names = 1)
tmax1 <- read.csv(paste0(ARC.store,"/Chelsa max Temperature Data 1901_1957.csv", sep = ""), row.names = 1)
tmax2 <- read.csv(paste0(ARC.store,"/Chelsa max Temperature Data 1958_2013.csv", sep = ""), row.names = 1)
tmin1 <- read.csv(paste0(ARC.store,"/Chelsa min Temperature Data 1901_1957.csv", sep = ""), row.names = 1)
tmin2 <- read.csv(paste0(ARC.store,"/Chelsa min Temperature Data 1958_2013.csv", sep = ""), row.names = 1)

# reshape the data to allow it to be combined together
prec2 <- cbind(prec1[,(dim(prec1)[2]-53):dim(prec1)[2]], prec2) # add the final 54 columns to the second data frame for each variable
tmax2 <- cbind(tmax1[,(dim(tmax1)[2]-53):dim(tmax1)[2]], tmax2)
tmin2 <- cbind(tmin1[,(dim(tmin1)[2]-53):dim(tmin1)[2]], tmin2)
prec1 <- select(prec1, -(c(dim(prec1)[2]-53):dim(prec1)[2])) # and snip them from the first
tmax1 <- select(tmax1, -(c(dim(tmax1)[2]-53):dim(tmax1)[2])) 
tmin1 <- select(tmin1, -(c(dim(tmin1)[2]-53):dim(tmin1)[2])) 

# Merge together the environmental datafiles.
# first make sure colnames match
colnames(prec2) <- colnames(tmax1) <- colnames(tmin1) <- colnames(tmax2) <- colnames(tmin2) <- names(prec1)
climate_all <- rbind(prec1,prec2,tmax1,tmax2,tmin1,tmin2)

# save the raw climate data as a checkpoint
write.csv(climate_all, file = "FILE_PATH/Raw climate data.csv", row.names = FALSE)

# save memory space 
rm(prec1,prec2,tmax1,tmax2,tmin1,tmin2)

# first split the datafile into a list of the different populations (some populations have different start and end years but are in the same locations)
climate_df <- split(climate_all, list(climate_all$Lon, climate_all$Lat, climate_all$SpeciesAccepted, climate_all$StartYear, climate_all$EndYear), drop = TRUE)

# create a storage output for the environmental variables.
climate_output <- data.frame("SpeciesAccepted" = rep(NA, length(climate_df)),
                             "StartYear" = rep(NA, length(climate_df)),
                             "EndYear" = rep(NA, length(climate_df)),
                             "Lat" = rep(NA, length(climate_df)),
                             "Lon" = rep(NA, length(climate_df)),
                             "Temp.Frequency.full" = rep(NA, length(climate_df)),
                             "Temp.Autocorrelation.full" = rep(NA, length(climate_df)),
                             "Temp.range.full" = rep(NA, length(climate_df)),
                             "Temp.Frequency.inter" = rep(NA, length(climate_df)),
                             "Temp.Autocorrelation.inter" = rep(NA, length(climate_df)),
                             "Temp.range.inter" = rep(NA, length(climate_df)),
                             "Temp.Frequency.intra" = rep(NA, length(climate_df)),
                             "Temp.Autocorrelation.intra" = rep(NA, length(climate_df)),
                             "Temp.range.intra" = rep(NA, length(climate_df)),
                             "Prec.Frequency.full" = rep(NA, length(climate_df)),
                             "Prec.Autocorrelation.full" = rep(NA, length(climate_df)),
                             "Prec.Frequency.inter" = rep(NA, length(climate_df)),
                             "Prec.Autocorrelation.inter" = rep(NA, length(climate_df)),
                             "Prec.Frequency.intra" = rep(NA, length(climate_df)),
                             "Prec.Autocorrelation.intra" = rep(NA, length(climate_df)))

# now work through each GPS point and determine the environmental variability
for (x in 1:length(climate_df)) {
  
  #Select GPS location
  pop.use <- climate_df[[x]]
  #store in output for ease later on
  climate_output$SpeciesAccepted[x] <- paste(pop.use$SpeciesAccepted[1])
  climate_output$StartYear[x] <- pop.use$StartYear[1]
  climate_output$EndYear[x] <- pop.use$EndYear[1]
  climate_output$Lat[x] <- pop.use$Lat[1]
  climate_output$Lon[x] <- pop.use$Lon[1] 
  
  #define start year and legacy year for the population
  start.norm <- pop.use$StartYear[1] - 50
  if(start.norm < 1901) {start.norm = 1901} # little fix because there is no data from before 1901
  end.norm <- pop.use$EndYear[1]
  #remove un-needed years
  pop.use <- pop.use[which(pop.use$year >= start.norm),]
  pop.use <- pop.use[which(pop.use$year <= end.norm),]
  #covert year values for grouping
  pop.use$year <- as.numeric(paste(pop.use$year))
  
  # Now separate out precipitation and temperature data as they need to be processed differently.
  pop.prec <- pop.use[which(pop.use$variable == "prec"),]
  pop <- pop.use[which(pop.use$variable != "prec"),] # min and max temperature records only
  
  # reformat focal data frames
  col.names <- names(pop.use[,1:9])
  # define indexing
  index.start <- seq(10,length(pop), by = 9)
  index.end <- seq(18,length(pop), by = 9)
  # define starting point
  pop.use <- pop[,1:9]
  pop.use.prec <- pop.prec[,1:9]
  
  # stack every nine columns on top of each other.
  for (ii in 1:length(index.start)) {
    # extract next chuck of data
    pop.extract <- pop[,index.start[ii]:index.end[ii]]
    pop.prec.extract <- pop.prec[,index.start[ii]:index.end[ii]]
    # reassign colnames
    names(pop.extract) <- names(pop.prec.extract) <- col.names
    # and stack it
    pop.use <- rbind(pop.use, pop.extract)
    pop.use.prec <- rbind(pop.use.prec, pop.prec.extract)
  }
  
  ######### Estimate variance metrics for the temperature data -------------------------------------------------
  
  # temperature in the records pre 1979 are recorded as 10ths of a degree, temperature recorded post 1979 are recorded 100ths of a degree - I therefore need to correct for that here
  pop.use[which(pop.use$year < 1979),]$value <- pop.use[which(pop.use$year < 1979),]$value/10
  pop.use[which(pop.use$year >= 1979),]$value <- pop.use[which(pop.use$year >= 1979),]$value/100
  
  # summarise annual variance (both inter and intra annual variance)
  try(pop.use.1 <- group_by(pop.use, month, year) %>%  
        summarise(MeanTemp = mean(value, na.rm = TRUE),
                  MaxTemp = max(value, na.rm = TRUE),
                  MinTemp = min(value, na.rm = TRUE),
                  TempRange = (MaxTemp - MinTemp)))
  # sometimes R orders the data by month not year so I need to correct for this before creating the temperature time series
  try(pop.use.1 <- pop.use.1[order(pop.use.1$year),])
  # generate the the environmental time series
  try(temp_seq <- ts(pop.use.1$MeanTemp, start = start.norm, frequency = 12))
  # Caclulate the measures of variability
  try(climate_output$Temp.Autocorrelation.full[x] <- autocorrelation(temp_seq, biasCorrection = TRUE))
  try(climate_output$Temp.range.full[x] <- mean(pop.use.1$TempRange, na.rm = TRUE))
  #####
  try(spectra <- spectrum(temp_seq, plot = F))
  try(spectra.x <- as.numeric(spectra[["freq"]]))
  try(spectra.y <- as.numeric(spectra[["spec"]]))
  try(spectra.mod <- lm(log(spectra.y)~log(spectra.x)))
  try(spectral.expo <- as.numeric(coef(spectra.mod)[2]))
  try(climate_output$Temp.Frequency.full[x] <- spectral.expo)
  
  # repeat for interannual variation
  try(pop.use.2 <- group_by(pop.use, year) %>%  
        summarise(MeanTemp = mean(value, na.rm = TRUE),
                  MaxTemp = max(value, na.rm = TRUE),
                  MinTemp = min(value, na.rm = TRUE),
                  TempRange = (MaxTemp - MinTemp)))
  try(pop.use.2 <- pop.use.2[order(pop.use.2$year),])
  # generate the the environmental time series
  try(temp_seq2 <- ts(pop.use.2$MeanTemp, start = start.norm, frequency = 1))
  # Caclulate the measures of variability
  try(climate_output$Temp.Autocorrelation.inter[x] <- autocorrelation(temp_seq2, biasCorrection = TRUE))
  try(climate_output$Temp.range.inter[x] <- mean(pop.use.2$TempRange, na.rm = TRUE))    
  #####
  try(spectra2 <- spectrum(temp_seq2, plot = F))
  try(spectra.x2 <- as.numeric(spectra2[["freq"]]))
  try(spectra.y2 <- as.numeric(spectra2[["spec"]]))
  try(spectra.mod2 <- lm(log(spectra.y2)~log(spectra.x2)))
  try(spectral.expo2 <- as.numeric(coef(spectra.mod2)[2]))
  try(climate_output$Temp.Frequency.inter[x] <- spectral.expo2)
  
  # repeat for intraannual variation
  try(pop.use.3 <- group_by(pop.use, month) %>%  
        summarise(MeanTemp = mean(value, na.rm = TRUE), 
                  MaxTemp = max(value, na.rm = TRUE),
                  MinTemp = min(value, na.rm = TRUE),
                  TempRange = (MaxTemp - MinTemp)))
  try(pop.use.3 <- pop.use.3[order(pop.use.3$month),])
  # generate the the environmental time series
  try(temp_seq3 <- ts(pop.use.3$MeanTemp, start = 1, frequency = 1))
  # Caclulate the measures of variability
  try(climate_output$Temp.Autocorrelation.intra[x] <- autocorrelation(temp_seq3, biasCorrection = TRUE))
  try(climate_output$Temp.range.intra[x] <- mean(pop.use.3$TempRange, na.rm = TRUE))
  #####
  try(spectra3 <- spectrum(temp_seq3, plot = F))
  try(spectra.x3 <- as.numeric(spectra3[["freq"]]))
  try(spectra.y3 <- as.numeric(spectra3[["spec"]]))
  try(spectra.mod3 <- lm(log(spectra.y3)~log(spectra.x3)))
  try(spectral.expo3 <- as.numeric(coef(spectra.mod3)[2]))
  try(climate_output$Temp.Frequency.intra[x] <- spectral.expo3)
  
  ######### Repeat for the precipitation data -------------------------------------------------
  
  # The precipitation data doesn't require reformatting and is consistent throughout the aruts and timeseries records
  # summarise annual variance (both inter and intra annual variance)
  try(pop.use.prec.1 <- group_by(pop.use.prec, month, year) %>%  
        summarise(MeanPrec = mean(value, na.rm = TRUE)))
  # sometimes R orders the data by month not year so I need to correct for this before creating the temperature time series
  try(pop.use.prec.1 <- pop.use.prec.1[order(pop.use.prec.1$year),])
  # generate the the environmental time series
  try(prec_seq <- ts(pop.use.prec.1$MeanPrec, start = start.norm, frequency = 12))
  # Caclulate the measures of variability
  try(climate_output$Prec.Autocorrelation.full[x] <- autocorrelation(prec_seq, biasCorrection = TRUE))
  #####
  try(pspectra <- spectrum(prec_seq, plot = F))
  try(pspectra.x <- as.numeric(pspectra[["freq"]]))
  try(pspectra.y <- as.numeric(pspectra[["spec"]]))
  try(pspectra.mod <- lm(log(pspectra.y)~log(pspectra.x)))
  try(pspectral.expo <- as.numeric(coef(pspectra.mod)[2]))
  try(climate_output$Prec.Frequency.full[x] <- pspectral.expo)
  
  # repeat for interannual variation
  try(pop.use.prec.2 <- group_by(pop.use.prec, year) %>%  
        summarise(MeanPrec = mean(value, na.rm = TRUE)))
  try(pop.use.prec.2 <- pop.use.prec.2[order(pop.use.prec.2$year),])
  # generate the the environmental time series
  try(prec_seq2 <- ts(pop.use.prec.2$MeanPrec, start = start.norm, frequency = 1))
  # Caclulate the measures of variability
  try(climate_output$Prec.Autocorrelation.inter[x] <- autocorrelation(prec_seq2, biasCorrection = TRUE))
  #####
  try(pspectra2 <- spectrum(prec_seq2, plot = F))
  try(pspectra.x2 <- as.numeric(pspectra2[["freq"]]))
  try(pspectra.y2 <- as.numeric(pspectra2[["spec"]]))
  try(pspectra.mod2 <- lm(log(pspectra.y2)~log(pspectra.x2)))
  try(pspectral.expo2 <- as.numeric(coef(pspectra.mod2)[2]))
  try(climate_output$Prec.Frequency.inter[x] <- pspectral.expo2)
  
  # repeat for intraannual variation
  try(pop.use.prec.3 <- group_by(pop.use.prec, month) %>%  
        summarise(MeanPrec = mean(value, na.rm = TRUE)))
  try(pop.use.prec.3 <- pop.use.prec.3[order(pop.use.prec.3$month),])
  # generate the the environmental time series
  try(prec_seq3 <- ts(pop.use.prec.3$MeanPrec, start = 1, frequency = 1))
  # Caclulate the measures of variability
  try(climate_output$Prec.Autocorrelation.intra[x] <- autocorrelation(prec_seq3, biasCorrection = TRUE))
  #####
  try(pspectra3 <- spectrum(prec_seq3, plot = F))
  try(pspectra.x3 <- as.numeric(pspectra3[["freq"]]))
  try(pspectra.y3 <- as.numeric(pspectra3[["spec"]]))
  try(pspectra.mod3 <- lm(log(pspectra.y3)~log(pspectra.x3)))
  try(pspectral.expo3 <- as.numeric(coef(pspectra.mod3)[2]))
  try(climate_output$Prec.Frequency.intra[x] <- pspectral.expo3)
  
} # end of extraction loop

# Combine the extracted terrestiral environmental data with the corresponding marine data
climate_combined <- rbind(climate_output, climate_marine)
# Then clean the data by removing the silly infs and NaNs 
climate_combined[sapply(climate_combined, is.infinite)] <- NA
climate_combined[sapply(climate_combined, is.nan)] <- NA #this converts all silly nan and inf to NAs.
climate_combined$SpeciesAccepted <- as.character(climate_combined$SpeciesAccepted)

# Now merge the environmental data into the main dataframe
full_data <- merge(data, unique(climate_combined), by = intersect(names(data), names(climate_combined)), all.x = TRUE)

# just in case any got through remove the silly infs and NaNs 
full_data[sapply(full_data, is.infinite)] <- NA
full_data[sapply(full_data, is.nan)] <- NA #this converts all silly nan and inf to NAs.

# and export file.
write.csv(full_data, file = "RawData2.csv", row.names = F) # Checkpoint

#***************************** End of Code ****************************** 
