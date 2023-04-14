# This script will work through all relevant demographic data to sort and remove outliers, before then transforming data as necessary. 
# Following this the script carries out a PCA displaying the environmental variability experienced by the selected populations.
# This PCA will only be used as supplementary material.
# Author: James Cant
# Last Modified: December 2020

# clear workspace
rm(list=ls(all=TRUE))

# load required packages
library(dplyr)
library(naniar)
library(RColorBrewer)
library(plotrix)
library(ggplot2)
library(ape)
library(geiger)
library(caret)
library(fuzzySim)

#set working directory
setwd("FILE_PATH")

# load in data file
Raw_data <- read.csv("RawData2.csv", stringsAsFactors = FALSE, row.names = 185)
# remove R's added rows.
Raw_data$X <- NULL
# load in the tree
phylo.tree <- read.tree("Phylogenetic subtree.tre")

# check that the tree still works as expected!
phylo.tree$tip.label <- sub(" ", "_", phylo.tree$tip.label) # re-loading the tree removes the under-score between genus and species names
# is it rooted?
is.rooted(phylo.tree)
# is it binary?
is.binary(phylo.tree)
# is it ultrametric
is.ultrametric(phylo.tree)
phylo.tree; plot(phylo.tree)
any(duplicated(phylo.tree$node.label))

# and check that the tree and datafile still match
name.check(phylo.tree, Raw_data)

# STEP 1: ---------------------------------------------------------------------------------------------
# Check and sort outliers, distributions and values that don't make sense, and subset the data based on correct matrices

# are the columns in the correct formats.
lapply(Raw_data,class) # this returns a list of the format for each column
# all the data columns that I am interested in are in the correct format.
summary(Raw_data) # data all looks good - no NaNs or Infs

# Remove any non-Ergodic, imprimitive and reducible matrices that have slipped through the net
Raw_data <- Raw_data[which(Raw_data$Ergodic == TRUE),]
Raw_data <- Raw_data[which(Raw_data$Primitive == TRUE),]
Raw_data <- Raw_data[which(Raw_data$Irreducible == TRUE),]
# it appears none slipped through but its good to check

# Compadre and Comadre contain both pooled, mean and individual matrices. These capture different types of data and so its inappropriate to include all.
# For my analysis I just want mean (with no corresponding individual matrices) and individual matrices.
Raw_data <- Raw_data[which(Raw_data$MatrixComposite %in% c("Mean", "Individual")),]
# again this subset was already carried out but its good to check no populations squeezed through.

# Now make sure that impossible and outlier values are excluded
# first start by checking the lambda values - if I take out the lambda values outside the range of 0-2 then I know that any matrices that definitely haven't worked have been removed.
# check the range of Lambda values.
hist(Raw_data$Lambda) #it looks like there are some very big lambda values - though there don't appear to be any negative entries
range(Raw_data$Lambda, na.rm = T) #that is confirmed!
# subset the dataset to remove all rows in which lambda is greater than 2 - as these matrices haven't worked for this analysis.
Raw_data <- Raw_data[which(Raw_data$Lambda <= 2),] #this removes 164 entries ( :( )
# recheck the data
hist(Raw_data$Lambda) #looks much better and its already normally distributed, and between the range of 0 and 2!! Sweet!

# I also need to remove (and identify!) the few populations that exhibit clonality.
# first identify the populations with clonality.
clonality <- Raw_data[which(Raw_data$SCloDRs != 0 & !is.na(Raw_data$SCloDRs) |
                            Raw_data$SCloPOs != 0 & !is.na(Raw_data$SCloPOs) |
                            Raw_data$SCloRs != 0 & !is.na(Raw_data$SCloRs) |
                            Raw_data$SCloAs != 0 & !is.na(Raw_data$SCloAs) |
                            Raw_data$SCloMAs != 0 & !is.na(Raw_data$SCloMAs) |
                            Raw_data$SCloMAtts != 0 & !is.na(Raw_data$SCloMAtts)),] # this will remove 142 populations, again not too bad. 
# what are the species names of the populations removed due to clonality.
clonality$SpeciesAccepted

# remove the populations exhibiting clonality - ie the remaining populations will only have zeros or NAs in the clonality sensitivity columns.
Raw_data <- Raw_data[which(Raw_data$SCloDRs == 0 | is.na(Raw_data$SCloDRs) &
                 Raw_data$SCloPOs == 0 | is.na(Raw_data$SCloPOs) &
                 Raw_data$SCloRs == 0 | is.na(Raw_data$SCloRs) &
                 Raw_data$SCloAs == 0 | is.na(Raw_data$SCloAs) &
                 Raw_data$SCloMAs == 0 | is.na(Raw_data$SCloMAs) &
                 Raw_data$SCloMAtts == 0 | is.na(Raw_data$SCloMAtts)),]

# Now that I have removed populations I need to ensure the species in the tree and the data set match.
check_species <- name.check(phylo.tree, Raw_data) #this returns a list of species names that are the your tree but not your data and visa versa
# now I need to make sure both tree and data match
phylo.tree <- drop.tip(phylo.tree, check_species$tree_not_data) #prune your tree so the two species lists match
# and subsample your data based on the pruned tree
Raw_data <- Raw_data[phylo.tree$tip.label,]
# just check that has worked!
name.check(phylo.tree, Raw_data) #this all now matches

# So following all this subsetting I am working with 2242 separate populations from
dim(Raw_data[unique(Raw_data$SpeciesAccepted),])[1] # 369 different species

# Now I can remove any outliers from the other parameters - note this doesn't remove the population entry.
# Function taken from Robs code.
removeOutliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.025, .975), na.rm = na.rm) #this determines the 25% and 75% limits of the data set based on the supplied vector (ie. where 95% of the data lies) 
  H <- 1.5 * IQR(x, na.rm = na.rm) # determines the interquartile range (variance of the dataset)
  y <- x
  y[x < (qnt[1] - H)] <- NA #any values below or above these ranges are then excluded.
  y[x > (qnt[2] + H)] <- NA
  y
}

# list the different numeric quantities that I am interested in.
numVars <- c("Rho1","SSurvDR","SGrowDR","SShriDR","SRepDR","SCloDR","ESurvDR","EGrowDR","EShriDR","ERepDR","ECloDR",
             "Rho2" ,"SSurvDRs","SGrowDRs","SShriDRs","SRepDRs","SCloDRs","ESurvDRs","EGrowDRs","EShriDRs","ERepDRs","ECloDRs",
             "Pi1","SSurvPO","SGrowPO","SShriPO","SRepPO","SCloPO","ESurvPO","EGrowPO","EShriPO","ERepPO","ECloPO",
             "Pi2","SSurvPOs","SGrowPOs","SShriPOs","SRepPOs","SCloPOs","ESurvPOs","EGrowPOs","EShriPOs","ERepPOs","ECloPOs",
             "Reactivity1","SSurvR","SGrowR","SShriR","SRepR","SCloR","ESurvR","EGrowR","EShriR","ERepR","ECloR",
             "Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs","SCloRs","ESurvRs","EGrowRs","EShriRs","ERepRs","ECloRs",
             "Attenuation1","SSurvA","SGrowA","SShriA","SRepA","SCloA","ESurvA","EGrowA","EShriA","ERepA","ECloA",
             "Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs","SCloAs","ESurvAs","EGrowAs","EShriAs","ERepAs","ECloAs",
             "MaxAmplification1" ,"SSurvMA","SGrowMA","SShriMA","SRepMA","SCloMA","ESurvMA","EGrowMA","EShriMA","ERepMA","ECloMA",
             "MaxAmplification2" ,"SSurvMAs","SGrowMAs","SShriMAs","SRepMAs","SCloMAs","ESurvMAs","EGrowMAs","EShriMAs","ERepMAs","ECloMAs",
             "MaxAttenuation1","SSurvMAtt","SGrowMAtt","SShriMAtt","SRepMAtt","SCloMAtt","ESurvMAtt","EGrowMAtt","EShriMAtt","ERepMAtt","ECloMAtt",
             "MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts","SCloMAtts","ESurvMAtts","EGrowMAtts","EShriMAtts","ERepMAtts","ECloMAtts",
             "Lambda","Inertia.upper","Inertia.lower", 
             "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full",
             "Temp.Frequency.inter", "Temp.Autocorrelation.inter", "Temp.range.inter", "Prec.Frequency.inter", "Prec.Autocorrelation.inter",
             "Temp.Frequency.intra", "Temp.Autocorrelation.intra", "Temp.range.intra", "Prec.Frequency.intra", "Prec.Autocorrelation.intra")
             
# loop through and remove outliers
for (i in numVars){
  Raw_data[,i]=removeOutliers(Raw_data[,i], na.rm = T)
}
# this changes outliers to NA but doesn't remove the population (this will need considering when constructing the PCAs)

# observe distributions and transform if necessary

# Damping ratio ----------------------------
par(mfrow = c(1,2))
hist(Raw_data[,"Rho2"], main = "Raw data", xlab = "", ylab = "")
# a BoxCox transformation will be used to achieve approximate normality.
# first check the range of the data (the type of transformation depends on the occurrence of zeros)
  range(Raw_data[,"Rho2"], na.rm = T) # no zeros so I single parameter transformation can be used
  BoxCoxTrans(Raw_data[,"Rho2"], na.rm = T) #lambda = -1.1
# this lambda value is equivalent to 1/x^0.8
  hist(1/(Raw_data[,"Rho2"]^1.1), xlab = "", ylab = "")
  Raw_data[,"Rho2"] <- 1/(Raw_data[,"Rho2"]^1.1)
# Sensitivities
par(mfrow = c(2,2))
hist(Raw_data[,"SSurvDRs"], main = "Raw data", xlab = "", ylab = "")
hist(Raw_data[,"SGrowDRs"], main = "Raw data", xlab = "", ylab = "")
hist(Raw_data[,"SShriDRs"], main = "Raw data", xlab = "", ylab = "")          
hist(Raw_data[,"SRepDRs"], main = "Raw data", xlab = "", ylab = "")          

# Period of Oscillation ----------------------
par(mfrow = c(1,2))
hist(Raw_data[,"Pi2"], main = "Raw data", xlab = "", ylab = "")
# apply the boxcox transformation again.
  range(Raw_data[,"Pi2"], na.rm = T) #still no zeros so only one parameter required
  BoxCoxTrans(Raw_data[,"Pi2"], na.rm = T) #lambda = -0.4
  hist(1/(Raw_data[,"Pi2"]^0.4), xlab = "", ylab = "")
  Raw_data[,"Pi2"] <- 1/(Raw_data[,"Pi2"]^0.4)
# Sensitivities
par(mfrow = c(2,2))
hist(Raw_data[,"SSurvPOs"], main = "Raw data", xlab = "", ylab = "") 
hist(Raw_data[,"SGrowPOs"], main = "Raw data", xlab = "", ylab = "")          
hist(Raw_data[,"SShriPOs"], main = "Raw data", xlab = "", ylab = "")          
hist(Raw_data[,"SRepPOs"], main = "Raw data", xlab = "", ylab = "")          

# Reactivity --------------------------------
par(mfrow = c(1,2))
hist(Raw_data[,"Reactivity2"], main = "Raw data", xlab = "", ylab = "")
# Box Cox
  range(Raw_data[,"Reactivity2"], na.rm = T) # no zeros.
  BoxCoxTrans(Raw_data[,"Reactivity2"], na.rm = T) #lambda = -0.6
  hist(1/(Raw_data[,"Reactivity2"]^0.6), xlab = "", ylab = "")
  Raw_data[,"Reactivity2"] <- 1/(Raw_data[,"Reactivity2"]^0.6)
# Sensitivities
par(mfrow = c(2,2))
hist(Raw_data[,"SSurvRs"], main = "Raw data", xlab = "", ylab = "")
# this data is heavily negatively skewed and will require a reflection transformation to centralise the mean.
  hist(1/(max(Raw_data[,"SSurvRs"], na.rm = T)-Raw_data[,"SSurvRs"]))
  hist(log(max(Raw_data[,"SSurvRs"], na.rm = T) - Raw_data[,"SSurvRs"]))
  hist(1/Raw_data[,"SSurvRs"])
  Raw_data[,"SSurvRs"] <- log(max(Raw_data[,"SSurvRs"], na.rm = T) - Raw_data[,"SSurvRs"])
hist(Raw_data[,"SGrowRs"], main = "Raw data", xlab = "", ylab = "") 
  hist(1/Raw_data[,"SGrowRs"], main = "inverse", xlab = "", ylab = "") #this centres the peak
  hist(log(max(Raw_data[,"SGrowRs"], na.rm = T) - Raw_data[,"SGrowRs"]))
  hist(1/(max(Raw_data[,"SGrowRs"], na.rm = T) - Raw_data[,"SGrowRs"]))
  Raw_data[,"SGrowRs"] <- log(max(Raw_data[,"SGrowRs"], na.rm = T) - Raw_data[,"SGrowRs"])
hist(Raw_data[,"SShriRs"], main = "Raw data", xlab = "", ylab = "")    
  hist(log(Raw_data[, "SShriRs"] + abs(min(Raw_data[, "SShriRs"], na.rm = T))))  
  Raw_data[, "SShriRs"] <- log(Raw_data[, "SShriRs"] + abs(min(Raw_data[, "SShriRs"], na.rm = T)))
hist(Raw_data[,"SRepRs"], main = "Raw data", xlab = "", ylab = "")
  hist(log(max(Raw_data[,"SRepRs"], na.rm = T) - Raw_data[,"SRepRs"]))
  Raw_data[, "SRepRs"] <- log(max(Raw_data[,"SRepRs"], na.rm = T) - Raw_data[,"SRepRs"])

# first-timestep Attenuation -------------------         
par(mfrow = c(1,2))
hist(Raw_data[,"Attenuation2"], main = "Raw data", xlab = "", ylab = "")
# Box Cox
  range(Raw_data[,"Attenuation2"], na.rm = T)
  BoxCoxTrans(Raw_data[,"Attenuation2"], na.rm = T) # lambda = 0.7
  hist(Raw_data[,"Attenuation2"]^0.7)
  Raw_data[,"Attenuation2"] <- Raw_data[,"Attenuation2"]^0.7
# Sensitivities           
par(mfrow = c(2,2))
hist(Raw_data[,"SSurvAs"], main = "Raw data", xlab = "", ylab = "")   
hist(Raw_data[,"SGrowAs"], main = "Raw data", xlab = "", ylab = "")
  hist(log(max(Raw_data[,"SGrowAs"], na.rm = T) - Raw_data[,"SGrowAs"]))
  Raw_data[,"SGrowAs"] <- log(max(Raw_data[,"SGrowAs"], na.rm = T) - Raw_data[,"SGrowAs"])
hist(Raw_data[,"SShriAs"], main = "Raw data", xlab = "", ylab = "")
hist(Raw_data[,"SRepAs"], main = "Raw data", xlab = "", ylab = "")           
  hist(log(max(Raw_data[,"SRepAs"], na.rm = T) - Raw_data[,"SRepAs"]))
  hist(1/(max(Raw_data[,"SRepAs"], na.rm = T) - Raw_data[,"SRepAs"]))
  Raw_data[,"SRepAs"] <- log(max(Raw_data[,"SRepAs"], na.rm = T) - Raw_data[,"SRepAs"])

# Maximum Amplification ------------------------         
par(mfrow = c(1,2))
hist(Raw_data[,"MaxAmplification2"], main = "Raw data", xlab = "", ylab = "")
# Box Cox
  range(Raw_data[,"MaxAmplification2"], na.rm = T)
  BoxCoxTrans(Raw_data[,"MaxAmplification2"], na.rm = T) #lambda = -0.5
  hist(1/(Raw_data[,"MaxAmplification2"]^0.5))
  Raw_data[,"MaxAmplification2"] <- 1/(Raw_data[,"MaxAmplification2"]^0.5)
# Sensitivities
par(mfrow = c(2,2))  
hist(Raw_data[,"SSurvMAs"], main = "Raw data", xlab = "", ylab = "")   
hist(Raw_data[,"SGrowMAs"], main = "Raw data", xlab = "", ylab = "")
  hist(1/Raw_data[,"SGrowMAs"], main = "Raw data", xlab = "", ylab = "") #centralises to a point.
  hist(log(max(Raw_data[,"SGrowMAs"], na.rm = T) - Raw_data[,"SGrowMAs"]))
  hist(1/(max(Raw_data[,"SGrowMAs"], na.rm = T) - Raw_data[,"SGrowMAs"]))
  hist((max(Raw_data[,"SGrowMAs"], na.rm = T) - Raw_data[,"SGrowMAs"]))
  Raw_data[, "SGrowMAs"] <- log(max(Raw_data[,"SGrowMAs"], na.rm = T) - Raw_data[,"SGrowMAs"])
hist(Raw_data[,"SShriMAs"], main = "Raw data", xlab = "", ylab = "")
hist(Raw_data[,"SRepMAs"], main = "Raw data", xlab = "", ylab = "")
  hist(1/Raw_data[,"SRepMAs"], main = "Raw data", xlab = "", ylab = "")
  hist(log(max(Raw_data[,"SRepMAs"], na.rm = T) - Raw_data[,"SRepMAs"]))
  hist(1/(max(Raw_data[,"SRepMAs"], na.rm = T) - Raw_data[,"SRepMAs"]))
  Raw_data[,"SRepMAs"] <- log(max(Raw_data[,"SRepMAs"], na.rm = T) - Raw_data[,"SRepMAs"])

# Maximal Attenuation ----------------------------
par(mfrow = c(1,2))
hist(Raw_data[,"MaxAttenuation2"], main = "Raw data", xlab = "", ylab = "")
# Box Cox
  range(Raw_data[,"MaxAttenuation2"], na.rm = T) #no zeros
  BoxCoxTrans(Raw_data[,"MaxAttenuation2"], na.rm = T) #lambda = 0.4
  hist(Raw_data[,"MaxAttenuation2"]^0.4)
  Raw_data[,"MaxAttenuation2"] <- Raw_data[,"MaxAttenuation2"]^0.4
#Sensitivities
par(mfrow = c(2,2))
hist(Raw_data[,"SSurvMAtts"], main = "Raw data", xlab = "", ylab = "")   
hist(Raw_data[,"SGrowMAtts"], main = "Raw data", xlab = "", ylab = "")
hist(Raw_data[,"SShriMAtts"], main = "Raw data", xlab = "", ylab = "")
  hist(log(Raw_data[,"SShriMAtts"]+1))
  hist(log(Raw_data[, "SShriMAtts"] + abs(min(Raw_data[, "SShriMAtts"], na.rm = T))))  
  Raw_data[,"SShriMAtts"] <- log(Raw_data[, "SShriMAtts"] + abs(min(Raw_data[, "SShriMAtts"], na.rm = T)))
hist(Raw_data[,"SRepMAtts"], main = "Raw data", xlab = "", ylab = "")           

# here BoxCox transformations were used to determine the lambda values nessecary to acheieve approximate normality in the transient measure distributions
# to achieve approximate normality sensitivity distributions were also log transformed where nessecary.
# in the case that a distribution was heavily negatively skewed then reflection was also used to centralise the mean (max(x) - x).

# Now check the distributions of the relevant environmental data and transform if nessecary
# Temperature
par(mfrow = c(1,1))
hist(Raw_data$Temp.Frequency.full)  # Normal
hist(Raw_data$Temp.Autocorrelation.full) # Normal
hist(Raw_data$Temp.range.full) # Normal
# Precipitation
hist(Raw_data$Prec.Frequency.full) # Normal
hist(Raw_data$Prec.Autocorrelation.full) # Normalish

# remove any NaNs that have been created by the transformations
summary(Raw_data) 
Raw_data[sapply(Raw_data,is.infinite)] <- NA
Raw_data[sapply(Raw_data,is.nan)] <- NA

# STEP 2: Evaluate the colinearity between environmental variables --------------------------------------------------------------
# In the environmental analysis I am investigating the influence of Precipitation and Temperature regimes on resilience.
# However there is no precipitation data for marine populations so I am going to exclude marine populations from the environmental analysis.
# How many populations will this exclude?
marine_subset <- Raw_data[which(Raw_data$Realm == "Marine"),] # this will remove 39 populations.
table(marine_subset$OrganismType) # 10 Algae populations, 29 Animal populations

# using the multicol function test evaluate the variance inflation factor of each of the three extracted measures of environmental variance.
multicol(Raw_data[which(Raw_data$Realm != "Marine"), c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")])
# There is no multicolinearity. VIF values of 1 indicate no corelation. Vlaues of between 1 and 5 indicate moderate colinearity (but not enough to be concerned).
# Values of over 5 are a problem.
# What about multicolinarity removing autocorrelation
multicol(Raw_data[which(Raw_data$Realm != "Marine"), c("Temp.Frequency.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")])
# prehaps excluding Thermal Autocorrelation is the most appropriate?

# STEP 3: Run the environmental PCA ---------------------------------------------------------------------------------------------
# Now I can condense the variables using a normal PCA (one PCA for full, inter and intra)
enviro.full <- Raw_data[which(Raw_data$Realm != "Marine"), c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")]    
#enviro.inter <- data[,c("Frequency.inter.store","Frequency.inter","Autocorrelation.inter","range.inter")]   
#enviro.intra <- data[,c("Frequency.intra.store", "Frequency.intra","Autocorrelation.intra","range.intra")]   

# How much data is missing?
vis_miss(enviro.full)
#vis_miss(enviro.inter)
#vis_miss(enviro.intra) #not a lot is the answer!

# seperate out the complete cases for the perfect PCA
enviro.full <- enviro.full[complete.cases(enviro.full),]
#enviro.inter <- enviro.inter[complete.cases(enviro.inter),]
#enviro.intra <- enviro.intra[complete.cases(enviro.intra),]

# Run the PCA
enviro_PCA1 <- prcomp(enviro.full[,c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")], scale = T, center = T)
#enviro_PCA2 <- prcomp(enviro.inter[,c("Frequency.inter","Autocorrelation.inter","range.inter")], scale = T, center = T)
#enviro_PCA3 <- prcomp(enviro.intra[,c("Frequency.intra","Autocorrelation.intra","range.intra")], scale = T, center = T)

# STEP 4: Plot and assess the PCAs ---------------------------------------------------------------------
# full time series
enviro_pc1 <- data.frame(PC1 = enviro_PCA1$x[,1], PC2 = enviro_PCA1$x[,2],
                         PC3 = enviro_PCA1$x[,3], PC4 = enviro_PCA1$x[,4],
                         PC5 = enviro_PCA1$x[,5], Freq = enviro.full[,1]) #store data for plotting
PCA1_loadings <- data.frame(enviro_PCA1$rotation)
hist(enviro_pc1$Freq) # this is still the untransformed data
range(enviro_pc1$Freq)

# define the colour scheme
pal = colorRampPalette(c("Blue", "White", "Red", "Brown"))

# plot the PCA
ggplot(data = enviro_pc1, aes(x = PC1, y = PC2)) + #base plot
  labs(y = paste0("PC2 (", round(summary(enviro_PCA1)$importance[2,2]*100, 1),"%)"), # add axis labels
       x = paste0("PC1 (", round(summary(enviro_PCA1)$importance[2,1]*100, 1),"%)")) +
  geom_point(size = 2, aes(col = Freq), data = enviro_pc1) + #add data points
  scale_color_gradientn(colors = pal(dim(enviro_pc1)[1]),
                        values = scales::rescale(c(2,0,-1,-2,-3), to = c(0,1)), limits = c(-3,3), #this sets the value limits of the colour scale - the rescaled values reflect the expected frequency values
                        # associated with the different coloured environments
                        guide = guide_colorbar(ticks = F, title = NULL, title.position = "top", # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = T,
                                               barwidth = 15, barheight = 1)) +
  # Add loading arrows (no transformations need to be accounted for)
  geom_segment(data = PCA1_loadings[1:5,], aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), arrow = arrow(length = unit(1, "picas")),
               color = "black", size = 1) +
  #annotate("text", x = (PCA1_loadings[1:5,]$PC1*4), y = (PCA1_loadings[1:5,]$PC2*4),
  #         label = c("Thermal Variance", "Thermal Autocorrelation", "Thermal Range", "Rainfall Variance", "Rainfall Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom", legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.title = element_text(size = 10))

# Extract PCA coordinates --------------------------------------------------------------
# good to store these in case they are needed later in the main analysis

# extract the principal components
# for this I am only keeping PC1 and PC2, with 1 describing the frequency and autocorrelation of environmental changed and PC2 describing the magnitude of change
Raw_data$e.full.PC1 <- NA
Raw_data$e.full.PC2 <- NA
Raw_data$e.full.PC3 <- NA
Raw_data$e.full.PC4 <- NA
Raw_data$e.full.PC5 <- NA#this won't be used but I'll extract it still

for (x in 1:dim(Raw_data)[1]) {
  for (k in 1:dim(enviro_pc1)[1]) {
    if (rownames(enviro_pc1)[k] == rownames(Raw_data)[x]) 
    {
    Raw_data$e.full.PC1[x] = enviro_pc1$PC1[k]
    Raw_data$e.full.PC2[x] = enviro_pc1$PC2[k]
    Raw_data$e.full.PC3[x] = enviro_pc1$PC3[k]
    Raw_data$e.full.PC4[x] = enviro_pc1$PC4[k]
    Raw_data$e.full.PC5[x] = enviro_pc1$PC5[k]
    }
  }
}

# and export file.
write.csv(Raw_data, file = "Formatted RawData.csv") # Checkpoint

# ************************************************************************** End of code *******************************