# This script is for carrying out a phylogenetically corrected partial least squares analysis evaluating the impact of environmental stochasticity on the transient dynamics of natural populations.
# Author: James Cant
# Date last modified: June 2020

# clear workspace 
rm(list=ls(all=TRUE))

# load required packages
library(dplyr)
library(naniar) 
library(phytools)
library(ggplot2)
library(ape)
library(RColorBrewer)
library(geiger)
library(Hmisc) 
library(plotrix)
library(caret)
library(colorspace)
library(Rphylopars)
library(compiler)
# Be sure to use the phylo.pls function
source("FILE_PATH/Phylogenetic Partial Least squares function.R")

# set working directory
setwd("FILE_PATH")

# load in data
data <- read.csv("Formatted CSV FILE.csv", stringsAsFactors = FALSE, row.names = 1)

# load in the tree
phylo.tree <- read.tree("Phylogenetic subtree.tre")

# check that the tree still works as expected!
phylo.tree$tip.label <- sub(" ", "_", phylo.tree$tip.label) # match up species names in the phylo tree with the rowname labels of the main dataframe
# is it rooted?
is.rooted(phylo.tree)
# is it binary?
is.binary(phylo.tree)
# is it ultrametric
is.ultrametric(phylo.tree)
phylo.tree; plot(phylo.tree)
any(duplicated(phylo.tree$node.label))

# and check that the tree and datafile still match
name.check(phylo.tree, data)
# a few names need dropping to reflect the subsetting of the main data set during the formatting of the data.
check.tree <- name.check(phylo.tree, data)
phylo.tree <- drop.tip(phylo.tree, check.tree$tree_not_data)
name.check(phylo.tree, data)
# that all works and matches!

# STEP 1: Run pls analyses -------------------------------------------------------------------------
# Now I have a formatted data set containing all the nessecary variables I can run a pls analysis evaluating the correlation of each of the 
# transient measures with both vital rate sensitivities and measures of environmental variables.

# first subset out the relevant data for this analysis.
key_data <- data[,c("Rho2", "SSurvDRs","SGrowDRs","SShriDRs","SRepDRs",
                    "Pi2","SSurvPOs","SGrowPOs","SShriPOs","SRepPOs",
                    "Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs",
                    "Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs",
                    "MaxAmplification2", "SSurvMAs","SGrowMAs","SShriMAs","SRepMAs",
                    "MaxAttenuation2", "SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts",
                    "Temp.Frequency.full","Temp.Autocorrelation.full", "Temp.range.full",
                    "Prec.Frequency.full", "Prec.Autocorrelation.full")]
                    #"e.full.PC1", "e.full.PC2", "e.full.PC3")]   # incase I want to use the environmental PCA componenents as my measure of environmental variability.
# clonality is being ignored as it is not a universal trait.

# For completeness reassign the seal populations as marine
data[data$SpeciesAccepted %in% c("Zalophus californianus", "Leptonychotes weddellii"), "Realm"] <- "Marine"

# Define a list of marine populations for use in subsetting the data set during the pls analysis part 2.
marine.list <- rownames(data[data$Realm == "Marine",]) # this is 39 populations as expected.

# view the missingness distribution
vis_miss(key_data)
# Most of the missing data relates the Period of Oscilations.
# My analysis will focus on only complete data sets but this adjustment will be done seperately for each transient measure data set (And is built into the pls function)
# However alongside this I am also going to evaluate the phylogenetic signal of the different variables and phylogenetically impute miising data to provide a supplementary analysis with the full dataset.



####### --------------------------------------- Analyses --------------------------------------



##### Start by sorting the data needed for each demographic variable (Imputation takes a long-time so its best to leave this running initially)
#### Damping ratio ----------------------------------------------------------------
damping_ratio_data <- as.matrix(key_data[, c("Rho2", "SSurvDRs","SGrowDRs","SShriDRs","SRepDRs",       # Transient measure and sensitivities
                                   "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full",
                                   "Prec.Frequency.full", "Prec.Autocorrelation.full")])  # Environmental data
names(damping_ratio_data) <- c("Rho2", "SSurvDRs","SGrowDRs","SShriDRs","SRepDRs", "Temp.Frequency.full",
                               "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full",
                               "Prec.Autocorrelation.full")

# check this data matches with the phylogenetic tree
name.check(phylo.tree, damping_ratio_data)

# Now determine the phylogenetic signal across the demographic variables relating to the damping ratio using Pagels lambda
# variables with a signal value greater than 0.65 are suitable for imputation
DR.phylo.signal <- apply(damping_ratio_data[,c("Rho2", "SSurvDRs", "SGrowDRs", "SShriDRs", "SRepDRs")],2, phylosig, tree = phylo.tree, method = "lambda", test = TRUE)

# Generate a phylogenetically imputed dataset to compare with the complete dataset. 
DR.impute <- data.frame(rownames(damping_ratio_data), damping_ratio_data[,1:5]); colnames(DR.impute)[1] <- "species"
imputed_list.DR <- phylopars(trait_data = DR.impute[,1:5], tree = phylo.tree, model = "BM")
imputed_list.DR_rep <- phylopars(trait_data = DR.impute[,c(1,6)], tree = phylo.tree, model = "BM") # for some weird reason when this column is included in the main imputation the system breaksdown
# stitch the relevant data back together
imputed_list.DR <- cbind(imputed_list.DR[["anc_recon"]][1:dim(data)[1],], imputed_list.DR_rep[["anc_recon"]][1:dim(data)[1],]); colnames(imputed_list.DR)[5] <- "SSRepDRs"

#### Period of Oscillation ----------------------------------------------------------------
PO_data <- as.matrix(key_data[, c("Pi2", "SSurvPOs","SGrowPOs","SShriPOs","SRepPOs",       # Transient measure and sensitivities
                                  "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full",
                                  "Prec.Frequency.full", "Prec.Autocorrelation.full")])  # Environmental data
names(PO_data) <- c("Pi2", "SSurvPOs","SGrowPOs","SShriPOs","SRepPOs", "Temp.Frequency.full",
                               "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full",
                               "Prec.Autocorrelation.full")

# check this data matches with the phylogenetic tree
name.check(phylo.tree, PO_data)

# Now determine the phylogenetic signal across the demographic variables relating to the damping ratio using Pagels lambda
# variables with a signal value greater than 0.65 are suitable for imputation
PO.phylo.signal <- apply(PO_data[, c("Pi2", "SSurvPOs","SGrowPOs","SShriPOs","SRepPOs")], 2, phylosig, tree = phylo.tree, method = "lambda", test = TRUE)

# Generate a phylogenetically imputed dataset to compare with the complete dataset.
# For some reason the phylopars function doesn't work here so I tried to call the phylo.impute function from phytools.
# However after 13 hours this showed no sign of completing - however due to the low phylogenetic signal in this set of variable I can actually only phylogenetically impute the Period of Oscilation value.
PO.impute <- data.frame(rownames(PO_data), PO_data[,1]); colnames(PO.impute)[1] <- "species"
imputed_list.PO <- phylopars(trait_data = PO.impute, tree = phylo.tree, model = "BM")


#### Reactivity  ------------------------------------------------------------------------
Reac_data <- as.matrix(key_data[, c("Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs",      # Transient measure and sensitivities
                                  "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full",
                                  "Prec.Frequency.full", "Prec.Autocorrelation.full")])  # Environmental data
names(Reac_data) <- c("Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs", "Temp.Frequency.full",
                    "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full",
                    "Prec.Autocorrelation.full")

# check this data matches with the phylogenetic tree
name.check(phylo.tree, Reac_data)

# Now determine the phylogenetic signal across the demographic variables relating to the damping ratio using Pagels lambda
# variables with a signal value greater than 0.65 are suitable for imputation
Reac.phylo.signal <- apply(Reac_data[, c("Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs")], 2, phylosig, tree = phylo.tree, method = "lambda", test = TRUE)

# Generate a phylogenetically imputed dataset to compare with the complete dataset. 
Reac.impute <- data.frame(rownames(Reac_data), Reac_data[,1:5]); colnames(Reac.impute)[1] <- "species"
imputed_list.Reac <- phylopars(trait_data = Reac.impute, tree = phylo.tree, model = "BM")


#### First-step Attenuation -----------------------------------------------------------
Att_data <- as.matrix(key_data[, c("Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs",     # Transient measure and sensitivities
                                    "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full",
                                    "Prec.Frequency.full", "Prec.Autocorrelation.full")])  # Environmental data
names(Att_data) <- c("Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs", "Temp.Frequency.full",
                      "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full",
                      "Prec.Autocorrelation.full")

# check this data matches with the phylogenetic tree
name.check(phylo.tree, Att_data)

# Now determine the phylogenetic signal across the demographic variables relating to the damping ratio using Pagels lambda
# variables with a signal value greater than 0.65 are suitable for imputation
Att.phylo.signal <- apply(Att_data[, c("Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs")], 2, phylosig, tree = phylo.tree, method = "lambda", test = TRUE)

# Generate a phylogenetically imputed dataset to compare with the complete dataset. 
Att.impute <- data.frame(rownames(Att_data), Att_data[,1:5]); colnames(Att.impute)[1] <- "species"
imputed_list.Att <- phylopars(trait_data = Att.impute, tree = phylo.tree, model = "BM")

      
#### Maximal Amplification ----------------------------------------------------------
MAmp_data <- as.matrix(key_data[, c("MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs",     # Transient measure and sensitivities
                                   "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full",
                                   "Prec.Frequency.full", "Prec.Autocorrelation.full")])  # Environmental data
names(MAmp_data) <- c("MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs", "Temp.Frequency.full",
                     "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full",
                     "Prec.Autocorrelation.full")

# check this data matches with the phylogenetic tree
name.check(phylo.tree, MAmp_data)

# Now determine the phylogenetic signal across the demographic variables relating to the damping ratio using Pagels lambda
# variables with a signal value greater than 0.65 are suitable for imputation
MAmp.phylo.signal <- apply(MAmp_data[, c("MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs")], 2, phylosig, tree = phylo.tree, method = "lambda", test = TRUE)

# Generate a phylogenetically imputed dataset to compare with the complete dataset. 
MAmp.impute <- data.frame(rownames(MAmp_data), MAmp_data[,1:5]); colnames(MAmp.impute)[1] <- "species"
imputed_list.MAmp <- phylopars(trait_data = MAmp.impute, tree = phylo.tree, model = "BM")


#### Maximal Attenuation ----------------------------------------------------------
MAtt_data <- as.matrix(key_data[, c("MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts",     # Transient measure and sensitivities
                                    "Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full",
                                    "Prec.Frequency.full", "Prec.Autocorrelation.full")])  # Environmental data
names(MAtt_data) <- c("MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts", "Temp.Frequency.full",
                      "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full",
                      "Prec.Autocorrelation.full")

# check this data matches with the phylogenetic tree
name.check(phylo.tree, MAtt_data)

# Now determine the phylogenetic signal across the demographic variables relating to the damping ratio using Pagels lambda
# variables with a signal value greater than 0.65 are suitable for imputation
MAtt.phylo.signal <- apply(MAtt_data[, c("MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts")], 2, phylosig, tree = phylo.tree, method = "lambda", test = TRUE)

# Generate a phylogenetically imputed dataset to compare with the complete dataset. 
MAtt.impute <- data.frame(rownames(MAtt_data), MAtt_data[,1:5]); colnames(MAtt.impute)[1] <- "species"
imputed_list.MAtt <- phylopars(trait_data = MAtt.impute, tree = phylo.tree, model = "BM")

##### store the outputs that take an age to run! -----------------------------------

# Firstly the phylogenetic signal of the variables
pagels.store <- data.frame(D.R = c(DR.phylo.signal[[1]][["lambda"]], DR.phylo.signal[[2]][["lambda"]], DR.phylo.signal[[3]][["lambda"]], DR.phylo.signal[[4]][["lambda"]], DR.phylo.signal[[5]][["lambda"]]),
                           P.O = c(PO.phylo.signal[[1]][["lambda"]], PO.phylo.signal[[2]][["lambda"]], PO.phylo.signal[[3]][["lambda"]], PO.phylo.signal[[4]][["lambda"]], PO.phylo.signal[[5]][["lambda"]]),
                           Rea = c(Reac.phylo.signal[[1]][["lambda"]], Reac.phylo.signal[[2]][["lambda"]], Reac.phylo.signal[[3]][["lambda"]], Reac.phylo.signal[[4]][["lambda"]], Reac.phylo.signal[[5]][["lambda"]]),
                           Att = c(Att.phylo.signal[[1]][["lambda"]], Att.phylo.signal[[2]][["lambda"]], Att.phylo.signal[[3]][["lambda"]], Att.phylo.signal[[4]][["lambda"]], Att.phylo.signal[[5]][["lambda"]]),
                           MAp = c(MAmp.phylo.signal[[1]][["lambda"]], MAmp.phylo.signal[[2]][["lambda"]], MAmp.phylo.signal[[3]][["lambda"]], MAmp.phylo.signal[[4]][["lambda"]], MAmp.phylo.signal[[5]][["lambda"]]),
                           MAt = c(MAtt.phylo.signal[[1]][["lambda"]], MAtt.phylo.signal[[2]][["lambda"]], MAtt.phylo.signal[[3]][["lambda"]], MAtt.phylo.signal[[4]][["lambda"]], MAtt.phylo.signal[[5]][["lambda"]])) 
# rename the rows for ease.
rownames(pagels.store) <- c("transient","survival","growth","shrinkage","reproduction")

# and the imputed datasets
imputed_data <- do.call(cbind, list(imputed_list.DR,
                                    #imputed_list.PO[["anc_recon"]][1:dim(data)[1],],
                                    imputed_list.Reac[["anc_recon"]][1:dim(data)[1],],
                                    imputed_list.Att[["anc_recon"]][1:dim(data)[1],],
                                    imputed_list.MAmp[["anc_recon"]][1:dim(data)[1],],
                                    imputed_list.MAtt[["anc_recon"]][1:dim(data)[1],]))

# and export to save
write.csv(pagels.store, file = "FILE_PATH/Phylogenetic signal.csv")
write.csv(imputed_data, file = "FILE_PATH/Imputed data.csv")


##### Run Pls analyses -------------------------------------------------------------

# load in imputed data
imputed_data <- read.csv("Imputed data.csv", row.names = 1)

# Just a quick formatting fix
imputed_list.DR <- imputed_data[,1:5]; #imputed_list.PO <- imputed_data[,6:10];
imputed_list.Reac <- imputed_data[6:10]; imputed_list.Att <- imputed_data[11:15]
imputed_list.MAmp <- imputed_data[,16:20]; imputed_list.MAtt <- imputed_data[21:25]

######### Damping ratio -------------------------------
# Combine together the correct variables, and extract the complete cases.
# Non-imputed dataset
perfect_DR_data <- damping_ratio_data[complete.cases(damping_ratio_data[,c("Rho2", "SSurvDRs","SGrowDRs","SShriDRs","SRepDRs")]),]  #remove any populations with NAs within the demographic variables.
# Imputed data
imputed_data.DR <- cbind(imputed_list.DR, damping_ratio_data[,6:10]) # this just confirms the dimensions of all the variables.
# Prune the tree to match the perfect data set (the imputed data will be able to use the full tree)
DR.check <- name.check(phylo.tree, perfect_DR_data)
DR.tree <- drop.tip(phylo.tree, DR.check$tree_not_data)

# define the variable sets for the analyses
# none-imputed dataset
DR.Y <- as.matrix(perfect_DR_data[, 1])
DR.X1 <- perfect_DR_data[,2:5]
DR.enviro1 <- perfect_DR_data[,6:10]
# imputed dataset
DR.Y2 <- as.matrix(imputed_data.DR[, 1]); rownames(DR.Y2) <- rownames(imputed_data.DR)
DR.X2 <- as.matrix(imputed_data.DR[,2:5])
DR.enviro2 <- as.matrix(imputed_data.DR[,6:10])
# ensure formating
rownames(DR.X1) <- rownames(DR.enviro1) <- rownames(DR.Y)
rownames(DR.X2) <- rownames(DR.enviro2) <- rownames(DR.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Damping.ratio.1 <- Phylo.PLS(DR.Y, DR.X1, DR.enviro1, tree = DR.tree, marine.list)
# Imputed data
Damping.ratio.2 <- Phylo.PLS(DR.Y2, DR.X2, DR.enviro2, tree = phylo.tree, marine.list)


######### Period of Oscillation -------------------------------
# Combine together the correct variables, and extract the complete cases.
# Non-imputed dataset
perfect_PO_data <- PO_data[complete.cases(PO_data[,c("Pi2", "SSurvPOs","SGrowPOs","SShriPOs","SRepPOs")]),]  #remove any populations with NAs within the demographic variables.
# Imputed data
imputed_data.PO <- as.data.frame(cbind(imputed_list.PO[["anc_recon"]][1:dim(data)[1],], PO_data[,2:10])) # this just confirms the dimensions of all the variables.
names(imputed_data.PO)[1] <- "Pi2"
# Prune the tree to match the perfect data set (the imputed data will be able to use the full tree)
PO.check <- name.check(phylo.tree, perfect_PO_data)
PO.tree <- drop.tip(phylo.tree, PO.check$tree_not_data)

# define the variable sets for the analyses
# none-imputed dataset
PO.Y <- as.matrix(perfect_PO_data[, 1])
PO.X1 <- perfect_PO_data[,2:5]
PO.enviro1 <- perfect_PO_data[,6:10]
# imputed dataset
PO.Y2 <- as.matrix(imputed_data.PO[, 1]); rownames(PO.Y2) <- rownames(imputed_data.PO)
PO.X2 <- as.matrix(imputed_data.PO[, 1]) # just a placeholder
PO.enviro2 <- as.matrix(imputed_data.PO[,6:10])
# ensure formating
rownames(PO.X1) <- rownames(PO.enviro1) <- rownames(PO.Y)
rownames(PO.X2) <- rownames(PO.enviro2) <- rownames(PO.Y2)

# Run the pls accounting for phylogeny.
# Normal data
PO.1 <- Phylo.PLS(PO.Y, PO.X1, PO.enviro1, tree = PO.tree, marine.list)
# Imputed data
PO.2 <- Phylo.PLS(PO.Y2, PO.X2, PO.enviro2, tree = phylo.tree, marine.list)


######### Reactivity ---------------------------------
# Combine together the correct variables, and extract the complete cases.
# Non-imputed dataset
perfect_reac_data <- Reac_data[complete.cases(Reac_data[,c("Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs")]),]  #remove any populations with NAs within the demographic variables.
# Imputed data
imputed_data.reac <- cbind(imputed_list.Reac, Reac_data[,6:10]) # this just confirms the dimensions of all the variables.
# Prune the tree to match the perfect data set (the imputed data will be able to use the full tree)
Reac.check <- name.check(phylo.tree, perfect_reac_data)
Reac.tree <- drop.tip(phylo.tree, Reac.check$tree_not_data)

# define the variable sets for the analyses
# none-imputed dataset
R.Y <- as.matrix(perfect_reac_data[, 1])
R.X1 <- perfect_reac_data[,2:5]
R.enviro1 <- perfect_reac_data[,6:10]
# imputed dataset
R.Y2 <- as.matrix(imputed_data.reac[, 1]); rownames(R.Y2) <- rownames(imputed_data.reac)
R.X2 <- as.matrix(imputed_data.reac[,2:5])
R.enviro2 <- as.matrix(imputed_data.reac[,6:10])
# ensure formating
rownames(R.X1) <- rownames(R.enviro1) <- rownames(R.Y)
rownames(R.X2) <- rownames(R.enviro2) <- rownames(R.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Reac.1 <- Phylo.PLS(R.Y, R.X1, R.enviro1, tree = Reac.tree, marine.list)
# Imputed data
Reac.2 <- Phylo.PLS(R.Y2, R.X2, R.enviro2, tree = phylo.tree, marine.list)


######### Attenuation ---------------------------------
# Combine together the correct variables, and extract the complete cases.
# Non-imputed dataset
perfect_Att_data <- Att_data[complete.cases(Att_data[,c("Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs")]),]  #remove any populations with NAs within the demographic variables.
# Imputed data
imputed_data.Att <- cbind(imputed_list.Att, Att_data[,6:10]) # this just confirms the dimensions of all the variables.
# Prune the tree to match the perfect data set (the imputed data will be able to use the full tree)
Att.check <- name.check(phylo.tree, perfect_Att_data)
Att.tree <- drop.tip(phylo.tree, Att.check$tree_not_data)

# define the variable sets for the analyses
# none-imputed dataset
Att.Y <- as.matrix(perfect_Att_data[, 1])
Att.X1 <- perfect_Att_data[,2:5]
Att.enviro1 <- perfect_Att_data[,6:10]
# imputed dataset
Att.Y2 <- as.matrix(imputed_data.Att[, 1]); rownames(Att.Y2) <- rownames(imputed_data.Att)
Att.X2 <- as.matrix(imputed_data.Att[,2:5])
Att.enviro2 <- as.matrix(imputed_data.Att[,6:10])
# ensure formating
rownames(Att.X1) <- rownames(Att.enviro1) <- rownames(Att.Y)
rownames(Att.X2) <- rownames(Att.enviro2) <- rownames(Att.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Att.1 <- Phylo.PLS(Att.Y, Att.X1, Att.enviro1, tree = Att.tree, marine.list)
# Imputed data
Att.2 <- Phylo.PLS(Att.Y2, Att.X2, Att.enviro2, tree = phylo.tree, marine.list)


######### Maximal Amplification ---------------------------------
# Combine together the correct variables, and extract the complete cases.
# Non-imputed dataset
perfect_MAmp_data <- MAmp_data[complete.cases(MAmp_data[,c("MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs")]),]  #remove any populations with NAs within the demographic variables.
# Imputed data
imputed_data.MAmp <- cbind(imputed_list.MAmp, MAmp_data[,6:10]) # this just confirms the dimensions of all the variables.
# Prune the tree to match the perfect data set (the imputed data will be able to use the full tree)
MAmp.check <- name.check(phylo.tree, perfect_MAmp_data)
MAmp.tree <- drop.tip(phylo.tree, MAmp.check$tree_not_data)

# define the variable sets for the analyses
# none-imputed dataset
MAmp.Y <- as.matrix(perfect_MAmp_data[, 1])
MAmp.X1 <- perfect_MAmp_data[,2:5]
MAmp.enviro1 <- perfect_MAmp_data[,6:10]
# imputed dataset
MAmp.Y2 <- as.matrix(imputed_data.MAmp[, 1]); rownames(MAmp.Y2) <- rownames(imputed_data.MAmp)
MAmp.X2 <- as.matrix(imputed_data.MAmp[,2:5])
MAmp.enviro2 <- as.matrix(imputed_data.MAmp[,6:10])
# ensure formating
rownames(MAmp.X1) <- rownames(MAmp.enviro1) <- rownames(MAmp.Y)
rownames(MAmp.X2) <- rownames(MAmp.enviro2) <- rownames(MAmp.Y2)

# Run the pls accounting for phylogeny.
# Normal data
MAmp.1 <- Phylo.PLS(MAmp.Y, MAmp.X1, MAmp.enviro1, tree = MAmp.tree, marine.list)
# Imputed data
MAmp.2 <- Phylo.PLS(MAmp.Y2, MAmp.X2, MAmp.enviro2, tree = phylo.tree, marine.list)


######### Maximal Attenuation ---------------------------------
# Combine together the correct variables, and extract the complete cases.
# Non-imputed dataset
perfect_MAtt_data <- MAtt_data[complete.cases(MAtt_data[,c("MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts")]),]  #remove any populations with NAs within the demographic variables.
# Imputed data
imputed_data.MAtt <- cbind(imputed_list.MAtt, MAtt_data[,6:10]) # this just confirms the dimensions of all the variables.
# Prune the tree to match the perfect data set (the imputed data will be able to use the full tree)
MAtt.check <- name.check(phylo.tree, perfect_MAtt_data)
MAtt.tree <- drop.tip(phylo.tree, MAtt.check$tree_not_data)

# define the variable sets for the analyses
# none-imputed dataset
MAtt.Y <- as.matrix(perfect_MAtt_data[, 1])
MAtt.X1 <- perfect_MAtt_data[,2:5]
MAtt.enviro1 <- perfect_MAtt_data[,6:10]
# imputed dataset
MAtt.Y2 <- as.matrix(imputed_data.MAtt[, 1]); rownames(MAtt.Y2) <- rownames(imputed_data.MAtt)
MAtt.X2 <- as.matrix(imputed_data.MAtt[,2:5])
MAtt.enviro2 <- as.matrix(imputed_data.MAtt[,6:10])
# ensure formating
rownames(MAtt.X1) <- rownames(MAtt.enviro1) <- rownames(MAtt.Y)
rownames(MAtt.X2) <- rownames(MAtt.enviro2) <- rownames(MAtt.Y2)

# Run the pls accounting for phylogeny.
# Normal data
MAtt.1 <- Phylo.PLS(MAtt.Y, MAtt.X1, MAtt.enviro1, tree = MAtt.tree, marine.list)
# Imputed data
MAtt.2 <- Phylo.PLS(MAtt.Y2, MAtt.X2, MAtt.enviro2, tree = phylo.tree, marine.list)


##### Extract and display key data ---------------------

### Hypothesis 1

### Damping ratio --------------------------
# Normal dataset
DR_1 <- data.frame(cbind(Damping.ratio.1[["Hypothesis 1 pls"]]$scores[,1],
                         Damping.ratio.1[["Hypothesis 1 pls"]]$scores[,2],
                         Damping.ratio.1[["Hypothesis 1 pls"]]$scores[,3],
                         Damping.ratio.1[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                         perfect_DR_data[,"Rho2"]))
colnames(DR_1) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"Rho") # puts the detail back in

# store PCA loadings.
DR_loadings_1 <- data.frame(cbind(Damping.ratio.1[["Hypothesis 1 pls"]]$loadings[,1],
                                  Damping.ratio.1[["Hypothesis 1 pls"]]$loadings[,2],
                                  Damping.ratio.1[["Hypothesis 1 pls"]]$loadings[,3],
                                  Damping.ratio.1[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(DR_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.DR <- explvar(Damping.ratio.1[["Hypothesis 1 pls"]]) # percentage variance explained
r2.DR <- R2(Damping.ratio.1[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
DR.coefficients = -coef(Damping.ratio.1[["Hypothesis 1 pls"]])
# standardise coefficients so that they sum to 1
sum.DR.coef = sum(sapply(DR.coefficients, abs))
DR.coefficients = DR.coefficients / sum.DR.coef
DR.coefficients = sort(DR.coefficients[, 1, 1])

# what range is required for the transient colour scale
range(DR_1$Rho, na.rm = T); median(DR_1$Rho, na.rm = T); mean(DR_1$Rho, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# create pls biplot
DR.perfect.plot <- ggplot(data = DR_1, aes(x = Comp1, y = Comp2)) + #base plot
                        labs(y = paste0("Component 2 (", round(p.var.DR[2], 1),"%)"), # add axis labels
                        x = paste0("Component 1 (", round(p.var.DR[1], 1),"%)")) +
                 geom_point(size = 3, aes(col = Rho), alpha = 0.7,  data = DR_1) + #add data points
                 scale_color_gradientn(colors = pal(dim(DR_1)[1]),
                                       values = scales::rescale(c(0.10,0.63,0.99999), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                                       guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                                              direction = "horizontal", reverse = T, label = F,
                                                              barwidth = 15, barheight = 1), na.value = "white") +
                 geom_segment(data = DR_loadings_1, aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
                              color = "black", size = 1) +
                 #annotate("text", x = (DR_loadings_1$Comp1*5), y = (DR_loadings_1$Comp2*5),
                        #label = c("Survival", "Growth", "Retrogression", "Reproduction")) +
                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
                       axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                       axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
                       legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot
level_order <- c('SSurvDRs', 'SGrowDRs', 'SShriDRs', 'SRepDRs')
DR.coefs.bar <- ggplot() + 
                geom_bar(aes(x = names(DR.coefficients), y = DR.coefficients), stat="identity", fill = "gold") +
                ylab("Proportional influence") +
                geom_hline(yintercept = 0) + 
                ylim(-1,1) +
                scale_x_discrete(limits = level_order, labels=c("", "", "", "")) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                      axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
                      axis.title.x = element_blank(), axis.title.y = element_blank())

### Not enough phylogenetic signal for using the imputed dataset



#### Period of Ocilation ---------------------------------------------------------------- 

# extract key data
# Perfect data
PO_1 <- data.frame(cbind(PO.1[["Hypothesis 1 pls"]]$scores[,1],
                         PO.1[["Hypothesis 1 pls"]]$scores[,2],
                         PO.1[["Hypothesis 1 pls"]]$scores[,3],
                         PO.1[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                         perfect_PO_data[,"Pi2"]))
colnames(PO_1) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"Pi") # puts the detail back in

# store PCA loadings.
PO_loadings_1 <- data.frame(cbind(PO.1[["Hypothesis 1 pls"]]$loadings[,1],
                                  PO.1[["Hypothesis 1 pls"]]$loadings[,2],
                                  PO.1[["Hypothesis 1 pls"]]$loadings[,3],
                                  PO.1[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(PO_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.PO <- explvar(PO.1[["Hypothesis 1 pls"]]) # percentage variance explained
r2.PO <- R2(PO.1[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
PO.coefficients = -coef(PO.1[["Hypothesis 1 pls"]])
# standardise coefficients so that they sum to 1
sum.PO.coef = sum(sapply(PO.coefficients, abs))
PO.coefficients = PO.coefficients / sum.PO.coef
PO.coefficients = sort(PO.coefficients[, 1 , 1])

# what range is required for the transient colour scale
range(PO_1$Pi, na.rm = T); median(PO_1$Pi, na.rm = T); mean(PO_1$Pi, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# create pls biplot
PO.perfect.plot <- ggplot(data = PO_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.PO[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.PO[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Pi), alpha = 0.7,  data = PO_1) + #add data points
  scale_color_gradientn(colors = pal(dim(PO_1)[1]),
                        values = scales::rescale(c(0.08,0.36,0.58), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = PO_loadings_1, aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (PO_loadings_1$Comp1*5), y = (PO_loadings_1$Comp2*5),
           #label = c("Survival", "Growth", "Retrogression", "Reproduction")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot
level_order <- c('SSurvPOs', 'SGrowPOs', 'SShriPOs', 'SRepPOs')
PO.coefs.bar <- ggplot() + 
  geom_bar(aes(x=names(PO.coefficients), y=PO.coefficients), stat="identity", fill = "gold") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())
# Again not enough phylogenetic signal here to use the imputed dataset.



#### Reactivity -------------------------------------------------------------------------

# extract key data
reac_1 <- data.frame(cbind(Reac.1[["Hypothesis 1 pls"]]$scores[,1],
                           Reac.1[["Hypothesis 1 pls"]]$scores[,2],
                           Reac.1[["Hypothesis 1 pls"]]$scores[,3],
                           Reac.1[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                           perfect_reac_data[,"Reactivity2"]))
colnames(reac_1) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"Reactivity") # puts the detail back in

# store PCA loadings.
reac_loadings_1 <- data.frame(cbind(Reac.1[["Hypothesis 1 pls"]]$loadings[,1],
                                    Reac.1[["Hypothesis 1 pls"]]$loadings[,2],
                                    Reac.1[["Hypothesis 1 pls"]]$loadings[,3],
                                    Reac.1[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(reac_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.R <- explvar(Reac.1[["Hypothesis 1 pls"]]) # percentage variance explained
r2.R <- R2(Reac.1[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where necessary)
reac.coefficients = c(coef(Reac.1[["Hypothesis 1 pls"]])[1:2], -coef(Reac.1[["Hypothesis 1 pls"]])[3], coef(Reac.1[["Hypothesis 1 pls"]])[4])
# standardise coefficients so that they sum to 1
sum.reac.coef = sum(sapply(reac.coefficients, abs))
reac.coefficients = reac.coefficients / sum.reac.coef
names(reac.coefficients) <- c("SSurvRs", "SGrowRs", "SShriRs", "SRepRs")
reac.coefficients = sort(reac.coefficients)

# what range is required for the transient colour scale
range(reac_1$Reactivity, na.rm = T); median(reac_1$Reactivity, na.rm = T); mean(reac_1$Reactivity, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# create pls biplot
R.perfect.plot <- ggplot(data = reac_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.R[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.R[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Reactivity), alpha = 0.7, data = reac_1) + #add data points
  scale_color_gradientn(colors = pal(dim(reac_1)[1]),
                        values = scales::rescale(c(0.01,0.63,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = reac_loadings_1[c(1:2,4),], aes(x = 0, y = 0, xend = (-Comp1*10), yend = (-Comp2*10)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  geom_segment(data = reac_loadings_1[c(3),], aes(x = 0, y = 0, xend = (Comp1*10), yend = (Comp2*10)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (-reac_loadings_1[c(1:2,4),]$Comp1*8), y = (-reac_loadings_1[c(1:2,4),]$Comp2*8),
          # label = c("Survival", "Growth", "Reproduction")) +
 # annotate("text", x = (reac_loadings_1[c(3),]$Comp1*8), y = (reac_loadings_1[c(3),]$Comp2*8),
          # label = c("Retrogression")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c('SSurvRs', 'SGrowRs', 'SShriRs', 'SRepRs')
R.coefs.bar <- ggplot() + 
  geom_bar(aes(x=names(reac.coefficients), y=reac.coefficients), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed dataset
reac_2 <- data.frame(cbind(Reac.2[["Hypothesis 1 pls"]]$scores[,1],
                           Reac.2[["Hypothesis 1 pls"]]$scores[,2],
                           Reac.2[["Hypothesis 1 pls"]]$scores[,3],
                           Reac.2[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                           imputed_data.reac[,"Reactivity2"]))
colnames(reac_2) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"Reactivity") # puts the detail back in

# store PCA loadings.
reac_loadings_2 <- data.frame(cbind(Reac.2[["Hypothesis 1 pls"]]$loadings[,1],
                                    Reac.2[["Hypothesis 1 pls"]]$loadings[,2],
                                    Reac.2[["Hypothesis 1 pls"]]$loadings[,3],
                                    Reac.2[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(reac_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.R2 <- explvar(Reac.2[["Hypothesis 1 pls"]]) # percentage variance explained
r2.R2 <- R2(Reac.2[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
reac.coefficients2 = c(coef(Reac.2[["Hypothesis 1 pls"]])[1:2], -coef(Reac.2[["Hypothesis 1 pls"]])[3], coef(Reac.2[["Hypothesis 1 pls"]])[4])
# standardise coefficients so that they sum to 1
sum.reac.coef2 = sum(sapply(reac.coefficients2, abs))
reac.coefficients2 = reac.coefficients2 / sum.reac.coef2
names(reac.coefficients2) <- c("SSurvRs", "SGrowRs", "SShriRs", "SRepRs")
reac.coefficients2 = sort(reac.coefficients2)

# what range is required for the transient colour scale
range(reac_2$Reactivity, na.rm = T); median(reac_2$Reactivity, na.rm = T); mean(reac_2$Reactivity, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# create pls biplot
R.imputed.plot <- ggplot(data = reac_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.R2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.R2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Reactivity), alpha = 0.7,  data = reac_2) + #add data points
  scale_color_gradientn(colors = pal(dim(reac_2)[1]),
                        values = scales::rescale(c(-0.73,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = reac_loadings_2[c(1:2,4),], aes(x = 0, y = 0, xend = (-Comp1*8), yend = (-Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  geom_segment(data = reac_loadings_2[c(3),], aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (-reac_loadings_2[c(1:2,4),]$Comp1*8), y = (-reac_loadings_2[c(1:2,4),]$Comp2*8),
         #  label = c("Survival", "Growth", "Reproduction")) +
  #annotate("text", x = (reac_loadings_2[c(3),]$Comp1*8), y = (reac_loadings_2[c(3),]$Comp2*8),
         #  label = c("Retrogression")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c('SSurvRs', 'SGrowRs', 'SShriRs', 'SRepRs')
R.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(reac.coefficients2), y=reac.coefficients2), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())



#### Attenuation ------------------------------------------------------------------------

# Perfect data
# extract key data
Att_1 <- data.frame(cbind(Att.1[["Hypothesis 1 pls"]]$scores[,1],
                          Att.1[["Hypothesis 1 pls"]]$scores[,2],
                          Att.1[["Hypothesis 1 pls"]]$scores[,3],
                          Att.1[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                          perfect_Att_data[,"Attenuation2"]))
colnames(Att_1) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"Attenuation") # puts the detail back in

# store PCA loadings.
Att_loadings_1 <- data.frame(cbind(Att.1[["Hypothesis 1 pls"]]$loadings[,1],
                                   Att.1[["Hypothesis 1 pls"]]$loadings[,2],
                                   Att.1[["Hypothesis 1 pls"]]$loadings[,3],
                                   Att.1[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(Att_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.Att <- explvar(Att.1[["Hypothesis 1 pls"]]) # percentage variance explained
r2.Att <- R2(Att.1[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
att.coefficients = c(coef(Att.1[["Hypothesis 1 pls"]])[1], -coef(Att.1[["Hypothesis 1 pls"]])[2], coef(Att.1[["Hypothesis 1 pls"]])[3], -coef(Att.1[["Hypothesis 1 pls"]])[4])
# standardise coefficients so that they sum to 1
sum.att.coef = sum(sapply(att.coefficients, abs))
att.coefficients = att.coefficients / sum.att.coef
names(att.coefficients) <- c("SSurvAs", "SGrowAs", "SShriAs", "SRepAs")
att.coefficients = sort(att.coefficients)

# what range is required for the transient colour scale
range(Att_1$Attenuation, na.rm = T); median(Att_1$Attenuation, na.rm = T); mean(Att_1$Attenuation, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# create pls biplot
Att.perfect.plot <- ggplot(data = Att_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Att[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Att[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Attenuation), alpha = 0.7, data = Att_1) + #add data points
  scale_color_gradientn(colors = pal(dim(Att_1)[1]),
                        values = scales::rescale(c(0.0,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Att_loadings_1[c(2,4),], aes(x = 0, y = 0, xend = (-Comp1*8), yend = (-Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  geom_segment(data = Att_loadings_1[c(1,3),], aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (-Att_loadings_1[c(2,4),]$Comp1*5), y = (-Att_loadings_1[c(2,4),]$Comp2*5),
     #      label = c("Growth", "Reproduction")) +
 # annotate("text", x = (Att_loadings_1[c(1,3),]$Comp1*5), y = (Att_loadings_1[c(1,3),]$Comp2*5),
          # label = c("Survival", "Retrogression")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c('SSurvAs', 'SGrowAs', 'SShriAs', 'SRepAs')
Att.coefs.bar <- ggplot() + 
  geom_bar(aes(x=names(att.coefficients), y=att.coefficients), stat="identity", fill = "darkgreen") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed dataset
Att_2 <- data.frame(cbind(Att.2[["Hypothesis 1 pls"]]$scores[,1],
                          Att.2[["Hypothesis 1 pls"]]$scores[,2],
                          Att.2[["Hypothesis 1 pls"]]$scores[,3],
                          Att.2[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                          imputed_data.Att[,"Attenuation2"]))
colnames(Att_2) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"Attenuation") # puts the detail back in

# store PCA loadings.
Att_loadings_2 <- data.frame(cbind(Att.2[["Hypothesis 1 pls"]]$loadings[,1],
                                   Att.2[["Hypothesis 1 pls"]]$loadings[,2],
                                   Att.2[["Hypothesis 1 pls"]]$loadings[,3],
                                   Att.2[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(Att_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.Att2 <- explvar(Att.2[["Hypothesis 1 pls"]]) # percentage variance explained
r2.Att2 <- R2(Att.2[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
att.coefficients2 = c(coef(Att.2[["Hypothesis 1 pls"]])[1], -coef(Att.2[["Hypothesis 1 pls"]])[2], coef(Att.2[["Hypothesis 1 pls"]])[3], -coef(Att.2[["Hypothesis 1 pls"]])[4])
# standardise coefficients so that they sum to 1
sum.att.coef2 = sum(sapply(att.coefficients2, abs))
att.coefficients2 = att.coefficients2 / sum.att.coef2
names(att.coefficients2) <- c("SSurvAs", "SGrowAs", "SShriAs", "SRepAs")
att.coefficients2 = sort(att.coefficients2)

# what range is required for the transient colour scale
range(Att_2$Attenuation, na.rm = T); median(Att_2$Attenuation, na.rm = T); mean(Att_2$Attenuation, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# create pls biplot
Att.imputed.plot <- ggplot(data = Att_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Att2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Att2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Attenuation), alpha = 0.7, data = Att_2) + #add data points
  scale_color_gradientn(colors = pal(dim(Att_2)[1]),
                        values = scales::rescale(c(0.0,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Att_loadings_2[c(2,4),], aes(x = 0, y = 0, xend = (-Comp1*8), yend = (-Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  geom_segment(data = Att_loadings_2[c(1,3),], aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (-Att_loadings_2[c(2,4),]$Comp1*5), y = (-Att_loadings_2[c(2,4),]$Comp2*5),
           #label = c("Growth", "Reproduction")) +
  #annotate("text", x = (Att_loadings_2[c(1,3),]$Comp1*5), y = (Att_loadings_2[c(1,3),]$Comp2*5),
           #label = c("Survival", "Retrogression")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c('SSurvAs', 'SGrowAs', 'SShriAs', 'SRepAs')
Att.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(att.coefficients2), y=att.coefficients2), stat="identity", fill = "darkgreen") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())



#### Max Amplification ----------------------------------------------------------------

# Perfect data
# extract key data
Amp_1 <- data.frame(cbind(MAmp.1[["Hypothesis 1 pls"]]$scores[,1],
                          MAmp.1[["Hypothesis 1 pls"]]$scores[,2],
                          MAmp.1[["Hypothesis 1 pls"]]$scores[,3],
                          MAmp.1[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                          perfect_MAmp_data[,"MaxAmplification2"]))
colnames(Amp_1) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"MaxAmplification") # puts the detail back in

# store PCA loadings.
Amp_loadings_1 <- data.frame(cbind(MAmp.1[["Hypothesis 1 pls"]]$loadings[,1],
                                   MAmp.1[["Hypothesis 1 pls"]]$loadings[,2],
                                   MAmp.1[["Hypothesis 1 pls"]]$loadings[,3],
                                   MAmp.1[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(Amp_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.Amp <- explvar(MAmp.1[["Hypothesis 1 pls"]]) # percentage variance explained
r2.Amp <- R2(MAmp.1[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Amp.coefficients = c(-coef(MAmp.1[["Hypothesis 1 pls"]])[1], coef(MAmp.1[["Hypothesis 1 pls"]])[2], -coef(MAmp.1[["Hypothesis 1 pls"]])[3], coef(MAmp.1[["Hypothesis 1 pls"]])[4])
# standardise coefficients so that they sum to 1
sum.Amp.coef = sum(sapply(Amp.coefficients, abs))
Amp.coefficients = Amp.coefficients / sum.Amp.coef
names(Amp.coefficients) <- c("SSurvMAs", "SGrowMAs", "SShriMAs", "SRepMAs")
Amp.coefficients = sort(Amp.coefficients)

# what range is required for the transient colour scale
range(Amp_1$MaxAmplification, na.rm = T); median(Amp_1$MaxAmplification, na.rm = T); mean(Amp_1$MaxAmplification, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# create pls biplot
Amp.perfect.plot <- ggplot(data = Amp_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Amp[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Amp[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAmplification), alpha = 0.7, data = Amp_1) + #add data points
  scale_color_gradientn(colors = pal(dim(Amp_1)[1]),
                        values = scales::rescale(c(0.0,0.60,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Amp_loadings_1[c(2,4),], aes(x = 0, y = 0, xend = (-Comp1*9), yend = (-Comp2*9)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  geom_segment(data = Amp_loadings_1[c(1,3),], aes(x = 0, y = 0, xend = (Comp1*9), yend = (Comp2*9)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (-Amp_loadings_1[c(2,4),]$Comp1*7), y = (-Amp_loadings_1[c(2,4),]$Comp2*7),
         #  label = c("Growth", "Reproduction")) +
  #annotate("text", x = (Amp_loadings_1[c(1,3),]$Comp1*7), y = (Amp_loadings_1[c(1,3),]$Comp2*7),
         #  label = c("Survival", "Retrogression")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c('SSurvMAs', 'SGrowMAs', 'SShriMAs', 'SRepMAs')
Amp.coefs.bar <- ggplot() + 
  geom_bar(aes(x=names(Amp.coefficients), y=Amp.coefficients), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed dataset
Amp_2 <- data.frame(cbind(MAmp.2[["Hypothesis 1 pls"]]$scores[,1],
                          MAmp.2[["Hypothesis 1 pls"]]$scores[,2],
                          MAmp.2[["Hypothesis 1 pls"]]$scores[,3],
                          MAmp.2[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                          imputed_data.MAmp[,"MaxAmplification2"]))
colnames(Amp_2) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"MaxAmplification") # puts the detail back in

# store PCA loadings.
Amp_loadings_2 <- data.frame(cbind(MAmp.2[["Hypothesis 1 pls"]]$loadings[,1],
                                   MAmp.2[["Hypothesis 1 pls"]]$loadings[,2],
                                   MAmp.2[["Hypothesis 1 pls"]]$loadings[,3],
                                   MAmp.2[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(Amp_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.Amp2 <- explvar(MAmp.2[["Hypothesis 1 pls"]]) # percentage variance explained
r2.Amp2 <- R2(MAmp.2[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Amp.coefficients2 = c(-coef(MAmp.2[["Hypothesis 1 pls"]])[1], coef(MAmp.2[["Hypothesis 1 pls"]])[2], -coef(MAmp.2[["Hypothesis 1 pls"]])[3], coef(MAmp.2[["Hypothesis 1 pls"]])[4])
# standardise coefficients so that they sum to 1
sum.Amp.coef2 = sum(sapply(Amp.coefficients2, abs))
Amp.coefficients2 = Amp.coefficients2 / sum.Amp.coef2
names(Amp.coefficients2) <- c("SSurvMAs", "SGrowMAs", "SShriMAs", "SRepMAs")
Amp.coefficients2 = sort(Amp.coefficients2)

# what range is required for the transient colour scale
range(Amp_2$MaxAmplification, na.rm = T); median(Amp_2$MaxAmplification, na.rm = T); mean(Amp_2$MaxAmplification, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# create pls biplot
Amp.imputed.plot <- ggplot(data = Amp_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Amp2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Amp2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAmplification), alpha = 0.7, data = Amp_2) + #add data points
  scale_color_gradientn(colors = pal(dim(Amp_2)[1]),
                        values = scales::rescale(c(-0.49,0.60,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Amp_loadings_2[c(2,4),], aes(x = 0, y = 0, xend = (-Comp1*7), yend = (-Comp2*7)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  geom_segment(data = Amp_loadings_2[c(1,3),], aes(x = 0, y = 0, xend = (Comp1*7), yend = (Comp2*7)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
 # annotate("text", x = (-Amp_loadings_2[c(2,4),]$Comp1*7), y = (-Amp_loadings_2[c(2,4),]$Comp2*7),
          # label = c("Growth", "Reproduction")) +
 # annotate("text", x = (Amp_loadings_2[c(1,3),]$Comp1*7), y = (Amp_loadings_2[c(1,3),]$Comp2*7),
          # label = c("Survival", "Retrogression")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c('SSurvMAs', 'SGrowMAs', 'SShriMAs', 'SRepMAs')
Amp.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(Amp.coefficients2), y=Amp.coefficients2), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())


### Maximal Attenuation ------------------------------------------

# Perfect data
# extract key data
MAtt_1 <- data.frame(cbind(MAtt.1[["Hypothesis 1 pls"]]$scores[,1],
                           MAtt.1[["Hypothesis 1 pls"]]$scores[,2],
                           MAtt.1[["Hypothesis 1 pls"]]$scores[,3],
                           MAtt.1[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                           perfect_MAtt_data[,"MaxAttenuation2"]))
colnames(MAtt_1) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"MaxAttenuation") # puts the detail back in

# store PCA loadings.
MAtt_loadings_1 <- data.frame(cbind(MAtt.1[["Hypothesis 1 pls"]]$loadings[,1],
                                    MAtt.1[["Hypothesis 1 pls"]]$loadings[,2],
                                    MAtt.1[["Hypothesis 1 pls"]]$loadings[,3],
                                    MAtt.1[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(MAtt_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.MAtt <- explvar(MAtt.1[["Hypothesis 1 pls"]]) # percentage variance explained
r2.MAtt <- R2(MAtt.1[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
MAtt.coefficients = c(coef(MAtt.1[["Hypothesis 1 pls"]])[1:4])
# standardise coefficients so that they sum to 1
sum.MAtt.coef = sum(sapply(MAtt.coefficients, abs))
MAtt.coefficients = MAtt.coefficients / sum.MAtt.coef
names(MAtt.coefficients) <- c("SSurvMAtts", "SGrowMAtts", "SShriMAtts", "SRepMAtts")
MAtt.coefficients = sort(MAtt.coefficients)

# what range is required for the transient colour scale
range(MAtt_1$MaxAttenuation, na.rm = T); median(MAtt_1$MaxAttenuation, na.rm = T); mean(MAtt_1$MaxAttenuation, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# create pls biplot
MAtt.perfect.plot <- ggplot(data = MAtt_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.MAtt[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.MAtt[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAttenuation), alpha = 0.7,  data = MAtt_1) + #add data points
  scale_color_gradientn(colors = pal(dim(MAtt_1)[1]),
                        values = scales::rescale(c(0.0,0.53,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = MAtt_loadings_1[c(1:4),], aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (MAtt_loadings_1[c(1:4),]$Comp1*5), y = (MAtt_loadings_1[c(1:4),]$Comp2*5),
          # label = c("Survival", "Growth", "Retrogression", "Reproduction")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot
level_order <- c('SSurvMAtts', 'SGrowMAtts', 'SShriMAtts', 'SRepMAtts')
MAtt.coefs.bar <- ggplot() + 
  geom_bar(aes(x=names(MAtt.coefficients), y=MAtt.coefficients), stat="identity", fill = "darkgreen") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed data
MAtt_2 <- data.frame(cbind(MAtt.2[["Hypothesis 1 pls"]]$scores[,1],
                           MAtt.2[["Hypothesis 1 pls"]]$scores[,2],
                           MAtt.2[["Hypothesis 1 pls"]]$scores[,3],
                           MAtt.2[["Hypothesis 1 pls"]]$scores[,4], #this accounts for the odd class used for storing pls score outputs
                           imputed_data.MAtt[,"MaxAttenuation2"]))
colnames(MAtt_2) <- c("Comp1", "Comp2", "Comp3", "Comp4" ,"MaxAttenuation") # puts the detail back in

# store PCA loadings.
MAtt_loadings_2 <- data.frame(cbind(MAtt.2[["Hypothesis 1 pls"]]$loadings[,1],
                                    MAtt.2[["Hypothesis 1 pls"]]$loadings[,2],
                                    MAtt.2[["Hypothesis 1 pls"]]$loadings[,3],
                                    MAtt.2[["Hypothesis 1 pls"]]$loadings[,4]))
colnames(MAtt_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4")

# extract the variance explained by each component
p.var.MAtt2 <- explvar(MAtt.2[["Hypothesis 1 pls"]]) # percentage variance explained
r2.MAtt2 <- R2(MAtt.2[["Hypothesis 1 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
MAtt.coefficients2 = c(coef(MAtt.2[["Hypothesis 1 pls"]])[1:4])
# standardise coefficients so that they sum to 1
sum.MAtt.coef2 = sum(sapply(MAtt.coefficients2, abs))
MAtt.coefficients2 = MAtt.coefficients2 / sum.MAtt.coef2
names(MAtt.coefficients2) <- c("SSurvMAtts", "SGrowMAtts", "SShriMAtts", "SRepMAtts")
MAtt.coefficients2 = sort(MAtt.coefficients2)

# what range is required for the transient colour scale
range(MAtt_2$MaxAttenuation, na.rm = T); median(MAtt_2$MaxAttenuation, na.rm = T); mean(MAtt_2$MaxAttenuation, na.rm = T)

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# create pls biplot
MAtt.imputed.plot <- ggplot(data = MAtt_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.MAtt2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.MAtt2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAttenuation), alpha = 0.7, data = MAtt_2) + #add data points
  scale_color_gradientn(colors = pal(dim(MAtt_2)[1]),
                        values = scales::rescale(c(0.0,0.55,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = MAtt_loadings_2[c(1:4),], aes(x = 0, y = 0, xend = (Comp1*8), yend = (Comp2*8)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
 # annotate("text", x = (MAtt_loadings_2[c(1:4),]$Comp1*5), y = (MAtt_loadings_2[c(1:4),]$Comp2*5),
           #label = c("Survival", "Growth", "Retrogression", "Reproduction")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot
level_order <- c('SSurvMAtts', 'SGrowMAtts', 'SShriMAtts', 'SRepMAtts')
MAtt.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(MAtt.coefficients2), y=MAtt.coefficients2), stat="identity", fill = "darkgreen") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())


####### ---------------------- Hypothesis 2 --------------------------------------

#### Damping ratio ----------------------------------------------------------------
# extract key data
DR_3 <- data.frame(cbind(Damping.ratio.1[["Hypothesis 2 pls"]]$scores[,1],
                         Damping.ratio.1[["Hypothesis 2 pls"]]$scores[,2],
                         Damping.ratio.1[["Hypothesis 2 pls"]]$scores[,3],
                         Damping.ratio.1[["Hypothesis 2 pls"]]$scores[,4],
                         Damping.ratio.1[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         Damping.ratio.1[["Transient variable"]]))
colnames(DR_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Rho") # puts the detail back in

# store PCA loadings.
DR_loadings_3 <- data.frame(cbind(Damping.ratio.1[["Hypothesis 2 pls"]]$loadings[,1],
                                  Damping.ratio.1[["Hypothesis 2 pls"]]$loadings[,2],
                                  Damping.ratio.1[["Hypothesis 2 pls"]]$loadings[,3],
                                  Damping.ratio.1[["Hypothesis 2 pls"]]$loadings[,4],
                                  Damping.ratio.1[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(DR_loadings_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.DR3 <- explvar(Damping.ratio.1[["Hypothesis 2 pls"]]) # percentage variance explained
r2.DR3 <- R2(Damping.ratio.1[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
DR.coefficients3 = -coef(Damping.ratio.1[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.DR.coef3 = sum(sapply(DR.coefficients3, abs))
DR.coefficients3 = DR.coefficients3 / sum.DR.coef3
DR.coefficients3 = sort(DR.coefficients3[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# descriptive stats for the damping ratio variable
range(DR_3$Rho, na.rm = T); median(DR_3$Rho, na.rm = T); mean(DR_3$Rho, na.rm = T)

# create pls biplot
DR.perfect.plot2 <- ggplot(data = DR_3, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.DR3[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.DR3[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Rho), alpha = 0.7, data = DR_3) + #add data points
  scale_color_gradientn(colors = pal(dim(DR_3)[1]),
                        values = scales::rescale(c(0.11,0.63,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = DR_loadings_3, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (DR_loadings_3$Comp1*4), y = (DR_loadings_3$Comp2*4),
           #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
DR.coefs.bar3 <- ggplot() + 
  geom_bar(aes(x=names(DR.coefficients3), y=DR.coefficients3), stat="identity", fill = "gold") +
  ylab("Proportional influence") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Again repeat for Imputed dataset
DR_4 <- data.frame(cbind(Damping.ratio.2[["Hypothesis 2 pls"]]$scores[,1],
                         Damping.ratio.2[["Hypothesis 2 pls"]]$scores[,2],
                         Damping.ratio.2[["Hypothesis 2 pls"]]$scores[,3],
                         Damping.ratio.2[["Hypothesis 2 pls"]]$scores[,4],
                         Damping.ratio.2[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         Damping.ratio.2[["Transient variable"]]))
colnames(DR_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Rho") # puts the detail back in

# store PCA loadings.
DR_loadings_4 <- data.frame(cbind(Damping.ratio.2[["Hypothesis 2 pls"]]$loadings[,1],
                                  Damping.ratio.2[["Hypothesis 2 pls"]]$loadings[,2],
                                  Damping.ratio.2[["Hypothesis 2 pls"]]$loadings[,3],
                                  Damping.ratio.2[["Hypothesis 2 pls"]]$loadings[,4],
                                  Damping.ratio.2[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(DR_loadings_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.DR4 <- explvar(Damping.ratio.2[["Hypothesis 2 pls"]]) # percentage variance explained
r2.DR4 <- R2(Damping.ratio.2[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
DR.coefficients4 = -coef(Damping.ratio.2[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.DR.coef4 = sum(sapply(DR.coefficients4, abs))
DR.coefficients4 = DR.coefficients4 / sum.DR.coef4
DR.coefficients4 = sort(DR.coefficients4[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# descriptive stats for the damping ratio variable
range(DR_4$Rho, na.rm = T); median(DR_4$Rho, na.rm = T); mean(DR_4$Rho, na.rm = T)

# create pls biplot
DR.imputed.plot2 <- ggplot(data = DR_4, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.DR4[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.DR4[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Rho), alpha = 0.7, data = DR_4) + #add data points
  scale_color_gradientn(colors = pal(dim(DR_4)[1]),
                        values = scales::rescale(c(-191.6047, 0.6006718 ,100), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = DR_loadings_4, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  annotate("text", x = (DR_loadings_4$Comp1*4), y = (DR_loadings_4$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
DR.coefs.bar4 <- ggplot() + 
  geom_bar(aes(x=names(DR.coefficients4), y=DR.coefficients4), stat="identity", fill = "gold") +
  ylab("Proportional influence") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())


#### Period of Ocilation ---------------------------------------------------------------- 
# extract key data
PO_3 <- data.frame(cbind(PO.1[["Hypothesis 2 pls"]]$scores[,1],
                         PO.1[["Hypothesis 2 pls"]]$scores[,2],
                         PO.1[["Hypothesis 2 pls"]]$scores[,3],
                         PO.1[["Hypothesis 2 pls"]]$scores[,4],
                         PO.1[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         PO.1[["Transient variable"]]))
colnames(PO_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5","Pi") # puts the detail back in

# store PCA loadings.
PO_loadings_3 <- data.frame(cbind(PO.1[["Hypothesis 2 pls"]]$loadings[,1],
                                  PO.1[["Hypothesis 2 pls"]]$loadings[,2],
                                  PO.1[["Hypothesis 2 pls"]]$loadings[,3],
                                  PO.1[["Hypothesis 2 pls"]]$loadings[,4],
                                  PO.1[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(PO_loadings_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.PO3 <- explvar(PO.1[["Hypothesis 2 pls"]]) # percentage variance explained
r2.PO3 <- R2(PO.1[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where necessary)
PO.coefficients3 = -coef(PO.1[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.PO.coef3 = sum(sapply(PO.coefficients3, abs))
PO.coefficients3 = PO.coefficients3 / sum.PO.coef3
PO.coefficients3 = sort(PO.coefficients3[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# descriptive stats for the damping ratio variable
range(PO_3$Pi, na.rm = T); median(PO_3$Pi, na.rm = T); mean(PO_3$Pi, na.rm = T)

# create pls biplot
PO.perfect.plot2 <- ggplot(data = PO_3, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.PO3[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.PO3[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Pi), alpha = 0.7, data = PO_3) + #add data points
  scale_color_gradientn(colors = pal(dim(PO_3)[1]),
                        values = scales::rescale(c(0.08,0.35,0.58), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = PO_loadings_3, aes(x = 0, y = 0, xend = (Comp1*3), yend = (Comp2*3)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (PO_loadings_3$Comp1*3), y = (PO_loadings_3$Comp2*3),
           #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
PO.coefs.bar3 <- ggplot() + 
  geom_bar(aes(x=names(PO.coefficients3), y=PO.coefficients3), stat="identity", fill = "gold") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed data
PO_4 <- data.frame(cbind(PO.2[["Hypothesis 2 pls"]]$scores[,1],
                         PO.2[["Hypothesis 2 pls"]]$scores[,2],
                         PO.2[["Hypothesis 2 pls"]]$scores[,3],
                         PO.2[["Hypothesis 2 pls"]]$scores[,4],
                         PO.2[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         PO.2[["Transient variable"]]))
colnames(PO_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5","Pi") # puts the detail back in

# store PCA loadings.
PO_loadings_4 <- data.frame(cbind(PO.2[["Hypothesis 2 pls"]]$loadings[,1],
                                  PO.2[["Hypothesis 2 pls"]]$loadings[,2],
                                  PO.2[["Hypothesis 2 pls"]]$loadings[,3],
                                  PO.2[["Hypothesis 2 pls"]]$loadings[,4],
                                  PO.2[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(PO_loadings_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.PO4 <- explvar(PO.2[["Hypothesis 2 pls"]]) # percentage variance explained
r2.PO4 <- R2(PO.2[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
PO.coefficients4 = -coef(PO.2[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.PO.coef4 = sum(sapply(PO.coefficients4, abs))
PO.coefficients4 = PO.coefficients4 / sum.PO.coef4
PO.coefficients4 = sort(PO.coefficients4[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# descriptive stats for the damping ratio variable
range(PO_4$Pi, na.rm = T); median(PO_4$Pi, na.rm = T); mean(PO_4$Pi, na.rm = T)

# create pls biplot
PO.imputed.plot2 <- ggplot(data = PO_4, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.PO4[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.PO4[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Pi), alpha = 0.7, data = PO_4) + #add data points
  scale_color_gradientn(colors = pal(dim(PO_4)[1]),
                        values = scales::rescale(c(0.08, 0.34, 0.58), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = NA) +
  geom_segment(data = PO_loadings_4, aes(x = 0, y = 0, xend = (Comp1*3), yend = (Comp2*3)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  annotate("text", x = (PO_loadings_4$Comp1*3), y = (PO_loadings_4$Comp2*3),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
PO.coefs.bar4 <- ggplot() + 
  geom_bar(aes(x=names(PO.coefficients4), y=PO.coefficients4), stat="identity", fill = "gold") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())



#### Reactivity -------------------------------------------------------------------------
# extract key data
reac_3 <- data.frame(cbind(Reac.1[["Hypothesis 2 pls"]]$scores[,1],
                           Reac.1[["Hypothesis 2 pls"]]$scores[,2],
                           Reac.1[["Hypothesis 2 pls"]]$scores[,3],
                           Reac.1[["Hypothesis 2 pls"]]$scores[,4],
                           Reac.1[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Reac.1[["Transient variable"]]))
colnames(reac_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Reactivity") # puts the detail back in

# store PCA loadings.
reac_loadings_3 <- data.frame(cbind(Reac.1[["Hypothesis 2 pls"]]$loadings[,1],
                                    Reac.1[["Hypothesis 2 pls"]]$loadings[,2],
                                    Reac.1[["Hypothesis 2 pls"]]$loadings[,3],
                                    Reac.1[["Hypothesis 2 pls"]]$loadings[,4],
                                    Reac.1[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(reac_loadings_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.R3 <- explvar(Reac.1[["Hypothesis 2 pls"]]) # percentage variance explained
r2.R3 <- R2(Reac.1[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
reac.coefficients3 = -coef(Reac.1[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.reac.coef3 = sum(sapply(reac.coefficients3, abs))
reac.coefficients3 = reac.coefficients3 / sum.reac.coef3
reac.coefficients3 = sort(reac.coefficients3[, 1, 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(reac_3$Reactivity, na.rm = T); median(reac_3$Reactivity, na.rm = T); mean(reac_3$Reactivity, na.rm = T)

# create pls biplot
R.perfect.plot2 <- ggplot(data = reac_3, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.R3[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.R3[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Reactivity), alpha = 0.7, data = reac_3) + #add data points
  scale_color_gradientn(colors = pal(dim(reac_3)[1]),
                        values = scales::rescale(c(0,0.62,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = reac_loadings_3, aes(x = 0, y = 0, xend = (Comp1*5), yend = (Comp2*5)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
 # annotate("text", x = (reac_loadings_3$Comp1*5), y = (reac_loadings_3$Comp2*5),
 #          label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
R.coefs.bar3 <- ggplot() + 
  geom_bar(aes(x=names(reac.coefficients3), y=reac.coefficients3), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed data
reac_4 <- data.frame(cbind(Reac.2[["Hypothesis 2 pls"]]$scores[,1],
                           Reac.2[["Hypothesis 2 pls"]]$scores[,2],
                           Reac.2[["Hypothesis 2 pls"]]$scores[,3],
                           Reac.2[["Hypothesis 2 pls"]]$scores[,4],
                           Reac.2[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Reac.2[["Transient variable"]]))
colnames(reac_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Reactivity") # puts the detail back in

# store PCA loadings.
reac_loadings_4 <- data.frame(cbind(Reac.2[["Hypothesis 2 pls"]]$loadings[,1],
                                    Reac.2[["Hypothesis 2 pls"]]$loadings[,2],
                                    Reac.2[["Hypothesis 2 pls"]]$loadings[,3],
                                    Reac.2[["Hypothesis 2 pls"]]$loadings[,4],
                                    Reac.2[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(reac_loadings_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.R4 <- explvar(Reac.2[["Hypothesis 2 pls"]]) # percentage variance explained
r2.R4 <- R2(Reac.2[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
reac.coefficients4 = -coef(Reac.2[["Hypothesis 2 pls"]])
# standardize coefficients so that they sum to 1
sum.reac.coef4 = sum(sapply(reac.coefficients4, abs))
reac.coefficients4 = reac.coefficients4 / sum.reac.coef4
reac.coefficients4 = sort(reac.coefficients4[, 1, 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(reac_4$Reactivity, na.rm = T); median(reac_4$Reactivity, na.rm = T); mean(reac_4$Reactivity, na.rm = T)

# create pls biplot
R.imputed.plot2 <- ggplot(data = reac_4, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.R4[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.R4[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Reactivity), alpha = 0.7, data = reac_4) + #add data points
  scale_color_gradientn(colors = pal(dim(reac_4)[1]),
                        values = scales::rescale(c(-0.72,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = reac_loadings_4, aes(x = 0, y = 0, xend = (Comp1*5), yend = (Comp2*5)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  annotate("text", x = (reac_loadings_4$Comp1*5), y = (reac_loadings_4$Comp2*5),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
R.coefs.bar4 <- ggplot() + 
  geom_bar(aes(x=names(reac.coefficients4), y=reac.coefficients4), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())




#### Attenuation ------------------------------------------------------------------------
# extract key data
Att_3 <-  data.frame(cbind(Att.1[["Hypothesis 2 pls"]]$scores[,1],
                           Att.1[["Hypothesis 2 pls"]]$scores[,2],
                           Att.1[["Hypothesis 2 pls"]]$scores[,3],
                           Att.1[["Hypothesis 2 pls"]]$scores[,4],
                           Att.1[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Att.1[["Transient variable"]]))
colnames(Att_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Attenuation") # puts the detail back in

# store PCA loadings.
Att_loadings_3  <- data.frame(cbind(Att.1[["Hypothesis 2 pls"]]$loadings[,1],
                                    Att.1[["Hypothesis 2 pls"]]$loadings[,2],
                                    Att.1[["Hypothesis 2 pls"]]$loadings[,3],
                                    Att.1[["Hypothesis 2 pls"]]$loadings[,4],
                                    Att.1[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Att_loadings_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Att3 <- explvar(Att.1[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Att3 <- R2(Att.1[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Att.coefficients3 = coef(Att.1[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Att.coef3 = sum(sapply(Att.coefficients3, abs))
Att.coefficients3 = Att.coefficients3 / sum.Att.coef3
Att.coefficients3 = sort(Att.coefficients3[, 1, 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(Att_3$Attenuation, na.rm = T); median(Att_3$Attenuation, na.rm = T); mean(Att_3$Attenuation, na.rm = T)

# create pls biplot
Att.perfect.plot2 <- ggplot(data = Att_3, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Att3[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Att3[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Attenuation), alpha = 0.7, data = Att_3) + #add data points
  scale_color_gradientn(colors = pal(dim(Att_3)[1]),
                        values = scales::rescale(c(0,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Att_loadings_3, aes(x = 0, y = 0, xend = (Comp1*5), yend = (Comp2*5)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (Att_loadings_3$Comp1*5), y = (Att_loadings_3$Comp2*5),
          # label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
Att.coefs.bar3 <- ggplot() + 
  geom_bar(aes(x=names(Att.coefficients3), y=Att.coefficients3), stat="identity", fill = "darkgreen") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed data
Att_4  <- data.frame(cbind(Att.2[["Hypothesis 2 pls"]]$scores[,1],
                           Att.2[["Hypothesis 2 pls"]]$scores[,2],
                           Att.2[["Hypothesis 2 pls"]]$scores[,3],
                           Att.2[["Hypothesis 2 pls"]]$scores[,4],
                           Att.2[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Att.2[["Transient variable"]]))
colnames(Att_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Attenuation") # puts the detail back in

# store PCA loadings.
Att_loadings_4  <- data.frame(cbind(Att.2[["Hypothesis 2 pls"]]$loadings[,1],
                                    Att.2[["Hypothesis 2 pls"]]$loadings[,2],
                                    Att.2[["Hypothesis 2 pls"]]$loadings[,3],
                                    Att.2[["Hypothesis 2 pls"]]$loadings[,4],
                                    Att.2[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Att_loadings_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Att4 <- explvar(Att.2[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Att4 <- R2(Att.2[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Att.coefficients4 = coef(Att.2[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Att.coef4 = sum(sapply(Att.coefficients4, abs))
Att.coefficients4 = Att.coefficients4 / sum.Att.coef4
Att.coefficients4 = sort(Att.coefficients4[, 1, 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(Att_4$Attenuation, na.rm = T); median(Att_4$Attenuation, na.rm = T); mean(Att_4$Attenuation, na.rm = T)

# create pls biplot
Att.imputed.plot2 <- ggplot(data = Att_4, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Att4[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Att4[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Attenuation), alpha = 0.7, data = Att_4) + #add data points
  scale_color_gradientn(colors = pal(dim(Att_4)[1]),
                        values = scales::rescale(c(0,0.66,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Att_loadings_4, aes(x = 0, y = 0, xend = (Comp1*5), yend = (Comp2*5)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  annotate("text", x = (Att_loadings_4$Comp1*5), y = (Att_loadings_4$Comp2*5),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
Att.coefs.bar4 <- ggplot() + 
  geom_bar(aes(x=names(Att.coefficients4), y=Att.coefficients4), stat="identity", fill = c("darkgreen")) +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())


#### Max Amplification ----------------------------------------------------------------
# extract key data
Amp_3 <- data.frame(cbind(MAmp.1[["Hypothesis 2 pls"]]$scores[,1],
                          MAmp.1[["Hypothesis 2 pls"]]$scores[,2],
                          MAmp.1[["Hypothesis 2 pls"]]$scores[,3],
                          MAmp.1[["Hypothesis 2 pls"]]$scores[,4],
                          MAmp.1[["Hypothesis 2 pls"]]$scores[,5],
                          #this accounts for the odd class used for storing pls score outputs
                          MAmp.1[["Transient variable"]]))
colnames(Amp_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "MaxAmplification") # puts the detail back in

# store PCA loadings.
Amp_loadings_3 <- data.frame(cbind(MAmp.1[["Hypothesis 2 pls"]]$loadings[,1],
                                   MAmp.1[["Hypothesis 2 pls"]]$loadings[,2],
                                   MAmp.1[["Hypothesis 2 pls"]]$loadings[,3],
                                   MAmp.1[["Hypothesis 2 pls"]]$loadings[,4],
                                   MAmp.1[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Amp_loadings_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Amp3 <- explvar(MAmp.1[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Amp3 <- R2(MAmp.1[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Amp.coefficients3 = -coef(MAmp.1[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Amp.coef3 = sum(sapply(Amp.coefficients3, abs))
Amp.coefficients3 = Amp.coefficients3 / sum.Amp.coef3
Amp.coefficients3 = sort(Amp.coefficients3[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(Amp_3$MaxAmplification, na.rm = T); median(Amp_3$MaxAmplification, na.rm = T); mean(Amp_3$MaxAmplification, na.rm = T)

# create pls biplot
Amp.perfect.plot2 <- ggplot(data = Amp_3, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Amp3[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Amp3[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAmplification), alpha = 0.7, data = Amp_3) + #add data points
  scale_color_gradientn(colors = pal(dim(Amp_3)[1]),
                        values = scales::rescale(c(0.0,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Amp_loadings_3, aes(x = 0, y = 0, xend = (Comp1*3), yend = (Comp2*3)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (Amp_loadings_3$Comp1*3), y = (Amp_loadings_3$Comp2*3),
           #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
Amp.coefs.bar3 <- ggplot() + 
  geom_bar(aes(x=names(Amp.coefficients3), y=Amp.coefficients3), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# And the imputed data set
Amp_4 <- data.frame(cbind(MAmp.2[["Hypothesis 2 pls"]]$scores[,1],
                          MAmp.2[["Hypothesis 2 pls"]]$scores[,2],
                          MAmp.2[["Hypothesis 2 pls"]]$scores[,3],
                          MAmp.2[["Hypothesis 2 pls"]]$scores[,4],
                          MAmp.2[["Hypothesis 2 pls"]]$scores[,5],
                          #this accounts for the odd class used for storing pls score outputs
                          MAmp.2[["Transient variable"]]))
colnames(Amp_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "MaxAmplification") # puts the detail back in

# store PCA loadings.
Amp_loadings_4 <- data.frame(cbind(MAmp.2[["Hypothesis 2 pls"]]$loadings[,1],
                                   MAmp.2[["Hypothesis 2 pls"]]$loadings[,2],
                                   MAmp.2[["Hypothesis 2 pls"]]$loadings[,3],
                                   MAmp.2[["Hypothesis 2 pls"]]$loadings[,4],
                                   MAmp.2[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Amp_loadings_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Amp4 <- explvar(MAmp.2[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Amp4 <- R2(MAmp.2[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Amp.coefficients4 = -coef(MAmp.2[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Amp.coef4 = sum(sapply(Amp.coefficients4, abs))
Amp.coefficients4 = Amp.coefficients4 / sum.Amp.coef4
Amp.coefficients4 = sort(Amp.coefficients4[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(Amp_4$MaxAmplification, na.rm = T); median(Amp_4$MaxAmplification, na.rm = T); mean(Amp_4$MaxAmplification, na.rm = T)

# create pls biplot
Amp.imputed.plot2 <- ggplot(data = Amp_4, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Amp4[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Amp4[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAmplification), alpha = 0.7, data = Amp_4) + #add data points
  scale_color_gradientn(colors = pal(dim(Amp_4)[1]),
                        values = scales::rescale(c(-0.49,0.6,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = Amp_loadings_4, aes(x = 0, y = 0, xend = (Comp1*3), yend = (Comp2*3)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  annotate("text", x = (Amp_loadings_4$Comp1*3), y = (Amp_loadings_4$Comp2*3),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
Amp.coefs.bar4 <- ggplot() + 
  geom_bar(aes(x=names(Amp.coefficients4), y=Amp.coefficients4), stat="identity", fill = "navyblue") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())



#### Maximal Attenuation -----------------------------------------------------------
# extract key data
MAtt_3 <- data.frame(cbind(MAtt.1[["Hypothesis 2 pls"]]$scores[,1],
                           MAtt.1[["Hypothesis 2 pls"]]$scores[,2],
                           MAtt.1[["Hypothesis 2 pls"]]$scores[,3],
                           MAtt.1[["Hypothesis 2 pls"]]$scores[,4],
                           MAtt.1[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           MAtt.1[["Transient variable"]]))
colnames(MAtt_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "MaxAttenuation") # puts the detail back in

# store PCA loadings.
MAtt_loadings_3 <- data.frame(cbind(MAtt.1[["Hypothesis 2 pls"]]$loadings[,1],
                                    MAtt.1[["Hypothesis 2 pls"]]$loadings[,2],
                                    MAtt.1[["Hypothesis 2 pls"]]$loadings[,3],
                                    MAtt.1[["Hypothesis 2 pls"]]$loadings[,4],
                                    MAtt.1[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(MAtt_loadings_3) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.MAtt3 <- explvar(MAtt.1[["Hypothesis 2 pls"]]) # percentage variance explained
r2.MAtt3 <- R2(MAtt.1[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
MAtt.coefficients3 = coef(MAtt.1[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.MAtt.coef3 = sum(sapply(MAtt.coefficients3, abs))
MAtt.coefficients3 = MAtt.coefficients3 / sum.MAtt.coef3
MAtt.coefficients3 = sort(MAtt.coefficients3[, 1, 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(MAtt_3$MaxAttenuation, na.rm = T); median(MAtt_3$MaxAttenuation, na.rm = T); mean(MAtt_3$MaxAttenuation, na.rm = T)

# create pls biplot
MAtt.perfect.plot2 <- ggplot(data = MAtt_3, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.MAtt3[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.MAtt3[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAttenuation), alpha = 0.7, data = MAtt_3) + #add data points
  scale_color_gradientn(colors = pal(dim(MAtt_3)[1]),
                        values = scales::rescale(c(0,0.59,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = MAtt_loadings_3, aes(x = 0, y = 0, xend = (Comp1*5), yend = (Comp2*5)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (MAtt_loadings_3$Comp1*5), y = (MAtt_loadings_3$Comp2*5),
          #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
MAtt.coefs.bar3 <- ggplot() + 
  geom_bar(aes(x=names(MAtt.coefficients3), y=MAtt.coefficients3), stat="identity", fill = "darkgreen") +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())

# Imputed data
MAtt_4 <- data.frame(cbind(MAtt.2[["Hypothesis 2 pls"]]$scores[,1],
                           MAtt.2[["Hypothesis 2 pls"]]$scores[,2],
                           MAtt.2[["Hypothesis 2 pls"]]$scores[,3],
                           MAtt.2[["Hypothesis 2 pls"]]$scores[,4],
                           MAtt.2[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           MAtt.2[["Transient variable"]]))
colnames(MAtt_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "MaxAttenuation") # puts the detail back in

# store PCA loadings.
MAtt_loadings_4  <- data.frame(cbind(MAtt.2[["Hypothesis 2 pls"]]$loadings[,1],
                                    MAtt.2[["Hypothesis 2 pls"]]$loadings[,2],
                                    MAtt.2[["Hypothesis 2 pls"]]$loadings[,3],
                                    MAtt.2[["Hypothesis 2 pls"]]$loadings[,4],
                                    MAtt.2[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(MAtt_loadings_4) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.MAtt4 <- explvar(MAtt.2[["Hypothesis 2 pls"]]) # percentage variance explained
r2.MAtt4 <- R2(MAtt.2[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where necessary)
MAtt.coefficients4 = coef(MAtt.2[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.MAtt.coef4 = sum(sapply(MAtt.coefficients4, abs))
MAtt.coefficients4 = MAtt.coefficients4 / sum.MAtt.coef4
MAtt.coefficients4 = sort(MAtt.coefficients4[, 1, 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(MAtt_4$MaxAttenuation, na.rm = T); median(MAtt_4$MaxAttenuation, na.rm = T); mean(MAtt_4$MaxAttenuation, na.rm = T)

# create pls biplot
MAtt.imputed.plot2 <- ggplot(data = MAtt_4, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.MAtt4[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.MAtt4[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MaxAttenuation), alpha = 0.7, data = MAtt_4) + #add data points
  scale_color_gradientn(colors = pal(dim(MAtt_4)[1]),
                        values = scales::rescale(c(0,0.55,1.0), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  # to account for the inverse transformations within the variables
  geom_segment(data = MAtt_loadings_4, aes(x = 0, y = 0, xend = (Comp1*5), yend = (Comp2*5)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  annotate("text", x = (MAtt_loadings_4$Comp1*5), y = (MAtt_loadings_4$Comp2*5),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())

# Create coefficients barplot 
MAtt.coefs.bar4 <- ggplot() + 
  geom_bar(aes(x=names(MAtt.coefficients4), y=MAtt.coefficients4), stat="identity", fill = c("darkgreen")) +
  geom_hline(yintercept = 0) + 
  ylim(-1,1) +
  scale_x_discrete(limits = level_order,
                   labels=c("", "", "", "")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank())


#**************************** End of Code *******************************************************
