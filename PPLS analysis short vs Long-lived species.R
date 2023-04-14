# This script is for exploring whether our observations of a lack of environmental influence on the transient characterisitics holds for short lived species
# and therefore whether it is merely an outcome of limited abiotic time series.
# To explore this question, this script re-runs a pPLS analysis for short lived and long-lived species.

# Date last modified: May 2021
# Author: James Cant
# ----------------------------------------------------------------------------------------

# load required packages
library(phytools)
library(ape)
library(geiger)
library(popbio)
library(MASS)
library(ggplot2)

# load converted pls SCript
source("FILE_PATH/PLS function for short-lived_long-lived analysis.R") # pPLS function

#################################################
# STEP 1: Determine how to split the populations.
#################################################

# Firstly tidy up the generation time and life expectancy data
# Remove outliers
removeOutliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.025, .975), na.rm = na.rm) #this determines the 25% and 75% limits of the data set based on the supplied vector (ie. where 95% of the data lies) 
  H <- 1.5 * IQR(x, na.rm = na.rm) # determines the interquartile range (variance of the dataset)
  y <- x
  y[x < (qnt[1] - H)] <- NA #any values below or above these ranges are then excluded.
  y[x > (qnt[2] + H)] <- NA
  y
}

numVars <- c("GT","LE")

for (i in numVars){
  data[,i]=removeOutliers(data[,i], na.rm = T)
}

# Now check how many populations would be omitted if the evaluation used generation time or life expectancy
summary(data$GT); summary(data$LE) # there is more data for life expectancy
# what is the distribution?
hist(data$GT)
hist(data$LE)

# Split the populations at a Generation time of 10 years. 
data$life_cat <- NA
data[which(data$LE <= 10), "life_cat"] <- "SL"
data[which(data$LE > 10), "life_cat"] <- "LL"
table(data$life_cat) # 577 long lived and 1606 short lived (2183 populations total)

############### Checkpoint
write.csv(data, file = "csv_file.csv")

#################################################
# STEP 2: Sort the data in preperation for running ppls analyses
#################################################
data <- read.csv("csv_file", row.names =  1) # re-loading the check points

# for this I am only interested in the actual transient characterisitics and the environmental measure of variability. 
# Normal data
key_data <- data[,c("Rho2", "Pi2", "Reactivity2", "Attenuation2", "MaxAmplification2", "MaxAttenuation2", 
                    "Temp.Frequency.full","Temp.Autocorrelation.full", "Temp.range.full",
                    "Prec.Frequency.full", "Prec.Autocorrelation.full", "life_cat")]

# Imputed data
key_data_I <- as.data.frame(cbind(imputed_list.PO[["anc_recon"]][1:dim(data)[1],], imputed_data[,c("Rho2", "Reactivity2", "Attenuation2", "MaxAmplification2", "MaxAttenuation2")]))
colnames(key_data_I) <- c("Imp_Pi2", "Imp_Rho2", "Imp_Reactivity2", "Imp_Attenuation2", "Imp_MaxAmplification2", "Imp_MaxAttenuation2")

# attach the imputed data before splitting the data by GT groups.
key_data <- cbind(key_data, key_data_I)
SL_data <- key_data[which(key_data$life_cat == "SL"),] # Short-lived
LL_data <- key_data[which(key_data$life_cat == "LL"),] # Long-lived
# and prune the tree
SL_name_check <- name.check(phylo.tree, SL_data)
SL.tree <- drop.tip(phylo.tree, SL_name_check$tree_not_data); name.check(SL.tree, SL_data)
LL_name_check <- name.check(phylo.tree, LL_data)
LL.tree <- drop.tip(phylo.tree, LL_name_check$tree_not_data); name.check(LL.tree, LL_data)

# For completeness reassign the seal populations as marine
data[data$SpeciesAccepted %in% c("Zalophus californianus", "Leptonychotes weddellii"), "Realm"] <- "Marine"

# Define a list of marine populations for use in subsetting the data set for the environmental pls analysis.
marine.list <- rownames(data[data$Realm == "Marine",]) # this is 39 populations as expected.

# define datasets of complete entries
### Short-lived populations
## Damping ratio
SL_DR_data <- SL_data[complete.cases(SL_data[,"Rho2"]),]; check <- name.check(SL.tree, SL_DR_data); SL.DR.tree <- drop.tip(SL.tree, check$tree_not_data)
## Period of Oscilation
SL_PO_data <- SL_data[complete.cases(SL_data[,"Pi2"]),]; check <- name.check(SL.tree, SL_PO_data); SL.PO.tree <- drop.tip(SL.tree, check$tree_not_data)
## Reactivity
SL_Reac_data <- SL_data[complete.cases(SL_data[,"Reactivity2"]),]; check <- name.check(SL.tree, SL_Reac_data); SL.Reac.tree <- drop.tip(SL.tree, check$tree_not_data)
## First-step Attenuation
SL_Att_data <- SL_data[complete.cases(SL_data[,"Attenuation2"]),]; check <- name.check(SL.tree, SL_Att_data)
## Maximal Amplification
SL_Amp_data <- SL_data[complete.cases(SL_data[,"MaxAmplification2"]),]; check <- name.check(SL.tree, SL_Amp_data); SL.Amp.tree <- drop.tip(SL.tree, check$tree_not_data)
## Maximal Attenuation
SL_MAtt_data <- SL_data[complete.cases(SL_data[,"MaxAttenuation2"]),]; check <- name.check(SL.tree, SL_MAtt_data)

### Long-lived populations
## Damping ratio
LL_DR_data <- LL_data[complete.cases(LL_data[,"Rho2"]),]; check <- name.check(LL.tree, LL_DR_data); LL.DR.tree <- drop.tip(LL.tree, check$tree_not_data)
## Period of Oscilation
LL_PO_data <- LL_data[complete.cases(LL_data[,"Pi2"]),]; check <- name.check(LL.tree, LL_PO_data); LL.PO.tree <- drop.tip(LL.tree, check$tree_not_data)
## Reactivity
LL_Reac_data <- LL_data[complete.cases(LL_data[,"Reactivity2"]),]; check <- name.check(LL.tree, LL_Reac_data)
## First-step Attenuation
LL_Att_data <- LL_data[complete.cases(LL_data[,"Attenuation2"]),]; check <- name.check(LL.tree, LL_Att_data)
## Maximal Amplification
LL_Amp_data <- LL_data[complete.cases(LL_data[,"MaxAmplification2"]),]; check <- name.check(LL.tree, LL_Amp_data)
## Maximal Attenuation
LL_MAtt_data <- LL_data[complete.cases(LL_data[,"MaxAttenuation2"]),]; check <- name.check(LL.tree, LL_MAtt_data)

#################################################
# STEP 3: Run analyses (only interested in the influence of environmental variables)
#################################################

############ Short-lived species ---------------------------------------------------

### Damping Ratio -------------------------
# define the variable sets for the analyses
# none-imputed dataset
DR.Y <- as.matrix(SL_DR_data[,1]); rownames(DR.Y) <- rownames(SL_DR_data) 
DR.enviro1 <- as.matrix(SL_DR_data[,7:11])
# imputed dataset
DR.Y2 <- as.matrix(SL_data[,14]); rownames(DR.Y2) <- rownames(SL_data)
DR.enviro2 <- as.matrix(SL_data[,7:11])
# ensure formating
rownames(DR.enviro1) <- rownames(DR.Y)
rownames(DR.enviro2) <- rownames(DR.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Damping.ratio.5 <- Phylo.PLS(DR.Y, DR.enviro1, tree = SL.DR.tree, marine.list)
# Imputed data
Damping.ratio.6 <- Phylo.PLS(DR.Y2, DR.enviro2, tree = SL.tree, marine.list)

### Period of Oscilation ------------------
# none-imputed dataset
PO.Y <- as.matrix(SL_PO_data[,2]); rownames(PO.Y) <- rownames(SL_PO_data) 
PO.enviro1 <- as.matrix(SL_PO_data[,7:11])
# imputed dataset
PO.Y2 <- as.matrix(SL_data[,13]); rownames(PO.Y2) <- rownames(SL_data)
PO.enviro2 <- as.matrix(SL_data[,7:11])
# ensure formating
rownames(PO.enviro1) <- rownames(PO.Y)
rownames(PO.enviro2) <- rownames(PO.Y2)

# Run the pls accounting for phylogeny.
# Normal data
PO.5 <- Phylo.PLS(PO.Y, PO.enviro1, tree = SL.PO.tree, marine.list)
# Imputed data
PO.6 <- Phylo.PLS(PO.Y2, PO.enviro2, tree = SL.tree, marine.list)

### Reactivity ----------------------------
# none-imputed dataset
Reac.Y <- as.matrix(SL_Reac_data[,3]); rownames(Reac.Y) <- rownames(SL_Reac_data) 
Reac.enviro1 <- as.matrix(SL_Reac_data[,7:11])
# imputed dataset
Reac.Y2 <- as.matrix(SL_data[,15]); rownames(Reac.Y2) <- rownames(SL_data)
Reac.enviro2 <- as.matrix(SL_data[,7:11])
# ensure formating
rownames(Reac.enviro1) <- rownames(Reac.Y)
rownames(Reac.enviro2) <- rownames(Reac.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Reac.5 <- Phylo.PLS(Reac.Y, Reac.enviro1, tree = SL.Reac.tree, marine.list)
# Imputed data
Reac.6 <- Phylo.PLS(Reac.Y2, Reac.enviro2, tree = SL.tree, marine.list)

### Attenuation ----------------------------
# none-imputed dataset
Att.Y <- as.matrix(SL_Att_data[,4]); rownames(Att.Y) <- rownames(SL_Att_data) 
Att.enviro1 <- as.matrix(SL_Att_data[,7:11])
# imputed dataset
Att.Y2 <- as.matrix(SL_data[,16]); rownames(Att.Y2) <- rownames(SL_data)
Att.enviro2 <- as.matrix(SL_data[,7:11])
# ensure formating
rownames(Att.enviro1) <- rownames(Att.Y)
rownames(Att.enviro2) <- rownames(Att.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Att.5 <- Phylo.PLS(Att.Y, Att.enviro1, tree = SL.tree, marine.list)
# Imputed data
Att.6 <- Phylo.PLS(Att.Y2, Att.enviro2, tree = SL.tree, marine.list)

### Maximal Amplification ----------------------------
# none-imputed dataset
Amp.Y <- as.matrix(SL_Amp_data[,5]); rownames(Amp.Y) <- rownames(SL_Amp_data) 
Amp.enviro1 <- as.matrix(SL_Amp_data[,7:11])
# imputed dataset
Amp.Y2 <- as.matrix(SL_data[,17]); rownames(Amp.Y2) <- rownames(SL_data)
Amp.enviro2 <- as.matrix(SL_data[,7:11])
# ensure formating
rownames(Amp.enviro1) <- rownames(Amp.Y)
rownames(Amp.enviro2) <- rownames(Amp.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Amp.5 <- Phylo.PLS(Amp.Y, Amp.enviro1, tree = SL.Amp.tree, marine.list)
# Imputed data
Amp.6 <- Phylo.PLS(Amp.Y2, Amp.enviro2, tree = SL.tree, marine.list)

### Maximal Attenuation ----------------------------
# none-imputed dataset
MAtt.Y <- as.matrix(SL_MAtt_data[,6]); rownames(MAtt.Y) <- rownames(SL_MAtt_data) 
MAtt.enviro1 <- as.matrix(SL_MAtt_data[,7:11])
# imputed dataset
MAtt.Y2 <- as.matrix(SL_data[,18]); rownames(MAtt.Y2) <- rownames(SL_data)
MAtt.enviro2 <- as.matrix(SL_data[,7:11])
# ensure formating
rownames(MAtt.enviro1) <- rownames(MAtt.Y)
rownames(MAtt.enviro2) <- rownames(MAtt.Y2)

# Run the pls accounting for phylogeny.
# Normal data
MAtt.5 <- Phylo.PLS(MAtt.Y, MAtt.enviro1, tree = SL.tree, marine.list)
# Imputed data
MAtt.6 <- Phylo.PLS(MAtt.Y2, MAtt.enviro2, tree = SL.tree, marine.list)

############ Long-lived species ---------------------------------------------------

### Damping Ratio -------------------------
# define the variable sets for the analyses
# none-imputed dataset
DR.Y <- as.matrix(LL_DR_data[,1]); rownames(DR.Y) <- rownames(LL_DR_data) 
DR.enviro1 <- as.matrix(LL_DR_data[,7:11])
# imputed dataset
DR.Y2 <- as.matrix(LL_data[,14]); rownames(DR.Y2) <- rownames(LL_data)
DR.enviro2 <- as.matrix(LL_data[,7:11])
# ensure formating
rownames(DR.enviro1) <- rownames(DR.Y)
rownames(DR.enviro2) <- rownames(DR.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Damping.ratio.3 <- Phylo.PLS(DR.Y, DR.enviro1, tree = LL.DR.tree, marine.list)
# Imputed data
Damping.ratio.4 <- Phylo.PLS(DR.Y2, DR.enviro2, tree = LL.tree, marine.list)

### Period of Oscilation ------------------
# none-imputed dataset
PO.Y <- as.matrix(LL_PO_data[,2]); rownames(PO.Y) <- rownames(LL_PO_data) 
PO.enviro1 <- as.matrix(LL_PO_data[,7:11])
# imputed dataset
PO.Y2 <- as.matrix(LL_data[,13]); rownames(PO.Y2) <- rownames(LL_data)
PO.enviro2 <- as.matrix(LL_data[,7:11])
# ensure formating
rownames(PO.enviro1) <- rownames(PO.Y)
rownames(PO.enviro2) <- rownames(PO.Y2)

# Run the pls accounting for phylogeny.
# Normal data
PO.3 <- Phylo.PLS(PO.Y, PO.enviro1, tree = LL.PO.tree, marine.list)
# Imputed data
PO.4 <- Phylo.PLS(PO.Y2, PO.enviro2, tree = LL.tree, marine.list)

### Reactivity ----------------------------
# none-imputed dataset
Reac.Y <- as.matrix(LL_Reac_data[,3]); rownames(Reac.Y) <- rownames(LL_Reac_data) 
Reac.enviro1 <- as.matrix(LL_Reac_data[,7:11])
# imputed dataset
Reac.Y2 <- as.matrix(LL_data[,15]); rownames(Reac.Y2) <- rownames(LL_data)
Reac.enviro2 <- as.matrix(LL_data[,7:11])
# ensure formating
rownames(Reac.enviro1) <- rownames(Reac.Y)
rownames(Reac.enviro2) <- rownames(Reac.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Reac.3 <- Phylo.PLS(Reac.Y, Reac.enviro1, tree = LL.Reac.tree, marine.list)
# Imputed data
Reac.4 <- Phylo.PLS(Reac.Y2, Reac.enviro2, tree = LL.tree, marine.list)

### Attenuation ----------------------------
# none-imputed dataset
Att.Y <- as.matrix(LL_Att_data[,4]); rownames(Att.Y) <- rownames(LL_Att_data) 
Att.enviro1 <- as.matrix(LL_Att_data[,7:11])
# imputed dataset
Att.Y2 <- as.matrix(LL_data[,16]); rownames(Att.Y2) <- rownames(LL_data)
Att.enviro2 <- as.matrix(LL_data[,7:11])
# ensure formating
rownames(Att.enviro1) <- rownames(Att.Y)
rownames(Att.enviro2) <- rownames(Att.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Att.3 <- Phylo.PLS(Att.Y, Att.enviro1, tree = LL.tree, marine.list)
# Imputed data
Att.4 <- Phylo.PLS(Att.Y2, Att.enviro2, tree = LL.tree, marine.list)

### Maximal Amplification ----------------------------
# none-imputed dataset
Amp.Y <- as.matrix(LL_Amp_data[,5]); rownames(Amp.Y) <- rownames(LL_Amp_data) 
Amp.enviro1 <- as.matrix(LL_Amp_data[,7:11])
# imputed dataset
Amp.Y2 <- as.matrix(LL_data[,17]); rownames(Amp.Y2) <- rownames(LL_data)
Amp.enviro2 <- as.matrix(LL_data[,7:11])
# ensure formating
rownames(Amp.enviro1) <- rownames(Amp.Y)
rownames(Amp.enviro2) <- rownames(Amp.Y2)

# Run the pls accounting for phylogeny.
# Normal data
Amp.3 <- Phylo.PLS(Amp.Y, Amp.enviro1, tree = LL.Amp.tree, marine.list)
# Imputed data
Amp.4 <- Phylo.PLS(Amp.Y2, Amp.enviro2, tree = LL.tree, marine.list)

### Maximal Attenuation ----------------------------
# none-imputed dataset
MAtt.Y <- as.matrix(LL_MAtt_data[,6]); rownames(MAtt.Y) <- rownames(LL_MAtt_data) 
MAtt.enviro1 <- as.matrix(LL_MAtt_data[,7:11])
# imputed dataset
MAtt.Y2 <- as.matrix(LL_data[,18]); rownames(MAtt.Y2) <- rownames(LL_data)
MAtt.enviro2 <- as.matrix(LL_data[,7:11])
# ensure formating
rownames(MAtt.enviro1) <- rownames(MAtt.Y)
rownames(MAtt.enviro2) <- rownames(MAtt.Y2)

# Run the pls accounting for phylogeny.
# Normal data
MAtt.3 <- Phylo.PLS(MAtt.Y, MAtt.enviro1, tree = LL.tree, marine.list)
# Imputed data
MAtt.4 <- Phylo.PLS(MAtt.Y2, MAtt.enviro2, tree = LL.tree, marine.list)

#################################################
# STEP 4: Plot outcomes
#################################################
# Given that it appears there is no relationship I am only interested in plotting the perfect datasets (this matches the data being used in the main manuscripts)

##### Short-lived populations ---------------------------------------------
#### Damping ratio --------------------
# extract key data
DR_1 <- data.frame(cbind(Damping.ratio.5[["Hypothesis 2 pls"]]$scores[,1],
                         Damping.ratio.5[["Hypothesis 2 pls"]]$scores[,2],
                         Damping.ratio.5[["Hypothesis 2 pls"]]$scores[,3],
                         Damping.ratio.5[["Hypothesis 2 pls"]]$scores[,4],
                         Damping.ratio.5[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         Damping.ratio.5[["Transient variable"]]))
colnames(DR_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Rho") # puts the detail back in

# store PCA loadings.
DR_loadings_1 <- data.frame(cbind(Damping.ratio.5[["Hypothesis 2 pls"]]$loadings[,1],
                                  Damping.ratio.5[["Hypothesis 2 pls"]]$loadings[,2],
                                  Damping.ratio.5[["Hypothesis 2 pls"]]$loadings[,3],
                                  Damping.ratio.5[["Hypothesis 2 pls"]]$loadings[,4],
                                  Damping.ratio.5[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(DR_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.DR1 <- explvar(Damping.ratio.5[["Hypothesis 2 pls"]]) # percentage variance explained
r2.DR1 <- R2(Damping.ratio.5[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
DR.coefficients1 = -coef(Damping.ratio.5[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.DR.coef1 = sum(sapply(DR.coefficients1, abs))
DR.coefficients1 = DR.coefficients1 / sum.DR.coef1
DR.coefficients1 = sort(DR.coefficients1[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# descriptive stats for the damping ratio variable
range(DR_1$Rho, na.rm = T); median(DR_1$Rho, na.rm = T); mean(DR_1$Rho, na.rm = T)

# create pls biplot
DR.SL.plot <- ggplot(data = DR_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.DR1[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.DR1[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Rho), alpha = 0.7, data = DR_1) + #add data points
  scale_color_gradientn(colors = pal(dim(DR_1)[1]),
                        values = scales::rescale(c(0.22,0.63,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = DR_loadings_1, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (DR_loadings_1$Comp1*4), y = (DR_loadings_1$Comp2*4),
  #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
DR.SL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
DR.coefs.bar1 <- ggplot() + 
  geom_bar(aes(x=names(DR.coefficients1), y=DR.coefficients1), stat="identity", fill = "gold") +
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
# display the plot
DR.coefs.bar1

#### Period of Oscillation --------------------
# extract key data
PO_1 <- data.frame(cbind(PO.5[["Hypothesis 2 pls"]]$scores[,1],
                         PO.5[["Hypothesis 2 pls"]]$scores[,2],
                         PO.5[["Hypothesis 2 pls"]]$scores[,3],
                         PO.5[["Hypothesis 2 pls"]]$scores[,4],
                         PO.5[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         PO.5[["Transient variable"]]))
colnames(PO_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Pi") # puts the detail back in

# store PCA loadings.
PO_loadings_1 <- data.frame(cbind(PO.5[["Hypothesis 2 pls"]]$loadings[,1],
                                  PO.5[["Hypothesis 2 pls"]]$loadings[,2],
                                  PO.5[["Hypothesis 2 pls"]]$loadings[,3],
                                  PO.5[["Hypothesis 2 pls"]]$loadings[,4],
                                  PO.5[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(PO_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.PO1 <- explvar(PO.5[["Hypothesis 2 pls"]]) # percentage variance explained
r2.PO1 <- R2(PO.5[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
PO.coefficients1 = -coef(PO.5[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.PO.coef1 = sum(sapply(PO.coefficients1, abs))
PO.coefficients1 = PO.coefficients1 / sum.PO.coef1
PO.coefficients1 = sort(PO.coefficients1[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("chocolate4", "gold", "lightgoldenrod1"))

# descriptive stats for the damping ratio variable
range(PO_1$Pi, na.rm = T); median(PO_1$Pi, na.rm = T); mean(PO_1$Pi, na.rm = T)

# create pls biplot
PO.SL.plot <- ggplot(data = PO_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.PO1[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.PO1[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Pi), alpha = 0.7, data = PO_1) + #add data points
  scale_color_gradientn(colors = pal(dim(PO_1)[1]),
                        values = scales::rescale(c(0.06,0.27,0.44), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = PO_loadings_1, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (PO_loadings_1$Comp1*4), y = (PO_loadings_1$Comp2*4),
           #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
PO.SL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
PO.coefs.bar1 <- ggplot() + 
  geom_bar(aes(x=names(PO.coefficients1), y=PO.coefficients1), stat="identity", fill = "gold") +
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
# display the plot
PO.coefs.bar1

#### Reactivity --------------------
# extract key data
Reac_1 <- data.frame(cbind(Reac.5[["Hypothesis 2 pls"]]$scores[,1],
                         Reac.5[["Hypothesis 2 pls"]]$scores[,2],
                         Reac.5[["Hypothesis 2 pls"]]$scores[,3],
                         Reac.5[["Hypothesis 2 pls"]]$scores[,4],
                         Reac.5[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         Reac.5[["Transient variable"]]))
colnames(Reac_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Reactivity") # puts the detail back in

# store PCA loadings.
Reac_loadings_1 <- data.frame(cbind(Reac.5[["Hypothesis 2 pls"]]$loadings[,1],
                                  Reac.5[["Hypothesis 2 pls"]]$loadings[,2],
                                  Reac.5[["Hypothesis 2 pls"]]$loadings[,3],
                                  Reac.5[["Hypothesis 2 pls"]]$loadings[,4],
                                  Reac.5[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Reac_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.R1 <- explvar(Reac.5[["Hypothesis 2 pls"]]) # percentage variance explained
r2.R1 <- R2(Reac.5[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
R.coefficients1 = -coef(Reac.5[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.R.coef1 = sum(sapply(R.coefficients1, abs))
R.coefficients1 = R.coefficients1 / sum.R.coef1
R.coefficients1 = sort(R.coefficients1[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(Reac_1$Reactivity, na.rm = T); median(Reac_1$Reactivity, na.rm = T); mean(Reac_1$Reactivity, na.rm = T)

# create pls biplot
R.SL.plot <- ggplot(data = Reac_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.R1[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.R1[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Reactivity), alpha = 0.7, data = Reac_1) + #add data points
  scale_color_gradientn(colors = pal(dim(Reac_1)[1]),
                        values = scales::rescale(c(0.03,0.55,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = Reac_loadings_1, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (Reac_loadings_1$Comp1*4), y = (Reac_loadings_1$Comp2*4),
          # label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
R.SL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
R.coefs.bar1 <- ggplot() + 
  geom_bar(aes(x=names(R.coefficients1), y=R.coefficients1), stat="identity", fill = "navyblue") +
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
# display the plot
R.coefs.bar1

#### Maximal Amplification --------------------
# extract key data
Amp_1 <- data.frame(cbind(Amp.5[["Hypothesis 2 pls"]]$scores[,1],
                           Amp.5[["Hypothesis 2 pls"]]$scores[,2],
                           Amp.5[["Hypothesis 2 pls"]]$scores[,3],
                           Amp.5[["Hypothesis 2 pls"]]$scores[,4],
                           Amp.5[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Amp.5[["Transient variable"]]))
colnames(Amp_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Amp") # puts the detail back in

# store PCA loadings.
Amp_loadings_1 <- data.frame(cbind(Amp.5[["Hypothesis 2 pls"]]$loadings[,1],
                                    Amp.5[["Hypothesis 2 pls"]]$loadings[,2],
                                    Amp.5[["Hypothesis 2 pls"]]$loadings[,3],
                                    Amp.5[["Hypothesis 2 pls"]]$loadings[,4],
                                    Amp.5[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Amp_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Amp1 <- explvar(Amp.5[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Amp1 <- R2(Amp.5[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Amp.coefficients1 = -coef(Amp.5[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Amp.coef1 = sum(sapply(Amp.coefficients1, abs))
Amp.coefficients1 = Amp.coefficients1 / sum.Amp.coef1
Amp.coefficients1 = sort(Amp.coefficients1[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(Amp_1$Amp, na.rm = T); median(Amp_1$Amp, na.rm = T); mean(Amp_1$Amp, na.rm = T)

# create pls biplot
Amp.SL.plot <- ggplot(data = Amp_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Amp1[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Amp1[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Amp), alpha = 0.7, data = Amp_1) + #add data points
  scale_color_gradientn(colors = pal(dim(Amp_1)[1]),
                        values = scales::rescale(c(0.02,0.53,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = Amp_loadings_1, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (Amp_loadings_1$Comp1*4), y = (Amp_loadings_1$Comp2*4),
           #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
Amp.SL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
Amp.coefs.bar1 <- ggplot() + 
  geom_bar(aes(x=names(Amp.coefficients1), y=Amp.coefficients1), stat="identity", fill = "navyblue") +
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
# display the plot
Amp.coefs.bar1

#### First-step Attenuation --------------------
# extract key data
Att_1 <- data.frame(cbind(Att.5[["Hypothesis 2 pls"]]$scores[,1],
                           Att.5[["Hypothesis 2 pls"]]$scores[,2],
                           Att.5[["Hypothesis 2 pls"]]$scores[,3],
                           Att.5[["Hypothesis 2 pls"]]$scores[,4],
                           Att.5[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Att.1[["Transient variable"]]))
colnames(Att_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Att") # puts the detail back in

# store PCA loadings.
Att_loadings_1 <- data.frame(cbind(Att.5[["Hypothesis 2 pls"]]$loadings[,1],
                                    Att.5[["Hypothesis 2 pls"]]$loadings[,2],
                                    Att.5[["Hypothesis 2 pls"]]$loadings[,3],
                                    Att.5[["Hypothesis 2 pls"]]$loadings[,4],
                                    Att.5[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Att_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Att1 <- explvar(Att.5[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Att1 <- R2(Att.5[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Att.coefficients1 = coef(Att.5[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Att.coef1 = sum(sapply(Att.coefficients1, abs))
Att.coefficients1 = Att.coefficients1 / sum.Att.coef1
Att.coefficients1 = sort(Att.coefficients1[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(Att_1$Att, na.rm = T); median(Att_1$Att, na.rm = T); mean(Att_1$Att, na.rm = T)

# create pls biplot
Att.SL.plot <- ggplot(data = Att_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Att1[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Att1[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Att), alpha = 0.7, data = Att_1) + #add data points
  scale_color_gradientn(colors = pal(dim(Att_1)[1]),
                        values = scales::rescale(c(0.0,0.51,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = Att_loadings_1, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (Att_loadings_1$Comp1*4), y = (Att_loadings_1$Comp2*4),
          # label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
Att.SL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
Att.coefs.bar1 <- ggplot() + 
  geom_bar(aes(x=names(Att.coefficients1), y=Att.coefficients1), stat="identity", fill = "darkgreen") +
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
# display the plot
Att.coefs.bar1

#### Maximal Attenuation --------------------
# extract key data
MAtt_1 <- data.frame(cbind(MAtt.5[["Hypothesis 2 pls"]]$scores[,1],
                          MAtt.5[["Hypothesis 2 pls"]]$scores[,2],
                          MAtt.5[["Hypothesis 2 pls"]]$scores[,3],
                          MAtt.5[["Hypothesis 2 pls"]]$scores[,4],
                          MAtt.5[["Hypothesis 2 pls"]]$scores[,5],
                          #this accounts for the odd class used for storing pls score outputs
                          MAtt.5[["Transient variable"]]))
colnames(MAtt_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "MAtt") # puts the detail back in

# store PCA loadings.
MAtt_loadings_1 <- data.frame(cbind(MAtt.5[["Hypothesis 2 pls"]]$loadings[,1],
                                   MAtt.5[["Hypothesis 2 pls"]]$loadings[,2],
                                   MAtt.5[["Hypothesis 2 pls"]]$loadings[,3],
                                   MAtt.5[["Hypothesis 2 pls"]]$loadings[,4],
                                   MAtt.5[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(MAtt_loadings_1) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.MAtt1 <- explvar(MAtt.5[["Hypothesis 2 pls"]]) # percentage variance explained
r2.MAtt1 <- R2(MAtt.5[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
MAtt.coefficients1 = coef(MAtt.5[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.MAtt.coef1 = sum(sapply(MAtt.coefficients1, abs))
MAtt.coefficients1 = MAtt.coefficients1 / sum.MAtt.coef1
MAtt.coefficients1 = sort(MAtt.coefficients1[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(MAtt_1$MAtt, na.rm = T); median(MAtt_1$MAtt, na.rm = T); mean(MAtt_1$MAtt, na.rm = T)

# create pls biplot
MAtt.SL.plot <- ggplot(data = MAtt_1, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.MAtt1[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.MAtt1[1], 1),"%)")) +
  geom_point(size = 3, aes(col = MAtt), alpha = 0.7, data = MAtt_1) + #add data points
  scale_color_gradientn(colors = pal(dim(MAtt_1)[1]),
                        values = scales::rescale(c(0.01,0.50,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = MAtt_loadings_1, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(1, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (MAtt_loadings_1$Comp1*4), y = (MAtt_loadings_1$Comp2*4),
           #label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
MAtt.SL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
MAtt.coefs.bar1 <- ggplot() + 
  geom_bar(aes(x=names(MAtt.coefficients1), y=MAtt.coefficients1), stat="identity", fill = "darkgreen") +
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
# display the plot
MAtt.coefs.bar1

##### Long-lived populations ---------------------------------------------
#### Damping ratio --------------------
# extract key data
DR_2 <- data.frame(cbind(Damping.ratio.3[["Hypothesis 2 pls"]]$scores[,1],
                         Damping.ratio.3[["Hypothesis 2 pls"]]$scores[,2],
                         Damping.ratio.3[["Hypothesis 2 pls"]]$scores[,3],
                         Damping.ratio.3[["Hypothesis 2 pls"]]$scores[,4],
                         Damping.ratio.3[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         Damping.ratio.3[["Transient variable"]]))
colnames(DR_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Rho") # puts the detail back in

# store PCA loadings.
DR_loadings_2 <- data.frame(cbind(Damping.ratio.3[["Hypothesis 2 pls"]]$loadings[,1],
                                  Damping.ratio.3[["Hypothesis 2 pls"]]$loadings[,2],
                                  Damping.ratio.3[["Hypothesis 2 pls"]]$loadings[,3],
                                  Damping.ratio.3[["Hypothesis 2 pls"]]$loadings[,4],
                                  Damping.ratio.3[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(DR_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.DR2 <- explvar(Damping.ratio.3[["Hypothesis 2 pls"]]) # percentage variance explained
r2.DR2 <- R2(Damping.ratio.3[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
DR.coefficients2 = -coef(Damping.ratio.3[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.DR.coef2 = sum(sapply(DR.coefficients2, abs))
DR.coefficients2 = DR.coefficients2 / sum.DR.coef2
DR.coefficients2 = sort(DR.coefficients2[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("dodgerblue3", "lightskyblue", "aliceblue"))

# descriptive stats for the damping ratio variable
range(DR_2$Rho, na.rm = T); median(DR_2$Rho, na.rm = T); mean(DR_2$Rho, na.rm = T)

# create pls biplot
DR.LL.plot <- ggplot(data = DR_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.DR2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.DR2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Rho), data = DR_2) + #add data points
  scale_color_gradientn(colors = pal(dim(DR_2)[1]),
                        values = scales::rescale(c(0.18,0.58,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = DR_loadings_2, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(0.5, "picas"), type = "open"),
               color = "black", size = 1) +
  annotate("text", x = (DR_loadings_2$Comp1*4), y = (DR_loadings_2$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
DR.LL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
DR.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(DR.coefficients2), y=DR.coefficients2), stat="identity", fill = "lightskyblue") +
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
# display the plot
DR.coefs.bar2

#### Period of Oscillation --------------------
# extract key data
PO_2 <- data.frame(cbind(PO.3[["Hypothesis 2 pls"]]$scores[,1],
                         PO.3[["Hypothesis 2 pls"]]$scores[,2],
                         PO.3[["Hypothesis 2 pls"]]$scores[,3],
                         PO.3[["Hypothesis 2 pls"]]$scores[,4],
                         PO.3[["Hypothesis 2 pls"]]$scores[,5],
                         #this accounts for the odd class used for storing pls score outputs
                         PO.3[["Transient variable"]]))
colnames(PO_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Pi") # puts the detail back in

# store PCA loadings.
PO_loadings_2 <- data.frame(cbind(PO.3[["Hypothesis 2 pls"]]$loadings[,1],
                                  PO.3[["Hypothesis 2 pls"]]$loadings[,2],
                                  PO.3[["Hypothesis 2 pls"]]$loadings[,3],
                                  PO.3[["Hypothesis 2 pls"]]$loadings[,4],
                                  PO.3[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(PO_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.PO2 <- explvar(PO.3[["Hypothesis 2 pls"]]) # percentage variance explained
r2.PO2 <- R2(PO.3[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
PO.coefficients2 = -coef(PO.3[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.PO.coef2 = sum(sapply(PO.coefficients2, abs))
PO.coefficients2 = PO.coefficients2 / sum.PO.coef2
PO.coefficients2 = sort(PO.coefficients2[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("dodgerblue3", "lightskyblue", "aliceblue"))

# descriptive stats for the damping ratio variable
range(PO_2$Pi, na.rm = T); median(PO_2$Pi, na.rm = T); mean(PO_2$Pi, na.rm = T)

# create pls biplot
PO.LL.plot <- ggplot(data = PO_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.PO2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.PO2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Pi), data = PO_2) + #add data points
  scale_color_gradientn(colors = pal(dim(PO_2)[1]),
                        values = scales::rescale(c(0.06,0.27,0.44), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = PO_loadings_2, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(0.5, "picas"), type = "open"),
               color = "black", size = 1) +
  annotate("text", x = (PO_loadings_2$Comp1*4), y = (PO_loadings_2$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
PO.LL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
PO.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(PO.coefficients2), y=PO.coefficients2), stat="identity", fill = "lightskyblue") +
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
# display the plot
PO.coefs.bar2

#### Reactivity --------------------
# extract key data
Reac_2 <- data.frame(cbind(Reac.3[["Hypothesis 2 pls"]]$scores[,1],
                           Reac.3[["Hypothesis 2 pls"]]$scores[,2],
                           Reac.3[["Hypothesis 2 pls"]]$scores[,3],
                           Reac.3[["Hypothesis 2 pls"]]$scores[,4],
                           Reac.3[["Hypothesis 2 pls"]]$scores[,5],
                           #this accounts for the odd class used for storing pls score outputs
                           Reac.3[["Transient variable"]]))
colnames(Reac_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Reactivity") # puts the detail back in

# store PCA loadings.
Reac_loadings_2 <- data.frame(cbind(Reac.3[["Hypothesis 2 pls"]]$loadings[,1],
                                    Reac.3[["Hypothesis 2 pls"]]$loadings[,2],
                                    Reac.3[["Hypothesis 2 pls"]]$loadings[,3],
                                    Reac.3[["Hypothesis 2 pls"]]$loadings[,4],
                                    Reac.3[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Reac_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.R2 <- explvar(Reac.3[["Hypothesis 2 pls"]]) # percentage variance explained
r2.R2 <- R2(Reac.3[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
R.coefficients2 = -coef(Reac.3[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.R.coef2 = sum(sapply(R.coefficients2, abs))
R.coefficients2 = R.coefficients2 / sum.R.coef2
R.coefficients2 = sort(R.coefficients2[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(Reac_2$Reactivity, na.rm = T); median(Reac_2$Reactivity, na.rm = T); mean(Reac_2$Reactivity, na.rm = T)

# create pls biplot
R.LL.plot <- ggplot(data = Reac_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.R2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.R2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Reactivity), data = Reac_2) + #add data points
  scale_color_gradientn(colors = pal(dim(Reac_2)[1]),
                        values = scales::rescale(c(0.03,0.55,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = Reac_loadings_2, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(0.5, "picas"), type = "open"),
               color = "black", size = 1) +
  annotate("text", x = (Reac_loadings_2$Comp1*4), y = (Reac_loadings_2$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
R.LL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
R.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(R.coefficients2), y=R.coefficients2), stat="identity", fill = "navyblue") +
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
# display the plot
R.coefs.bar2

#### Maximal Amplification --------------------
# extract key data
Amp_2 <- data.frame(cbind(Amp.3[["Hypothesis 2 pls"]]$scores[,1],
                          Amp.3[["Hypothesis 2 pls"]]$scores[,2],
                          Amp.3[["Hypothesis 2 pls"]]$scores[,3],
                          Amp.3[["Hypothesis 2 pls"]]$scores[,4],
                          Amp.3[["Hypothesis 2 pls"]]$scores[,5],
                          #this accounts for the odd class used for storing pls score outputs
                          Amp.3[["Transient variable"]]))
colnames(Amp_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Amp") # puts the detail back in

# store PCA loadings.
Amp_loadings_2 <- data.frame(cbind(Amp.3[["Hypothesis 2 pls"]]$loadings[,1],
                                   Amp.3[["Hypothesis 2 pls"]]$loadings[,2],
                                   Amp.3[["Hypothesis 2 pls"]]$loadings[,3],
                                   Amp.3[["Hypothesis 2 pls"]]$loadings[,4],
                                   Amp.3[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Amp_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Amp2 <- explvar(Amp.3[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Amp2 <- R2(Amp.3[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Amp.coefficients2 = -coef(Amp.3[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Amp.coef2 = sum(sapply(Amp.coefficients2, abs))
Amp.coefficients2 = Amp.coefficients2 / sum.Amp.coef2
Amp.coefficients2 = sort(Amp.coefficients2[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("navyblue", "dodgerblue", "lightblue"))

# descriptive stats for the damping ratio variable
range(Amp_2$Amp, na.rm = T); median(Amp_2$Amp, na.rm = T); mean(Amp_2$Amp, na.rm = T)

# create pls biplot
Amp.LL.plot <- ggplot(data = Amp_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Amp2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Amp2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Amp), data = Amp_2) + #add data points
  scale_color_gradientn(colors = pal(dim(Amp_2)[1]),
                        values = scales::rescale(c(0.02,0.59,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = T, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = Amp_loadings_2, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(0.5, "picas"), type = "open"),
               color = "black", size = 1) +
  annotate("text", x = (Amp_loadings_2$Comp1*4), y = (Amp_loadings_2$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
Amp.LL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
Amp.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(Amp.coefficients2), y=Amp.coefficients2), stat="identity", fill = "navyblue") +
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
# display the plot
Amp.coefs.bar2

#### First-step Attenuation --------------------
# extract key data
Att_2 <- data.frame(cbind(Att.3[["Hypothesis 2 pls"]]$scores[,1],
                          Att.3[["Hypothesis 2 pls"]]$scores[,2],
                          Att.3[["Hypothesis 2 pls"]]$scores[,3],
                          Att.3[["Hypothesis 2 pls"]]$scores[,4],
                          Att.3[["Hypothesis 2 pls"]]$scores[,5],
                          #this accounts for the odd class used for storing pls score outputs
                          Att.3[["Transient variable"]]))
colnames(Att_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Att") # puts the detail back in

# store PCA loadings.
Att_loadings_2 <- data.frame(cbind(Att.3[["Hypothesis 2 pls"]]$loadings[,1],
                                   Att.3[["Hypothesis 2 pls"]]$loadings[,2],
                                   Att.3[["Hypothesis 2 pls"]]$loadings[,3],
                                   Att.3[["Hypothesis 2 pls"]]$loadings[,4],
                                   Att.3[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(Att_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.Att2 <- explvar(Att.3[["Hypothesis 2 pls"]]) # percentage variance explained
r2.Att2 <- R2(Att.3[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
Att.coefficients2 = coef(Att.3[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.Att.coef2 = sum(sapply(Att.coefficients2, abs))
Att.coefficients2 = Att.coefficients2 / sum.Att.coef2
Att.coefficients2 = sort(Att.coefficients2[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(Att_2$Att, na.rm = T); median(Att_2$Att, na.rm = T); mean(Att_2$Att, na.rm = T)

# create pls biplot
Att.LL.plot <- ggplot(data = Att_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.Att2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.Att2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Att), data = Att_2) + #add data points
  scale_color_gradientn(colors = pal(dim(Att_2)[1]),
                        values = scales::rescale(c(0.0,0.56,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = Att_loadings_2, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(0.5, "picas"), type = "open"),
               color = "black", size = 1) +
  annotate("text", x = (Att_loadings_2$Comp1*4), y = (Att_loadings_2$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
Att.LL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
Att.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(Att.coefficients2), y=Att.coefficients2), stat="identity", fill = "darkgreen") +
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
# display the plot
Att.coefs.bar2

#### Maximal Attenuation --------------------
# extract key data
MAtt_2 <- data.frame(cbind(MAtt.3[["Hypothesis 2 pls"]]$scores[,1],
                          MAtt.3[["Hypothesis 2 pls"]]$scores[,2],
                          MAtt.3[["Hypothesis 2 pls"]]$scores[,3],
                          MAtt.3[["Hypothesis 2 pls"]]$scores[,4],
                          MAtt.3[["Hypothesis 2 pls"]]$scores[,5],
                          #this accounts for the odd class used for storing pls score outputs
                          MAtt.3[["Transient variable"]]))
colnames(MAtt_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Att") # puts the detail back in

# store PCA loadings.
MAtt_loadings_2 <- data.frame(cbind(MAtt.3[["Hypothesis 2 pls"]]$loadings[,1],
                                   MAtt.3[["Hypothesis 2 pls"]]$loadings[,2],
                                   MAtt.3[["Hypothesis 2 pls"]]$loadings[,3],
                                   MAtt.3[["Hypothesis 2 pls"]]$loadings[,4],
                                   MAtt.3[["Hypothesis 2 pls"]]$loadings[,5]))
colnames(MAtt_loadings_2) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5")

# extract the variance explained by each component
p.var.MAtt2 <- explvar(MAtt.3[["Hypothesis 2 pls"]]) # percentage variance explained
r2.MAtt2 <- R2(MAtt.3[["Hypothesis 2 pls"]]) # relative effect on R squared.

# Extract coefficients (accounting for inverse transformations where nessecary)
MAtt.coefficients2 = coef(MAtt.3[["Hypothesis 2 pls"]])
# standardise coefficients so that they sum to 1
sum.MAtt.coef2 = sum(sapply(MAtt.coefficients2, abs))
MAtt.coefficients2 = MAtt.coefficients2 / sum.MAtt.coef2
MAtt.coefficients2 = sort(MAtt.coefficients2[, 1 , 1])

# define colour scheme for overlaying transient measure over the pls weightings plot
pal = colorRampPalette(c("khaki","darkolivegreen3","darkgreen"))

# descriptive stats for the damping ratio variable
range(MAtt_2$Att, na.rm = T); median(MAtt_2$Att, na.rm = T); mean(MAtt_2$Att, na.rm = T)

# create pls biplot
MAtt.LL.plot <- ggplot(data = MAtt_2, aes(x = Comp1, y = Comp2)) + #base plot
  labs(y = paste0("Component 2 (", round(p.var.MAtt2[2], 1),"%)"), # add axis labels
       x = paste0("Component 1 (", round(p.var.MAtt2[1], 1),"%)")) +
  geom_point(size = 3, aes(col = Att), data = MAtt_2) + #add data points
  scale_color_gradientn(colors = pal(dim(MAtt_2)[1]),
                        values = scales::rescale(c(0.0,0.56,1), to = c(0,1)), limits = c(0,1), #this sets the value limits of the colour scale
                        guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                               direction = "horizontal", reverse = F, label = F,
                                               barwidth = 15, barheight = 1), na.value = "white") +
  geom_segment(data = MAtt_loadings_2, aes(x = 0, y = 0, xend = (Comp1*4), yend = (Comp2*4)), arrow = arrow(length = unit(0.5, "picas"), type = "open"),
               color = "black", size = 1) +
  annotate("text", x = (MAtt_loadings_2$Comp1*4), y = (MAtt_loadings_2$Comp2*4),
           label = c("Temp Frequency", "Temp Autocorrelation", "Thermal Range", "Prec Frequency", "Prec Autocorrelation")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom",
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), plot.margin = unit(c(2,3,1,1), "cm"),
        legend.text = element_text(), legend.title = element_blank())
# display the plot
MAtt.LL.plot

# Create coefficients barplot
level_order <- c("Temp.Frequency.full", "Temp.Autocorrelation.full", "Temp.range.full", "Prec.Frequency.full", "Prec.Autocorrelation.full")
MAtt.coefs.bar2 <- ggplot() + 
  geom_bar(aes(x=names(MAtt.coefficients2), y=MAtt.coefficients2), stat="identity", fill = "darkgreen") +
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
# display the plot
MAtt.coefs.bar2

# So importantly the speed of life has no effect.
# ******************************************************************** End of Code *******************************************************