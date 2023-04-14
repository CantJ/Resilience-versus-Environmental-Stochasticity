# ---------------------------------------------------
# This function is for carrying out a phylogenetically corrected partial least squares analysis.
# Canonical means to be ordered by some pre defined weighting and so allows the PLS analysis to be adjusted for phylogeny
# This function is an adaptation arising from code provided in: 
#                       1. D.C. Adams & R.N. Felice (2014) Plos One 9(4): e94335. doi: 10.1371/journal.pone.0094335
#                       2. The phytools package (https://github.com/liamrevell/phytools/tree/master/R)
#                       3. The cppls function in the pls package 
# Date Last modified: 15th June 2020
# Author: James Cant
# --------------------------------------------------

# This function carries out two partial squares regression analysis to investigate
# the selection pressures on transient dynamics
# and the correlation between environmental measures and these selection pressures.

Phylo.PLS <-function(BlockY, BlockX1, BlockX2, tree, marine.list){ 
  # this function takes three variable sets and the associated phylogenetic tree.
  
  # BlockY is a vector/matrix of numerical values describing a transient property 
  # measured from various populations (i.e. Damping ratio, Period of Oscillation, 
  # Reactivity, Maximal Amplification, First-step Attenuation, Maximal Attenuation).
  
  # BlockX1 is a matrix describing the population specific sensitivities of the 
  # transient parameter (BlockY) to each of the vital rates Survival, Progression, 
  # Retrogression and Reproduction.
  
  # BlockX2 is a matrix describing the environmental variability to which each 
  # population is exposed. This variable set will consist of two abiotic variables: 
  # Variance frequency and thermal range (variance magnitude).   
  
  # BlockY and BlockX1 will be used in the analysis of Hypothesis 1: Evaluating 
  # whether the components of resilience align with the fast-slow continuum of 
  # life history strategies, with the expectation that a trade-off exists 
  # between recovery and stability (resistance or resonance), with stability 
  # increased by greater investments in survival, and recovery correlating with 
  # reproductive investment. 
  
  # BlockY and BlockX2 will be used in the analysis of Hypothesis 2: Investigating 
  # the role played by the intensity and frequency of environmental variability in 
  # defining a population's resilience charatersitics, with exposure to broader 
  # scales of abiotic variance selecting for enhanced resistance/resonance potential, 
  # whereas environments characterised by higher frequency variation 
  # would promote recovery capcities.
  
  # The function uses the phylogenetic tree to ensure each analysis is weighted 
  # by phylogenetic covariance.
  
  # The marine.list is then a list of species (in the same format as the row names of the Block matrices) that are marine.
  # This list will be used to remove marine populations from the analysis of hypothesis 2.
  
  # 1. Function loads required packages
  library(ape) # phylogenetic functions
  library(phytools) # phylogenetic functions
  library(pls) # partial least squares analysis
  
  # 2. Confirm that the input phylogenetic tree is of class 'phylo'
  if (class(tree) != "phylo"){
    stop("phy must be of class 'phylo.")}   
  
  # 3. Count the number of taxa and create a vector of the species names:
  # Block 1
  num.taxa.Y<-nrow(BlockY)  
  namesY<-rownames(BlockY)
  if (is.null(namesY)){ # a little break if R cannot find species names
    stop("No specimen names in BlockY.")}
  # Block 2
  num.taxa.X1<-nrow(BlockX1)
  namesX1<-rownames(BlockX1)
  if (is.null(namesX1)){
    stop("No specimen names in Block X1.")} #confirms species names in block 2
  # Block 3
  num.taxa.X2<-nrow(BlockX2)
  namesX2<-rownames(BlockX2)
  if (is.null(namesX2)){
    stop("No specimen names in Block X2.")} #confirms species names in block 3
  
  # 4. Confirm species lists match between the phylo tree and between each variable block. 
  # first do the dimensions (species numbers) match.
  if (length(match(tree$tip.label, namesY)) != num.taxa.Y && length(tree$tip.label) < num.taxa.Y){ #confirm that there are equal numbers of taxa in the tree and the dataset
    stop("Tree is missing some taxa present in the data matrix") }
  if (length(match(tree$tip.label, namesY)) != num.taxa.Y && num.taxa.Y < length(tree$tip.label)){ #confirm that there are equal numbers of taxa in the tree and the dataset
    stop("Tree contains some taxa not present in present in the data matrix") } 
  if (length(BlockY[which(is.na(BlockY)),]) != 0) {
    stop("Transient data contains missing values. Estimate these first.")  } #no missing values are allowed
  if (length(BlockY[which(is.na(BlockX1)),]) != 0) {
    stop("Sensitivity matrix contains missing values. Estimate these first.") }
  if (length(BlockX2[which(is.na(BlockX2))]) != 0) {
    cat("Environmental data matrix contains missing values.") }
  
  # This section checks that the species IDs match between the blocks and the phylo tree.
  if (is.null(namesY) == FALSE && is.null(namesX1) == FALSE && is.null(namesX2) == FALSE) {
      mtch.A <- namesY[is.na(match(namesY, namesX1))]
    if (length(mtch.A) > 0) {
      stop("Specimen names in data sets are not the same.")} #confirms that the species names match in both blocks of data
      mtch.B <- namesX1[is.na(match(namesX1, namesX2))]
    if (length(mtch.B) > 0) {
      stop("Taxa labels on tree and taxa matrix are not the same.")} #confirms that the species names match in the data and the phylogeny
      mtch.C <- namesX2[is.na(match(namesY, namesX2))]
    if (length(mtch.C) > 0) {
        stop("Taxa labels on tree and taxa matrix are not the same.")}
  } 
  
  # 5. Format variable blocks ready for use in the PPLS analysis ensuring they are in correct matrix format.
  # and then brings variables together that require phylogenetic correction
  # Block Y
  BlockY <- as.matrix(BlockY) #transient variable
  # Block X1
  BlockX1 <- as.matrix(BlockX1) # vital rate sensitivties
  #demographic info - Transient measures and sensitivities (requires phylogenetic adjustment).
  demog.data <- cbind(BlockY, BlockX1)
  # Block X2
  enviro.data <- as.matrix(BlockX2) # environmental variance variables 

  # 6. Create phylogenetic covariance matrix - this allows for the weighting in the pls analysis below.
  # this will be the same matrix regardless of which hypothesis is being worked on
  C <- vcv.phylo(tree, anc.nodes = FALSE) # finds the phylogenetic variance-covariance matrix for the input phylogeny (this function will assume a Brownian motion correlation)
  Nspec <- nrow(C) # calculate the number of species in the vcv matrix
  # The covariance matrix C is used to compute other phylogenetic metrics using the tree and specific dataset, again this is done assuming a brownian motion distribution (denoted by lambda = 1)
  temp <- phyl.vcv(demog.data, C, lambda = 1.0)
  C <- temp$C # this output is the same as the vcv.phylo function but is just code to ensure all further workings are singing to the same tune.
  a <- t(temp$alpha) #estimation of phylogenetic mean - the character values at the root of the phylogeny (estimate the common ancestral trait values)
  
  # 7. Transform data to adjust for phylogenetic relationship
  C <- C[rownames(demog.data),rownames(demog.data)] #sorts the VCV matrix to be in the same order as the data matrix requiring adjustment
  eigC <- eigen(C) # eigenanalysis of covariance matrix
  one <- matrix(1, nrow = Nspec, ncol = 1)  #generates a vector of 1's with length = number of taxa in phylo tree
  D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors)) #transformation matrix D
  Phy.data <- D.mat %*% (demog.data - one %*% a) # this is the phylogenetically transformed data (equation 4 in the manuscript listed as Reference 1 above)
  # Split the demographic data for use below
  transientY <- as.matrix(Phy.data[,1:dim(BlockY)[2]]) # phylogenetically adjusted transient measure
  SensX <- Phy.data[,(dim(BlockY)[2]+1):dim(Phy.data)[2]] # phylogenetically adjusted vital rate sensitivities
  
  # In hypothesis 2 the analysis is only interested in terrestrial populations. Here the provided list of marine species will be used to remove these populations.
  # from the second part of the analysis.
  rownames(transientY) <- rownames(BlockY)
  transientY2 <- transientY[!(rownames(transientY) %in% marine.list)]
  enviro.data <- enviro.data[!(rownames(enviro.data) %in% marine.list),]
  
  ##### Hypothesis 1 #####
  # Using vital rate sensitivities to test for selection pressures or trade-offs between measures of the transient properties of natural populations.
   
  # 8a. Carryout PLS analysis.

  # Run the analyses
  # There is no need to include weighting factors here as phylogeny has already been adjusted for
  # estimate the correlation between Y and each X variable (Y = transient variable, X = Vital rate sensitivity matrix)
  cor1 <- cor(transientY, SensX) 
  # run the PLS (using an unweighted format)
  pls.H1 <- plsr(transientY~SensX, scale = TRUE, centre = TRUE)
  #########################
  
  ##### Hypothesis 2 ######
  
  # First bring together the relevant matrix blocks to remove incomplete cases (within the environmental variables)
  blocks.combined <- cbind(transientY2, enviro.data)
  blocks.combined <- blocks.combined[complete.cases(blocks.combined),]
  # and then split back into there two original matrix blocks
  transientY2 <- as.matrix(blocks.combined[,1:dim(BlockY)[2]])
  enviro.data <- blocks.combined[,(dim(BlockY)[2]+1):dim(blocks.combined)[2]]
  
  # Run the analyses
  # Testing for correlations between measures of environmental variability and the transient properties of natural populations
  # estimate the correlation between Y and each X variable (Y = transient variable, X = Environmental variables matrix)
  cor2 <- cor(transientY2, enviro.data) 
  # 8b. Carry out second PLS analysis
  pls.H2 <- plsr(transientY2 ~ enviro.data, scale = TRUE, centre = TRUE)
  
  # N.B. For this analysis replacing transientY2 with SensX (after subsetting) could also work for investigating how vital rate sensitivities correlate with environmental variability.
  # This would provide an alternative angle for investigating the role of the environment in constraining resilience charactersitics.
  #########################
  
  # A quick peice of indexing to subset the transient variable to match the data used in the second pls analysis
  BlockY <- BlockY[which(rownames(BlockY) %in% rownames(transientY2)),]
  
  # 9.Return keyout outputs
  return(list("Vital Rate Correlation" = cor1, "Hypothesis 1 pls" = pls.H1, 
              "Environmental variance correlation" = cor2, "Hypothesis 2 pls" = pls.H2,
              "Transient variable" = BlockY))
}
