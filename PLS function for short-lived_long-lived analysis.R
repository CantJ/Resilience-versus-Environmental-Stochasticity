# ---------------------------------------------------
# This function is for carrying out a phylogenetically corrected Canonical partial least squares analysis.
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

Phylo.PLS <-function(BlockY, BlockX, tree, marine.list){ 
  # this function takes three variable sets and the associated phylogenetic tree.
  
  # BlockY is a vector/matrix of numerical values describing a transient property 
  # measured from various populations (i.e. Damping ratio, Period of Oscillation, 
  # Reactivity, Maximal Amplification, First-step Attenuation, Maximal Attenuation).
  
  # BlockX is a matrix describing the environmental variability to which each 
  # population is exposed. This variable set will consist of two abiotic variables: 
  # Variance frequency and thermal range (variance magnitude).   
  
  # BlockY and BlockX will be used in checking the analysis of Hypothesis 2, and whether it holds for long or short lived species.
  
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
  num.taxa.X<-nrow(BlockX)
  namesX<-rownames(BlockX)
  if (is.null(namesX)){
    stop("No specimen names in BlockX.")} #confirms species names in block 2
  
  # 4. Format variable blocks ready for use in the PPLS analysis ensuring they are in correct matrix format.
  # and then brings variables together that require phylogenetic correction
  # Block Y
  BlockY <- as.matrix(BlockY) #transient variable
  # Demographic info - Transient measures (requires phylogenetic adjustment).
  demog.data <- BlockY
  # Block X
  enviro.data <- as.matrix(BlockX) # environmental variance variables 
  
  # 5. Create phylogenetic covariance matrix - this allows for the weighting in the pls analysis below.
  # this will be the same matrix regardless of which hypothesis is being worked on
  C <- vcv.phylo(tree, anc.nodes = FALSE) # finds the phylogenetic variance-covariance matrix for the input phylogeny (this function will assume a Brownian motion correlation)
  Nspec <- nrow(C) # calculate the number of species in the vcv matrix
  # The covariance matrix C is used to compute other phylogenetic metrics using the tree and specific dataset, again this is done assuming a brownian motion distribution (denoted by lambda = 1)
  temp <- phyl.vcv(demog.data, C, lambda = 1.0)
  C <- temp$C # this output is the same as the vcv.phylo function but is just code to ensure all further workings are singing to the same tune.
  a <- t(temp$alpha) #estimation of phylogenetic mean - the character values at the root of the phylogeny (estimate the common ancestral trait values)
  
  # 6. Transform data to adjust for phylogenetic relationship
  C <- C[rownames(demog.data),rownames(demog.data)] #sorts the VCV matrix to be in the same order as the data matrix requiring adjustment
  eigC <- eigen(C) # eigenanalysis of covariance matrix
  one <- matrix(1, nrow = Nspec, ncol = 1)  #generates a vector of 1's with length = number of taxa in phylo tree
  D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors)) #transformation matrix D
  Phy.data <- D.mat %*% (demog.data - one %*% a) # this is the phylogenetically transformed data (equation 4 in the manuscript listed as Reference 1 above)
  # Split the demographic data for use below
  transientY <- as.matrix(Phy.data) # phylogenetically adjusted transient measure
  
  # In hypothesis 2 the analysis is only interested in terrestrial populations. Here the provided list of marine species will be used to remove these populations.
  # from the second part of the analysis.
  rownames(transientY) <- rownames(BlockY)
  transientY <- transientY[!(rownames(transientY) %in% marine.list)]
  enviro.data <- enviro.data[!(rownames(enviro.data) %in% marine.list),]
  
  ##### Hypothesis 2 ######
  
  # First bring together the relevant matrix blocks to remove incomplete cases (within the environmental variables)
  blocks.combined <- cbind(transientY, enviro.data)
  blocks.combined <- blocks.combined[complete.cases(blocks.combined),]
  # and then split back into the two original matrix blocks
  transientY2 <- as.matrix(blocks.combined[,1:dim(BlockY)[2]])
  enviro.data <- blocks.combined[,(dim(BlockY)[2]+1):dim(blocks.combined)[2]]
  
  # 7. Run the analyses
  # Testing for correlations between measures of environmental variability and the transient properties of natural populations
  # estimate the correlation between Y and each X variable (Y = transient variable, X = Environmental variables matrix)
  cor1 <- cor(transientY2, enviro.data) 
  # 7b. Carry out second PLS analysis
  pls.H1 <- plsr(transientY2 ~ enviro.data, scale = TRUE, centre = TRUE)
  
  # A quick peice of indexing to subset the transient variable to match the data used in the second pls analysis
  BlockY <- BlockY[which(rownames(BlockY) %in% rownames(transientY2)),]
  
  # 9.Return keyout outputs
  return(list("Environmental variance correlation" = cor1, "Hypothesis 2 pls" = pls.H1,
              "Transient variable" = BlockY))
}
