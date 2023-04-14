##############################################
# Functions for extracting transient characteristics and their vital rate sensitivities from a series of Matrix Population Models
##############################################

### 1. Base functions
# These functions will be built into the more complex functions below.

# Damping ratio sensitivity
SensDR <- function(matA){
  DR <- eigen.analysis(matA)$damping.ratio #calculate damping ratio
  S1 <- sens(matA, eval = "max", all= FALSE) #determine sensitivity of dominant eigenvalue for only observed transitions
  S2 <- sens(matA, eval = 2, all = FALSE) #determine sensitivity of subdomiant eigenvalue for observed transitions (this accounts for the fact the subdomiant eigenvalue is complex, consiting of a real and imaginary component)
  Lambda2 = abs(eigen(matA)$values)[2] #determine the absolute value of the subdominant eigenvalue
  # this produces the parameters needed for equation 9.117 (Caswell 2001 pg 245)
  Srho <- (1/Lambda2)*(S1-((DR/Lambda2)*(Re(S2)+Im(S2)))) #this returns the element level sensitivities of the damping ratio 
  return(Srho)
}

# Damping ratio elasticity
ElasDR <- function(Srho, matA){
  DR <- eigen.analysis(matA)$damping.ratio #calculate damping ratio
  # this function now corresponds to equation 9.119 (Caswell 2001 pg 245)
  Erho <- hadamard.prod((1/DR)*Srho, matA) #the hadamard product is element by element multiplication of two matrices of equal dimensions.
  return(Erho)
}

# period of oscillation (Caswell 2001, pg 100, equations 4.98 & 4.99)
Period.O <- function(matA){
  Lambda2_real <- Re(eigen(matA)$value)[2] #the real component of the complex number
  Lambda2_imag <- Im(eigen(matA)$value)[2] #the imaginary component of the complex number
  theta <- atan(Lambda2_imag/Lambda2_real) #equivalent to the angle formed by the subdominent eigenvalue in the complex plane (equation 4.99 Caswell 2001 pg 100)
  PO <- abs((2*pi)/theta) # the period of oscillation
  return(PO)
}

# Period of oscillation sensitivity
SensPO <- function(matA){
  Lambda2_real <- Re(eigen(matA)$values)[2] #the real component of the complex number
  Lambda2_imag <- Im(eigen(matA)$values)[2] #the imaginary component of the complex number
  Lambda2 <- abs(eigen(matA)$values)[2] #the absolute value of the complex subdomiant eigenvalue
  theta <- atan(Lambda2_imag/Lambda2_real) #equivalent to the angle formed by the subdominent eigenvalue in the complex plane (equation 4.99 Caswell 2001 pg 100)
  # atan is the inverse tangent function (tan^-1)
  S2 <- sens(matA, eval = 2, all = FALSE)
  # this produces the parameters needed for equation 9.121 (Caswell 2001 pg 247)
  Speriod <- ((-2*pi)/((theta^2)*(Lambda2^2)))*(Im(S2)-Re(S2))
  return(Speriod)
}

# Period of oscillation elasticitiy.
ElasPO <- function(Speriod, PO, matA){
  # this function corresponds to equation 9.122 (Caswell 2001 pg 247)
  Eperiod <- hadamard.prod((1/PO)*Speriod, matA)
  return(Eperiod)
}

# Sensitivity & Elasticity of matrix reactivity - brute force method (Morris & Doak 2002)
SensReact <- function(matA){
  react1 <- reac(matA, bound = "upper") #determine the reactivity of the original model
  # and the dimensions of the model
  drow <- dim(matA)[1] #number of rows
  dcol <- dim(matA)[2] #number of columns
  #create blank sensitivity matrix
  SensA <- matrix(NA, drow, dcol)
  ElasA <- matrix(NA, drow,dcol)
  # now loop through each element, changing it slightly and recalculating the reactivity of the model. This becomes that elements reactivity sensitivity value.
  for (i in 1:drow) {
    for(j in 1:dcol) {
      matNew <- matA #resave the original matrix to avoid changing it accidently (this also ensures the model is only changing one element at a time)
      perturb <- 0.01 #makes it simple to change the size of the perturbation
      matNew[i,j] <- matNew[i,j] + (matNew[i,j]*perturb) #make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)
      elementdiff <- matNew[i,j]*perturb #store the actuall change in the element
      react2 <- reac(matNew, bound = "upper") #recalculate the matrix's reactivity
      SensA[i,j] <- (react2 - react1)/elementdiff # using the change in reactivity and the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
      if(SensA[i,j] %in% c(Inf, NA)){SensA[i,j] = 0}
      ElasA[i,j] <- ((react2-react1)/react1)/(elementdiff/matA[i,j])
    }
  }
  SensA[is.nan(SensA)]=0
  ElasA[is.nan(ElasA)]=0
  return(list(SensA,ElasA)) #this now returns the sensitivity and elasticity matrices for reactivity of matA.
}

# Sensitivity & Elasticity of the first time step attenuation of the matrix - brute force method (Morris & Doak 2002)
SensAtten <- function(matA){
  Atten1 <- reac(matA, bound = "lower") #determine the reactivity of the original model
  # and the dimensions of the model
  drow <- dim(matA)[1] #number of rows
  dcol <- dim(matA)[2] #number of columns
  #create blank sensitivity matrix
  SensA <- matrix(NA, drow, dcol)
  ElasA <- matrix(NA, drow,dcol)
  # now loop through each element, changing it slightly and recalculating the reactivity of the model. This becomes that elements reactivity sensitivity value.
  for (i in 1:drow) {
    for(j in 1:dcol) {
      matNew <- matA #resave the original matrix to avoid changing it accidently (this also ensures the model is only changing one element at a time)
      perturb <- 0.01 #makes it simple to change the size of the perturbation
      matNew[i,j] <- matNew[i,j] + (matNew[i,j]*perturb) #make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)
      elementdiff <- matNew[i,j]*perturb
      Atten2 <- reac(matNew, bound = "lower") #recalculate the matrix's reactivity
      SensA[i,j] <- (Atten2 - Atten1)/elementdiff # using the change in reactivity and the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
      if(SensA[i,j] %in% c(Inf, NA)){SensA[i,j] = 0}
      ElasA[i,j] <- ((Atten2-Atten1)/Atten1)/(elementdiff/matA[i,j])
    }
  }
  SensA[is.nan(SensA)]=0
  ElasA[is.nan(ElasA)]=0
  return(list(SensA,ElasA)) #this now returns the sensitivity and elasticity matrices for reactivity of matA.
}

# Sensitivity & Elasticity of the maximum amplification of the matrix - brute force method (Morris & Doak 2002)
SensAmp <- function(matA){
  Amp1 <- maxamp(matA, vector = "n") #determine the reactivity of the original model
  # and the dimensions of the model
  drow <- dim(matA)[1] #number of rows
  dcol <- dim(matA)[2] #number of columns
  #create blank sensitivity matrix
  SensA <- matrix(NA, drow, dcol)
  ElasA <- matrix(NA, drow, dcol)
  # now loop through each element, changing it slightly and recalculating the reactivity of the model. This becomes that elements reactivity sensitivity value.
  for (i in 1:drow) {
    for(j in 1:dcol) {
      matNew <- matA #resave the original matrix to avoid changing it accidently (this also ensures the model is only changing one element at a time)
      perturb <- 0.01 #makes it simple to change the size of the perturbation
      matNew[i,j] <- matNew[i,j] + (matNew[i,j]*perturb) #make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)
      elementdiff <- matNew[i,j]*perturb
      Amp2 <- maxamp(matNew, vector = "n") #recalculate the matrix's reactivity
      SensA[i,j] <- (Amp2 - Amp1)/elementdiff # using the change in reactivity and the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
      if(SensA[i,j] %in% c(Inf, NA)){SensA[i,j] = 0}
      ElasA[i,j] <- ((Amp2-Amp1)/Amp1)/(elementdiff/matA[i,j])
    }
  }
  SensA[is.nan(SensA)]=0
  ElasA[is.nan(ElasA)]=0
  return(list(SensA,ElasA)) #this now returns the sensitivity and elasticity matrices for reactivity of matA.
}

# Sensitivity & Elasticity of the maximum attenuation of the matrix - brute force method (Morris & Doak 2002)
SensAtt <- function(matA){ #careful this is a similar function name
  Maxatt1 <- maxatt(matA, vector = "n") #determine the max attenuation of the original model
  # and the dimensions of the model
  drow <- dim(matA)[1] #number of rows
  dcol <- dim(matA)[2] #number of columns
  #create blank sensitivity matrix
  SensA <- matrix(NA, drow, dcol)
  ElasA <- matrix(NA, drow, dcol)
  # now loop through each element, changing it slightly and recalculating the max attenuation of the model. This becomes that elements reactivity sensitivity value.
  for (i in 1:drow) {
    for(j in 1:dcol) {
      matNew <- matA #resave the original matrix to avoid changing it accidently (this also ensures the model is only changing one element at a time)
      perturb <- 0.01 #makes it simple to change the size of the perturbation
      matNew[i,j] <- matNew[i,j] + (matNew[i,j]*perturb) #make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)
      elementdiff <- matNew[i,j]*perturb
      Maxatt2 <- maxatt(matNew, vector = "n") #recalculate the matrix's attenuation
      SensA[i,j] <- (Maxatt2 - Maxatt1)/elementdiff # using the change in attenuation the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
      if(SensA[i,j] %in% c(Inf, NA)){SensA[i,j] = 0}
      ElasA[i,j] <- ((Maxatt2-Maxatt1)/Maxatt1)/(elementdiff/matA[i,j])
    }
  }
  SensA[is.nan(SensA)]=0
  ElasA[is.nan(ElasA)]=0
  return(list(SensA,ElasA)) #this now returns the sensitivity and elasticity matrices for attenuation of matA.
}

### 2. Vital-Rate level functions 

#Function to calculate vital rate level sensitivities and elasticities of the Damping ratio (for unadjusted matrices)
vitalRatePerturbationDR <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate the damping ratio
  DR = eigen.analysis(matA)$damping.ratio
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the sensitivities of the damping ratio to each element 
  sensA=SensDR(matA)
  # now the associated elasticities
  ElasA=ElasDR(Srho = sensA, matA = matA)
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/DR #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/DR #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/DR #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/DR #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out1 = data.frame("rho" = DR, #Damping ratio
                    
                    "SSurvival"=sum(sensSurv,na.rm=T),
                    "SGrowth"=sum(sensGrow,na.rm=T),
                    "SShrinkage"=sum(sensShri,na.rm=T),
                    "SReproduction"=sum(sensFec,na.rm=T),
                    "SClonality"=sum(sensClo,na.rm=T),
                    
                    "ESurvival"=sum(elasSurv,na.rm=T),
                    "EGrowth"=sum(elasGrow,na.rm=T),
                    "EShrinkage"=sum(elasShri,na.rm=T),
                    "EReproduction"=sum(elasFec,na.rm=T),
                    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out1$rho)){out1$rho <- NA}
  if (is.nan(out1$rho)) {out1$rho <- NA}
  if (is.na(out1$rho)) {out1[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out1)
}

#Function to calculate vital rate level sensitivities and elasticities of Period of Occilation (none adjusted matrices)
vitalRatePerturbationPO <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate the period of oscillation
  PO = Period.O(matA)
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the sensitivities of the period of oscilation to each element 
  sensA=SensPO(matA)
  # now the associated elasticities
  ElasA=ElasPO(Speriod = sensA, PO = PO, matA = matA)
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/PO #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/PO #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/PO #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/PO #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out2 = data.frame(
    "Pi"= PO, #period of oscilation
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out2$Pi)){out2$Pi <- NA}
  if (is.nan(out2$Pi)) {out2$Pi <- NA}
  if (is.na(out2$Pi)) {out2[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out2)
}

#Function to calculate vital rate level sensitivities and elasticities of Reactivity (none adjusted matrices)
vitalRatePerturbationReac <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate the reactivity of the matrix
  ReacA = reac(matA,bound="upper")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the reactivity to each element 
  SensR <- SensReact(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensR[[1]]
  ElasA <- SensR[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/ReacA #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/ReacA #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/ReacA #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/ReacA #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out3 = data.frame(
    "reactivity"= ReacA, #reactivity
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out3$reactivity)){out3$reactivity <- NA}
  if (is.nan(out3$reactivity)) {out3$reactivity <- NA}
  if (is.na(out3$reactivity)) {out3[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out3)
}

#Function to calculate vital rate level sensitivities and elasticities of first time step attenuation (none adjusted matrices)
vitalRatePerturbationAtten <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate the reactivity of the matrix
  AttenA = reac(matA,bound="lower")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the reactivity to each element 
  SensR <- SensAtten(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensR[[1]]
  ElasA <- SensR[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/AttenA #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/AttenA #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/AttenA #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/AttenA #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out7 = data.frame(
    "attenuation"= AttenA, #reactivity
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out7$attenuation)){out7$attenuation <- NA}
  if (is.nan(out7$attenuation)) {out7$attenuation <- NA}
  if (is.na(out7$attenuation)) {out7[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out7)
}

#Function to calculate vital rate level sensitivities and elasticities of Maximum Amplification (non-standardised matrices)
vitalRatePerturbationMaxAmp <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  
  #calculate the maximum amplification of the matrix
  MaxA = maxamp(matA,vector = "n")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the maximum amplification to each element 
  SensMaxAmp <- SensAmp(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensMaxAmp[[1]]
  ElasA <- SensMaxAmp[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/MaxA #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/MaxA #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/MaxA #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/MaxA #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out9 = data.frame(
    "MaxAmp"= MaxA, #Maximum amplification
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out9$MaxAmp)){out9$MaxAmp <- NA}
  if (is.nan(out9$MaxAmp)) {out9$MaxAmp <- NA}
  if (is.na(out9$MaxAmp)) {out9[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out9)
}

#Function to calculate vital rate level sensitivities and elasticities of Maximum Attenuation (none standardised matrices)
vitalRatePerturbationMaxAtt <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  
  #calculate the maximum attenuation of the matrix
  MaxAtt = maxatt(matA,vector = "n")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the maximum attenuation to each element 
  SensMaxAtt <- SensAtt(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensMaxAtt[[1]]
  ElasA <- SensMaxAtt[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/MaxAtt #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/MaxAtt #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/MaxAtt #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/MaxAtt #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out11 = data.frame(
    "MaxAtt"= MaxAtt, #Maximum amplification
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out11$MaxAtt)){out11$Maxtt <- NA}
  if (is.nan(out11$MaxAtt)) {out11$MaxAtt <- NA}
  if (is.na(out11$MaxAtt)) {out11[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out11)
}

#Function to calculate vital rate level sensitivities and elasticities of the Damping ratio (standardised matrices)
S.vitalRatePerturbationDR <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  
  #calculate lambda
  lambda1 <- Re(eigen(matA)$values)[1]
  
  #standardise submatrices
  matU = matU/lambda1
  matF = matF/lambda1
  matC = matC/lambda1
  matA = matA/lambda1
  
  #calculate the damping ratio
  DR = eigen.analysis(matA)$damping.ratio
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the sensitivities of the damping ratio to each element 
  sensA=SensDR(matA)
  # now the associated elasticities
  ElasA=ElasDR(Srho = sensA, matA = matA)
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/DR #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/DR #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/DR #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/DR #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out4 = data.frame("rho" = DR, #Damping ratio
                    
                    "SSurvival"=sum(sensSurv,na.rm=T),
                    "SGrowth"=sum(sensGrow,na.rm=T),
                    "SShrinkage"=sum(sensShri,na.rm=T),
                    "SReproduction"=sum(sensFec,na.rm=T),
                    "SClonality"=sum(sensClo,na.rm=T),
                    
                    "ESurvival"=sum(elasSurv,na.rm=T),
                    "EGrowth"=sum(elasGrow,na.rm=T),
                    "EShrinkage"=sum(elasShri,na.rm=T),
                    "EReproduction"=sum(elasFec,na.rm=T),
                    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out4$rho)){out4$rho <- NA}
  if (is.nan(out4$rho)) {out4$rho <- NA}
  if (is.na(out4$rho)) {out4[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out4)
}

#Function to calculate vital rate level sensitivities and elasticities of Period of Occilation (standardised matrices)
S.vitalRatePerturbationPO <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate lambda
  lambda1 <- Re(eigen(matA)$values)[1]
  
  #standardise submatrices
  matU = matU/lambda1
  matF = matF/lambda1
  matC = matC/lambda1
  matA = matA/lambda1
  
  #calculate the period of oscillation
  PO = Period.O(matA)
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the sensitivities of the period of oscilation to each element 
  sensA=SensPO(matA)
  # now the associated elasticities
  ElasA=ElasPO(Speriod = sensA, PO = PO, matA = matA)
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/PO #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/PO #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/PO #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/PO #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out5 = data.frame(
    "Pi"= PO, #period of oscilation
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out5$Pi)){out5$Pi <- NA}
  if (is.nan(out5$Pi)) {out5$Pi <- NA}
  if (is.na(out5$Pi)) {out5[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out5)
}

#Function to calculate vital rate level sensitivities and elasticities of Reactivity (standardised matrices)
S.vitalRatePerturbationReac <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate lambda
  lambda1 <- Re(eigen(matA)$values)[1]
  
  #standardised submatrices
  matU = matU/lambda1
  matF = matF/lambda1
  matC = matC/lambda1
  matA = matA/lambda1
  
  #calculate the reactivity of the matrix
  ReacA = reac(matA,bound="upper")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the reactivity to each element 
  SensR <- SensReact(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensR[[1]]
  ElasA <- SensR[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/ReacA #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/ReacA #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/ReacA #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/ReacA #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out6 = data.frame(
    "reactivity"= ReacA, #reactivity
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out6$reactivity)){out6$reactivity <- NA}
  if (is.nan(out6$reactivity)) {out6$reactivity <- NA}
  if (is.na(out6$reactivity)) {out6[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out6)
}

#Function to calculate vital rate level sensitivities and elasticities of first time step attenuation (standardised matrix)
S.vitalRatePerturbationAtten <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate lambda
  lambda1 <- Re(eigen(matA)$values)[1]
  
  #standardised submatrices
  matU = matU/lambda1
  matF = matF/lambda1
  matC = matC/lambda1
  matA = matA/lambda1
  
  #calculate the reactivity of the matrix
  AttenA = reac(matA,bound="lower")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the reactivity to each element 
  SensR <- SensAtten(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensR[[1]]
  ElasA <- SensR[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/AttenA #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/AttenA #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/AttenA #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/AttenA #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out8 = data.frame(
    "attenuation"= AttenA, #reactivity
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out8$attenuation)){out8$attenuation <- NA}
  if (is.nan(out8$attenuation)) {out8$attenuation <- NA}
  if (is.na(out8$attenuation)) {out8[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out8)
}

#Function to calculate vital rate level sensitivities and elasticities of Maximum Amplification (standardised matrices)
S.vitalRatePerturbationMaxAmp <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate lambda
  lambda1 <- Re(eigen(matA)$values)[1]
  
  #standardised submatrices
  matU = matU/lambda1
  matF = matF/lambda1
  matC = matC/lambda1
  matA = matA/lambda1
  
  #calculate the maximum amplification of the matrix
  MaxA = maxamp(matA,vector = "n")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the maximum amplification to each element 
  SensMaxAmp <- SensAmp(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensMaxAmp[[1]]
  ElasA <- SensMaxAmp[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/MaxA #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/MaxA #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/MaxA #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/MaxA #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out10 = data.frame(
    "MaxAmp"= MaxA, #Maximum amplification
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out10$MaxAmp)){out10$MaxAmp <- NA}
  if (is.nan(out10$MaxAmp)) {out10$MaxAmp <- NA}
  if (is.na(out10$MaxAmp)) {out10[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out10)
}

#Function to calculate vital rate level sensitivities and elasticities of Maximum Attenuation (standardised matrices)
S.vitalRatePerturbationMaxAtt <- function(matU, matF, matC){
  
  # fix any elements entered as NA to zeros
  matU[is.na(matU)]=0
  matF[is.na(matF)]=0
  matC[is.na(matC)]=0
  
  # combine submatrices to generate overal population matrix consisting of all vital rates
  matA=matU+matF+matC
  #calculate lambda
  lambda1 <- Re(eigen(matA)$values)[1]
  
  #standardised submatrices
  matU = matU/lambda1
  matF = matF/lambda1
  matC = matC/lambda1
  matA = matA/lambda1
  
  #calculate the maximum attenuation of the matrix
  MaxAtt = maxatt(matA,vector = "n")
  # determine matrix dimensions
  matDim=dim(matA)[1]
  # create a stage-based survival vector
  sigma=colSums(matU)
  # determine the elasticities and sensitivities of the maximum attenuation to each element 
  SensMaxAtt <- SensAtt(matA) #this produces a list object containing both sensitivities and elasticities
  # these two lines of code seperate them out.
  sensA <- SensMaxAtt[[1]]
  ElasA <- SensMaxAtt[[2]]
  # create a blank survival independant matrix
  noSurvA=matrix(NA,matDim,matDim)
  # fill it by dividing each element by its associated stage-based survival to return a survival independant matrix. 
  for (i in 1:matDim){noSurvA[,i]=matA[,i]/sigma[i]}
  
  ### survival ###
  # create adjusted survival only sensitivitiy matrix
  matSensSurv=noSurvA*sensA
  sensSurv=colSums(matSensSurv) #calculate column sums (stage based survival sensitivities vector)
  elasSurv=sigma*sensSurv/MaxAtt #adjust to calculate elasticities 
  
  # create blank matrices for growth and shrinkage
  matSensGrowShri=matSensGrow=matSensShri=matrix(NA,matDim,matDim) #blank sensitivity matrices
  matElasGrowShri=matElasGrow=matElasShri=matElasFec=matElasClo=matrix(NA,matDim,matDim) #blank elasticitiy matrices
  
  ### Growth & Retrogression ###
  for (i in 1:matDim){matSensGrowShri[,i]=sensA[i,i]*(-sigma[i])+sensA[,i]*(+sigma[i])} #fill blank matrices following equations 8 & 9 Franco and Silvertown 2004
  matSensGrowShri[which(matU==0)]=0 #if survival is 0 then so are growth and shrinkage
  matSensGrow[lower.tri(matSensGrowShri,diag=T)]=matSensGrowShri[lower.tri(matSensGrowShri,diag=T)] #growth sensitivites only correspond to elements in the lower half of the matrix.
  sensGrow=colSums(matSensGrow,na.rm=T) #calculate the stage-based growth sensitivities
  matSensShri[upper.tri(matSensGrowShri,diag=T)]=matSensGrowShri[upper.tri(matSensGrowShri,diag=T)] #retrogression sensitivities corespond to the upper half of the matrix
  sensShri=colSums(matSensShri,na.rm=T) #calculate the stage-based retrogression sensitivities
  matElasGrowShri=noSurvA*matSensGrowShri/MaxAtt #create and growth and shrinkage elasticitiy matrix
  matElasGrow[lower.tri(matElasGrowShri,diag=T)]=matElasGrowShri[lower.tri(matElasGrowShri,diag=T)] #lower half elasticities correspond to growth
  elasGrow=colSums(matElasGrow,na.rm=T) 
  matElasShri[upper.tri(matElasGrowShri,diag=T)]=matElasGrowShri[upper.tri(matElasGrowShri,diag=T)] #upper half elasticities correspond to shrinkage.
  elasShri=colSums(matElasShri,na.rm=T)
  
  ### Fecundity ###
  matSensFec=matrix(0,matDim,matDim) #create a blank fecundity matrix - with zeros instead of NAs
  matSensFec[which(matF>0)]=sensA[which(matF>0)]*sigma[which(colSums(matF)>0)] #the senstivity elements involving reproductive info multiplied by the associated stage survival to give a sensitivity matrix for just reproduction (equation 10 Franco & Silvertown 2004)
  sensFec=colSums(matSensFec) #calculate stage-based fecundity sensitivities
  matElasFec=noSurvA*matSensFec/MaxAtt #calculate fecundity elasticities 
  elasFec=colSums(matElasFec)
  
  ### Clonality ###
  matSensClo=matrix(0,matDim,matDim) # blank clonality matrix
  matSensClo[which(matC>0)]=sensA[which(matC>0)]*sigma[which(colSums(matC)>0)] #the senstivity elements involving clonal info multiplied by the associated stage survival to give a sensitivity matrix for just clonality (equation 11 Franco & Silvertown 2004)
  sensClo=colSums(matSensClo) #calculate stage-based clonality sensitivities
  matElasClo=noSurvA*matSensClo/MaxAtt #calculate clonality elasticities
  elasClo=colSums(matElasClo)
  
  # store the desired outputs from the function
  out12 = data.frame(
    "MaxAtt"= MaxAtt, #Maximum amplification
    
    "SSurvival"=sum(sensSurv,na.rm=T),
    "SGrowth"=sum(sensGrow,na.rm=T),
    "SShrinkage"=sum(sensShri,na.rm=T),
    "SReproduction"=sum(sensFec,na.rm=T),
    "SClonality"=sum(sensClo,na.rm=T),
    
    "ESurvival"=sum(elasSurv,na.rm=T),
    "EGrowth"=sum(elasGrow,na.rm=T),
    "EShrinkage"=sum(elasShri,na.rm=T),
    "EReproduction"=sum(elasFec,na.rm=T),
    "EClonality"=sum(elasClo,na.rm=T))
  
  # a little fix to prevent the function returning zeros when in fact it cannot calculate values.
  if (is.infinite(out12$MaxAtt)){out12$Maxtt <- NA}
  if (is.nan(out12$MaxAtt)) {out12$MaxAtt <- NA}
  if (is.na(out12$MaxAtt)) {out12[,c("SSurvival","SGrowth","SShrinkage","SReproduction","SClonality","ESurvival","EGrowth","EShrinkage","EReproduction","EClonality")] <- NA}
  
  return(out12)
}

# 3. Function for converting post-reproductive matrices into pre-reproductive matrices
convert2pre <- function(matA, matF, matC, matU, post){
  # create output storage
  mat_corrected_pre <- list()
  # if the matrix is not post-reproductive then simply return the original matrices
  if(!isTRUE(post)){
    mat_corrected_pre[[1]] <- matA
    mat_corrected_pre[[2]] <- matU
    mat_corrected_pre[[3]] <- matF
    mat_corrected_pre[[4]] <- matC
    } 
  # otherwise reproduction needs fixing
  else{
    mat <- matA
    matF <- matF
    matC <- matC
    matS <- matU

    # determine matrix dimensions
    matrix_size <- nrow(mat)
    
    # determine stage survival
    surv_vec<-apply(matS,2,sum)
    surv_mat<-matrix(surv_vec,nrow=matrix_size,ncol=matrix_size,byrow=T)
    # divide seedbank survival by stage based survival (i.e the chance of recruitment is now no longer contingent on the survival of adults)
    newmat1<-matF/surv_mat
    newmat1[is.nan(newmat1)] = 0 # correction for NaNs
    
    # generate corrected pre-reproductive reproductive matrix
    newmatF <- matS%*%newmat1
    newmatF[is.nan(newmatF)] = 0
    
    # and use it to recalculate the overall population matrix
    newmat <- newmatF+matC+matS
    newmat[is.nan(newmat)] = 0
    
    # apply a small little fix
      if(mat[1,1]==0){
        mat_corrected_pre[[1]] <- newmat[2:matrix_size,2:matrix_size]
        mat_corrected_pre[[2]] <- matS[2:matrix_size,2:matrix_size]
        mat_corrected_pre[[3]] <- newmatF[2:matrix_size,2:matrix_size]
        mat_corrected_pre[[4]] <- matC[2:matrix_size,2:matrix_size]
      }else{
        mat_corrected_pre[[1]] <- newmat
        mat_corrected_pre[[2]] <- matS
        mat_corrected_pre[[3]] <- newmatF
        mat_corrected_pre[[4]] <- matC
        }
      if(any(is.infinite(newmat))){
        mat_corrected_pre[1] <- NA
        mat_corrected_pre[2] <- NA
        mat_corrected_pre[3] <- NA
        mat_corrected_pre[4] <- NA}
    }
  
  return(mat_corrected_pre)
}
