# This script is for extracting a bunch of demographic data relating to the transient nature of matrix models sourced from COMPADRE and COMADRE with associated GPS coordinates.

# Authors: James Cant
# Last Modified: April 2021
# --------------------------------------------------------------------------------------------------------------------------------------------------

# clear workspace
rm(list=ls(all=TRUE))

# set working directory
setwd("FILE_PATH")

# load required programs
library(popbio)
library(popdemo)
library(Rcompadre)
library(matrixcalc)

# load transient sensitivity functions
source("FILE_PATH/Transient functions.R")

##############################################
# STEP 2: Load and Extract demographic information
##############################################

#Upload COMPADRE and COMADRE
load("COMPADRE_v.X.X.X.RData")
compadre <- as_cdb(compadre) # convert to CompadreDB format
load("COMADRE_v.X.X.X.RData")
comadre <- as_cdb(comadre)

#Take subsets for species with GPS info, and for matrices from wild, unmanipulated populations.
indexCOMADRE=which(!is.na(comadre$Lat) & comadre$MatrixTreatment == "Unmanipulated" & comadre$MatrixCaptivity == "W" & comadre$MatrixDimension >= 3)
indexCOMPADRE=which(!is.na(compadre$Lat) & compadre$MatrixTreatment == "Unmanipulated" & compadre$MatrixCaptivity == "W" & compadre$MatrixDimension >= 3)
# for now I have left in populations that where surveyed at varying frequencies (perodicities) so as to not limit my sample size too much. 
# But mean matrices are being selected because individual matrices only account for the dynamics in the year in which they where measured, whereas the mean dynamics can be correlated better with legacy effects of the environment.

#run the subsetting across all levels of the data files.
comadre <- comadre[indexCOMADRE]
compadre <- compadre[indexCOMPADRE]

#dataframe to store outputs from the extraction
long = dim(comadre)[1] + dim(compadre)[1]
output=data.frame(
  #metadata
  "SpeciesAuthor"=rep(NA,long),
  "SpeciesAccepted"=rep(NA,long),
  "OrganismType"=rep(NA,long),
  "Family"=rep(NA,long),
  "Class"=rep(NA,long),
  "Kingdom"=rep(NA,long),
  "Authors"=rep(NA,long),
  "Journal"=rep(NA,long),
  "YearPublication"=rep(NA,long),
  "Country"=rep(NA,long),
  "Continent"=rep(NA,long),
  "Ecoregion"=rep(NA,long),
  "MatrixDimension"=rep(NA,long),
  
  "Population"=rep(NA,long),
  "StartYear"=rep(NA,long),
  "StartMonth"=rep(NA,long),
  "StartSeason"=rep(NA,long),
  "EndYear"=rep(NA,long),
  "EndMonth"=rep(NA,long),
  "EndSeason"=rep(NA,long),
  "Periodicity"=rep(NA,long),
  "Treatment"=rep(NA,long),
  "MatrixComposite"=rep(NA,long),
  "Lat"=rep(NA,long),
  "Lon"=rep(NA,long),
  "Alt"=rep(NA,long),
  
  #damping ratio analyses
  "Rho1" = rep(NA,long),
  "SSurvDR"=rep(NA,long),
  "SGrowDR"=rep(NA,long),
  "SShriDR"=rep(NA,long),
  "SRepDR"=rep(NA,long),
  "SCloDR"=rep(NA,long),
  
  "ESurvDR"=rep(NA,long),
  "EGrowDR"=rep(NA,long),
  "EShriDR"=rep(NA,long),
  "ERepDR"=rep(NA,long),
  "ECloDR"=rep(NA,long),
  
  "Rho2" = rep(NA,long),
  "SSurvDRs"=rep(NA,long),
  "SGrowDRs"=rep(NA,long),
  "SShriDRs"=rep(NA,long),
  "SRepDRs"=rep(NA,long),
  "SCloDRs"=rep(NA,long),
  
  "ESurvDRs"=rep(NA,long),
  "EGrowDRs"=rep(NA,long),
  "EShriDRs"=rep(NA,long),
  "ERepDRs"=rep(NA,long),
  "ECloDRs"=rep(NA,long),
  
  #Period of Occilations analyses
  "Pi1" = rep(NA, long),
  "SSurvPO"=rep(NA,long),
  "SGrowPO"=rep(NA,long),
  "SShriPO"=rep(NA,long),
  "SRepPO"=rep(NA,long),
  "SCloPO"=rep(NA,long),
  
  "ESurvPO"=rep(NA,long),
  "EGrowPO"=rep(NA,long),
  "EShriPO"=rep(NA,long),
  "ERepPO"=rep(NA,long),
  "ECloPO"=rep(NA,long),
  
  "Pi2" = rep(NA, long),
  "SSurvPOs"=rep(NA,long),
  "SGrowPOs"=rep(NA,long),
  "SShriPOs"=rep(NA,long),
  "SRepPOs"=rep(NA,long),
  "SCloPOs"=rep(NA,long),
  
  "ESurvPOs"=rep(NA,long),
  "EGrowPOs"=rep(NA,long),
  "EShriPOs"=rep(NA,long),
  "ERepPOs"=rep(NA,long),
  "ECloPOs"=rep(NA,long),
  
  #reactivity analyses
  "Reactivity1" = rep(NA, long),
  "SSurvR"=rep(NA,long),
  "SGrowR"=rep(NA,long),
  "SShriR"=rep(NA,long),
  "SRepR"=rep(NA,long),
  "SCloR"=rep(NA,long),
  
  "ESurvR"=rep(NA,long),
  "EGrowR"=rep(NA,long),
  "EShriR"=rep(NA,long),
  "ERepR"=rep(NA,long),
  "ECloR"=rep(NA,long),
  
  "Reactivity2" = rep(NA, long),
  "SSurvRs"=rep(NA,long),
  "SGrowRs"=rep(NA,long),
  "SShriRs"=rep(NA,long),
  "SRepRs"=rep(NA,long),
  "SCloRs"=rep(NA,long),
  
  "ESurvRs"=rep(NA,long),
  "EGrowRs"=rep(NA,long),
  "EShriRs"=rep(NA,long),
  "ERepRs"=rep(NA,long),
  "ECloRs"=rep(NA,long),
  
  #firsttimestepattenuation analyses
  "Attenuation1" = rep(NA, long),
  "SSurvA"=rep(NA,long),
  "SGrowA"=rep(NA,long),
  "SShriA"=rep(NA,long),
  "SRepA"=rep(NA,long),
  "SCloA"=rep(NA,long),
  
  "ESurvA"=rep(NA,long),
  "EGrowA"=rep(NA,long),
  "EShriA"=rep(NA,long),
  "ERepA"=rep(NA,long),
  "ECloA"=rep(NA,long),
  
  "Attenuation2" = rep(NA, long),
  "SSurvAs"=rep(NA,long),
  "SGrowAs"=rep(NA,long),
  "SShriAs"=rep(NA,long),
  "SRepAs"=rep(NA,long),
  "SCloAs"=rep(NA,long),
  
  "ESurvAs"=rep(NA,long),
  "EGrowAs"=rep(NA,long),
  "EShriAs"=rep(NA,long),
  "ERepAs"=rep(NA,long),
  "ECloAs"=rep(NA,long),
  
  #Maximum Amplification analyses
  "MaxAmplification1" = rep(NA, long),
  "SSurvMA"=rep(NA,long),
  "SGrowMA"=rep(NA,long),
  "SShriMA"=rep(NA,long),
  "SRepMA"=rep(NA,long),
  "SCloMA"=rep(NA,long),
  
  "ESurvMA"=rep(NA,long),
  "EGrowMA"=rep(NA,long),
  "EShriMA"=rep(NA,long),
  "ERepMA"=rep(NA,long),
  "ECloMA"=rep(NA,long),
  
  "MaxAmplification2" = rep(NA, long),
  "SSurvMAs"=rep(NA,long),
  "SGrowMAs"=rep(NA,long),
  "SShriMAs"=rep(NA,long),
  "SRepMAs"=rep(NA,long),
  "SCloMAs"=rep(NA,long),
  
  "ESurvMAs"=rep(NA,long),
  "EGrowMAs"=rep(NA,long),
  "EShriMAs"=rep(NA,long),
  "ERepMAs"=rep(NA,long),
  "ECloMAs"=rep(NA,long),
  
  #Maximum Attenuation analyses
  "MaxAttenuation1" = rep(NA, long),
  "SSurvMAtt"=rep(NA,long),
  "SGrowMAtt"=rep(NA,long),
  "SShriMAtt"=rep(NA,long),
  "SRepMAtt"=rep(NA,long),
  "SCloMAtt"=rep(NA,long),
  
  "ESurvMAtt"=rep(NA,long),
  "EGrowMAtt"=rep(NA,long),
  "EShriMAtt"=rep(NA,long),
  "ERepMAtt"=rep(NA,long),
  "ECloMAtt"=rep(NA,long),
  
  "MaxAttenuation2" = rep(NA, long),
  "SSurvMAtts"=rep(NA,long),
  "SGrowMAtts"=rep(NA,long),
  "SShriMAtts"=rep(NA,long),
  "SRepMAtts"=rep(NA,long),
  "SCloMAtts"=rep(NA,long),
  
  "ESurvMAtts"=rep(NA,long),
  "EGrowMAtts"=rep(NA,long),
  "EShriMAtts"=rep(NA,long),
  "ERepMAtts"=rep(NA,long),
  "ECloMAtts"=rep(NA,long),
  
  # other demographic parameters
  "Lambda"=rep(NA,long),
  "Inertia.upper"=rep(NA,long),
  "Inertia.lower"=rep(NA,long),
  "GT" =rep(NA,long),
  "LE" =rep(NA,long),
  "Primitive"=rep(NA,long),
  "Ergodic"=rep(NA,long),
  "Irreducible"=rep(NA,long)) #these will allow for the testing of matrix suitability.

# Now to run a loop filling this output dataframe with extracted values for each matrix.
count = 0
for (i in 1:long){
  if (i <= dim(compadre)[1]) {d=compadre} else {d=comadre} #this allows the loop to switch between compadre and comadre data files depending on where it is in the run through
  if (i == (dim(compadre)[1]+1)) {count = 0}
  count=count+1
  
  #extract the metadata
  try(output[i,c("SpeciesAuthor","SpeciesAccepted","OrganismType","Family","Class","Kingdom","Authors","Journal",
                 "YearPublication","Country","Continent","Ecoregion","Lat","Lon","Alt","StartYear","StartSeason","StartMonth","EndYear",
                 "EndSeason","EndMonth","Population","Treatment","MatrixComposite","Periodicity")]<-
        data.frame(d[count,c("SpeciesAuthor","SpeciesAccepted","OrganismType","Family","Class","Kingdom","Authors","Journal","YearPublication",
                                "Country","Continent","Ecoregion","Lat","Lon","Altitude","MatrixStartYear","MatrixStartSeason","MatrixStartMonth",
                                "MatrixEndYear","MatrixEndSeason","MatrixEndMonth","MatrixPopulation","MatrixTreatment","MatrixComposite","ProjectionInterval")])[,-1])
  
  # determine if the matrices corresponding with the selected population is post- (True) or pre-reproductive (False).
  post <- mpm_has_prop(d$mat[[count]]) 
  
  #Extract the matrices for use in the analyses below.
  try(matU <- matU(d$mat[[count]]))
  try(matU[is.na(matU)] <- 0) #this is also part of each of the function but is good to make sure its being done as a fix.
  try(matF <- matF(d$mat[[count]]))
  try(matF[is.na(matF)] <- 0)
  try(matC <- matC(d$mat[[count]]))
  try(matC[is.na(matC)] <- 0)
  try(matA <- matU + matF + matC)
  
  # Now convert post-reproductive matrices as necessary, and re extract.
  try(mat.store <- convert2pre(matA, matF, matC, matU, post = post))
  try(matU <- mat.store[[2]])
  try(matF <- mat.store[[3]])
  try(matC <- mat.store[[4]])
  try(matA <- mat.store[[1]])
  
  # Now Extract demographic data.
  # matrix dimension
  if (length(mat.store[[1]]) >= 2) {
    try(output$MatrixDimension[i] <- dim(matA)[1])
    
    # Lambda
    try(output$Lambda[i]<-Re(eigen(matA)$values)[1])
    
    # these are to determine the suitability of the matrices to further analysis
    # Primitivity
    try(output$Primitive[i]<-isPrimitive(matA)) #these are used to test the suitablity of demographic properties.
    
    # Ergodicity
    try(output$Ergodic[i]<-isErgodic(matA))
    
    # Irreducibility
    try(output$Irreducible[i]<-isIrreducible(matA))
    
    # the next few lines of code are to correct any inappropriately calculated values. These lines have been adapted from Robs code.
    try(if (output$Irreducible[i]==TRUE) { #the following analyses are only possible on Irreducible matrices
      
      #Population inertia
      try(output$Inertia.upper[i]<-inertia(matA, bound = "upper"))
      try(output$Inertia.lower[i]<-inertia(matA, bound = "lower"))
      
      # Damping ratio
      try(output[i,c("Rho1","SSurvDR","SGrowDR","SShriDR","SRepDR","SCloDR","ESurvDR","EGrowDR","EShriDR","ERepDR","ECloDR")] <- 
            vitalRatePerturbationDR(matU,matF,matC))
      
      try(output[i,c("Rho2","SSurvDRs","SGrowDRs","SShriDRs","SRepDRs","SCloDRs","ESurvDRs","EGrowDRs","EShriDRs","ERepDRs","ECloDRs")] <- 
            S.vitalRatePerturbationDR(matU,matF,matC))
      
      # Period of oscillation
      try(output[i,c("Pi1","SSurvPO","SGrowPO","SShriPO","SRepPO","SCloPO","ESurvPO","EGrowPO","EShriPO","ERepPO","ECloPO")] <- 
            vitalRatePerturbationPO(matU,matF,matC))
      
      try(output[i,c("Pi2","SSurvPOs","SGrowPOs","SShriPOs","SRepPOs","SCloPOs","ESurvPOs","EGrowPOs","EShriPOs","ERepPOs","ECloPOs")] <- 
            S.vitalRatePerturbationPO(matU,matF,matC))   
      
      # Reactivity
      try(output[i,c("Reactivity1","SSurvR","SGrowR","SShriR","SRepR","SCloR","ESurvR","EGrowR","EShriR","ERepR","ECloR")] <- 
            vitalRatePerturbationReac(matU,matF,matC))
      
      try(output[i,c("Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs","SCloRs","ESurvRs","EGrowRs","EShriRs","ERepRs","ECloRs")] <- 
            S.vitalRatePerturbationReac(matU,matF,matC))
      
      # First time step attenuation
      try(output[i,c("Attenuation1","SSurvA","SGrowA","SShriA","SRepA","SCloA","ESurvA","EGrowA","EShriA","ERepA","ECloA")] <- 
            vitalRatePerturbationAtten(matU,matF,matC))
      
      try(output[i,c("Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs","SCloAs","ESurvAs","EGrowAs","EShriAs","ERepAs","ECloAs")] <- 
            S.vitalRatePerturbationAtten(matU,matF,matC))
      
      # Maximum Amplification attenuation
      try(output[i,c("MaxAmplification1","SSurvMA","SGrowMA","SShriMA","SRepMA","SCloMA","ESurvMA","EGrowMA","EShriMA","ERepMA","ECloMA")] <- 
            vitalRatePerturbationMaxAmp(matU,matF,matC))
      
      try(output[i,c("MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs","SCloMAs","ESurvMAs","EGrowMAs","EShriMAs","ERepMAs","ECloMAs")] <- 
            S.vitalRatePerturbationMaxAmp(matU,matF,matC))
      
      # Maximum Attenuation
      try(output[i,c("MaxAttenuation1","SSurvMAtt","SGrowMAtt","SShriMAtt","SRepMAtt","SCloMAtt","ESurvMAtt","EGrowMAtt","EShriMAtt","ERepMAtt","ECloMAtt")] <- 
            vitalRatePerturbationMaxAtt(matU,matF,matC))
      
      try(output[i,c("MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts","SCloMAtts","ESurvMAtts","EGrowMAtts","EShriMAtts","ERepMAtts","ECloMAtts")] <- 
            S.vitalRatePerturbationMaxAtt(matU,matF,matC))
      
    } else {
      output[i,c("Reactivity1","SSurvR","SGrowR","SShriR","SRepR","SCloR","ESurvR","EGrowR","EShriR","ERepR","ECloR",
                 "Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs","SCloRs","ESurvRs","EGrowRs","EShriRs","ERepRs","ECloRs",
                 "Attenuation1","SSurvA","SGrowA","SShriA","SRepA","SCloA","ESurvA","EGrowA","EShriA","ERepA","ECloA",
                 "Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs","SCloAs","ESurvAs","EGrowAs","EShriAs","ERepAs","ECloAs",
                 "MaxAmplification1","SSurvMA","SGrowMA","SShriMA","SRepMA","SCloMA","ESurvMA","EGrowMA","EShriMA","ERepMA","ECloMA",
                 "MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs","SCloMAs","ESurvMAs","EGrowMAs","EShriMAs","ERepMAs","ECloMAs",
                 "MaxAttenuation1","SSurvMAtt","SGrowMAtt","SShriMAtt","SRepMAtt","SCloMAtt","ESurvMAtt","EGrowMAtt","EShriMAtt","ERepMAtt","ECloMAtt",
                 "MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts","SCloMAtts","ESurvMAtts","EGrowMAtts","EShriMAtts","ERepMAtts","ECloMAtts"
      )] <- NA
    })
    
    # just a little tidy up.
    try(if (output$Irreducible[i]==FALSE | output$Primitive[i]==FALSE) {output[i,c("Reactivity1","SSurvR","SGrowR","SShriR","SRepR","SCloR","ESurvR","EGrowR","EShriR","ERepR","ECloR",
                                                                                   "Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs","SCloRs","ESurvRs","EGrowRs","EShriRs","ERepRs","ECloRs",
                                                                                   "Attenuation1","SSurvA","SGrowA","SShriA","SRepA","SCloA","ESurvA","EGrowA","EShriA","ERepA","ECloA",
                                                                                   "Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs","SCloAs","ESurvAs","EGrowAs","EShriAs","ERepAs","ECloAs",
                                                                                   "MaxAmplification1","SSurvMA","SGrowMA","SShriMA","SRepMA","SCloMA","ESurvMA","EGrowMA","EShriMA","ERepMA","ECloMA", 
                                                                                   "MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs","SCloMAs","ESurvMAs","EGrowMAs","EShriMAs","ERepMAs","ECloMAs", 
                                                                                   "MaxAttenuation1","SSurvMAtt","SGrowMAtt","SShriMAtt","SRepMAtt","SCloMAtt","ESurvMAtt","EGrowMAtt","EShriMAtt","ERepMAtt","ECloMAtt", 
                                                                                   "MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts","SCloMAtts","ESurvMAtts","EGrowMAtts","EShriMAtts","ERepMAtts","ECloMAtts"
    )] <- NA})
    
    #Zero-ing / NAing values that should not have been calculated
    if (sum(matF)==0){ #if no recruitment was observed the matrix model is not complete and cannot be used.
      output[i,c( "Rho1","SSurvDR","SGrowDR","SShriDR","SRepDR","SCloDR","ESurvDR","EGrowDR","EShriDR","ERepDR","ECloDR",
                  "Rho2","SSurvDRs","SGrowDRs","SShriDRs","SRepDRs","SCloDRs","ESurvDRs","EGrowDRs","EShriDRs","ERepDRs","ECloDRs",
                  "Pi1","SSurvPO","SGrowPO","SShriPO","SRepPO","SCloPO","ESurvPO","EGrowPO","EShriPO","ERepPO","ECloPO",
                  "Pi2","SSurvPOs","SGrowPOs","SShriPOs","SRepPOs","SCloPOs","ESurvPOs","EGrowPOs","EShriPOs","ERepPOs","ECloPOs",
                  "Reactivity1","SSurvR","SGrowR","SShriR","SRepR","SCloR","ESurvR","EGrowR","EShriR","ERepR","ECloR",
                  "Reactivity2","SSurvRs","SGrowRs","SShriRs","SRepRs","SCloRs","ESurvRs","EGrowRs","EShriRs","ERepRs","ECloRs",
                  "Attenuation1","SSurvA","SGrowA","SShriA","SRepA","SCloA","ESurvA","EGrowA","EShriA","ERepA","ECloA",
                  "Attenuation2","SSurvAs","SGrowAs","SShriAs","SRepAs","SCloAs","ESurvAs","EGrowAs","EShriAs","ERepAs","ECloAs",
                  "MaxAmplification1","SSurvMA","SGrowMA","SShriMA","SRepMA","SCloMA","ESurvMA","EGrowMA","EShriMA","ERepMA","ECloMA", 
                  "MaxAmplification2","SSurvMAs","SGrowMAs","SShriMAs","SRepMAs","SCloMAs","ESurvMAs","EGrowMAs","EShriMAs","ERepMAs","ECloMAs", 
                  "MaxAttenuation1","SSurvMAtt","SGrowMAtt","SShriMAtt","SRepMAtt","SCloMAtt","ESurvMAtt","EGrowMAtt","EShriMAtt","ERepMAtt","ECloMAtt", 
                  "MaxAttenuation2","SSurvMAtts","SGrowMAtts","SShriMAtts","SRepMAtts","SCloMAtts","ESurvMAtts","EGrowMAtts","EShriMAtts","ERepMAtts","ECloMAtts", 
                  "Lambda")]<-NA}
  }
  #allows you to observe progress
  print(output[i,c("SpeciesAuthor","Treatment","StartYear","EndYear")])
}

#save the extracted data as a csv file
write.csv(output, file = "RawData.csv")