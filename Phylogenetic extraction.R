# This script is for extracting phylogenetic information for a series of selected populations
# Date last modified: December 2020
# This will need running before any further analysis comparing between populations.

# clear workspace
rm(list=ls(all=TRUE))

#set working directory
setwd("FILE_PATH")

# load in data file
data <- read.csv("CSVFILE.csv")

# load required packages
library(ape)
library(caper)
library(diversitree)
library(geiger)
library(phytools)
library(RColorBrewer)
library(rotl) 
library(taxize)
library(stringr)
devtools::install_github("GuangchuangYu/ggtree")
library(ggtree)
library(ggplot2)
library(ggimage)
# load function for populating a tree with multiple populations of the same species
source("FILE_PATH/Tree tricker function.R")

# Sort the data for extraction ------------------------------------------------------------------------------------
# first I need to fix the species names in the overall dataframe
data$SpeciesAccepted <- sub("_", " ", data$SpeciesAccepted) #this will remove any pesky underscores in species names

# Now I can take a list of every species that I am working with for this I will only need unique species
SpeciesList <- data.frame(Type = data$OrganismType[which(!duplicated(data$SpeciesAccepted))], 
                          Species = unique(data$SpeciesAccepted))
SpeciesList$Species <- as.character(SpeciesList$Species) # <- taxise requires species lists in character format

# Now I want to follow the Phylogenetic analysis procedure carried out by Pol so I will now check species names in my data frame are taxonmically accepted
# this is done using Taxize.
# first check that there are no spelling/formatting errors in the species list
name.check <- gnr_resolve(names = SpeciesList$Species, canonical = TRUE, best_match_only = TRUE)
name.check$user_supplied_name == name.check$matched_name2
# for the species for which the matched name differs from the supplied name make the edit in the species list.
SpeciesList$Species <- name.check$matched_name2 

# main data frame name edits
for (i in 1:dim(data)[1]){
  for (j in 1:dim(name.check)[1]){
    if (data$SpeciesAccepted[i] == name.check$user_supplied_name[j]) {data$SpeciesAccepted[i] <- name.check$matched_name2[j]}
  }
}

# Step 1: Match the species to its ott_id ----------------------------------------------------------------------------

# firstly the species names need to match their ott_id values
# the matching names function only works on small strings of data so I need break my species list down into more managable chunks.

#define organism groups to split the species list into.
taxa_types <- unique(SpeciesList$Type); print(taxa_types)
#Animals
Animals <- taxa_types[c(3,4,5,7,13,14,15,18,19,21)]
#Plants
Plants <- taxa_types[c(1,6,9,10,11,12,16,17,22)] #combined the plant group is still too big!
Herbs <- taxa_types[c(2,20)]
#Algae
Algae <- taxa_types[8]
# and split the species list
Animals <- SpeciesList[which(SpeciesList$Type %in% Animals),]
Plants <- SpeciesList[which(SpeciesList$Type %in% Plants),]
Herbs <- SpeciesList[which(SpeciesList$Type %in% Herbs),]
Algae <- SpeciesList[which(SpeciesList$Type %in% Algae),]
# just a quick check to make sure it has worked and no species have been missed!
dim(Algae)[1]+dim(Animals)[1]+dim(Plants)[1]+dim(Herbs)[1] == dim(SpeciesList)[1]

# Now I can determine the ott_id values
resolved_names_animals <- tnrs_match_names(Animals$Species) # this fails to find Hystrix refossa - so I will research atr a higher level
Animals[which(Animals$Species == "Hystrix refossa"),]$Species <- "Hystrix"
resolved_names_animals <- tnrs_match_names(Animals$Species) # all works now
resolved_names_plants <- tnrs_match_names(Plants$Species) # error
resolved_names_herbs <- tnrs_match_names(Herbs$Species) # error
resolved_names_algae <- tnrs_match_names(Algae$Species)
# trying to determine the ott_id values here returns an error for both plants and herbs due to the inability of the function to find some of the subspecies names.
# therefore I need to correct these names and re-run the function.
Plants[which(Plants$Species %in% c("Geonoma pohliana weddelliana",
                                   "Vella pseudocytisus paui",
                                   "Magnolia macrophylla dealbata",
                                   "Quercus mongolica crispula" )),]$Species <- c("Geonoma pohliana", # weddelliana
                                                                                  "Vella pseudocytisus", #paui
                                                                                  "Magnolia macrophylla", # dealbata
                                                                                  "Quercus mongolica") # crispula 
Herbs[which(Herbs$Species %in% c("Anthyllis vulneraria alpicola",
                                 "Arenaria grandiflora bolosii",
                                 "Echinospartum ibericum algibicum",
                                 "Eriogonum longifolium gnaphalifolium",
                                 "Pityopsis aspera aspera",
                                 "Silene douglasii oraria",
                                 "Viola sagittata ovata")),]$Species <- c("Anthyllis vulneraria", #alpicola
                                                                                  "Arenaria grandiflora", #bolosii
                                                                                  "Echinospartum ibericum", #algibicum
                                                                                  "Eriogonum longifolium", #gnaphalifolium
                                                                                  "Pityopsis aspera", #aspera
                                                                                  "Silene douglasii", #oraria
                                                                                  "Viola sagittata") #ovata

# and rerun the resolved name search for Plants and herbs
resolved_names_plants <- tnrs_match_names(Plants$Species) 
resolved_names_herbs <- tnrs_match_names(Herbs$Species)

# I will now insert found ott_id values into the taxa seperated species lists
Animals$ott_id <- NA; Plants$ott_id <- NA; Herbs$ott_id <- NA; Algae$ott_id <- NA
# but before I can do this I need to account for the removal of one species (due to duplication) from both data sets
# in the plants dataframe it is row 424 that is to be removed and from the herbs it is 31.
Plants <- Plants[which(rownames(Plants) != 424),]
Herbs <- Herbs[which(rownames(Herbs) != 31),]
# now add the ott_id values to the main species data frames
Animals$ott_id <- resolved_names_animals$ott_id; Plants$ott_id <- resolved_names_plants$ott_id; Herbs$ott_id <- resolved_names_herbs$ott_id; Algae$ott_id <- resolved_names_algae$ott_id

# on previous efforts the species Lagopus muta japonica in the animal subset throws back an error so the following code a fix for that if needed.
resolved_extra_animals <- tnrs_match_names("Lagopus muta")
Animals[which(Animals$Species == "Lagopus muta japonica"),]$ott_id <- resolved_extra_animals$ott_id 

# I also need to correct for any changes in species id that occured when downloading ott_id values
Animals$Species_updated <- resolved_names_animals$unique_name; Plants$Species_updated <- resolved_names_plants$unique_name
Herbs$Species_updated <- resolved_names_herbs$unique_name; Algae$Species_updated <- resolved_names_algae$unique_name;

# and add in the missing entries for plants, animals and herbs if needed
Animals[which(Animals$Species == "Lagopus muta japonica"),]$Species_updated <- "Lagopus muta"

# Now all species names need to be corrected in the orginal file.
full.list <- rbind(Animals, Plants, Herbs, Algae)
# loop through each entry and change if nessecary
for (i in 1:dim(data)[1]){
  for (j in 1:dim(full.list)[1]){
    if (data$SpeciesAccepted[i] == full.list$Species[j]) {data$SpeciesAccepted[i] <- full.list$Species_updated[j]}
  }
}
# this loop misses the species where I had to manually remove the subspecies. So I will fix that here.
#plants
data[which(data$SpeciesAccepted == "Geonoma pohliana weddelliana"), "SpeciesAccepted"] <- "Geonoma pohliana"
data[which(data$SpeciesAccepted == "Vella pseudocytisus paui"), "SpeciesAccepted"] <- "Vella pseudocytisus"
data[which(data$SpeciesAccepted == "Magnolia macrophylla dealbata"), "SpeciesAccepted"] <- "Magnolia macrophylla"
data[which(data$SpeciesAccepted == "Quercus mongolica crispula"), "SpeciesAccepted"] <- "Quercus mongolica"
# herbs                                                                           
data[which(data$SpeciesAccepted == "Anthyllis vulneraria alpicola"), "SpeciesAccepted"] <- "Anthyllis vulneraria"
data[which(data$SpeciesAccepted == "Arenaria grandiflora bolosii"), "SpeciesAccepted"] <- "Arenaria grandiflora"
data[which(data$SpeciesAccepted == "Echinospartum ibericum algibicum"), "SpeciesAccepted"] <- "Echinospartum ibericum"
data[which(data$SpeciesAccepted == "Eriogonum longifolium gnaphalifolium"), "SpeciesAccepted"] <- "Eriogonum longifolium"
data[which(data$SpeciesAccepted == "Pityopsis aspera aspera"), "SpeciesAccepted"] <- "Pityopsis aspera"
data[which(data$SpeciesAccepted == "Silene douglasii oraria"), "SpeciesAccepted"] <- "Silene douglasii"
data[which(data$SpeciesAccepted == "Viola sagittata ovata"), "SpeciesAccepted"] <- "Viola sagittata"

# animals (if needed)
data[which(data$SpeciesAccepted == "Lagopus muta japonica"), "SpeciesAccepted"] <- "Lagopus muta"

# Step 2: Get the tree corresponding to my taxa -------------------------------------------------------------------

# for each of the resolved tree groups I can extract subtrees (this prevents overloading the program)
# Animal subtree
animal_tree <- tol_induced_subtree(ott_ids = Animals$ott_id)
# Algae subtree
algae_tree <- tol_induced_subtree(ott_ids = Algae$ott_id)
# Plant subtree
plants_overall <- rbind(Plants, Herbs) #combined together plants and herbs as these would all be in the same subtree
plant_tree <- tol_induced_subtree(ott_ids = plants_overall$ott_id)
# five of the ott_ids found for a plant population are not present within the synthetic tree list and so is throwing an error.
print(plants_overall[which(plants_overall$ott_id %in% c("16952","3904575","5144270","5525303","5738006")),]$Species_updated)
# these need dropping from the analysis
plants_overall <- plants_overall[which(plants_overall$Species_updated != "Asplenium x adulterinum"),]
plants_overall <- plants_overall[which(plants_overall$Species_updated != "Horkelia congesta"),]
plants_overall <- plants_overall[which(plants_overall$Species_updated != "Pilosella x floribunda"),]
plants_overall <- plants_overall[which(plants_overall$Species_updated != "Polemonium van-bruntiae"),]
plants_overall <- plants_overall[which(plants_overall$Species_updated != "Rubus praecox"),]
# re-create the tree
plant_tree <- tol_induced_subtree(ott_ids = plants_overall$ott_id) #this works now! 
# also remove the error species from the overall data file
data <- data[which(data$SpeciesAccepted != c("Asplenium x adulterinum", "Horkelia congesta", "Pilosella x floribunda",
                                             "Polemonium van-bruntiae", "Rubus praecox")),]

# Now I need to check if these trees are rooted.
is.rooted(plant_tree); plot(plant_tree)
is.rooted(animal_tree); plot(animal_tree)
is.rooted(algae_tree); plot(algae_tree)
# to improve the accepted phylogeny of this tree I need to ensure the demospongia as the outgroup
# which species are in the sponge family?
print(Animals[which(Animals$Type == "Demospongiae"),]$Species)
# from looking at the structure of the animal tree it is rooted with sponges as the outgroup.

# Now I can bind the subtrees.
combine_tree <- bind.tree(x = algae_tree, y = plant_tree, where = "root")
# determine the node on the animals tree for where to bind the subtrees.
plot(animal_tree, show.node.label = TRUE)
combine_tree <- bind.tree(x = animal_tree, y = combine_tree, where = "root") #this will ensure the plant and animal trees out bound with the sponge being to closest living relative

# Now I need to check the details of the tree?
# is it rooted?
is.rooted(combine_tree); plotTree(combine_tree, fsize = 0.5)
# try rooting it with the sponges as the outgroup as Pol did with his.
combine_tree <- root(combine_tree, outgroup = c("Xestospongia_muta_ott880230","Amphimedon_compressa_ott981298"))
# now re-check if it is rooted
is.rooted(combine_tree); plotTree(combine_tree, fsize = 0.5)
# is it binary?
is.binary(combine_tree) # there is lots of polytomies. The function multi2di can fix those
combine_tree <- multi2di(combine_tree, random = TRUE)
is.binary(combine_tree) # that has worked!
# is the the tree Ultrameteric
is.ultrametric(combine_tree) # no because the tree has no branch lengths.
# add the branch lengths and a time dimension to the tree
# firstly are there any duplicated node labels as these need removing 
any(duplicated(combine_tree$node.label))
combine_tree <- makeNodeLabel(combine_tree)
# now I can calculate branch lengths
phylo_tree <- compute.brlen(combine_tree, method = "Grafen")
phylo_tree <- chronos(phylo_tree, model = "correlated")
# Now check that it is all good!
is.rooted(phylo_tree); is.binary(phylo_tree); is.ultrametric(phylo_tree)
phylo_tree; plotTree(phylo_tree, fsize = 0.2)
any(duplicated(phylo_tree$node.label))
# everything is working!

# Now I can clean the Species list based on the tree I have
#first format the names in the tree
phylo_tree$tip.label[79] <- "Phyllanthus_emblica_ott399352" # a little fix to align the notation of all species names
phylo_tree$tip.label <- sub("_ott", " ", phylo_tree$tip.label)
phylo_tree$tip.label <- word(phylo_tree$tip.label, 1)
phylo_tree$tip.label <- sub("_", " ", phylo_tree$tip.label)
phylo_tree$tip.label <- sub("_", " ", phylo_tree$tip.label) # this needs repeating for any subspecies

# ---------------------------------------------------------------
# Step 3 - now for the moment of truth - adding to the tree any replicate populations for each species
# first add a column to store the populations 'new' names (those that match what the tree is using)
data$SpeciesTree <- NA

# check the species compatability between the tree and the data set
name.check(phylo_tree, data, data.names = data$SpeciesAccepted)
# Just need to address the odd formatting of the species Phyllanthus_emblica and Hystrix refossa
data[which(data$SpeciesAccepted == "Phyllanthus emblica (species in kingdom Archaeplastida)"), "SpeciesAccepted"] <- "Phyllanthus emblica"
phylo_tree$tip.label[2] <- "Hystrix"
# re-run the check
name.check(phylo_tree, data, data.names = data$SpeciesAccepted)
# these last discrepencies need dropping from the data as clearly the tree of life cannot find them
# first tidy the tree
check_species <- name.check(phylo_tree, data, data.names = data$SpeciesAccepted)
phylo_tree <- drop.tip(phylo_tree, check_species$tree_not_data)
# and tidy the dataset
data <- data[which(data$SpeciesAccepted %in% phylo_tree$tip.label),]
# and run one final check
name.check(phylo_tree, data, data.names = data$SpeciesAccepted) # all good to add in duplicate populations

# and run the data and tree through the function
tricked_tree <- tricking_tree(data, phylo_tree)

# split the output back to a data frame and phylogenetic tree
data <- tricked_tree[[2]]
phylo_tree <- tricked_tree[[1]]
# did it work???
name.check(phylo_tree, data, data.names = data$SpeciesTree)
# just a couple of data names need changing to match the tree tip labels
data$SpeciesTree <- sub(" ", "_", data$SpeciesTree)
# and recheck
name.check(phylo_tree, data, data.names = data$SpeciesTree)
# all fixed

# Now just check the tree still works
is.rooted(phylo_tree); is.binary(phylo_tree); is.ultrametric(phylo_tree)
phylo_tree; plotTree(phylo_tree, fsize = 0.1)
any(duplicated(phylo_tree$node.label)) # this needs correcting
phylo_tree <- makeNodeLabel(phylo_tree)
# everything is working!

# and export adjusted data file and phylogenetic tree --------------------------------------------
write.csv(data, file = "CSVFILE2.csv") # Checkpoint 
# this is now a datafile of immediate and legacy environmental variability along with the transient dynamics of associated populations.
# I have also incorperated information regarding the phylogenies of each population (this is the final set up data set and is ready for analysis)
write.tree(phylo_tree, file = "Phylogenetic subtree.tre")

#Plotting time -----------------------------------------------------------------------------------

# Now I can plot the phylogeny
# first determine kingdom taxonomy
SpeciesList <- data[,c("SpeciesTree", "OrganismType")]
SpeciesList$kingdom <- NA
SpeciesList[which(SpeciesList$OrganismType %in% taxa_types[c(3,4,5,7,13,14,15,18,19,21)]),]$kingdom <- "Animalia"
SpeciesList[which(SpeciesList$OrganismType %in% taxa_types[c(1,2,6,9,10,11,12,16,17,20,22)]),]$kingdom <- "Plantae"
SpeciesList[which(SpeciesList$OrganismType %in% taxa_types[8]),]$kingdom <- "Protista"
SpeciesList$kingdom <- as.factor(SpeciesList$kingdom)

# just check everything matches species wise!
rownames(SpeciesList) <- SpeciesList$SpeciesTree
check.tree <- name.check(phylo_tree, SpeciesList)
phylo_tree <- drop.tip(phylo_tree, check.tree$tree_not_data)
name.check(phylo_tree, SpeciesList)
# add define the tip.label categories
tip.label.categories <- SpeciesList[,c(1,3)]

# now create the tree
tree <- ggtree(phylo_tree, layout = "circular", size = 0.7)
# extract animal silhouettes
phylopic_info <- ggimage::phylopic_uid(sub("_", " ", phylo_tree$tip.label[which(phylo_tree$tip.label %in% c("Giraffa_camelopardalis", "Geranium_sylvaticum", "Cirsium_acaule", "Ascophyllum_nodosum"))]))
# reformat the names to match the tree tip labels
rownames(phylopic_info) <- phylopic_info$name <- c("Geranium_sylvaticum", "Cirsium_acaule", "Giraffa_camelopardalis", "Ascophyllum_nodosum") 

# add the tip label categories
tree <- tree %<+% tip.label.categories +
  geom_tippoint(aes(col = kingdom), size = 6, shape = 20) +
  scale_color_manual(values = c("ivory2","darkred","grey63")) +
  labs(colour = NULL) #Removes legend title.

# add animal silhouettes
tree <- tree %<+% phylopic_info +
  geom_tiplab(aes(image = uid), geom = "phylopic", offset = 0.15)

# ******************************************************** End of Code *************************
