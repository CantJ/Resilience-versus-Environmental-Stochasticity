# This script sets up a function for incorporating multiple populations of the same species into a phylogenetic tree. 
# Original author: Thomas Merrien (thomas.merrien@agroparistech.fr)
# Script modified by James cant to fit personal set up and coding style.
# Date last modified: June 2020


# ------------------------------------------------------------------------------

# This function allows for the incorporation of multiple populations of a same species in to a tree, separated by a small distance
# called epsilon. This distance has to be very small but it should not be too small in order to avoid rounding values
# leading to an ultrametric tree.
# In this version, the new population is always added on top of the population already in the tree or on top of
# the one inserted before if they are coming from the same species.
# The inputs of the function are a data frame with multiple replicates of the same species with species names in a 
# column labelled SpeciesAccepted and a classic phylogenetic tree involving only unique species. The output of the function is a phylogenetic tree
# with species replicates

# -------------------------------------------------------------------------------

# create function
tricking_tree<-function(df, tree){
  # set the distance that will be used between populations of the same species
  epsilon=0.0000001 
  
  for (i in 1:length(unique(df$SpeciesAccepted))){ #For each unique species
    
    nb_pop <- length(grep(unique(df$SpeciesAccepted)[i], df$SpeciesAccepted)) #find the number of replicates
    
    cat(" Working on", unique(df$SpeciesAccepted)[i], "which has", nb_pop,"population replicates.\n") # print a progress read out
    
    if ((nb_pop > 1) == TRUE){ #if there is more than one population
      
      if ((nb_pop > 2) == TRUE){
        
        order <- sample(2:(nb_pop)) #create a random vector with the order of the population
        
      } else {
        
        order <- 2
        
      }
      
      pop_names <- paste(unique(df$SpeciesAccepted)[i], order, sep = "_") #create a list of populations along with their population identifier
      
      pop_names <- c(paste(unique(df$SpeciesAccepted)[i]), pop_names) #add to the list the original species already in the tree.
      
      # use the population list to assign population names to the raw data file.
      pop_index <- rownames(df[which(df$SpeciesAccepted == unique(df$SpeciesAccepted)[i]),]) #identify where in the data file the populations in question are being used.
      df[pop_index, "SpeciesTree"] <- pop_names
      
      # format the names in the population vector to match the reformatted tree
      pop_names <- sub(" ", "_", pop_names); tree$tip.label <- sub(" ", "_", tree$tip.label)
      # a little fix to ensure that populations identified to the sub-species level stay in the same format.
      # As the sub function will only change one space at a time.
      pop_names <- sub(" ", "_", pop_names); tree$tip.label <- sub(" ", "_", tree$tip.label)
      
      for (j in 1:length(order)){ #incorporate the replicated populations into the phylogenetic tree.
          
          tree<-bind.tip(tree, pop_names[j+1], edge.length = order[j] * epsilon, 
                         where=which(tree$tip.label==pop_names[j])) #Add each new population into the tree
          
          # a little fix to ensure that populations identified to the sub-species level stay in the same format. 
          
          index <- mrca(tree)[pop_names[j+1], pop_names[j]] #last nodes between two populations - provide the function with two species (in this case populations) and it will return the last common ansceter (node value)
          
          tree$edge.length[which(tree$edge[,1] == index)] <- (nb_pop - j) * epsilon #Set the right length for the new branch (all sub-populations will be equally diverged)
      
      }
    
    } else {
        
      pop_index <- rownames(df[which(df$SpeciesAccepted == unique(df$SpeciesAccepted)[i]),]) #identify where in the data file the populations in question are being used.
      df[pop_index, "SpeciesTree"] <- unique(df$SpeciesAccepted)[i]
      
    }
    
    cat("Done!\n") # completion read out
    
  }
  tree<-force.ultrametric(tree) #To ensure that the tree remains ultrametric after all the additions
  df[, "SpeciesTree"] <- sub(" ", "_", df[, "SpeciesTree"]) # reformat the species names to be used as rownames to match those in the phylogenetic tree.
  
  # End function
  return(c(tree, df))
}

# --------------------------------------------------- end of function --------------------------------------------