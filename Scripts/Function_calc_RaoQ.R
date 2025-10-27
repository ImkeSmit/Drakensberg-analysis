#calculate RaoQ, one trait at a time
#This is unweighted RaoQ based on euclidean distance. 
#It uses the mean traits of species
#Figure out: does trait filling make sense?
#Should we use raw traits when we have them in cells and only mean traits to fill gaps?

library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)
library(FD)

calc_RaoQ <- function(mean_traits, abun_matrix) {
  
  traitlist <- c(colnames(mean_traits))
  abun_matrix2 <- as.data.frame(abun_matrix)
  
  abun_long <- abun_matrix2 |> 
    rownames_to_column(var = "cellref") |> 
    pivot_longer(cols = !cellref, names_to = "taxon", values_to = "cover")
  
  traits_long <- mean_traits |> 
    rownames_to_column(var = "taxon") |> 
    pivot_longer(!taxon, names_to = "trait", values_to = "value")
  
  trait_sum <- abun_long |> 
    full_join(traits_long, by = "taxon") |> 
    filter(cover > 0) |>
    mutate(value = if_else(is.na(value), 0, value))  |>
    group_by(cellref, trait) |> 
    mutate(trait_sum = sum(value)) |> 
    ungroup()
  
  
  for(t in 1:length(traitlist)) {
    
    chosen_trait <- mean_traits[, t] # preserve matrix format
    names(chosen_trait)<- row.names(mean_traits)
    
    #Check for communities with zero-sum abundances
    #Get cells that do not have any measurements of chosen_trait
    
    abun_temp <- abun_matrix2
    
    problems <- trait_sum |> 
      filter(trait == traitlist[t]) |> 
      filter(trait_sum == 0)
    
    
    if(nrow(problems) > 0) {
      #remove these cells from the abundance matrix
      abun_temp <- abun_temp[-which(rownames(abun_temp) %in% c(problems$cellref)) , , drop = F]
    }
    
    #identify species that do not occur in any of the remaining cells, and remove them
    abundance_sums <- colSums(abun_temp)
    empty_names <- names(abundance_sums[which(abundance_sums == 0)])
    
    if(length(empty_names) > 0) {
      abun_temp <- abun_temp[, - which(colnames(abun_temp) %in% c(empty_names)), drop = F]
      #also remove these sp from the trait matrix
      chosen_trait <- chosen_trait[-which(names(chosen_trait) %in% c(empty_names)), drop = F]
    } 
    
    #final check that all communities have at least one sp
    zero_sum_comm <- which(specnumber(abun_temp) == 0)
    if(length(zero_sum_comm) > 0) {
      abun_temp <- abun_temp[-zero_sum_comm, ]
    }
    
    
    #calculate RaoQ 
    FD_cells <- dbFD(chosen_trait, abun_temp,
                     w.abun = F, #do not weight RaoQ by abundances
                     corr = "cailliez", 
                     calc.FRic = F, 
                     scale.RaoQ = F, 
                     calc.FGR = F, 
                     calc.FDiv = F, 
                     calc.CWM = F, 
                     messages = F)
    #RAoQ = 0 if there is only one distinct trait value in a cell
    
    if(t==1) {
      
      RaoQ_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                                 RaoQ = FD_cells$RaoQ, 
                                 trait = traitlist[t])
      rownames(RaoQ_results) <- NULL
      
    }else {
      more_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                                 RaoQ = FD_cells$RaoQ,
                                 trait = traitlist[t])
      rownames(more_results) <- NULL
      
      RaoQ_results <- rbind(RaoQ_results, more_results)
    }
    
  } 
  return(RaoQ_results) } 


test <- calc_RaoQ(mean_traits, nullcomm_cells[[1]])
