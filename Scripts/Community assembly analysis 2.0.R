#Community assembly analysis 2.0###
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)
library(traitstrap)
library(FD)

drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) |> 
  pivot_longer(cols = c(Wet_mass_mg, Dry_mass_mg, Chlorophyll_mg_per_m2, Ft, Height_cm, 
                        Thickness_mm, Leaf_area_mm2, SLA, LDMC), names_to = "trait", values_to = "value")
  

#combine occurrence and trait data
FT_join <- drak |> 
  inner_join(FT, by = c("site", "grid", "cell", "taxon")) |> #do full join so that we can trait fill
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

#Trait filling, leave for now
#comb <- trait_fill(
#  comm = drak,
#  traits = FT,
#  scale_hierarchy = c("site", "grid"),
#  global = F,
#  taxon_col = "taxon",
#  trait_col = "trait",
#  value_col = "value",
#  abundance_col = "cover",
#  keep_all = FALSE,
#  min_n_in_sample = 5,
#  complete_only = FALSE
#)


###Now we can compute RaoQ at the cell, grid, and site levels###
##Get mean traits for species
mean_traits <- FT_join |> 
  group_by(taxon, trait) |> 
  summarise(mean_trait = mean(value, na.rm = T)) |> 
  pivot_wider(names_from = trait, values_from = mean_trait) |> 
  ungroup() |> 
  arrange(taxon)

mean_traits <- as.data.frame(mean_traits)
row.names(mean_traits) <- mean_traits$taxon
mean_traits <- mean_traits[, -1]

#create abundance matrix
abun_matrix <- FT_join |> 
  select(cellref, taxon, cover) |> 
  distinct(cellref, taxon, .keep_all = T) |> 
  ungroup() |> 
  arrange(taxon) |> 
  mutate(cover = ceiling(cover)) |> #change cover values to integer to use in null models
  pivot_wider(names_from = taxon, values_from = cover) 

abun_matrix <- as.data.frame(abun_matrix)
row.names(abun_matrix) <- abun_matrix$cellref
abun_matrix <- abun_matrix[, -1]

#replace NA values with 0
for(r in 1:nrow(abun_matrix)) {
  for(c in 1:ncol(abun_matrix)) {
    
    if(is.na(abun_matrix[r,c])) {
      abun_matrix[r,c] <- 0
    }
  }
}



#calculate RaoQ for cells, one trait at a time
#This is unweighted RaoQ based on euclidean distance. 
#It uses the mean traits of species
#Figure out: does trait filling make sense?
#Should we use raw traits when we have them in cells and only mean traits to fill gaps?

calc_RaoQ <- function(mean_traits, abun_matrix, FT_join) {

traitlist <- c(colnames(mean_traits))

for(t in 1:length(traitlist)) {
  
  chosen_trait <- mean_traits[, t]
  names(chosen_trait)<- row.names(mean_traits)
  
  abun_matrix2 <- abun_matrix
  
  #Check for communities with zero-sum abundances
  #Get cells from FT_join that do not have any measurements of chosen_trait
  problems <- FT_join |> 
    dplyr::filter(trait == traitlist[t]) |> 
    mutate(value = if_else(is.na(value), 0, value))  |> 
    group_by(cellref) |> 
    mutate(sum_chlor = sum(value)) |> 
    ungroup() |> 
    filter(sum_chlor == 0) |> 
    distinct(cellref)
  
  if(nrow(problems) > 0) {
  #remove these cells from the abundance matrix
  abun_matrix2 <- abun_matrix2[-which(rownames(abun_matrix2) %in% c(problems$cellref)) , ]
  }
  
  #identify species that do not occur in any of the remaining cells, and remove them
  abundance_sums <- colSums(abun_matrix2)
  empty_names <- names(abundance_sums[which(abundance_sums == 0)])
  
  if(length(empty_names) > 0) {
  abun_matrix2 <- abun_matrix2[, - which(colnames(abun_matrix2) %in% c(empty_names))]
  #also remove these sp from the trait matrix
  chosen_trait <- chosen_trait[-which(names(chosen_trait) %in% c(empty_names))]
  } 
 
  #calculate RaoQ 
  FD_cells <- dbFD(chosen_trait, abun_matrix2,
                   w.abun = F, #do not weight RaoQ by abundances
                   corr = "cailliez", 
                   calc.FRic = F, 
                   scale.RaoQ = F, 
                   calc.FGR = F, 
                   calc.FDiv = F, 
                   calc.CWM = F)
  
  if(t==1) {
    
    RaoQ_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                               RaoQ = FD_cells$RaoQ, 
                               trait = traitlist[t])
  
  }else {
    more_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                               RaoQ = FD_cells$RaoQ,
                               trait = traitlist[t])
    
    RaoQ_results <- rbind(RaoQ_results, more_results)
  }

} 
return(RaoQ_results) }


RQ_obs <- calc_RaoQ(mean_traits, abun_matrix, FT_join)


###Nullmodels####
generate_C3_null <- function(comm, traits) {
  null_comm <- comm * 0  # initialize matrix
  
  for (i in 1:nrow(comm)) {
    site <- comm[i, ]
    
    # Species present at site (richness)
    richness <- sum(site > 0)
    
    # Randomly choose species (without replacement) to occupy this site
    chosen_species <- sample(colnames(comm), richness, replace = FALSE)
    
    #For C3, we can pick abundances from sites in the matrix where the sp occurs
    for (s in 1:length(chosen_species)) {
      possible_abundances <-  comm |> 
        select(all_of(chosen_species)) |> 
        pull(chosen_species[s]) |> 
        discard(~ .x <= 0)
      
      #abundances that are more frequent are more likely to be sampled
      chosen_abundance <- sample(possible_abundances, 1)
      
      #assign abundance to species and site
      null_comm[i, chosen_species[s]] <- chosen_abundance
    } #end loop through species
  }#end loop through sites
  
  # Now shuffle traits within abundance classes
  traits_null <- traits
  for (cl in unique(traits$ab_class)) {
    trait_values <- traits$height[traits$ab_class == cl]
    traits_null$height[traits$ab_class == cl] <- sample(trait_values)
  }
  
  return(list(comm = null_comm, traits = traits_null))
}





