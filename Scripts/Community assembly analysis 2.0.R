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

traitlist <- c(colnames(mean_traits))
RaoQ_results <- data.frame(cellref = NA, 
                           RaoQ = NA, 
                           trait = NA)

for(trait in 1:length(traitlist)) {
  
  chosen_trait <- mean_traits[, trait]
  names(chosen_trait)<- row.names(mean_traits)

  
  FD_cells <- dbFD(chosen_trait, abun_matrix,
                   w.abun = F, #do not weight RaoQ by abundances
                   corr = "cailliez", 
                   calc.FRic = F, 
                   scale.RaoQ = F, 
                   calc.FGR = F, 
                   calc.FDiv = F, 
                   calc.CWM = F)
  #fails because some communities have zero sum abundances
  #I guess there are sommunities in which none of the sp had chlorphyll measured
  
  if(trait==1) {
    RaoQ_results$cellref <- names(FD_cells$RaoQ)
    RaoQ_results$RaoQ <- FD_cells$RaoQ
    RaoQ_results$trait <- traitlist[trait]
  }else {
    more_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                               RaoQ = FD_cells$RaoQ,
                               trait = traitlist[trait])
    
    RaoQ_results <- rbind(RaoQ_results, more_results)
  }
  
}


FD_cells <- dbFD(mean_traits, abun_matrix,
             w.abun = F, #do not weight RaoQ by abundances
             stand.x = T, #standardise traits to mean 0 and unit variance before doing calc
             corr = "cailliez", 
             calc.FRic = F, 
             scale.RaoQ = F, 
             calc.FGR = F, 
             calc.FDiv = F, 
             calc.CWM = F)
FD_cells <- FD_cells$RaoQ



