###Assess trait coverage of Micro-climb data####

library(tidyverse)
library(openxlsx)
library(tidylog)

FT <- read.csv("All_data/clean_data/micro-climb_traits.csv", row.names = 1) |> 
  rename(taxon = Taxon)

occ <- read.csv("All_data/clean_data/micro_climb_occurrence.csv", row.names = 1)

comb <- occ |>
  left_join(FT, by = c("Cell_ID", "taxon"))

#how many sp in occ do we have at least one trait value for?
nsp_occ <- occ |> 
  distinct(taxon) |> 
  summarise(n = n()) #401 total sp in occ

nsp_traits <- comb |> 
  filter(!is.na(Sample_ID)) |> #keep only sp for which we have traits
  distinct(taxon) |> 
  summarise(n = n()) #159 sp have at least one trait val

p_traits <- nsp_traits/nsp_occ *100 #we have at least one trait val for 39% of species


###How much cover is accounted for by these 146 species?
trait_sp <-  comb |> #names of sp for which we have traits
  filter(!is.na(Sample_ID)) |> #keep only sp for which we have traits
  distinct(taxon)

total_cover <- sum(occ$cover, na.rm = TRUE) #total veg cover in the occ data

trait_cover <- comb |> 
  filter(taxon %in% c(trait_sp$taxon)) |> #remember these are the SPECIES for which we have at least one trait
  summarise(trait_cover = sum(cover))

percent_trait_cover <- trait_cover/total_cover *100 #89%

####Get trait coverage per grid
coverage <- cell_trait_coverage(level = "grid")

