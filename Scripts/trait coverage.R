###Assess trait coverage of Micro-climb data####

library(tidyverse)
library(openxlsx)
library(tidylog)

FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  mutate(cellref = paste0(Site,Grid,Cell)) |> 
  rename(taxon = Taxon)

occ <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx")

comb <- occ |>
  left_join(FT, by = c("cellref", "taxon"))

#how many sp in occ do we have at least one trait value for?
nsp_occ <- occ |> 
  distinct(taxon) |> 
  summarise(n = n()) #414 total sp in occ

nsp_traits <- comb |> 
  filter(!is.na(Sample_ID)) |> #keep only sp for which we have traits
  distinct(taxon) |> 
  summarise(n = n())

p_traits <- nsp_traits/nsp_occ *100 #we have at least one trait val for 35% of species


