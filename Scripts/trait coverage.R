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

#