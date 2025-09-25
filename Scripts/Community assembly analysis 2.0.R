#Community assembly analysis 2.0###
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)

drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) 
  

#combine occurrence and trait data
comb <- drak |> 
  full_join(FT, by = c("site", "grid", "cell", "taxon")) |> #do full join so that we can trait fill
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

#Trait filling###



