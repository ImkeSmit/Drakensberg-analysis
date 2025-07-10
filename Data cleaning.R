###Script to clean MIcro-Climb occurrence data####
library(tidyverse)
library(tidylog)
library(readxl)
library(janitor)

####Import all occurrence data####
gg_spring <- read_excel("All_data/raw_occurrence_data/GoldenGate/GoldenGate_grids_2019_09_27.xlsx", 
                        sheet = "Species_data_spring") |> 
  select(!Richness) |> 
  pivot_longer(cols = 4:146, names_to = "taxon", values_to = "cover") |> 
  rename(grid = Grid, column = Column, row = Row) |> 
  mutate(site = "GG", 
         cellref = paste0(site, grid, column, row))

#check that all grids and cells are there:
length(unique(gg_spring$grid)) #only 4
length(unique(gg_spring$cellref)) #640 all cells from the 4 grids present


gg_summer <- read_excel("All_data/raw_occurrence_data/GoldenGate/GoldenGate_grids_2019_09_27.xlsx", 
                        sheet = "Species_data_summer", col_names = as.character(c(1:212))) |> 
  row_to_names(row_number = 2) |> 
  clean_names() |> #makes all colnames lowercase and with underscores
  select(!c(richness, lichen, moss, Vascular_cover, grass_richness, nongrass_richness, grass_cover)) |>
  pivot_longer(cols = 4:205, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "GG", 
         grid = as.numeric(grid), 
         row = as.numeric(row),
         cellref = paste0(site, grid, column, row), 
         taxon = )

#check that all grids and cells are there:
length(unique(gg_summer$grid)) #all 8
length(unique(gg_summer$cellref)) #1280 all cells from the 8 grids present
  


wh <- read_excel("All_data/raw_occurrence_data/Witsieshoek/Data entry 22 Mar 2023.xlsx", sheet = 1)

bk <- read_excel("All_data/raw_occurrence_data/Bokong/BNR_vegetation_survey_data_updatedMarch2023.xlsx", sheet = "Veg_data", 
                 col_names = as.character(c(1:1677)))
