###Script to clean MIcro-Climb occurrence data####
library(tidyverse)
library(tidylog)
library(readxl)
library(janitor)
library(ggplot2)
library(openxlsx)

####GOLDEN GATE####
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
  select(!c(richness, lichen, moss, vascular_cover, grass_richness, nongrass_richness, grass_cover)) |>
  pivot_longer(cols = 4:205, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "GG", 
         grid = as.numeric(grid), 
         row = as.numeric(row),
         cellref = paste0(site, grid, column, row), 
         cover = as.numeric(cover)) |> 
  filter(!cover == 0)

#check that all grids and cells are there:
length(unique(gg_summer$grid)) #all 8
length(unique(gg_summer$cellref)) #1280 all cells from the 8 grids present

#check that cover values make sense
ggplot(gg_summer) +
  geom_histogram(aes(x = cover))
min(gg_summer$cover)
max(gg_summer$cover)

####WITSIESHOEK####
wh <- read_excel("All_data/raw_occurrence_data/Witsieshoek/Data entry 22 Mar 2023.xlsx", sheet = 1) |> 
  clean_names() |> 
  select(!species_richness) |> 
  pivot_longer(cols= 4:187, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "WH", 
         cellref = paste0(site, grid, column, row)) |> 
  filter(!cover == 0)

#check that all grids and cells are there:
length(unique(wh$grid)) #all 7
length(unique(wh$cellref)) #1120 all cells from the 8 grids present

#check that cover values make sense
ggplot(wh) +
  geom_histogram(aes(x = cover))
min(wh$cover)
max(wh$cover)


####BOKONG####
bk <- read_excel("All_data/raw_occurrence_data/Bokong/BNR_vegetation_survey_data_updatedMarch2023.xlsx", sheet = "Veg_data", 
                 col_names = as.character(c(1:1677))) |>
  select(!2) |> 
  t() |> 
  row_to_names(row_number = 1) |> 
  clean_names() |> 
  as_tibble() |> 
  select(!date) |> 
  pivot_longer(cols= 4:103, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "BK", 
         cellref = paste0(site, grid, column, row), 
         cover = str_replace_all(cover, ",", ".")) |> 
  filter(!is.na(cover)) |> 
  mutate(cover2 = as.numeric(cover), 
         grid = as.numeric(grid), 
         row = as.numeric(row))

#there are still NA's in the cover column, let's sort them out
bk[which(is.na(bk$cover2)) , ] #get Na rows
#let's just give all these a value of 0.5
bk[which(is.na(bk$cover2)) , which(colnames(bk) == "cover2")] <- 0.5

#delete the character cover column and rename cover2
bk <- bk |> 
  select(!cover) |> 
  rename(cover = cover2)

#check that cover values make sense
ggplot(bk) +
  geom_histogram(aes(x = cover))
min(bk$cover)
max(bk$cover)
unique(bk$cover) #there are a few problems here : 0, 0.2 and 1510.0

bk <- bk[-which(bk$cover == 0), ] #I assume cover = 0 means the species was absent, remove it
bk[which(bk$cover == 0.2), which(colnames(bk) == "cover")] <- 0.5 #felicia rosulata, give cover of 0.5
bk[which(bk$cover == 1510.0), which(colnames(bk) == "cover")] <-15  #tenaxia disticha, can be big grass, give cover of 15

####Bind sites together and extract speciesnames####
allsites <- gg_summer |> 
  bind_rows(wh) |> 
  bind_rows(bk) |> 
  mutate(dataset = "veg_survey")

all_sp <- allsites |> 
  distinct(dataset, site, taxon) #get all unique site and taxon combinations
#export to excel to create name key
write.xlsx(all_sp, "All_data/clean_data/micro_climb_survey_names.xlsx", sheetName = "veg_survey_names")

#checking some species:
gg_summer[grep("argyrolobium", gg_summer$taxon), ]
test <- gg_summer[grep("asclep", gg_summer$taxon), ]
gg_summer[grep("gladiolus", gg_summer$taxon), ]


####Import trait data####
#for now, we'll just get the names from the trait data, and clean the rest later

gg_trait_names <- read_excel("All_data/raw_trait_data/GG_dataset_Functional_traits.xlsx") |> 
  rename(taxon = Species) |> 
  mutate(taxon = str_to_lower(taxon),
         taxon = str_squish(taxon),
         taxon = str_replace_all(taxon, " ", "_")) |> 
  select(taxon) |> 
  mutate(site = "GG", dataset = "FT") |> 
  distinct(taxon, site, dataset)

wh_trait_names <-read_excel("All_data/raw_trait_data/WH_Functional_traits_dataset.xlsx", sheet = "Data_entry") |> 
  select(Date:Notes) |> 
  rename(taxon = Species) |> 
  mutate(taxon = str_to_lower(taxon),
         taxon = str_squish(taxon),
         taxon = str_replace_all(taxon, " ", "_")) |> 
  select(taxon) |> 
  mutate(site = "WH", dataset = "FT") |> 
  distinct(taxon, site, dataset)

bk_trait_names <- read_excel("All_data/raw_trait_data/FT_Bokong_19Nov.xlsx", sheet = "FT measurements") |> 
  rename(taxon = Species) |> 
  mutate(taxon = str_to_lower(taxon),
         taxon = str_squish(taxon),
         taxon = str_replace_all(taxon, " ", "_")) |> 
  select(taxon) |> 
  mutate(site = "BK", dataset = "FT") |> 
  distinct(taxon, site, dataset)

all_trait_sp <- gg_trait_names |> 
  bind_rows(wh_trait_names) |> 
  bind_rows(bk_trait_names) 

#create a trait names file
write.xlsx(all_trait_sp, "All_data/clean_data/micro_climb_trait_names.xlsx", sheetName = "trait_names", overwrite = F)

