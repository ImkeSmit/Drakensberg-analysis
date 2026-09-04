###Script to clean MIcro-Climb occurrence data####
library(tidyverse)
library(tidylog)
library(readxl)
library(janitor)
library(ggplot2)
library(openxlsx)
library(conflicted)
conflict_prefer_all("tidylog", quiet = TRUE)

####GOLDEN GATE####
gg_summer <- read_excel("All_data/raw_data/raw_occurrence_data/GoldenGate/GoldenGate_grids_2019_09_27.xlsx", 
                        sheet = "Species_data_summer", col_names = as.character(c(1:212))) |> 
  row_to_names(row_number = 2) |> 
  clean_names() |> #makes all colnames lowercase and with underscores
  select(!c(richness, lichen, moss, vascular_cover, grass_richness, nongrass_richness, grass_cover)) |>
  pivot_longer(cols = 4:205, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "GG", 
         grid = as.numeric(grid), 
         row = as.numeric(row),
         cellref = paste0(site, "_", "G", grid, "_", column, row), 
         cover = as.numeric(cover))


##Cells with zero species richness
gg_summer |> 
  group_by(cellref) |> 
  summarise(sum_cover = sum(cover)) |> 
  filter(sum_cover == 0) #19 cells with zero species richness
#these seem to be cells with 100% rock cover

#remove species with no cover from long format data
gg_summer2 <- gg_summer |> 
  filter(!cover == 0)

#check that all grids and cells are there:
length(unique(gg_summer2$grid)) #all 8
length(unique(gg_summer2$cellref)) #19 cells with no species missing

#check that cover values make sense
ggplot(gg_summer2) +
  geom_histogram(aes(x = cover))
min(gg_summer$cover)
max(gg_summer$cover)
class(gg_summer$cover)

####WITSIESHOEK####
wh <- read_excel("All_data/raw_data/raw_occurrence_data/Witsieshoek/Data entry 22 Mar 2023.xlsx", sheet = 1) |> 
  clean_names() |> 
  select(!species_richness) |> 
  pivot_longer(cols= 4:187, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "WH", 
         cellref = paste0(site, "_", "G", grid, "_", column, row)) 

##Cells with zero species richness
wh |> 
  group_by(cellref) |> 
  summarise(sum_cover = sum(cover)) |> 
  filter(sum_cover == 0) #4 cells with zero species richness
#these seem to be cells with 100% rock cover

#remove species with no cover from long format data
wh2 <- wh |> 
  filter(!cover == 0)

#check that all grids and cells are there:
length(unique(wh2$grid)) #all 7
length(unique(wh2$cellref)) #1116 

#check that cover values make sense
ggplot(wh2) +
  geom_histogram(aes(x = cover))
min(wh2$cover)
max(wh2$cover)


####BOKONG####
bk <- read_excel("All_data/raw_data/raw_occurrence_data/Bokong/BNR_vegetation_survey_data_updatedMarch2023.xlsx", sheet = "Veg_data", 
                 col_names = as.character(c(1:1677))) |>
  select(!2) |> 
  t() |> 
  row_to_names(row_number = 1) |> 
  clean_names() |> 
  as_tibble() |> 
  select(!date) |> 
  pivot_longer(cols= 4:103, names_to = "taxon", values_to = "cover") |> 
  mutate(site = "BK", 
         cellref = paste0(site, "_", "G", grid, "_", column, row), 
         cover = str_replace_all(cover, ",", ".")) |> 
  filter(!is.na(cover)) |> 
  mutate(cover2 = as.numeric(cover), 
         grid = as.numeric(grid), 
         row = as.numeric(row))

#there are still NA's in the cover column, let's sort them out
bk[which(is.na(bk$cover2)) , ] #get Na rows
#Go back to paper data sheets and check these values
bk[which(bk$cellref == "BK_G1_A13" & bk$taxon == "new_alternate_anthospermum") , which(colnames(bk) == "cover2")] <- 0.5
bk[which(bk$cellref == "BK_G2_A5" & bk$taxon == "festuca_caprina") , which(colnames(bk) == "cover2")] <- 25
bk <- bk[-which(bk$cellref == "BK_G6_G12" & bk$taxon == "afroaster_erucifolius") , ] #sp not present on sheet
bk[which(bk$cellref == "BK_G7_B20" & bk$taxon == "gazania_krebsiana") , which(colnames(bk) == "cover2")] <- 0.5


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

bk <- bk[-which(bk$cover == 0), ] #sp not present on sheet, remove
bk[which(bk$cover == 0.2), which(colnames(bk) == "cover")] <- 0.5 #checked on paper sheet
bk[which(bk$cover == 1510.0), which(colnames(bk) == "cover")] <-15  #checked on paper sheet


##Cells with zero species richness
bk |> 
  group_by(cellref) |> 
  summarise(sum_cover = sum(cover)) |> 
  filter(sum_cover == 0) #none

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


####get names from trait data####
#for now, we'll just get the names from the trait data, and clean the rest of the FT data later

gg_trait_names <- read_excel("All_data/raw_data/raw_trait_data/GG_dataset_Functional_traits.xlsx") |> 
  rename(taxon = Species) |> 
  mutate(taxon = str_to_lower(taxon),
         taxon = str_squish(taxon),
         taxon = str_replace_all(taxon, " ", "_")) |> 
  select(taxon) |> 
  mutate(site = "GG", dataset = "FT") |> 
  distinct(taxon, site, dataset)

wh_trait_names <-read_excel("All_data/raw_data/raw_trait_data/WH_Functional_traits_dataset.xlsx", sheet = "Data_entry") |> 
  select(Date:Notes) |> 
  rename(taxon = Species) |> 
  mutate(taxon = str_to_lower(taxon),
         taxon = str_squish(taxon),
         taxon = str_replace_all(taxon, " ", "_")) |> 
  select(taxon) |> 
  mutate(site = "WH", dataset = "FT") |> 
  distinct(taxon, site, dataset)

bk_trait_names <- read_excel("All_data/raw_data/raw_trait_data/FT_Bokong_Edited.xlsx", sheet = "FT measurements") |> 
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



###Now clean and standardise names####
#Load standardise names function
#run Function_standardise_names.R


#import naming system 
name_trail <- read.xlsx("All_data/clean_data/Species names/micro_climb_ALL_names_editing.xlsx", sheet = "editing")

#clean names of veg surveys on elevation gradient
gg_summer_clean_names <- standardise_names(gg_summer2, "taxon", naming_system = name_trail, 
                                    "taxon", c("synonym1", "synonym2", "synonym3", "synonym4"))

wh_clean_names <- standardise_names(wh2, "taxon", naming_system = name_trail, 
                                    "taxon", c("synonym1", "synonym2", "synonym3", "synonym4"))

bk_clean_names <- standardise_names(bk, "taxon", naming_system = name_trail, 
                                    "taxon", c("synonym1", "synonym2", "synonym3", "synonym4"))


#Quality control the names changes
unique(gg_summer_clean_names$change_tracker)
unique(wh_clean_names$change_tracker)
unique(bk_clean_names$change_tracker)


#Bind all three sites together and export
micro_climb_veg_survey <- gg_summer_clean_names |> 
  bind_rows(wh_clean_names) |> 
  bind_rows(bk_clean_names) |> 
  select(!change_tracker) |> 
  rename(Cell_ID = cellref) |> 
  select(Cell_ID, site, grid, column, row, taxon, cover)

write.csv(micro_climb_veg_survey, "All_data/clean_data/micro_climb_occurrence.csv")
