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
#select the trait we want to work with 
height <- FT_join |> 
  filter(trait == "Height_cm") |> 
  filter(!is.na(value)) 

chosen_one <- height$cellref[1]

FT_subset <- height |> 
  filter(cellref == chosen_one) |>  #don't know why we have 2 trait measurements for erica alopecurus
  select(taxon, trait, value) |> 
  group_by(taxon) |> 
  mutate(value = mean(value)) |> 
  ungroup() |> 
  distinct(taxon, .keep_all = T) |> 
  pivot_wider(names_from = "trait", values_from = "value")
trait <- c(FT_subset$Height_cm)
names(trait) <- c(FT_subset$taxon)

abun_subset <- FT_join |> 
  filter(trait == "Height_cm", 
         cellref == chosen_one) |> 
  filter(!is.na(value)) |> 
  select(taxon, cover) |> 
  group_by(taxon) |> 
  mutate(cover = mean(cover)) |> 
  ungroup() |> 
  distinct(taxon, .keep_all = T) |> 
  pivot_wider(names_from = "taxon", values_from = "cover")
abun <- c(abun_subset[1,])
  
RaoQ <- dbFD(x = trait, 
             a = abun)
#something wrong with species labels here... how to transform to matrix with rownames??




