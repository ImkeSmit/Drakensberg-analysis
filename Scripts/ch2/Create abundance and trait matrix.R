###This is the script used to create the abundance and trait matrix used in further analyses
#we only do this once to avoid inconsistencies in downstream analyses
library(openxlsx)
library(tidyverse)
library(tidylog)

####Import community and trait data####
#occurrence data
drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))

#trait data
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) |> 
  pivot_longer(cols = c(Wet_mass_mg, Dry_mass_mg, Chlorophyll_mg_per_m2, Ft, Height_cm, 
                        Thickness_mm, Leaf_area_mm2, SLA, LDMC), names_to = "trait", values_to = "value")


#Get mean traits for species
mean_traits <- FT |> 
  filter(trait %in% c("Height_cm", "Leaf_area_mm2","Thickness_mm", "SLA", "LDMC")) |> 
  group_by(taxon, trait) |> 
  summarise(mean_trait = mean(value, na.rm = T)) |> 
  pivot_wider(names_from = trait, values_from = mean_trait) |> 
  ungroup() |> 
  arrange(taxon)

mean_traits <- as.data.frame(mean_traits)
row.names(mean_traits) <- mean_traits$taxon
mean_traits <- mean_traits[, -1]

#Do inner join between trait and cover data to get sp that match ebtween the two
FT_join <- drak |> 
  inner_join(mean_traits, by = "taxon") |> #inner join to only work with taxa that have trait data
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

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

#remove sites that have only one species 
no <- specnumber(abun_matrix)
abun_matrix <- abun_matrix[-which(no == 1), ]

#Save abundance matrix:
write.csv(abun_matrix, "All_data/comm_assembly_results/abun_matrix.csv")

#Save trait matrix: 
#first get row and columns names right
mean_traits <- as.data.frame(mean_traits)
row.names(mean_traits) <- mean_traits$taxon
mean_traits <- mean_traits[, -1]

write.csv(mean_traits, "All_data/comm_assembly_results/mean_traits.csv")


mean_traits2 <- FT |> 
  filter(trait %in% c("Height_cm", "Leaf_area_mm2","Thickness_mm", "SLA", "LDMC")) |> 
  group_by(taxon, trait) |> 
  summarise(mean_trait = mean(value, na.rm = T)) |> 
  pivot_wider(names_from = trait, values_from = mean_trait) |> 
  ungroup() |> 
  arrange(taxon)

FT_join2 <- drak |> 
  inner_join(mean_traits2, by = "taxon") |> #inner join to only work with taxa that have trait data
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

#create abundance matrix
abun_matrix2 <- FT_join2 |> 
  select(cellref, taxon, cover) |> 
  distinct(cellref, taxon, .keep_all = T) |> 
  ungroup() |> 
  arrange(taxon) |> 
  mutate(cover = ceiling(cover)) |> #change cover values to integer to use in null models
  pivot_wider(names_from = taxon, values_from = cover) 

abun_matrix2 <- as.data.frame(abun_matrix2)
row.names(abun_matrix2) <- abun_matrix2$cellref
abun_matrix2 <- abun_matrix2[, -1]

#replace NA values with 0
for(r in 1:nrow(abun_matrix2)) {
  for(c in 1:ncol(abun_matrix2)) {
    
    if(is.na(abun_matrix2[r,c])) {
      abun_matrix2[r,c] <- 0
    }
  }
}

#remove sites that have only one species 
no <- specnumber(abun_matrix2)
abun_matrix2 <- abun_matrix2[-which(no == 1), ]

