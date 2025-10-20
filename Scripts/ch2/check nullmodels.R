#test null model
#occurrence data
drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))

#trait data
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) |> 
  pivot_longer(cols = c(Wet_mass_mg, Dry_mass_mg, Chlorophyll_mg_per_m2, Ft, Height_cm, 
                        Thickness_mm, Leaf_area_mm2, SLA, LDMC), names_to = "trait", values_to = "value")

#combine occurrence and trait data
FT_join <- drak |> 
  inner_join(FT, by = c("site", "grid", "cell", "taxon")) |> #inner join to only work with taxa that have trait data
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

#Get mean traits for species
mean_traits <- FT_join |> 
  filter(trait %in% c("Height_cm", "Leaf_area_mm2", "SLA", "LDMC")) |> 
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


#create null models
abun_mat_small <- abun_matrix[c(1:10), c(1:20)]


set.seed(99)
nullcomm_cells <- generate_C5_null_imp(abun_mat_small, 1, pool = "entire")

#let's check colsums and rowsums
colSums(abun_mat_small)
colSums(nullcomm_cells[[1]]) #colsums (sum of sp abundances) stay the same

rowSums(abun_mat_small)
rowSums(nullcomm_cells[[1]]) #rowsums (sum of abundances in a site) not staying the same
#that is ok

sum(colSums(abun_matrix))
sum(colSums(nullcomm_cells[[1]])) #total sum of matrix also not staying the same

sum(rowSums(abun_mat_small))
sum(rowSums(nullcomm_cells[[1]]))
#when the matrix is big, these are not equel. Problem???


test <- matrix(c(1,1,2,2,2,1),nrow = 3, ncol = 2)
swapped <- swap_matrix(test, n_swaps = 1e5)
