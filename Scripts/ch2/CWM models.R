####How does CWM traits vary with elevation?####
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)
library(traitstrap)
library(FD)
library(ggridges)

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

#combine occurrence and trait data
FT_join <- drak |> 
  inner_join(FT, by = c("site", "grid", "cell", "taxon")) |> #inner join to only work with taxa that have trait data
  mutate(cellref = paste0(site, grid, cell)) |> 
  dplyr::select(!c(column, row))

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
  dplyr::select(cellref, taxon, cover) |> 
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

no <- specnumber(abun_matrix)
abun_matrix <- abun_matrix[which(no == 1), ]


#compute CWM of each trait for each cell
cwm <- functcomp(x = mean_traits, a = as.matrix(abun_matrix))
cwm <- cwm |>
  rownames_to_column(var = "cellref") |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  pivot_longer(cols = c("Height_cm", "LDMC", "Leaf_area_mm2", "SLA"), names_to = "trait", values_to = "cwm_value")
cwm$elevation <- as.factor(cwm$elevation)  

cwm_ridges <- cwm |> 
  ggplot(aes(x = cwm_value, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~trait, scales = "free_x") +
  theme_classic() 

ggsave(filename = "cwm_elevation.png", plot = cwm_ridges, path = "Figures")


###Which species lie where on the cwm trait spectrum? 
cwm_xt <- cwm %>%
  group_by(elevation, trait) %>%
  summarise(
    highest_cell = cellref[which.max(cwm_value)],
    highest_val  = max(cwm_value),
    lowest_cell  = cellref[which.min(cwm_value)],
    lowest_val   = min(cwm_value),
    median_cell  = cellref[which.min(abs(cwm_value - median(cwm_value)))],
    median_val   = median(cwm_value)
  )

abun_matrix[which(row.names(abun_matrix) == "GG2H15"), ]
abun_matrix["GG2H15", abun_matrix[, "GG2H15"] > 0, drop = FALSE]

site <- "WH7B5"

present_df <- data.frame(
  species   = colnames(abun_matrix)[abun_matrix[site, ] > 0],
  abundance = abun_matrix[site, abun_matrix[site, ] > 0]
) #something is wrong.. a lot of these extreme values come from cells with one sp, I thought they were discarded...

