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
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1)


#compute CWM of each trait for each cell
cwm <- functcomp(x = mean_traits, a = as.matrix(abun_matrix))
cwm <- cwm |>
  rownames_to_column(var = "cellref") |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  pivot_longer(cols = c("Height_cm", "LDMC", "Leaf_area_mm2", "SLA", "Thickness_mm"), names_to = "trait", values_to = "cwm_value")
cwm$elevation <- as.factor(cwm$elevation)  

cwm_ridges <- cwm |> 
  ggplot(aes(x = cwm_value, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~trait, scales = "free_x", ncol = 2, nrow = 3) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(filename = "cwm_elevation.png", plot = cwm_ridges, path = "Figures", 
       width = 1200, height = 1500, units = "px")


###Which species lie where on the cwm trait spectrum? 
#get cells which have extremem or median cwm values
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

#get the composition of each of these cells
for(i in c(cwm_xt$highest_cell, cwm_xt$lowest_cell, cwm_xt$median_cell)) {
  site <- i

if (i == "GG4H14")  {
present_df <- data.frame(
  site = site,
  species   = colnames(abun_matrix)[abun_matrix[site, ] > 0],
  abundance = unlist(abun_matrix[site, abun_matrix[site, ] > 0]))
rownames(present_df) <- NULL
} else {
  temp_present_df <- data.frame(
    site = site,
    species   = colnames(abun_matrix)[abun_matrix[site, ] > 0],
    abundance = unlist(abun_matrix[site, abun_matrix[site, ] > 0]))
  rownames(temp_present_df) <- NULL
  
  present_df <- rbind(present_df, temp_present_df)
}}

#join the composition to the cwm_xt table
cell_composition <- present_df |> 
  rename(cell = site) |> 
  group_by(cell) |> 
  summarise(comp_list = list(species))
  
cwm_xt_join <- cwm_xt |> 
  left_join(cell_composition, by = join_by(highest_cell == cell)) |> 
  rename(higest_comp = comp_list) |> 
  left_join(cell_composition, by = join_by(lowest_cell == cell)) |> 
  rename(lowest_comp = comp_list) |> 
  left_join(cell_composition, by = join_by(median_cell == cell)) |> 
  rename(median_comp = comp_list)

#export
write.xlsx(cwm_xt_join, "All_data/comm_assembly_results/composition_extreme_cwm.xlsx")
