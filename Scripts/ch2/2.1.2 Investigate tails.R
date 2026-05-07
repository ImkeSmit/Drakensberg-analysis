####Investigate long tails of SES###

cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA), 
         site = case_when(grepl("BK", cellref) == T ~ "BK", #add elevation variable
                          grepl("WH", cellref) == T ~ "WH",
                          grepl("GG", cellref) == T ~ "GG",.default = NA),
         grid = str_sub(cellref, 1, 3), 
         column = str_sub(cellref, 4,4), 
         row = as.numeric(str_sub(cellref, 5,6)), 
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(Cell_ID = paste0(site, "_G", str_sub(cellref, 3, 3), "_", column, row), 
         elevation = as.numeric(elevation)) |> 
  select(-c(cellref, site, grid, column, row)) 

ggplot(cell_ses, aes(x = SES)) +
  geom_histogram()+
  facet_wrap(~trait, scales = "free") +
  theme_classic()

#long tails height
height_tails <- cell_ses |> 
  filter(trait == "Height_cm", SES > 2.5) 

#lets look at their composition
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1) |> 
  rownames_to_column(var = "Cell_ID") |> 
  pivot_longer(cols = !Cell_ID, names_to = "taxon", values_to = "cover") |> 
  filter(cover>0) |> 
  mutate(
  site = case_when(grepl("BK", Cell_ID) == T ~ "BK", #add elevation variable
                   grepl("WH", Cell_ID) == T ~ "WH",
                   grepl("GG", Cell_ID) == T ~ "GG",.default = NA),
  grid = str_sub(Cell_ID, 1, 3), 
  column = str_sub(Cell_ID, 4,4), 
  row = as.numeric(str_sub(Cell_ID, 5,6)),
  ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(Cell_ID = paste0(site, "_G", str_sub(Cell_ID, 3, 3), "_", column, row)) |> 
  select(-c(site, grid, column, row)) 

composition_height_tails <- height_tails |> 
  inner_join(abun_matrix)

spnumber_height_tails <- composition_height_tails |> 
  group_by(Cell_ID) |> 
  summarise(sprichness = n())
hist(spnumber_height_tails$sprichness)


#long tails leaf area
la_tails <- cell_ses |> 
  filter(trait == "Leaf_area_mm2", SES > 3) 

composition_la_tails <- la_tails |> 
  inner_join(abun_matrix)

spnumber_la_tails <- composition_la_tails |> 
  group_by(Cell_ID) |> 
  summarise(sprichness = n())
hist(spnumber_la_tails$sprichness)


#let slook at the null model
nullcomm_cells <- readRDS("All_data/comm_assembly_results/nullmodel_C5_cells.rds")
nullcomm1 <- as.data.frame(nullcomm_cells[[1]]) |> 
  rownames_to_column(var = "Cell_ID") |> 
  pivot_longer(cols = !Cell_ID, names_to = "taxon", values_to = "cover") |> 
  filter(cover>0, Cell_ID == "WH5A8")

observed_comm <- abun_matrix |> 
  filter(Cell_ID == "WH_G5_A8")

