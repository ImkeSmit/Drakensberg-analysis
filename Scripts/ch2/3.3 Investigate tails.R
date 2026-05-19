####Investigate long tails of SES###

cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  rename(Cell_ID = cellref)

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
  filter(cover>0) 

composition_height_tails <- height_tails |> 
  inner_join(abun_matrix, by = "Cell_ID")

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

