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
