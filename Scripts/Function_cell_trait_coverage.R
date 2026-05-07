#=========================#
#====Cell_trait_coverage==#
#=========================#
####Function to identify cells that have trait measurements for less than 80% of the cover###
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
  select(-c(site, grid, column, row, ncolumn))

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1) |> 
  rownames_to_column(var = "taxon")

cell_trait_coverage <- function() {
  abundances <- abun_matrix |> 
    filter(cover>0) #table with covers
  mean_traits <- mean_traits #table with mean traits for species
  
  Cell_IDlist <- c(unique(abundances$Cell_ID))
  
  for(c in Cell_IDlist) {
    cell_abun <- abundances[abundances$Cell_ID ==c, ]
    
    merge <- cell_abun |> 
      left_join(mean_traits, by = "taxon")
    
    total_cov <- sum(merge$cover)
  }
  
  
  
  
  
  
}