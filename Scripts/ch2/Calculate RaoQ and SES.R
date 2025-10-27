#Calculate RaoQ and SES###
#at different scales and with different pools
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)
library(traitstrap)
library(FD)
library(future.apply)

#load required functions first

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


####SES at cell scale, C5, pool = entire####
#observed RaoQ, not scaled
RQ_obs_cells <- calc_RaoQ(mean_traits, abun_matrix)

#Create null models
set.seed(123)
nullcomm_cells <- generate_C5_null_imp(abun_matrix, 999, pool = "entire")
saveRDS(nullcomm_cells, file = "All_data/comm_assembly_results/nullmodel_C5_cells.rds")

#Calculate SES, with unscaled RaoQ#
#we need to calculate RaoQ for each of the observed null communities

# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)  # or plan(multicore) on Linux/mac

chunks <- split(nullcomm_cells, ceiling(seq_along(nullcomm_cells) / 50))
#each chunk has 50 null matrices
RaoQ_results <- list()

for (i in seq_along(chunks)) { #run chunks sequentially
  message("Processing chunk ", i, " of ", length(chunks))
  
  chunk_list <- chunks[[i]]  # smaller subset only
  
  #each core receives one null matrix to run calc_RaoQ on. this happens until all matrices in the chunk are finished
  sub_RQ <- future_lapply(seq_along(chunk_list), function(idx) {
    chosen_null <- chunk_list[[idx]]
    RQ_result <- calc_RaoQ(mean_traits, chosen_null)
    RQ_result$counter <- paste0("null matrix ", (i - 1) * 100 + idx)
    RQ_result
  }, future.seed = TRUE) #generates a unique reproducible sub seed for each worker
  #ensures that results are reproducible, and that there is no overlap in random processes for each core
  
  RaoQ_results[[i]] <- bind_rows(sub_RQ) #results from the chunk are merged
  rm(sub_RQ, chunk_list); gc()
}

null_RQ <- bind_rows(RaoQ_results)

plan(sequential)


#SES of each cell
RQ_cells_summary <- null_RQ |> 
  group_by(trait, cellref) |> 
  summarise(sd_null = sd(RaoQ), 
            mean_null = mean(RaoQ)) |> 
  filter(sd_null > 0) |> #cannot divide by zero in SES calculation
  inner_join(RQ_obs_cells, by = c("trait", "cellref")) |> 
  mutate(SES = (RaoQ - mean_null)/sd_null)

write.csv(RQ_cells_summary, "All_data/comm_assembly_results/RQ_cells_C5_entire.csv")



#####
####SES at cell scale, C2, pool = site####
#observed RaoQ, not scaled
RQ_obs_cells <- calc_RaoQ(mean_traits, abun_matrix)

#Create null models
set.seed(123)
nullcomm_cells <- generate_C2_null(abun_matrix, 999, pool = "site")
saveRDS(nullcomm_cells, file = "All_data/comm_assembly_results/nullmodel_C2_cells.rds")

#Calculate SES, with unscaled RaoQ#
#we need to calculate RaoQ for each of the observed null communities

# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)  # or plan(multicore) on Linux/mac

chunks <- split(nullcomm_cells, ceiling(seq_along(nullcomm_cells) / 50))
#each chunk has 50 null matrices
RaoQ_results <- list()

for (i in seq_along(chunks)) { #run chunks sequentially
  message("Processing chunk ", i, " of ", length(chunks))
  
  chunk_list <- chunks[[i]]  # smaller subset only
  
  #each core receives one null matrix to run calc_RaoQ on. this happens until all matrices in the chunk are finished
  sub_RQ <- future_lapply(seq_along(chunk_list), function(idx) {
    chosen_null <- chunk_list[[idx]]
    RQ_result <- calc_RaoQ(mean_traits, chosen_null)
    RQ_result$counter <- paste0("null matrix ", (i - 1) * 100 + idx)
    RQ_result
  }, future.seed = TRUE) #generates a unique reproducible sub seed for each worker
  #ensures that results are reproducible, and that there is no overlap in random processes for each core
  
  RaoQ_results[[i]] <- bind_rows(sub_RQ) #results from the chunk are merged
  rm(sub_RQ, chunk_list); gc()
}

null_RQ <- bind_rows(RaoQ_results)

plan(sequential)


#SES of each cell
RQ_cells_summary <- null_RQ |> 
  group_by(trait, cellref) |> 
  summarise(sd_null = sd(RaoQ), 
            mean_null = mean(RaoQ)) |> 
  filter(sd_null > 0) |> #cannot divide by zero in SES calculation
  inner_join(RQ_obs_cells, by = c("trait", "cellref")) |> 
  mutate(SES = (RaoQ - mean_null)/sd_null)

write.csv(RQ_cells_summary, "All_data/comm_assembly_results/RQ_cells_C5_entire.csv")


####SES at the grid scale, pool = entire####
##create abundance matrix at the grid level
abun_matrix_grid <- FT_join |> 
  mutate(gridref = paste0(site, grid)) |> 
  select(gridref, taxon, cover) |> 
  distinct(gridref, taxon, .keep_all = T) |>
  mutate(cover = ceiling(cover)) |> #change cover values to integer to use in null models
  group_by(gridref, taxon) |> 
  summarise(gridlvl_cover = sum(cover)) |> #gridlevel cover is just the sum of cell level covers, ok?
  ungroup() |> 
  arrange(taxon) |> 
  pivot_wider(names_from = taxon, values_from = gridlvl_cover) 

abun_matrix_grid <- as.data.frame(abun_matrix_grid)
row.names(abun_matrix_grid) <- abun_matrix_grid$gridref
abun_matrix_grid <- abun_matrix_grid[, -1]

#replace NA values with 0
for(r in 1:nrow(abun_matrix_grid)) {
  for(c in 1:ncol(abun_matrix_grid)) {
    
    if(is.na(abun_matrix_grid[r,c])) {
      abun_matrix_grid[r,c] <- 0
    }
  }
}

#observed RaoQ
RQ_obs_grid <- calc_RaoQ(mean_traits, abun_matrix_grid)

#Create null models
set.seed(123) #set the same seed here, I guess so?
nullcomm_grids <- generate_C5_null(abun_matrix_grid, 999, pool = "entire")

#Calculate RaoQ for each of the null grids
# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)  # or plan(multicore) on Linux/mac

chunks <- split(nullcomm_grids, ceiling(seq_along(nullcomm_grids) / 50))
#each chunk has 50 null matrices
RaoQ_grids_results <- list()
error_log <- list()  # store any errors

for (i in seq_along(chunks)) { #run chunks sequentially
  message("Processing chunk ", i, " of ", length(chunks))
  
  chunk_list <- chunks[[i]]  # smaller subset only
  
  #each core receives one null matrix to run calc_RaoQ on. this happens until all matrices in the chunk are finished
  sub_RQ <- future_lapply(seq_along(chunk_list), function(idx) {
    chosen_null <- chunk_list[[idx]]
  
    # Safely run calc_RaoQ without errors stopping it
  result <- tryCatch({
    RQ_result <- calc_RaoQ(mean_traits, chosen_null)
    RQ_result$counter <- paste0("null matrix ", (i - 1) * 50 + idx)
    list(success = TRUE, result = RQ_result)
  }, error = function(e) {
    list(success = FALSE, 
         error = e$message, 
         chunk = i, 
         index_in_chunk = idx,
         null_id = (i - 1) * 50 + idx)
  })
  
  return(result)
}, future.seed = TRUE)

# Separate successes and errors
successes <- lapply(sub_RQ, function(x) if (x$success) x$result else NULL)
errors <- lapply(sub_RQ, function(x) if (!x$success) x else NULL)

# Store results
RaoQ_grids_results[[i]] <- bind_rows(Filter(Negate(is.null), successes))
error_log[[i]] <- Filter(Negate(is.null), errors)

rm(sub_RQ, chunk_list); gc() #garbage collection, returns memory to OS
}

null_RQ_grids <- bind_rows(RaoQ_grids_results)

plan(sequential)


###Look at errors
failed <- unlist(error_log, recursive = FALSE)
length(failed)  # how many failed
failed[[1]]  # details for the first failure


#SES of each grid
RQ_grids_summary <- null_RQ_grids |> 
  group_by(trait, cellref) |> 
  summarise(sd_null = sd(RaoQ), 
            mean_null = mean(RaoQ)) |> 
  filter(sd_null > 0) |> #cannot divide by zero in SES calculation
  inner_join(RQ_obs_grid, by = c("trait", "cellref")) |> 
  mutate(SES = (RaoQ - mean_null)/sd_null)

write.csv(RQ_grids_summary, "All_data/comm_assembly_results/RQ_grids_C5_entire.csv")







