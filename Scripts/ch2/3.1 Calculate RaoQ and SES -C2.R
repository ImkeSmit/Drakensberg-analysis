#Calculate RaoQ and SES###
#at different scales and with different pools
#Using the C2 null model to test for limiting similarity
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
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1)


####All species - SES at cell scale, C2, pool = site, weighted RaoQ####
#observed RaoQ, weighted by abundance
RQ_obs_cells <- calc_RaoQ_weighted(mean_traits, abun_matrix)

#Create null models
set.seed(150)
nullcomm_cells <- generate_C2_null(abun_matrix, 999, pool = "site")

saveRDS(nullcomm_cells, file = "All_data/comm_assembly_results/nullmodel_C2_cells.rds")
nullcomm_cells <- readRDS("All_data/comm_assembly_results/nullmodel_C2_cells.rds")
#Calculate SES, with unscaled RaoQ#
#we need to calculate RaoQ for each of the observed null communities

# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)  

chunks <- split(nullcomm_cells, ceiling(seq_along(nullcomm_cells) / 50))
#each chunk has 50 null matrices
RaoQ_results <- list()

for (i in seq_along(chunks)) { #run chunks sequentially
  message("Processing chunk ", i, " of ", length(chunks))
  
  chunk_list <- chunks[[i]]  # smaller subset only
  
  #each core receives one null matrix to run calc_RaoQ on. this happens until all matrices in the chunk are finished
  sub_RQ <- future_lapply(seq_along(chunk_list), function(idx) {
    chosen_null <- chunk_list[[idx]]
    RQ_result <- calc_RaoQ_weighted(mean_traits, chosen_null)
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

write.csv(RQ_cells_summary, "All_data/comm_assembly_results/RQ_weighted_cells_C2_entire.csv")

##NB! Using this null model, RQnull may have been calculated with less sites and less species than RQobs
#This is because the C2 nullmodel may create empty species columns.
#these species and any subsequent empty sites are removed in the RQ calculation


####Forbs only - SES at cell scale, C5, pool = entire, weighted RaoQ####

#import matrix with only forbs occurrences
forb_abun_matrix <- read.csv("All_data/comm_assembly_results/abun_forbs.csv", row.names = 1)

forb_mean_traits <- read.csv("All_data/comm_assembly_results/traits_forbs.csv", row.names = 1)

#observed RaoQ, weighted by abundance
RQ_obs_cells <- calc_RaoQ_weighted(forb_mean_traits, forb_abun_matrix)

#Create null models
set.seed(150)
nullcomm_cells <- generate_C5_null_imp(forb_abun_matrix, 999, pool = "entire")

saveRDS(nullcomm_cells, file = "All_data/comm_assembly_results/forbs_only_nullmodel_C5_cells.rds")
nullcomm_cells <- readRDS("All_data/comm_assembly_results/forbs_only_nullmodel_C5_cells.rds")
#Calculate SES, with unscaled RaoQ#
#we need to calculate RaoQ for each of the observed null communities

# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)  

chunks <- split(nullcomm_cells, ceiling(seq_along(nullcomm_cells) / 50))
#each chunk has 50 null matrices
RaoQ_results <- list()

for (i in seq_along(chunks)) { #run chunks sequentially
  message("Processing chunk ", i, " of ", length(chunks))
  
  chunk_list <- chunks[[i]]  # smaller subset only
  
  #each core receives one null matrix to run calc_RaoQ on. this happens until all matrices in the chunk are finished
  sub_RQ <- future_lapply(seq_along(chunk_list), function(idx) {
    chosen_null <- chunk_list[[idx]]
    RQ_result <- calc_RaoQ_weighted(forb_mean_traits, chosen_null)
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

write.csv(RQ_cells_summary, "All_data/comm_assembly_results/forbs_only_RQ_weighted_cells_C5_entire.csv")


####Graminoids only - SES at cell scale, C5, pool = entire, weighted RaoQ####

#import matrix with only forbs occurrences
gram_abun_matrix <- read.csv("All_data/comm_assembly_results/abun_graminoids.csv", row.names = 1)

gram_mean_traits <- read.csv("All_data/comm_assembly_results/traits_graminoids.csv", row.names = 1)

#observed RaoQ, weighted by abundance
RQ_obs_cells <- calc_RaoQ_weighted(gram_mean_traits, gram_abun_matrix)

#Create null models
set.seed(150)
nullcomm_cells <- generate_C5_null_imp(gram_abun_matrix, 999, pool = "entire")

saveRDS(nullcomm_cells, file = "All_data/comm_assembly_results/graminoids_only_nullmodel_C5_cells.rds")
nullcomm_cells <- readRDS("All_data/comm_assembly_results/graminoids_only_nullmodel_C5_cells.rds")
#Calculate SES, with unscaled RaoQ#
#we need to calculate RaoQ for each of the observed null communities

# Set up parallel backend
plan(multisession, workers = parallel::detectCores() - 2)  

chunks <- split(nullcomm_cells, ceiling(seq_along(nullcomm_cells) / 50))
#each chunk has 50 null matrices
RaoQ_results <- list()

for (i in seq_along(chunks)) { #run chunks sequentially
  message("Processing chunk ", i, " of ", length(chunks))
  
  chunk_list <- chunks[[i]]  # smaller subset only
  
  #each core receives one null matrix to run calc_RaoQ on. this happens until all matrices in the chunk are finished
  sub_RQ <- future_lapply(seq_along(chunk_list), function(idx) {
    chosen_null <- chunk_list[[idx]]
    RQ_result <- calc_RaoQ_weighted(gram_mean_traits, chosen_null)
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

write.csv(RQ_cells_summary, "All_data/comm_assembly_results/graminoids_only_RQ_weighted_cells_C5_entire.csv")



