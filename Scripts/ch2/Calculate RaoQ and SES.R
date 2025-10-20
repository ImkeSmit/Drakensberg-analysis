#Community assembly analysis 2.0###
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
#library(vegan)
library(traitstrap)
library(FD)
library(future.apply)

####Create required functions####

#~~~~~~~~~~~~~~~~~~
#calculate RaoQ, one trait at a time
#This is unweighted RaoQ based on euclidean distance. 
#It uses the mean traits of species
#Figure out: does trait filling make sense?
#Should we use raw traits when we have them in cells and only mean traits to fill gaps?

calc_RaoQ <- function(mean_traits, abun_matrix) {
  
  traitlist <- c(colnames(mean_traits))
  abun_matrix2 <- as.data.frame(abun_matrix)
  
  abun_long <- abun_matrix2 |> 
    rownames_to_column(var = "cellref") |> 
    pivot_longer(cols = !cellref, names_to = "taxon", values_to = "cover")
  
  traits_long <- mean_traits |> 
    rownames_to_column(var = "taxon") |> 
    pivot_longer(!taxon, names_to = "trait", values_to = "value")
  
  trait_sum <- abun_long |> 
    full_join(traits_long, by = "taxon") |> 
    filter(cover > 0) |>
    mutate(value = if_else(is.na(value), 0, value))  |>
    group_by(cellref, trait) |> 
    mutate(trait_sum = sum(value)) |> 
    ungroup()

  
  for(t in 1:length(traitlist)) {
    
    chosen_trait <- mean_traits[, t] # preserve matrix format
    names(chosen_trait)<- row.names(mean_traits)
    
    #Check for communities with zero-sum abundances
    #Get cells that do not have any measurements of chosen_trait
    
    abun_temp <- abun_matrix2
    
    problems <- trait_sum |> 
      filter(trait == traitlist[t]) |> 
      filter(trait_sum == 0)

    
    if(nrow(problems) > 0) {
      #remove these cells from the abundance matrix
      abun_temp <- abun_temp[-which(rownames(abun_temp) %in% c(problems$cellref)) , , drop = F]
    }
    
    #identify species that do not occur in any of the remaining cells, and remove them
    abundance_sums <- colSums(abun_temp)
    empty_names <- names(abundance_sums[which(abundance_sums == 0)])
    
    if(length(empty_names) > 0) {
      abun_temp <- abun_temp[, - which(colnames(abun_temp) %in% c(empty_names)), drop = F]
      #also remove these sp from the trait matrix
      chosen_trait <- chosen_trait[-which(names(chosen_trait) %in% c(empty_names)), drop = F]
    } 
    
    #calculate RaoQ 
    FD_cells <- dbFD(chosen_trait, abun_temp,
                     w.abun = F, #do not weight RaoQ by abundances
                     corr = "cailliez", 
                     calc.FRic = F, 
                     scale.RaoQ = F, 
                     calc.FGR = F, 
                     calc.FDiv = F, 
                     calc.CWM = F, 
                     messages = F)
    #RAoQ = 0 if there is only one distinct trait value in a cell
    
    if(t==1) {
      
      RaoQ_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                                 RaoQ = FD_cells$RaoQ, 
                                 trait = traitlist[t])
      rownames(RaoQ_results) <- NULL
      
    }else {
      more_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                                 RaoQ = FD_cells$RaoQ,
                                 trait = traitlist[t])
      rownames(more_results) <- NULL
      
      RaoQ_results <- rbind(RaoQ_results, more_results)
    }
    
  } 
  return(RaoQ_results) } 



#~~~~~~~~~~~~~~~~~~~~~~~
###Nullmodels###

###Null model C5###
#good at detecting EF
#first randomise presences accross sites. 
#then assign abundances to those sites, pulling from the abundances where the sp occurs in the observed comm

generate_C5_null <- function(comm, iterations, pool) {
  # Ensure comm is a matrix for speed
  comm <- as.matrix(comm)
  rownames_comm <- rownames(comm)
  colnames_comm <- colnames(comm)
  
  # Precompute site groups
  three_sites <- c("GG", "WH", "BK")
  grids_22 <- c("BK1","BK2","BK3","BK4","BK5","BK6","BK7",
                "WH1","WH2","WH3","WH4","WH5","WH6","WH7",
                "GG1","GG2","GG3","GG4","GG5","GG6","GG7","GG8")
  
  # Precompute per-species nonzero abundances for faster sampling
  abundance_lists <- lapply(1:ncol(comm), function(j) comm[comm[, j] > 0, j])
  names(abundance_lists) <- colnames_comm
  
  # Function to run a single iteration
  single_iter <- function(iter) {
    null_comm <- matrix(0, nrow = nrow(comm), ncol = ncol(comm),
                        dimnames = list(rownames_comm, colnames_comm))
    
    for (i in seq_len(nrow(comm))) {
      site <- comm[i, ]
      richness <- sum(site)
      if (richness == 0) next
      
      # Choose species pool
      if (pool == "entire") {
        chosen_species <- sample(colnames_comm, richness, replace = FALSE)
      } else if (pool == "site") {
        pool_site <- three_sites[grep(paste(three_sites, collapse = "|"), rownames_comm[i])]
        sp_pool <- comm[grep(pool_site, rownames_comm), , drop = FALSE]
        sp_pool <- sp_pool[, colSums(sp_pool) > 0, drop = FALSE]
        chosen_species <- sample(colnames(sp_pool), richness, replace = FALSE)
      } else if (pool == "grid") {
        pool_site <- grids_22[grep(paste(grids_22, collapse = "|"), rownames_comm[i])]
        sp_pool <- comm[grep(pool_site, rownames_comm), , drop = FALSE]
        sp_pool <- sp_pool[, colSums(sp_pool) > 0, drop = FALSE]
        chosen_species <- sample(colnames(sp_pool), richness, replace = FALSE)
      }
      
      # Assign abundances
      for (sp in chosen_species) {
        abunds <- abundance_lists[[sp]]
        if (length(abunds) > 0) {
          null_comm[i, sp] <- sample(abunds, 1)
        }
      }
    }
    null_comm
  }
  

  null_list <- lapply(seq_len(iterations), single_iter)
  
  message(sprintf(
    "C5 null model completed (%d iterations, pool = %s).",
    iterations, pool
  ))
  
  return(null_list)
}


###improved function
generate_C5_null_imp <- function(comm, iterations = 10, pool = "entire") {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Please install the 'vegan' package: install.packages('vegan')")
  }
  
  comm <- as.matrix(comm)
  rownames_comm <- rownames(comm)
  colnames_comm <- colnames(comm)
  
  # Presence–absence and abundance data
  pa <- (comm > 0) * 1
  species_freq <- colSums(pa)
  site_richness <- rowSums(pa)
  
  # Per-species abundance distributions
  abundance_lists <- lapply(colnames_comm, function(sp) comm[comm[, sp] > 0, sp])
  names(abundance_lists) <- colnames_comm
  
  # Define site groups
  three_sites <- c("GG", "WH", "BK")
  grids_22 <- c("BK1","BK2","BK3","BK4","BK5","BK6","BK7",
                "WH1","WH2","WH3","WH4","WH5","WH6","WH7",
                "GG1","GG2","GG3","GG4","GG5","GG6","GG7","GG8")
  
  get_sites_in_pool <- function(site_name) {
    if (pool == "entire") return(rownames_comm)
    if (pool == "site") {
      site_prefix <- three_sites[grepl(paste(three_sites, collapse="|"), site_name)]
      return(rownames_comm[grepl(site_prefix, rownames_comm)])
    }
    if (pool == "grid") {
      grid_prefix <- grids_22[grepl(paste(grids_22, collapse="|"), site_name)]
      return(rownames_comm[grepl(grid_prefix, rownames_comm)])
    }
    stop("Invalid pool type")
  }
  
  # ---- Function for one iteration ----
  single_iter <- function(iter) {
    pres_matrix <- pa
    
    if (pool == "entire") {
      perm <- vegan::permatswap(pres_matrix, fixedmar = "both", mtype = "prab", times = 1)
      pres_matrix <- perm$perm[[1]]
    } else {
      pres_matrix[,] <- 0
      site_groups <- unique(sapply(rownames_comm, function(x)
        get_sites_in_pool(x)[1]))
      for (sg in site_groups) {
        group_sites <- get_sites_in_pool(sg)
        group_rows <- which(rownames_comm %in% group_sites)
        sp_in_group <- colnames_comm[colSums(pa[group_rows, , drop = FALSE]) > 0]
        submat <- pa[group_rows, sp_in_group, drop = FALSE]
        
        if (sum(submat) > 0) {
          perm_sub <- vegan::permatswap(submat, fixedmar = "both", mtype = "prab", times = 1)
          pres_matrix[group_rows, sp_in_group] <- perm_sub$perm[[1]]
        }
      }
    }
    
    # Optional: check that species frequencies were preserved
    freq_after <- colSums(pres_matrix)
    if (!all(freq_after == species_freq)) {
      warning("Species frequencies changed slightly after swap; abundance assignment uses replacement.")
    }
    
    # ---- Assign abundances ----
    null_comm <- matrix(0, nrow = nrow(comm), ncol = ncol(comm),
                        dimnames = list(rownames_comm, colnames_comm))
    
    for (sp in colnames_comm) {
      sp_sites <- which(pres_matrix[, sp] == 1)
      abunds <- abundance_lists[[sp]]
      
      n_sites <- length(sp_sites)
      n_abunds <- length(abunds)
      
      if (n_sites > 0 && n_abunds > 0) {
        # If mismatch between number of sites and abundances, use replace = TRUE
        replace_flag <- n_sites > n_abunds
        null_comm[sp_sites, sp] <- sample(abunds, n_sites, replace = replace_flag)
      }
    }
    
    null_comm
  }
  
  # ---- Run iterations ----
  null_list <- lapply(seq_len(iterations), single_iter)
  
  message(sprintf(
    "C5 null model completed (%d iterations, pool = %s).",
    iterations, pool
  ))
  
  return(null_list)
}


    




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


####SES at cell scale, pool = entire####
#observed RaoQ
RQ_obs_cells <- calc_RaoQ(mean_traits, abun_matrix)

#Create null models
set.seed(123)
nullcomm_cells <- generate_C5_null(abun_matrix, 999, pool = "entire")

#Calculate SES#
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

#some graphs
RQ_ele <- drak |> 
  select(cellref, site) |> 
  distinct() |> 
  right_join(RQ_cells_summary, by = "cellref")
RQ_ele$site <- factor(RQ_ele$site, levels = c("GG", "WH", "BK"))

RQ_ele |> 
  dplyr::filter(SES < 5) |> #lots of outlier SES values, logical??
  ggplot(aes(x = site, y = SES)) +
  geom_boxplot() +
  facet_wrap(~trait) 


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







