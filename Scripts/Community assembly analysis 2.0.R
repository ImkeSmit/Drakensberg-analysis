#Community assembly analysis 2.0###
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
#library(vegan)
library(traitstrap)
library(FD)

####Create required functions####

#~~~~~~~~~~~~~~~~~~
#calculate RaoQ, one trait at a time
#This is unweighted RaoQ based on euclidean distance. 
#It uses the mean traits of species
#Figure out: does trait filling make sense?
#Should we use raw traits when we have them in cells and only mean traits to fill gaps?

calc_RaoQ <- function(mean_traits, abun_matrix) {
  
  traitlist <- c(colnames(mean_traits))
  
  for(t in 1:length(traitlist)) {
    
    chosen_trait <- mean_traits[, t]
    names(chosen_trait)<- row.names(mean_traits)
    
    abun_matrix2 <- abun_matrix
    
    #Check for communities with zero-sum abundances
    #Get cells from that do not have any measurements of chosen_trait
    abun_long <- abun_matrix2 |> 
      rownames_to_column( var = "cellref") |> 
      pivot_longer(cols = !cellref, names_to = "taxon", values_to = "cover")
    
    traits_long <- data.frame(value = mean_traits[, t], taxon = row.names(mean_traits))
    
    problems <- abun_long |> 
      full_join(traits_long, by = "taxon") |> 
      filter(cover > 0) |> 
      mutate(value = if_else(is.na(value), 0, value))  |> 
      group_by(cellref) |> 
      mutate(sum_chlor = sum(value)) |> 
      ungroup() |> 
      filter(sum_chlor == 0) |> 
      distinct(cellref)
    
    
    if(nrow(problems) > 0) {
      #remove these cells from the abundance matrix
      abun_matrix2 <- abun_matrix2[-which(rownames(abun_matrix2) %in% c(problems$cellref)) , ]
    }
    
    #identify species that do not occur in any of the remaining cells, and remove them
    abundance_sums <- colSums(abun_matrix2)
    empty_names <- names(abundance_sums[which(abundance_sums == 0)])
    
    if(length(empty_names) > 0) {
      abun_matrix2 <- abun_matrix2[, - which(colnames(abun_matrix2) %in% c(empty_names))]
      #also remove these sp from the trait matrix
      chosen_trait <- chosen_trait[-which(names(chosen_trait) %in% c(empty_names))]
    } 
    
    #calculate RaoQ 
    FD_cells <- dbFD(chosen_trait, abun_matrix2,
                     w.abun = F, #do not weight RaoQ by abundances
                     corr = "cailliez", 
                     calc.FRic = F, 
                     scale.RaoQ = F, 
                     calc.FGR = F, 
                     calc.FDiv = F, 
                     calc.CWM = F)
    #RAoQ = 0 if there is only one distinct trait value in a cell
    
    if(t==1) {
      
      RaoQ_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                                 RaoQ = FD_cells$RaoQ, 
                                 trait = traitlist[t])
      
    }else {
      more_results <- data.frame(cellref = names(FD_cells$RaoQ), 
                                 RaoQ = FD_cells$RaoQ,
                                 trait = traitlist[t])
      
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
generate_C5_null <- function(comm, iterations, pool) { #either entire, site, grid
  
  null_list <- vector(mode = "list", length = iterations)
  
  for(n in 1:iterations) {
    null_comm <- comm * 0  # initialize matrix
    
    for (i in 1:nrow(comm)) {
      site <- comm[i, ]
      
      # Species present at site (richness)
      richness <- sum(site > 0)
      
      if(pool == "entire") {
      # Randomly choose species (without replacement) to occupy this site
      #Choose from all species in the matrix
      chosen_species <- sample(colnames(comm), richness, replace = FALSE)
      
      } else { if(pool == "site"){
        
        #which site to draw sp from
        three_sites <- c("GG", "WH", "BK")
        pool_site <- three_sites[str_detect(rownames(site), three_sites)]
        
        #isolate the species pool
        sp_pool <- comm[grep(pool_site, rownames(comm)) , ]
        sp_pool <- sp_pool[, which(colSums(sp_pool) > 0)]
        
        chosen_species <- sample(colnames(sp_pool), richness, replace = FALSE)
        
      } else { #pool == grid
        
        #which grid to draw sp from
        grids_22 <- c("BK1","BK2","BK3", "BK4", "BK5","BK6", "BK7", 
                      "WH1","WH2", "WH3","WH4","WH5", "WH6","WH7", 
                       "GG1", "GG2", "GG3","GG4","GG5", "GG6","GG7", "GG8")
        pool_site <- grids_22[str_detect(rownames(site), grids_22)]
        
        #isolate the species pool
        sp_pool <- comm[grep(pool_site, rownames(comm)) , ]
        sp_pool <- sp_pool[, which(colSums(sp_pool) > 0)]
        
        chosen_species <- sample(colnames(sp_pool), richness, replace = FALSE)
      
        }}
      
      #For C3, we can pick abundances from sites in the matrix where the sp occurs
      for (s in 1:length(chosen_species)) {
        possible_abundances <-  comm |> 
          select(all_of(chosen_species)) |> 
          pull(chosen_species[s]) |> 
          discard(~ .x <= 0)
        
        #abundances that are more frequent are more likely to be sampled
        chosen_abundance <- sample(possible_abundances, 1)
        
        #assign abundance to species and site
        null_comm[i, chosen_species[s]] <- chosen_abundance
      } #end loop through species
    }#end loop through sites
    
    null_list[[n]] <- null_comm
    message(paste0("iteration ", n, " completed"))
  }#finish iteration
  
  message(paste0(
    "Null model C5 completed (",
    iterations, " iterations). Species pool = ", pool,
    ". Presences randomised across sites and abundances chosen from species abundances in the observed community."
  ))
  
  return(null_list)
}

generate_C5_null_fast <- function(comm, iterations, pool) {
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
      richness <- sum(site > 0)
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
nullcomm_cells <- generate_C5_null_fast(abun_matrix, 999, pool = "entire")

#Calculate SES#
#we need to calculate RaoQ for each of the observed null communities
for(l in 1:length(nullcomm_cells)) {
  chosen_null <- nullcomm_cells[[l]]
  
  RQ_result <- calc_RaoQ(mean_traits = mean_traits, abun_matrix = chosen_null)
  RQ_result$counter <- paste0("null matrix ", l)
  
  if(l == 1) {
    null_RQ <-  RQ_result
  } else {
    null_RQ <- rbind(null_RQ, RQ_result)
  }}

#SES of each cell
RQ_cells_summary <- null_RQ |> 
  group_by(trait, cellref) |> 
  summarise(sd_null = sd(RaoQ), 
            mean_null = mean(RaoQ)) |> 
  filter(sd_null > 0) |> #cannot divide by zero in SES calculation
  inner_join(RQ_obs_cells, by = c("trait", "cellref")) |> 
  mutate(SES = (RaoQ - mean_null)/sd_null)

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


####SES at the grid scale####
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
nullcomm_grids <- generate_C5_null(abun_matrix_grid, 3, pool = "site")

#Calculate RaoQ for each of the observed null communities
for(l in 1:length(nullcomm_grids)) {
  chosen_null <- nullcomm_grids[[l]]
  
  RQ_result_grid <- calc_RaoQ(mean_traits = mean_traits, abun_matrix = chosen_null)
  RQ_result_grid$counter <- paste0("null matrix ", l)
  
  if(l == 1) {
    null_RQ_grid <-  RQ_result_grid
  } else {
    null_RQ_grid <- rbind(null_RQ_grid, RQ_result_grid)
  }}
  
#SES of each grid
RQ_grids_summary <- null_RQ_grid |> 
  group_by(trait, cellref) |> 
  summarise(sd_null = sd(RaoQ), 
            mean_null = mean(RaoQ)) |> 
  ungroup() |> 
  filter(sd_null > 0) |> #cannot divide by zero in SES calculation
  inner_join(RQ_obs_grid, by = c("trait", "cellref")) |> 
  mutate(SES = (RaoQ - mean_null)/sd_null)


