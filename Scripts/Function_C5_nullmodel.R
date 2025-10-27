#This is a script to generate the C5 null model
#Species abundances(including 0) are shuffled within the columns. 
#Therefore, the frequencies of species, and their total abundances stay the same
#The species occurring in each site change, but the richness of each site stays the same

library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)
library(FD)


###improved function
generate_C5_null_imp <- function(comm, iterations = 10, pool = "entire") {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Please install the 'vegan' package: install.packages('vegan')")
  }
  
  comm <- as.matrix(comm)
  rownames_comm <- rownames(comm)
  colnames_comm <- colnames(comm)
  
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
  
  # ---- Single iteration ----
  single_iter <- function(iter) {
    null_comm <- comm * 0
    
    if (pool == "entire") {
      #swap abundances among the sites within a column 
      perm <- vegan::permatswap(comm, method = "swsh", shuffle ="samp", fixedmar = "columns", mtype = "count", times = 1)
      null_comm <- perm$perm[[1]]
      
    } else {
      site_groups <- unique(sapply(rownames_comm, function(x)
        get_sites_in_pool(x)[1]))
      for (sg in site_groups) {
        group_sites <- get_sites_in_pool(sg)
        group_rows <- which(rownames_comm %in% group_sites)
        submat <- comm[group_rows, , drop = FALSE]
        
        perm_sub <- vegan::permatswap(submat, method = "swsh", shuffle ="samp", fixedmar = "columns", mtype = "count", times = 1)
        null_comm[group_rows, ] <- perm_sub$perm[[1]]
      }
    }
    
    # Replace NAs (if any) with zeros
    null_comm[is.na(null_comm)] <- 0
    
    null_comm
  }
  
  # ---- Run iterations ----
  null_list <- vector("list", iterations)
  for (n in seq_len(iterations)) {
    null_list[[n]] <- single_iter(n)
  }
  
  message("✅ Finished generating ", iterations, " C5 null model(s).")
  return(null_list)
}

###test null model####
##Import community and trait data
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

testnull <- generate_C5_null_imp(abun_matrix, 10, pool = "entire")

colSums(abun_matrix)
colSums(testnull[[1]]) #colsums are the same, that is correct

rowSums(abun_matrix)
rowSums(testnull[[1]]) #rowsums differ, that is correct

sum(colSums(abun_matrix))
sum(colSums(testnull[[1]])) #total matrix sum is the same, that is correct

specnumber(abun_matrix)
specnumber(testnull[[1]]) #richness of plots are the same, that is correct
