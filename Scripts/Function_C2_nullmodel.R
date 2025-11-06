#This is a script to generate the C2 null model
#Species abundances(including 0) are shuffled within the rows. 
#Therefore, the species richness and total abundances of samples stay the same, 
#but the frequency of species and species total abundances do not

library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(vegan)
library(FD)

###C2 null model function
generate_C2_null <- function(comm, iterations = 10, pool = "entire") {
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
  
  #function to get the sites (plots) that are in teh specified pool
  get_sites_in_pool <- function(site_name) {
    if (pool == "entire") return(rownames_comm)
    if (pool == "site") {
      site_prefix <- three_sites[str_detect(site_name, three_sites)]
      return(rownames_comm[grepl(site_prefix, rownames_comm)])
    }
    if (pool == "grid") {
      grid_prefix <- grids_22[str_detect(site_name, grids_22)]
      return(rownames_comm[grepl(grid_prefix, rownames_comm)])
    }
    stop("Invalid pool type")
  }
  
  # ---- Single iteration ----
  single_iter <- function(iter) {
    null_comm <- comm * 0
    
    if (pool == "entire") {
      #swap abundances among the sites within a column 
      perm <- vegan::permatswap(comm, method = "swsh", shuffle ="samp", fixedmar = "rows", mtype = "count", times = 1)
      null_comm <- perm$perm[[1]]
      
    } else {
      site_groups <- unique(sapply(rownames_comm, function(x)
        get_sites_in_pool(x)[1]))
      #if pool not entire, get the sites that are in the pool for each plot
      for (sg in site_groups) {
        group_sites <- get_sites_in_pool(sg)
        group_rows <- which(rownames_comm %in% group_sites)
        submat <- comm[group_rows, , drop = FALSE]
        
        perm_sub <- vegan::permatswap(submat, method = "swsh", shuffle ="samp", fixedmar = "rows", mtype = "count", times = 1)
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
  
  message("✅ Finished generating ", iterations, " C2 null models")
  return(null_list)
}

