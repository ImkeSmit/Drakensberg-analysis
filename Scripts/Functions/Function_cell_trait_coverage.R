#=========================#
#====Cell_trait_coverage==#
#=========================#
####Function to identify cells that have trait measurements for less than 80% of the cover###
#level = "grid" or "cell"

cell_trait_coverage <- function(level) {
  #import occurrence data without any sites or sp removed
  drak <-read.csv("All_data/clean_data/micro_climb_occurrence.csv", row.names = 1) |> 
    mutate(grid = paste0(site, grid))
  
  #import mean traits
  mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv") |>
    rename(taxon = X)
    
  
  
  abundances <- drak |> 
    filter(cover>0) #table with covers

  Cell_IDlist <- c(unique(abundances$Cell_ID))
  gridlist <- c(unique(abundances$grid))
  
  if(level == "cell"){
  
  for(c in Cell_IDlist) {
    cell_abun <- abundances[abundances$Cell_ID == c, ]
    
    merge <- cell_abun |> 
      left_join(mean_traits, by = "taxon")
    
    total_cov <- sum(merge$cover)
    
    trait_cov <- merge |> 
      filter(if_all(Height_cm:Thickness_mm, ~ !is.na(.)))
    trait_cov<- sum(trait_cov$cover)
    
    trait_coverage <- trait_cov/total_cov
    
    #What percent of species did we measure traits for?
    total_sp <- length(merge$taxon)
    
    trait_sp <- merge |> 
      filter(if_all(Height_cm:Thickness_mm, ~ !is.na(.))) |> 
      summarise(nsp = n())
    
    sp_coverage <- trait_sp$nsp/total_sp
    
    #put results in a table
    if(c == Cell_IDlist[1]) {
    result <- data.frame(Cell_ID = c, Trait_Coverage = as.numeric(trait_coverage), sp_coverage = as.numeric(sp_coverage))
      }else {
        temp<- data.frame(Cell_ID = c, Trait_Coverage = as.numeric(trait_coverage), sp_coverage = as.numeric(sp_coverage))
        result<- rbind(result, temp)
      }
    
  }}#finish loop through cells
  
  
  else{ #if level = "grid"
    for(g in gridlist) {
      grid_abun <- abundances[abundances$grid == g, ]
      
      merge <- grid_abun |> 
        left_join(mean_traits, by = "taxon")
      
      #what percent of vascular cover did we measure traits for?
      total_cov <- sum(merge$cover)
      
      trait_cov <- merge |> 
        filter(if_all(Height_cm:Thickness_mm, ~ !is.na(.)))
      trait_cov<- sum(trait_cov$cover)
      
      trait_coverage <- trait_cov/total_cov
      
      #What percent of species did we measure traits for?
      total_sp <- length(merge$taxon)
      
      trait_sp <- merge |> 
        filter(if_all(Height_cm:Thickness_mm, ~ !is.na(.))) |> 
        summarise(nsp = n())
      
      sp_coverage <- trait_sp$nsp/total_sp
      
      
      #put results in a table
      if(g == gridlist[1]) {
        result <- data.frame(grid = g, Trait_Coverage = as.numeric(trait_coverage), sp_coverage = as.numeric(sp_coverage))
      }else {
        temp<- data.frame(grid = g, Trait_Coverage = as.numeric(trait_coverage), sp_coverage = as.numeric(sp_coverage))
        result<- rbind(result, temp)
      }
    }#finish loop through grids
    }
  
  return(result)
}

#example use
#test<- cell_trait_coverage(level = "cell")
