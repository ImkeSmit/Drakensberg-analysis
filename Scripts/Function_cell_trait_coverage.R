#=========================#
#====Cell_trait_coverage==#
#=========================#
####Function to identify cells that have trait measurements for less than 80% of the cover###

cell_trait_coverage <- function() {
  #import occurrence data without any sites or sp removed
  drak <-read.csv("All_data/clean_data/micro_climb_occurrence.csv", row.names = 1) 
  
  #import mean traits
  mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1) |> 
    rownames_to_column(var = "taxon")
  
  
  abundances <- drak |> 
    filter(cover>0) #table with covers
  mean_traits <- mean_traits #table with mean traits for species
  
  Cell_IDlist <- c(unique(abundances$Cell_ID))
  
  for(c in Cell_IDlist) {
    cell_abun <- abundances[abundances$Cell_ID ==c, ]
    
    merge <- cell_abun |> 
      left_join(mean_traits, by = "taxon")
    
    total_cov <- sum(merge$cover)
    
    trait_cov <- merge |> 
      filter(if_all(Height_cm:Thickness_mm, ~ !is.na(.)))
    trait_cov<- sum(trait_cov$cover)
    
    trait_coverage <- trait_cov/total_cov
  }
  
  if(c == Cell_IDlist[1]) {
    result <- cbind(c, trait_coverage)
  }else {
    temp<- cbind(c, trait_coverage)
    result<- rbind(result, temp)
  }
  return(result)
}