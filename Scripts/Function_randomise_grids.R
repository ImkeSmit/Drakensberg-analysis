#==========================================#
###########RANDOMISE GRID###################
#==========================================#
#This function shuffles the values of one column within grids

randomise_grids <- function(data, var) {

  grids_22 <- c(unique(data$grid)
  
  for(g in grids_22) {
    #select grid to randomise
    data <- as.data.frame(data)
    values <- data[which(data$grid == g),  which(colnames(data) == var)]
    rnames <- data[which(data$grid == g),  which(colnames(data) == "Cell_ID")]
    subs <- data.frame(Cell_ID = rnames, SES = values)
    
    set.seed(39)
    
    #order of SES values
    order <- c(1:nrow(subs))
    #new order
    new_order <- sample(order)
    #shuffle subs according to new order
    shuffled_var <- subs[order(new_order), which(colnames(subs) == var)]
    
    if(g == "GG1") {
      new_var <- data.frame(Cell_ID = rnames, randomised_SES = shuffled_var)
    }
    
    temp <- data.frame(Cell_ID = rnames, randomised_SES = shuffled_var)
    new_var <- rbind(new_var, temp)
    
    randomised_dat <- data |> 
      left_join(new_var, by = "Cell_ID")
  }
}#end function

