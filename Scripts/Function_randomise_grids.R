#==========================================#
###########RANDOMISE GRID###################
#==========================================#
#This function shuffles the values of one column within grids

randomise_grids <- function(data, var, iterations) {

  grids_22 <- c(unique(data$grid))
  
  #----single iteration-------#
  single_iter <- function(iter) {
  
    for(g in grids_22) {
      #select grid to randomise
      data <- as.data.frame(data)
      values <- data[which(data$grid == g),  which(colnames(data) %in% var)]
      Cell_ID <- data[which(data$grid == g),  which(colnames(data) == "Cell_ID")]
      subs <- cbind(Cell_ID, values)
      
      #order of values
      order <- c(1:nrow(subs))
      #new order
      new_order <- sample(order)
      #shuffle subs according to new order
      shuffled_var <- subs[order(new_order), which(colnames(subs) %in% var)]
      
      if(g == "GG1") {
        new_var <- cbind(Cell_ID, shuffled_var)
      } else{
      temp <- cbind(Cell_ID, shuffled_var)
      new_var <- rbind(new_var, temp)
      }
      #remove unshuffled columns
      data_min_var <- data[, -which(colnames(data) %in% var)] 
      #join shuffled columns
      randomised_dat <- data_min_var|> 
        left_join(new_var, by = "Cell_ID")
    }#end loop through grids
    return(randomised_dat)
  }#end single_iter
  
  # ---- Run iterations ----#
  randomised_list <- vector("list", iterations)
  for (n in seq_len(iterations)) {
    randomised_list[[n]] <- single_iter(n)
  }
  
  return(randomised_list)
}#end function

#example use
#test <- randomise_grids(data = Hdat_filled, 
#                        var = "SES",
#                        iterations = 10)
