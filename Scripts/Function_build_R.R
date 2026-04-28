#============================================#
#===================Build R==================#
#============================================#

#This function build R, the correlation matrix, from the decay_df created by Function_grid_correlation_structure


build_R <- function(data,        #dataframe of cell scale variables containing columns grid, x_coord, row 
                    grid_params, #dataframe returned by Function_grid_correlation_structure
                    cutoff       #column in grid_params after which spatial autocorrelation is zero.
                                 #use "first_nonsig_lag" or "range_dist"
                    ){
#==========================================================
#Function: build one 160×160 correlation block for each grid
# Coordinates are integer grid steps; first_nonsig_lag is in the same units
make_grid_corr <- function(b, first_nonsig_lag, coords) {
  
  # Euclidean distance matrix in grid-step units
  dmat <- as.matrix(dist(coords[, c("x_coord", "row")])) #row is in steps of 1-20
  
  # Exponential correlation, zeroed beyond first_nonsig_lag steps
  corr <- exp(-b * dmat)
  corr[dmat > first_nonsig_lag] <- 0 #replace with range_dist to make it more conservative
  diag(corr) <- 1
  
  ###Fix negative eigenvalues###
  #For a matrix to be valid the variance explained along any axis should be positive, zero or negative variance is impossible
  #However, because of the truncation when building R (autocorrelation beyon nonsig lag = 0), negative eigenvalues can arise
  #E.g. Point A correlates with point B (close together)
  #Point B correlates with point C (close together)
  #But A→C distance exceeds the lag cutoff, so their correlation is forced to 0
  #This can produce negative eigenvalues
  
  
  # Iterative correction: alternate between
  #   (1) flooring negative off-diagonals to 0
  #   (2) flooring negative eigenvalues to small positive
  # until both conditions are satisfied simultaneously
  max_iter <- 20
  for (iter in seq_len(max_iter)) {
    
    # Step 1: fix any negative off-diagonal values
    corr[corr < 0] <- 0
    diag(corr)     <- 1
    
    # Step 2: check and fix positive definiteness
    eig <- eigen(corr, symmetric = TRUE)
    if (all(eig$values >= 1e-8)) break  # both conditions met, done
    
    eig$values <- pmax(eig$values, 1e-8)
    corr       <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    
    # Rescale diagonal back to 1
    d    <- 1 / sqrt(diag(corr))
    corr <- diag(d) %*% corr %*% diag(d)
  }
  
  if (iter == max_iter) warning("Convergence not reached for a grid block — inspect manually")
  
  corr
}#end make grid coor
#========================================================

# Run make_grid_corr for each grid
# Preserve the grid order as they appear in Hdat_filled
grid_order <- unique(data$grid)
stopifnot(all(grid_order %in% grid_params$grid))

block_list <- vector("list", length(grid_order))

for (i in seq_along(grid_order)) {
  g        <- grid_order[i]
  params   <- grid_params[grid_params$grid == g, ]
  coords_g <- data[data$grid == g, c("x_coord", "row")] #row is the y coord in steps of 1-20
  
  # Safety check: each grid must have exactly 160 points
  stopifnot(nrow(coords_g) == 160)
  
  block_list[[i]] <- make_grid_corr(
    b                = params$b,
    first_nonsig_lag = params[, cutoff], 
    #replace with range_dist to make it more conservative
    #using the first nonsignificant lag produces negative eigenvalues, especially in grids with slower decay
    coords           = coords_g
  )
}


##Build a block diagonal matrix from the list of matrices. 
#This places the autoccrelation structure of each grid along the diagonal of a large matrix, with xeroes everywhere else
#GG1-block  GG2-block  ...  BK7-block
#GG1  [ 160×160  |    0      |  0  |    0    ]
#GG2  [    0     | 160×160   |  0  |    0    ]
 #.   [    0     |    0      |  .  |    0    ]
#BK7  [    0     |    0      |  0  | 160×160 ]
R2 <- as.matrix(Matrix::bdiag(block_list))

#Sanity checks 
cat("R dimensions       :", dim(R2), "\n")          # 3520 × 3520
cat("All diagonal = 1    :", all(abs(diag(R2) - 1) < 1e-8), "\n")
cat("Symmetric           :", isSymmetric(R2), "\n")
cat("Min eigenvalue      :", min(eigen(R2, only.values = TRUE)$values), "\n")  # must be > 0
cat("Off-diagonal range  :", range(R2[row(R2) != col(R2)]), "\n")  # should be [0, <1)

return(R2)

} #end function
