#=========================================================#
#Get the spatial autocorrelation structure for each grid#
#=========================================================#

#This function:
#Calculates Moran's I at differet lag distances
#Fits a negative exponential curve to the correlogram and calculates the parameters of the function:
#a = the value of MoranI at lag = 0, i.e. between adjacent cells
#b = decay rate. How rapidly spatial autocorrelation breaks down with increasing lag distance. 
#large b -> autocorrelation disappears quickly over short distances
#c = asymptote. The value that Moran's I tends to at large lag distances
#range dist = the lag distance at which autocorrelation has decayed to 5% of its initial value. 
#Beyond this distance, cells can be treated as independent

grid_correlation_structure <- function(grid_vector, #vector of grid names c("GG1", "GG2")
                                       data         #Gridded data with response and predictor variables
                                       formula,     #formula to use in lm() to get residuals, a character string "y~x1 +x2"
                                       k_specified) {   #how many neighbours to use in the calculation of Moran's I, an integer
#create list to store correlograms
grid_correlograms <- vector(mode = "list", length = 22)
names(grid_correlograms) <- grid_vector
#create decay df to store parameters of exponential function
decay_df <- data.frame(grid = grid_vector, b = NA, a = NA, c = NA, range_dist = NA, first_nonsig_lag = NA)

for(g in grid_vector) {
  cat("Processing grid:", g, "\n")
  
  #subset one grid
  one_grid_dat <- data |>  filter(grid == g)
  
  #linear model
  lm <- lm(as.formula(formula), 
           data = one_grid_dat)
  one_grid_resid <- c(resid(lm))
  
  
  #Build neighbour list within grid
  grid_coords <- cbind(one_grid_dat$x_coord, one_grid_dat$y_coord)
  k_local    <- min(k_specified, nrow(one_grid_dat) - 1) #make sure grid has enough cells to compute 4 nearest neighbours
  grid_neighbours <- knn2nb(knearneigh(grid_coords, k = k_local))
  
  #Generate correlogram
  max_order <- min(15, floor(nrow(one_grid_dat) / 5))
  one_grid_cor <- tryCatch(
    sp.correlogram(grid_neighbours, one_grid_resid, method = "I", order = max_order, 
                   zero.policy = T), #tolerate zero neighbour sets (e.g., cells with less than 4 neighbours)
    error = function(e) { cat("  Correlogram failed for grid", g, ":", e$message, "\n"); NULL }
  )
  
  if (is.null(one_grid_cor)) next
  grid_correlograms[[which(names(grid_correlograms) == g)]]<- one_grid_cor
  
  
  #Extract Moran's I values
  morans_df <- data.frame(
    lag   = 1:max_order,
    I     = one_grid_cor$res[, 1],
    lower = one_grid_cor$res[, 1] - 1.96 * sqrt(one_grid_cor$res[, 3]),
    upper = one_grid_cor$res[, 1] + 1.96 * sqrt(one_grid_cor$res[, 3]) 
  )
  #evaluate whether moranI not different from zero
  morans_df$nonsig_diff_from_zero = morans_df$lower < 0 & morans_df$upper > 0 
  
  # --- Fit negative exponential ---
  fit_g <- tryCatch(
    nls(
      I ~ a * exp(-b * lag) + c,
      data    = morans_df,
      start   = list(a = max(morans_df$I), b = 0.3, c = 0),
      control = nls.control(maxiter = 200)
    ),
    error = function(e) { cat("  NLS failed for grid", g, ":", e$message, "\n"); NULL }
  )
  
  if (is.null(fit_g)) next
  
  # --- Store results ---
  coefs_g <- coef(fit_g)
  decay_df[decay_df$grid == g, c("a", "b", "c")] <- coefs_g[c("a", "b", "c")]
  decay_df[decay_df$grid == g, "range_dist"]      <- -log(0.05) / coefs_g["b"]
  decay_df[decay_df$grid == g, "first_nonsig_lag"] <- which(morans_df$nonsig_diff_from_zero == TRUE)[1]
  #a = the value of MoranI at lag = 0, i.e. between adjacent cells
  #b = decay rate. How rapidly spatial autocorrelation breaks down with increasing lag distance. 
  #large b -> autocorrelation disappears quickly over short distances
  #c = asymptote. The value that Moran's I tends to at large lag distances
  #range dist = the lag distance at which autocorrelation has decayed to 5% of its initial value. 
  #Beyond this distance, cells can be treated as independent



# --- Plot all 22 correlograms in one figure ---
par(mfrow = c(4, 6), mar = c(2, 2, 2, 1))  # adjust layout to taste

  if (is.null(one_grid_cor)) next
  
  max_order <- nrow(one_grid_cor$res)
  
  # Base plot
  plot(morans_df$lag, morans_df$I,
       pch  = 16, cex = 0.7,
       xlab = "Lag", ylab = "Moran's I",
       ylim = c(min(morans_df$lower, 0), max(morans_df$upper)),
       main = paste("Grid", g))
  
  arrows(morans_df$lag, morans_df$lower,
         morans_df$lag, morans_df$upper,
         length = 0.03, angle = 90, code = 3, col = "grey60")
  
  abline(h = 0, lty = 3, col = "grey40")
  
  # Overlay fitted curve if available
  row_g <- decay_df[decay_df$grid == g, ]
  if (!is.na(row_g$b)) {
    lag_seq  <- seq(1, max_order, length.out = 200)
    I_pred   <- row_g$a * exp(-row_g$b * lag_seq) + row_g$c
    lines(lag_seq, I_pred, col = "firebrick", lwd = 1.5)
    abline(v = row_g$range_dist, lty = 2, col = "steelblue")
    
    p <- recordPlot()
  }
}

return(decay_df) #return decay df
return(p) #return plot
}#close function