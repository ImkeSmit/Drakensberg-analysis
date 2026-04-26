####Function to get the neighbours of cells belonging to the same grid####
library(spdep)

# Assuming Hdat_filled has a column called 'grid' identifying grid membership

# Build a neighbour list restricted to same-grid cells
same_grid_neighbours <- function(data, #df with SES and environmental data
                         grid_col, #column identifying grid membership, text
                         coord_cols, #columns with the coordinates, text vector
                         k = 4) {    #How many cells to use in moranI calculation
  
  grids <- unique(data[[grid_col]])
  n <- nrow(data)
  nb_list <- vector("list", n)
  
  for (g in grids) {
    # Row indices belonging to this grid
    idx <- which(data[[grid_col]] == g)
    
    if (length(idx) < 2) {
      # Isolated cell — no neighbours
      for (i in idx) nb_list[[i]] <- integer(0)
      next
    }
    
    # Coordinates of cells in this grid only
    sub_coords <- as.matrix(data[idx, coord_cols])
    
    # k cannot exceed number of available neighbours in the grid
    k_local <- min(k, length(idx) - 1)
    
    # Find k nearest neighbours within the grid
    knn <- knearneigh(sub_coords, k = k_local)$nn  # matrix: rows = cells, cols = neighbour ranks
    
    # Map local indices back to global row indices and populate nb_list
    for (local_i in seq_along(idx)) {
      global_i <- idx[local_i]
      nb_list[[global_i]] <- idx[knn[local_i, ]]
    }
  }
  
  # Convert to nb class — nb_list entries must be sorted integer vectors
  nb_list <- lapply(nb_list, function(x) sort(as.integer(x)))
  class(nb_list) <- "nb"
  attr(nb_list, "region.id") <- as.character(seq_len(n))
  attr(nb_list, "call") <- match.call()
  nb_list
}

# Usage
coords_cols <- c("x_coord", "y_coord")

neibs <- same_grid_nb(Hdat_filled, grid_col = "grid", coord_cols = coords_cols, k = 8)

# Validate (optional but recommended)
summary(neibs)

# Compute Moran's I correlogram as before
spcor <- sp.correlogram(neibs, lm_resid, method = "I", order = 15)
plot(spcor, ylim = c(0, 1))