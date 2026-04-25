# =============================================================================
# Ecological Grid: Fill Missing Values Using 8 Nearest Neighbours
# =============================================================================
# For each NA value in a numeric column, the function finds the cell's eight
# neighbouring cells (Moore neighbourhood) using x_coord and y_coord, then
# replaces the NA with the mean of all available (non-NA) neighbour values.
#
# Neighbours are identified within the same grid group only (e.g. GG_G1),
# so cells from different grids are never used as neighbours.
#
# Assumptions:
#   - Grid position is given by columns `x_coord` and `y_coord`
#   - The grid identifier is stored in a column called `grid`
#     (character, of the form GG1, GG2 ... GG7)
#   - Edge/corner cells are filled using however many real neighbours exist
#     (5 for edges, 3 for corners); cells where ALL neighbours are also NA
#     are left unfilled
# =============================================================================

library(dplyr)

impute_cells <- function(df, cols_to_impute) {
  
  # --- Input checks -----------------------------------------------------------
  missing_cols <- setdiff(cols_to_impute, names(df))
  if (length(missing_cols) > 0) {
    stop("The following columns are not in the data frame: ",
         paste(missing_cols, collapse = ", "))
  }
  
  non_numeric <- cols_to_impute[!sapply(df[cols_to_impute], is.numeric)]
  if (length(non_numeric) > 0) {
    stop("The following columns are not numeric: ",
         paste(non_numeric, collapse = ", "))
  }
  
  required_cols <- c("x_coord", "y_coord", "grid")
  missing_required <- setdiff(required_cols, names(df))
  if (length(missing_required) > 0) {
    stop("The following required columns are missing from the data frame: ",
         paste(missing_required, collapse = ", "))
  }
  
  # --- Helper: mean of 8-neighbour values for one cell -----------------------
  neighbour_mean <- function(target_x, target_y, grid_df, value_col) {
    
    offsets <- expand.grid(dx = -1:1, dy = -1:1) |>
      filter(!(dx == 0 & dy == 0))
    
    neighbour_vals <- grid_df |>
      inner_join(
        offsets |>
          mutate(x_coord = target_x + dx,
                 y_coord = target_y + dy),
        by = c("x_coord", "y_coord")
      ) |>
      pull(all_of(value_col))
    
    if (all(is.na(neighbour_vals))) return(NA_real_)
    
    mean(neighbour_vals, na.rm = TRUE)
  }
  
  # --- Imputation: loop over grids, then columns, then missing rows ----------
  df |>
    group_by(grid) |>
    group_modify(function(grid_df, grid_key) {
      
      for (col in cols_to_impute) {
        if (!any(is.na(grid_df[[col]]))) next   # skip if no NAs in this col
        
        missing_rows <- which(is.na(grid_df[[col]]))
        
        for (i in missing_rows) {
          grid_df[[col]][i] <- neighbour_mean(
            target_x  = grid_df$x_coord[i],
            target_y  = grid_df$y_coord[i],
            grid_df   = grid_df,
            value_col = col
          )
        }
      }
      grid_df
    }) |>
    ungroup()
}
