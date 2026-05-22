# =============================================================================
# Spatial Sampling: Select Non-Autocorrelated Cells from an Ecological Grid
# =============================================================================
# Grid: 20 cols (x: 1-20) x 8 rows (y: 1-8), 1x1m cells
# Spatial independence threshold: cells must be > 4 lags apart (i.e., >4m)
# NA cells cannot be selected; adjacent cells are substituted instead.
# =============================================================================


select_independent_cells <- function(data,
                                     grid_var,        # character: grouping column (e.g. "grid_id")
                                     x,               # character: x-coordinate column (e.g. "x")
                                     y,               # character: y-coordinate column (e.g. "y")
                                     value_col,       # character: measurement column (e.g. "ndvi")
#                                     lag_threshold    = 4,   # min distance between selected cells
                                     max_search_radius = 3   # max search radius for NA substitutes
) {
  
  set.seed(42)  # For reproducibility
  
  # HELPER FUNCTIONS

  # Euclidean distance between two cells.
  # Accepts scalars or equal-length numeric vectors.
  cell_distance <- function(x1, y1, x2, y2) {
    x1 <- as.numeric(unlist(x1))
    y1 <- as.numeric(unlist(y1))
    x2 <- as.numeric(unlist(x2))
    y2 <- as.numeric(unlist(y2))
    sqrt((x1 - x2)^2 + (y1 - y2)^2)
  }
  
  # TRUE if the candidate cell is > lag_threshold away from every
  # already-selected cell.
  is_independent <- function(cx, cy, selected_df, lag_threshold) {
    if (nrow(selected_df) == 0) return(TRUE)
    dists <- cell_distance(cx, cy,
                           selected_df[[x]],   # [[ ]] with character name
                           selected_df[[y]])
    all(dists > lag_threshold)
  }
  
  # Find the nearest non-NA neighbour within max_radius of (cx, cy).
  # Returns a single-row data frame, or NULL if none found.
  find_adjacent_non_na <- function(cx, cy, grid_df, max_radius) {
    candidates <- grid_df |>
      filter(!is.na(.data[[value_col]])) |>          # .data[[]] for tidy eval
      mutate(dist = cell_distance(cx, cy,
                                  .data[[x]],
                                  .data[[y]])) |>
      filter(dist > 0, dist <= max_radius) |>
      arrange(dist)
    if (nrow(candidates) == 0) return(NULL)
    candidates[1, ]
  }
  
  ####Get lag distances for each grid
  ##First we have to impute cells otherwise correlation doesn't work
  data_filled <- impute_cells(df = data, 
                              cols_to_impute = colnames(data)[c(25, 31:35)])
  
  #get autocorrelation structure of each grid
  decay_df <- grid_correlation_structure(grid_vector = c(unique(data_filled$grid)), 
                                         data = data_filled, 
                                         formula = "SES ~ zrock_cover + znorthness + zsoil_moist + zsoil_depth + zslope_height", 
                                         k_specified = 4)
  
  mean_lag_distance = round(mean(decay_df$first_nonsig_lag, na.rm = T))
  
  #Replace lag distances >6 with 6 because that is as big as we can go
  decay_df<- decay_df |> 
    mutate(first_nonsig_lag = case_when(first_nonsig_lag > 5 ~ 6, 
                                        is.na(first_nonsig_lag) ~ mean_lag_distance, #replace NA's with the mean
                                        .default = first_nonsig_lag))
  
  
  #### MAIN LOOP — iterate over every unique grid

  grid_list <- unique(data[[grid_var]])
  all_results <- vector("list", length(grid_list))   # pre-allocate for speed
  
  for (g_idx in seq_along(grid_list)) {
    
    #lag threshold per grid
    g      <- grid_list[g_idx]
    lag_threshold <- decay_df[which(decay_df$grid == g), which(colnames(decay_df) == "first_nonsig_lag")]

    grid_df <- data |> filter(.data[[grid_var]] == g)
    
    # 1. Shuffle candidate cells 
    candidate_pool <- grid_df |> slice_sample(prop = 1)
    
    # 2. Greedy spatial sampling 
    selected_cells <- tibble()
    
    for (i in seq_len(nrow(candidate_pool))) {
      
      cell <- candidate_pool[i, ]
      
      # Case A: cell has a value and is spatially independent
      if (!is.na(cell[[value_col]]) &&
          is_independent(cell[[x]], cell[[y]], selected_cells, lag_threshold)) {
        selected_cells <- bind_rows(
          selected_cells,
          mutate(cell,
                 substituted = FALSE,
                 original_x  = NA_real_,
                 original_y  = NA_real_)
        )
        
      } else if (is.na(cell[[value_col]])) {
        # Case B: cell is NA — find nearest non-NA neighbour
        neighbour <- find_adjacent_non_na(cell[[x]], cell[[y]],
                                          grid_df, max_search_radius)
        
        if (!is.null(neighbour)) {
          # Ensure the neighbour hasn't already been selected
          already_picked <- nrow(selected_cells) > 0 &&
            any(selected_cells[[x]] == neighbour[[x]] &
                  selected_cells[[y]] == neighbour[[y]])
          
          if (!already_picked &&
              is_independent(neighbour[[x]], neighbour[[y]],
                             selected_cells, lag_threshold)) {
            selected_cells <- bind_rows(
              selected_cells,
              mutate(neighbour,
                     substituted = TRUE,
                     original_x  = as.numeric(cell[[x]]),
                     original_y  = as.numeric(cell[[y]]))
            )
          }
        }
      } # end if/else
    } # end cell loop
    
    
    # 3. Verify independence for this grid
    if (nrow(selected_cells) > 1) {
      min_dist <- min(combn(nrow(selected_cells), 2, function(idx) {
        cell_distance(selected_cells[[x]][idx[1]],
                      selected_cells[[y]][idx[1]],
                      selected_cells[[x]][idx[2]],
                      selected_cells[[y]][idx[2]])
      }))
      all_ok <- min_dist > lag_threshold
    }
    
    # Tag each row with its grid identifier, then store
    selected_cells[[grid_var]] <- g
    all_results[[g_idx]]       <- selected_cells
    
  } # end grid loop
  
  # COMBINE AND RETURN
  result <- bind_rows(all_results)
  #Check how many cells in each grid
  check<- result |> 
    group_by(grid) |> 
    summarise(n = n())
  print(check, n = 22)

    return(result)
} # end select_independent_cells()




