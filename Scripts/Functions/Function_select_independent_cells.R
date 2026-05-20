# =============================================================================
# Spatial Sampling: Select Non-Autocorrelated Cells from an Ecological Grid
# =============================================================================
# Grid: 20 cols (x: 1-20) x 8 rows (y: 1-8), 1x1m cells
# Spatial independence threshold: cells must be > 4 lags apart (i.e., >4m)
# NA cells cannot be selected; adjacent cells are substituted instead.
# =============================================================================

library(tidyverse)


select_independent_cells <- function(data,
                                     grid_var,        # character: grouping column (e.g. "grid_id")
                                     x,               # character: x-coordinate column (e.g. "x")
                                     y,               # character: y-coordinate column (e.g. "y")
                                     value_col,       # character: measurement column (e.g. "ndvi")
                                     lag_threshold    = 4,   # min distance between selected cells
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
  
  # MAIN LOOP — iterate over every unique grid

  grid_list <- unique(data[[grid_var]])
  all_results <- vector("list", length(grid_list))   # pre-allocate for speed
  
  for (g_idx in seq_along(grid_list)) {
    
    g      <- grid_list[g_idx]
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
  return(result)
} # end select_independent_cells()

Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA) #also select cells which have no SES measurement. This is necessary to make the grid complete
         ) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))

test<- select_independent_cells2(Hdat, grid_var = "grid", x = "x_coord", y = "y_coord", value_col = "SES",
                          max_search_radius = 3, lag_threshold = 4)



# -----------------------------------------------------------------------------
# 7. VISUALISE THE GRID AND SELECTED CELLS
# -----------------------------------------------------------------------------

# Build a plot-ready version of the full grid
plot_df <- grid_df %>%
  mutate(status = case_when(
    is.na(value) ~ "NA cell",
    TRUE         ~ "Available"
  ))

# Mark selected and substituted cells
selected_plot <- selected_cells %>%
  mutate(status = if_else(substituted, "Substituted (was NA)", "Selected"))

# Replace matching rows in plot_df
plot_df <- plot_df %>%
  rows_update(selected_plot %>% select(x, y, status), by = c("x", "y"))

ggplot() +
  # Full grid as background tiles
  geom_tile(data = plot_df,
            aes(x = x, y = y, fill = status),
            colour = "white", linewidth = 0.4) +
  # Arrows showing NA → substitute substitutions
  geom_segment(data = selected_cells %>% filter(substituted),
               aes(x = original_x, y = original_y,
                   xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "Available"             = "#d0e8d0",
    "NA cell"               = "#f0b8b8",
    "Selected"              = "#2a7d4f",
    "Substituted (was NA)"  = "#e8a000"
  )) +
  scale_x_continuous(breaks = 1:20) +
  scale_y_continuous(breaks = 1:8) +
  coord_equal() +
  labs(
    title    = "Spatially Independent Cell Selection",
    subtitle = paste0("Lag threshold: >4 m | Selected: ", nrow(selected_cells),
                      " cells (", sum(selected_cells$substituted), " substituted)"),
    x        = "x (m)",
    y        = "y (m)",
    fill     = "Cell status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("spatial_sampling_plot.png", width = 12, height = 5, dpi = 150)
cat("\nPlot saved to: spatial_sampling_plot.png\n")

