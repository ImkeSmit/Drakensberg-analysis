# =============================================================================
# Spatial Sampling: Select Non-Autocorrelated Cells from an Ecological Grid
# =============================================================================
# Grid: 20 cols (x: 1-20) x 8 rows (y: 1-8), 1x1m cells
# Spatial independence threshold: cells must be > 4 lags apart (i.e., >4m)
# NA cells cannot be selected; adjacent cells are substituted instead.
# =============================================================================

library(tidyverse)

set.seed(42)  # For reproducibility

# -----------------------------------------------------------------------------
# 1. SIMULATE GRID DATA
#    Replace this section with your actual data import, e.g.:
#    grid_df <- read.csv("your_data.csv")   # must have columns: x, y, value
# -----------------------------------------------------------------------------

grid_df <- expand.grid(x = 1:20, y = 1:8) %>%
  mutate(value = rnorm(n()))

# Introduce some NA cells to demonstrate NA-handling logic
na_indices <- sample(nrow(grid_df), size = 20)
grid_df$value[na_indices] <- NA

cat("Grid dimensions:", max(grid_df$x), "x", max(grid_df$y), "\n")
cat("Total cells:", nrow(grid_df), "\n")
cat("NA cells:", sum(is.na(grid_df$value)), "\n\n")

Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA) #also select cells which have no SES measurement. This is necessary to make the grid complete
 ) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))

check <- Hdat |> group_by(site, grid) |> 
  summarise(n = n())
 data = Hdat
 grid_var = "grid"
 x = "x_coord"
 y = "y_coord"
####Start function
select_independent_cells <- function(data, grid_var, x, y, lag_threshold, max_search_radius, 
                                     value_col){

  
  grid_list<- c(unique(working_dat[, grid_var]))[[1]]
  
  for(g in grid_list) {
    grid_df <- data |> filter(grid == g)

# -----------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Euclidean distance between two cells (in lag units = metres here)
cell_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

#' Check whether a candidate cell is spatially independent of all
#' already-selected cells (distance must be > lag_threshold).
is_independent <- function(cx, cy, selected_df, lag_threshold = 4) {
  if (nrow(selected_df) == 0) return(TRUE)
  dists <- cell_distance(cx, cy, selected_df[, x], selected_df[, y])
  all(dists > lag_threshold)
}

#' Return the valid (non-NA) neighbours of a cell, ordered by proximity.
#' Neighbours are searched in an expanding radius until one is found or
#' the search radius exceeds the grid bounds.
find_adjacent_non_na <- function(cx, cy, grid_df, max_radius = 3) {
  candidates <- grid_df %>%
    filter(!is.na(value)) %>%
    mutate(dist = cell_distance(cx, cy, x, y)) %>%
    filter(dist > 0, dist <= max_radius) %>%         # exclude the cell itself
    arrange(dist)                                      # nearest first
  if (nrow(candidates) == 0) return(NULL)
  candidates[1, ]                                      # return closest neighbour
}

# -----------------------------------------------------------------------------
# 3. GENERATE CANDIDATE CELLS IN RANDOM ORDER
#    Randomise the visit order so the selected set is not biased toward any
#    corner of the grid.
# -----------------------------------------------------------------------------

candidate_pool <- grid_df %>%
  slice_sample(prop = 1)   # shuffle all cells

# -----------------------------------------------------------------------------
# 4. GREEDY SPATIAL SAMPLING
#    Walk through the shuffled candidate pool. For each cell:
#      a) If it has a value and is independent of selected cells → select it.
#      b) If it is NA → try to substitute a nearby non-NA neighbour that is
#         also independent of already-selected cells.
#      c) Otherwise → skip.
# -----------------------------------------------------------------------------

LAG_THRESHOLD <- lag_threshold          # cells must be > 4 m apart
MAX_SEARCH_RADIUS <- max_search_radius  # max metres to search for a NA substitute

selected_cells <- data.frame()   # accumulates the final selection

for (i in seq_len(nrow(candidate_pool))) {
  
  cell <- candidate_pool[i, ]
  
  # ---- Case A: cell has a value ----
  if (!is.na(cell[, value_col])) {
    if (is_independent(cell[, x], cell[, y], selected_cells, LAG_THRESHOLD)) {
      selected_cells <- bind_rows(selected_cells,
                                  mutate(cell, substituted = FALSE,
                                         original_x = NA, original_y = NA))
    }
    next #move to next i
  }
  
  # ---- Case B: cell is NA – look for an adjacent substitute ----
  neighbour <- find_adjacent_non_na(cell[, x], cell[, y], grid_df, MAX_SEARCH_RADIUS)
  
  if (is.null(neighbour)) next   # no non-NA neighbour found within radius
  
  # Check the substitute is not already selected and is spatially independent
  already_picked <- nrow(selected_cells) > 0 &&
    any(selected_cells$x == neighbour$x & selected_cells$y == neighbour$y)
  
  if (!already_picked &&
      is_independent(neighbour$x, neighbour$y, selected_cells, LAG_THRESHOLD)) {
    selected_cells <- bind_rows(selected_cells,
                                mutate(neighbour,
                                       substituted   = TRUE,
                                       original_x    = cell[, x],
                                       original_y    = cell[, y]))
  }
}

# -----------------------------------------------------------------------------
# 5. RESULTS SUMMARY
# -----------------------------------------------------------------------------

cat("=== Sampling Results ===\n")
cat("Selected cells (total)  :", nrow(selected_cells), "\n")
cat("Directly selected       :", sum(!selected_cells$substituted), "\n")
cat("NA-substituted          :", sum( selected_cells$substituted), "\n\n")

print(selected_cells %>% select(x, y, value, substituted, original_x, original_y))

# -----------------------------------------------------------------------------
# 6. VERIFY SPATIAL INDEPENDENCE
#    Confirm that every pair of selected cells exceeds the lag threshold.
# -----------------------------------------------------------------------------

if (nrow(selected_cells) > 1) {
  pairwise_ok <- combn(nrow(selected_cells), 2, function(idx) {
    cell_distance(selected_cells$x[idx[1]], selected_cells$y[idx[1]],
                  selected_cells$x[idx[2]], selected_cells$y[idx[2]]) > LAG_THRESHOLD
  })
  
  cat("\n=== Independence Check ===\n")
  cat("All pairs > 4 lags apart:", all(pairwise_ok), "\n")
  cat("Min pairwise distance   :", round(
    min(combn(nrow(selected_cells), 2, function(idx)
      cell_distance(selected_cells$x[idx[1]], selected_cells$y[idx[1]],
                    selected_cells$x[idx[2]], selected_cells$y[idx[2]]))), 3), "m\n")
}

}#end loop through grids
  return(result)
}#end function

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

# -----------------------------------------------------------------------------
# 8. EXPORT SELECTED CELLS
# -----------------------------------------------------------------------------

write.csv(selected_cells, "selected_cells.csv", row.names = FALSE)
cat("Selected cells saved to: selected_cells.csv\n")