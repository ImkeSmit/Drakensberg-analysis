# =============================================================================
# fill_missing_cells()
# =============================================================================
# Identifies missing x/y coordinate combinations in an ecological grid dataset
# and adds new rows with appropriate Cell_ID, site, and grid values.
# Measurement columns (SES, rock_cover, mean_soil_depth) are filled with NA.
#
# Usage:
#   df_complete <- fill_missing_cells(df)
# =============================================================================

library(tidyverse)

fill_missing_cells <- function(df) {
  
  # ---------------------------------------------------------------------------
  # Site/grid lookup: y_coord ranges → site and grid number
  # ---------------------------------------------------------------------------
  # Site GG: y 1–160,   grids 1–8, each spanning 20 y-units
  # Site WH: y 161–300, grids 1–7, each spanning 20 y-units
  # Site BK: y 301–440, grids 1–7, each spanning 20 y-units
  site_grid_lookup <- bind_rows(
    tibble(site = "GG", grid = 1:8, y_min = seq(1,   141, by = 20), y_max = seq(20,  160, by = 20)),
    tibble(site = "WH", grid = 1:7, y_min = seq(161, 281, by = 20), y_max = seq(180, 300, by = 20)),
    tibble(site = "BK", grid = 1:7, y_min = seq(301, 421, by = 20), y_max = seq(320, 440, by = 20))
  )
  
  # ---------------------------------------------------------------------------
  # Find missing x/y combinations
  # ---------------------------------------------------------------------------
  missing_coords <- expand_grid(x_coord = 1:8, y_coord = 1:440) %>%
    anti_join(df %>% select(x_coord, y_coord), by = c("x_coord", "y_coord"))
  
  cat("Missing coordinate combinations found:", nrow(missing_coords), "\n")
  
  if (nrow(missing_coords) == 0) {
    cat("Dataset is already complete — no rows added.\n")
    return(df)
  }
  
  # ---------------------------------------------------------------------------
  # Build new rows: assign site/grid and construct Cell_ID
  # ---------------------------------------------------------------------------
  new_rows <- missing_coords %>%
    left_join(site_grid_lookup, by = join_by(y_coord >= y_min, y_coord <= y_max)) %>%
    mutate(
      Cell_ID         = paste0(site, "_G", grid, "_", LETTERS[x_coord], y_coord),
      SES             = NA_real_,
      rock_cover      = NA_real_,
      mean_soil_depth = NA_real_
    ) %>%
    select(-y_min, -y_max)
  
  # ---------------------------------------------------------------------------
  # Combine, sort, and return
  # ---------------------------------------------------------------------------
  df_complete <- bind_rows(df, new_rows) %>%
    arrange(y_coord, x_coord)
  
  cat("Original rows:  ", nrow(df), "\n")
  cat("New rows added: ", nrow(new_rows), "\n")
  cat("Total rows:     ", nrow(df_complete), "\n")
  
  return(df_complete)
}


# =============================================================================
# Example usage
# =============================================================================

df <- tribble(
  ~Cell_ID,    ~site, ~grid, ~x_coord, ~y_coord, ~SES,   ~rock_cover, ~mean_soil_depth,
  "BK_G3_A1",  "BK",  3,     1,        341,      -0.373, 0,           38.3,
  "BK_G3_B1",  "BK",  3,     2,        341,      -0.386, 0,           37,
  "BK_G3_C1",  "BK",  3,     3,        341,      -0.417, 0,           37,
  "BK_G3_D1",  "BK",  3,     4,        341,      -0.398, 10,          19,
  "BK_G3_E1",  "BK",  3,     5,        341,      -0.522, 10,          23,
  "BK_G3_F1",  "BK",  3,     6,        341,      -0.731, 50,          9.33,
  "BK_G3_G1",  "BK",  3,     7,        341,       0.075, 35,          7.67,
  "BK_G3_H1",  "BK",  3,     8,        341,      -0.064, 0,           15.7,
  "BK_G3_A2",  "BK",  3,     1,        342,      -0.474, 15,          20.3,
  "BK_G3_B2",  "BK",  3,     2,        342,      -0.406, 2,           36
)

df_complete <- fill_missing_cells(df)
