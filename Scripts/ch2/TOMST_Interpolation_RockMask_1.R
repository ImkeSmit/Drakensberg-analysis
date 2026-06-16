# =============================================================================
# TOMST MICROCLIMATE INTERPOLATION PIPELINE with rock cover masking
#
# Change from previous version:
#   Cells with rock_cover >= ROCK_THRESHOLD are excluded from both:
#     (a) IDW target points receive NA, not an interpolated value
#     (b) Logger source points if a logger sits on bare rock, it is
#         excluded so it doesn't pull the interpolation surface toward
#         conditions unrepresentative of vegetated cells
#
#   This is the same logic used for the campaign moisture NAs bare rock
#   cells are biologically uninformative for soil variables and their
#   inclusion was dragging down cross-validation R� for moisture.
#
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(gstat)
library(sf)
library(ggplot2)

setwd("C:/Users/User/Documents/Academics 5/Env data merge")  
# =============================================================================
# SETTINGS
# =============================================================================

SITE_MAP       <- c("BOK"="BK", "GOL"="GG", "WIT"="WH")
IMP_THRESHOLD  <- 0.20
ROCK_THRESHOLD <- 40    # cells with rock_cover >= this are masked out
                         # consistent with the campaign moisture NA logic
R2_THRESHOLD   <- 0.10  # poor grid fallback
IDP            <- 2      # IDW power parameter: p=2 is standard 
                          # higher p = more local influence, lower p = smoother surface
NMAX           <- Inf    # use all loggers per grid 

# TOMST_VARS and CELL_LEVEL_VARS are defined in Section 2 after extraction
# All variables are attempted at cell level 
# =============================================================================
# SECTION 1 Load data 
# =============================================================================

tomst_raw  <- read.csv("microclimate_indices.csv", stringsAsFactors=FALSE)
centroids  <- read.csv("All_Sites_Cell_Centroids.csv", stringsAsFactors=FALSE)
env_data   <- read.csv("All_Sites_Environmental_Data.csv", stringsAsFactors=FALSE)

# Rock cover lookup used for masking throughout
rock_lookup <- env_data %>%
  select(Cell_ID, rock_cover) %>%
  mutate(high_rock = rock_cover >= ROCK_THRESHOLD)

cat("High-rock cells (rock_cover >=", ROCK_THRESHOLD, "%):",
    sum(rock_lookup$high_rock, na.rm=TRUE), "of", nrow(rock_lookup), "\n")

cat("By site:\n")
rock_lookup %>%
  mutate(site = sub("_G.*", "", Cell_ID)) %>%
  group_by(site) %>%
  summarise(n_high_rock = sum(high_rock, na.rm=TRUE),
            pct = round(100*mean(high_rock, na.rm=TRUE), 1),
            .groups="drop") %>%
  print()

# =============================================================================
# SECTION 2 Extract and pivot TOMST data
# Raw file columns: site, Site, Grid, Column, Row, sensor, period, stat,
#                   value, impprop
# Site codes in raw file: BOK / GOL / WIT  ->  mapped to BK / GG / WH
# =============================================================================

TOMST_VARS <- tribble(
  ~sensor,  ~period,   ~stat,   ~output_name,
  "T1",     "annual",  "avg",   "T1_mean_annual",
  "T1",     "winter",  "min",   "T1_winter_min",
  "T1",     "annual",  "gdh5",  "T1_gdh5",
  "T1",     "annual",  "fdh",   "T1_fdh",
  "T3",     "annual",  "avg",   "T3_mean_annual",
  "T3",     "winter",  "min",   "T3_winter_min",
  "T3",     "annual",  "gdh5",  "T3_gdh5",
  "T3",     "annual",  "fdh",   "T3_fdh",
  "moist",  "annual",  "avg",   "moist_mean_annual",
  "moist",  "annual",  "cv",    "moist_cv",
  "moist",  "summer",  "avg",   "moist_summer_avg"
)

# Extract each variable, applying the imputation filter per logger
tomst_extracted <- TOMST_VARS %>%
  pmap_dfr(function(sensor, period, stat, output_name) {
    tomst_raw %>%
      filter(
        .data$sensor == !!sensor,
        .data$period == !!period,
        .data$stat   == !!stat
      ) %>%
      mutate(
        output_name = output_name,
        value_clean = ifelse(impprop > IMP_THRESHOLD, NA_real_, value)
      ) %>%
      select(Site, Grid, Column, Row, output_name, value = value_clean)
  })

# Pivot to wide format and build Cell_ID
tomst_wide <- tomst_extracted %>%
  pivot_wider(names_from = output_name, values_from = value) %>%
  mutate(
    site_code = SITE_MAP[Site],
    Grid      = as.integer(Grid),
    Column    = toupper(trimws(Column)),
    Row       = as.integer(Row),
    Cell_ID   = paste0(site_code, "_G", Grid, "_", Column, Row)
  )

# All variables attempted at cell level 
CELL_LEVEL_VARS <- TOMST_VARS$output_name

cat("\nLoggers after extraction:", nrow(tomst_wide), "\n")
cat("Sample Cell_IDs:\n")
print(head(tomst_wide %>% select(site_code, Grid, Column, Row, Cell_ID), 5))

# =============================================================================
# SECTION 3 Join logger positions to coordinates
# =============================================================================

tomst_spatial <- tomst_wide %>%
  left_join(centroids %>% select(Cell_ID, lat, lon), by="Cell_ID") %>%
  # Flag loggers sitting on high-rock cells
  left_join(rock_lookup %>% select(Cell_ID, high_rock), by="Cell_ID")

n_rock_loggers <- sum(tomst_spatial$high_rock, na.rm=TRUE)
cat("Loggers on high-rock cells (excluded from interpolation):",
    n_rock_loggers, "\n\n")

# =============================================================================
# SECTION 4 IDW interpolation function WITH rock masking
# =============================================================================

interpolate_grid_idw <- function(site_code, grid_num, tomst_sp,
                                 all_centroids, rock_lkp, var_name, idp=IDP) {

  # Source points: logger cells with valid values AND not high-rock
  source_pts <- tomst_sp %>%
    filter(site_code == !!site_code,
           Grid      == grid_num,
           (is.na(high_rock) | !high_rock)) %>%  # exclude confirmed rock cells; NA rock_cover treated as non-rock
    select(Cell_ID, lat, lon, value=all_of(var_name)) %>%
    filter(!is.na(value), !is.na(lat), !is.na(lon))

  # Minimum 2 loggers to fit IDW (gstat requires at least 2 source points).
  # Validation uses min 3 so one can be held out while 2 remain for fitting.
  if (nrow(source_pts) < 2) {
    message("  Skipping ", site_code, " G", grid_num, " ", var_name,
            " — fewer than 2 valid non-rock loggers")
    return(NULL)
  }

  # Target points: all cells in this grid that are NOT high-rock
  target_pts <- all_centroids %>%
    filter(site == site_code, grid == grid_num) %>%
    select(Cell_ID, lat, lon) %>%
    left_join(rock_lkp %>% select(Cell_ID, high_rock), by="Cell_ID") %>%
    filter(is.na(high_rock) | !high_rock) %>%  # exclude confirmed rock cells
    select(Cell_ID, lat, lon)

  if (nrow(target_pts) == 0) return(NULL)

  # Project to UTM 36S for accurate distance calculation
  source_sf <- st_as_sf(source_pts, coords=c("lon","lat"), crs=4326) %>%
    st_transform(32736)
  target_sf <- st_as_sf(target_pts, coords=c("lon","lat"), crs=4326) %>%
    st_transform(32736)

  source_coords <- st_coordinates(source_sf)
  target_coords <- st_coordinates(target_sf)

  # st_transform preserves row order, so source_coords rows align with source_pts rows
  source_df <- data.frame(X=source_coords[,1], Y=source_coords[,2],
                           value=source_pts$value)
  target_df <- data.frame(X=target_coords[,1], Y=target_coords[,2])

  # Cap predictions to observed logger range (prevent extrapolation)
  obs_min <- min(source_df$value)
  obs_max <- max(source_df$value)

  idw_model <- gstat(formula=value~1, locations=~X+Y,
                     data=source_df, nmax=nrow(source_df), set=list(idp=IDP))
  predicted <- predict(idw_model, newdata=target_df)

  result <- target_pts %>%
    mutate(!!var_name := pmax(pmin(predicted$var1.pred, obs_max), obs_min))

  return(result)
}

# =============================================================================
# SECTION 5 Run interpolation
# =============================================================================

cat("=== Running IDW interpolation with rock masking ===\n")

site_grids <- tomst_spatial %>%
  distinct(site_code, Grid) %>%
  arrange(site_code, Grid)

all_interpolated <- map_dfr(CELL_LEVEL_VARS, function(var) {
  cat("Interpolating:", var, "\n")
  map2_dfr(site_grids$site_code, site_grids$Grid, function(sc, gn) {
    interpolate_grid_idw(sc, gn, tomst_spatial, centroids,
                         rock_lookup, var)
  })
})

tomst_interpolated <- all_interpolated %>%
  group_by(Cell_ID) %>%
  summarise(across(all_of(CELL_LEVEL_VARS),
                   ~first(.[!is.na(.)])),
            .groups="drop")

cat("\nInterpolated surface:", nrow(tomst_interpolated), "cells\n")

# High-rock cells will be absent from tomst_interpolated — add them back as NA
# so the output has all cells and rock cells simply have NA for TOMST variables
all_cell_ids <- centroids %>%
  filter(site %in% unique(tomst_spatial$site_code)) %>%
  select(Cell_ID)

tomst_interpolated <- all_cell_ids %>%
  left_join(tomst_interpolated, by="Cell_ID")

cat("After adding back rock cells as NA:", nrow(tomst_interpolated), "cells\n")
cat("\nNA counts per variable:\n")
tomst_interpolated %>%
  summarise(across(all_of(CELL_LEVEL_VARS), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to="variable", values_to="n_NA") %>%
  print()

# =============================================================================
# SECTION 6 Leave-one-logger-out validation (rock masking applied)
# =============================================================================

cat("\n=== Leave-one-logger-out validation (non-rock loggers only) ===\n")

run_loloo_validation <- function(site_code, grid_num, tomst_sp,
                                  var_name, idp=IDP) {

  # Only use non-rock loggers in validation
  grid_loggers <- tomst_sp %>%
    filter(site_code == !!site_code,
           Grid      == grid_num,
           (is.na(high_rock) | !high_rock)) %>%
    select(Cell_ID, lat, lon, value=all_of(var_name)) %>%
    filter(!is.na(value), !is.na(lat), !is.na(lon))

  if (nrow(grid_loggers) < 3) return(NULL)

  map_dfr(seq_len(nrow(grid_loggers)), function(i) {
    held_out <- grid_loggers[i, ]
    training <- grid_loggers[-i, ]

    train_sf <- st_as_sf(training, coords=c("lon","lat"), crs=4326) %>%
      st_transform(32736)
    held_sf  <- st_as_sf(held_out, coords=c("lon","lat"), crs=4326) %>%
      st_transform(32736)

    # st_transform preserves row order — coordinates align with original data rows
    train_df <- data.frame(X=st_coordinates(train_sf)[,1],
                           Y=st_coordinates(train_sf)[,2],
                           value=training$value)
    held_df  <- data.frame(X=st_coordinates(held_sf)[,1],
                           Y=st_coordinates(held_sf)[,2])

    obs_min <- min(training$value)
    obs_max <- max(training$value)

    idw_mod <- gstat(formula=value~1, locations=~X+Y,
                     data=train_df, nmax=nrow(train_df), set=list(idp=IDP))
    pred    <- predict(idw_mod, newdata=held_df)
    pred_val <- pmax(pmin(pred$var1.pred, obs_max), obs_min)

    tibble(site_code=site_code, Grid=grid_num, variable=var_name,
           Cell_ID=held_out$Cell_ID, observed=held_out$value,
           predicted=pred_val, error=pred_val - held_out$value)
  })
}

validation_results <- map_dfr(CELL_LEVEL_VARS, function(var) {
  cat("Validating:", var, "\n")
  map2_dfr(site_grids$site_code, site_grids$Grid, function(sc, gn) {
    run_loloo_validation(sc, gn, tomst_spatial, var)
  })
})

# =============================================================================
# SECTION 7 Validation summary
# =============================================================================

validation_grid <- validation_results %>%
  group_by(site_code, Grid, variable) %>%
  summarise(
    n        = n(),
    R2       = ifelse(var(observed) > 0,
                      cor(observed, predicted)^2, NA_real_),
    RMSE     = sqrt(mean(error^2)),
    bias     = mean(error),
    .groups  = "drop"
  )

cat("\n=== Validation R² by variable (non-rock loggers only) ===\n")
validation_grid %>%
  group_by(variable) %>%
  summarise(
    mean_R2  = round(mean(R2, na.rm=TRUE), 3),
    n_poor   = sum(R2 < R2_THRESHOLD, na.rm=TRUE),
    n_grids  = n(),
    .groups  = "drop"
  ) %>%
  arrange(desc(mean_R2)) %>%
  print()

# Compare moisture R² with vs without masking for reporting
cat("\nMoisture grid-level R² (all grids):\n")
validation_grid %>%
  filter(grepl("moist", variable)) %>%
  select(site_code, Grid, variable, R2) %>%
  arrange(variable, site_code, Grid) %>%
  print(n=Inf)

# =============================================================================
# SECTION 8 Poor grid fallback 
# =============================================================================

poor_grids <- validation_grid %>%
  filter(R2 < R2_THRESHOLD | is.na(R2)) %>%
  select(site_code, Grid, variable)

cat("\nPoor grids (R² <", R2_THRESHOLD, ") — using grid mean fallback:\n")
print(poor_grids)

# Grid means from non-rock loggers only
grid_means <- tomst_spatial %>%
  filter(is.na(high_rock) | !high_rock) %>%
  group_by(site_code, Grid) %>%
  summarise(across(all_of(CELL_LEVEL_VARS), ~mean(., na.rm=TRUE)),
            .groups="drop") %>%
  pivot_longer(cols=all_of(CELL_LEVEL_VARS),
               names_to="variable", values_to="grid_mean")

# Apply fallback: poor grids get grid mean, good grids keep IDW values
final_output <- tomst_interpolated %>%
  mutate(
    site_code = sub("_G.*", "", Cell_ID),
    Grid      = as.integer(sub(".*_G(\\d+)_.*", "\\1", Cell_ID))
  ) %>%
  pivot_longer(cols=all_of(CELL_LEVEL_VARS),
               names_to="variable", values_to="idw_value") %>%
  left_join(poor_grids %>% mutate(is_poor=TRUE),
            by=c("site_code","Grid","variable")) %>%
  left_join(grid_means, by=c("site_code","Grid","variable")) %>%
  mutate(final_value = ifelse(!is.na(is_poor), grid_mean, idw_value)) %>%
  select(-idw_value, -is_poor, -grid_mean) %>%
  pivot_wider(names_from=variable, values_from=final_value) %>%
  select(-site_code, -Grid)

# =============================================================================
# SECTION 9 Save output
# =============================================================================

write.csv(final_output,
          "All_Sites_TOMST_Interpolated_RockMasked.csv",
          row.names=FALSE)

cat("\nSaved: All_Sites_TOMST_Interpolated_RockMasked.csv\n")
cat("Rows:", nrow(final_output), "\n")

# Quick summary of what changed for moisture
cat("\n=== Summary: moisture cells affected by rock masking ===\n")
rock_lookup %>%
  mutate(site = sub("_G.*", "", Cell_ID)) %>%
  filter(high_rock) %>%
  group_by(site) %>%
  summarise(n_masked=n(), .groups="drop") %>%
  print()

cat("\nThese cells have NA for all TOMST variables in the output.\n")
cat("Rock threshold used: rock_cover >=", ROCK_THRESHOLD, "%\n")
cat("Adjust ROCK_THRESHOLD at top of script if needed.\n")
