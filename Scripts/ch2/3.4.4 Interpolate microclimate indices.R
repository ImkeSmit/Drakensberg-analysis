#Script for interpolating TOMSt indices####



imke_indices <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices.csv", row.names = 1)
env_data   <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv")
centroids  <-read.xlsx("All_data/clean_data/Environmental data/all_cell_coords.xlsx") |> 
  mutate(Cell_ID = paste0(Site, "_G", GRID, "_", COLUMN, ROW)) |> 
         rename(site = Site, 
                grid = GRID)

meter_centroids <- env_data |> #centroids I made myself instead of the geographic coordinates
  select(Cell_ID, site, grid, row, column) |> 
  mutate(y_new = case_when(site == "GG" ~ row, 
                           site == "WH" ~ row+160, 
                           site == "BK" ~ row+160+140, .default = NA)) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                             grepl("2", grid) ~ y_new+20, 
                             grepl("3", grid) ~ y_new+20*2,
                             grepl("4", grid) ~ y_new+20*3, 
                             grepl("5", grid) ~ y_new+20*4, 
                             grepl("6", grid) ~ y_new+20*5, 
                             grepl("7", grid) ~ y_new+20*6, 
                             grepl("8", grid) ~ y_new+20*7, .default = NA)) |> 
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         x_coord = ncolumn) |> 
  select(Cell_ID, x_coord, y_coord, site, grid)


###Identify cells which have too high rock cover####

ROCK_THRESHOLD <- 40

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


####Join loggers to cell centroids
indices_spatial <- imke_indices %>%
  left_join(meter_centroids %>% select(Cell_ID, x_coord, y_coord), by="Cell_ID") %>%
  # Flag loggers sitting on high-rock cells
  left_join(rock_lookup %>% select(Cell_ID, high_rock), by="Cell_ID") |> 
  mutate(grid_unique = str_split_i(Cell_ID, "_", 2), 
         grid_number = substr(grid_unique, 2, 2))

n_rock_loggers <- sum(indices_spatial$high_rock, na.rm=TRUE)
cat("Loggers on high-rock cells (excluded from interpolation):",
    n_rock_loggers, "\n\n")



# =============================================================================
# SECTION 4 IDW interpolation function WITH rock masking
# =============================================================================
IMP_THRESHOLD  <- 0.20
ROCK_THRESHOLD <- 40    # cells with rock_cover >= this are masked out
# consistent with the campaign moisture NA logic
R2_THRESHOLD   <- 0.10  # poor grid fallback
IDP            <- 2      # IDW power parameter: p=2 is standard 
# higher p = more local influence, lower p = smoother surface
NMAX           <- Inf    # use all loggers per grid 


interpolate_grid_idw <- function(site_code, grid_num, tomst_sp,
                                 all_centroids, rock_lkp, var_name, idp=IDP) {
  
  # Source points: logger cells with valid values AND not high-rock
  source_pts <- tomst_sp %>%
    filter(site == !!site_code,
           grid_number      == grid_num,
           (is.na(high_rock) | !high_rock)) %>%  # exclude confirmed rock cells; NA rock_cover treated as non-rock
    select(Cell_ID, x_coord, y_coord, value=all_of(var_name)) %>%
    filter(!is.na(value), !is.na(x_coord), !is.na(y_coord))
  
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
    select(Cell_ID, x_coord, y_coord) %>%
    left_join(rock_lkp %>% select(Cell_ID, high_rock), by="Cell_ID") %>%
    filter(is.na(high_rock) | !high_rock) %>%  # exclude confirmed rock cells
    select(Cell_ID, x_coord, y_coord)
  
  if (nrow(target_pts) == 0) return(NULL)
  
  # Project to UTM 36S for accurate distance calculation
  #source_sf <- st_as_sf(source_pts, coords=c("CellY","CellX"), crs=4326) %>%
   # st_transform(32736)
  #target_sf <- st_as_sf(target_pts, coords=c("CellY","CellX"), crs=4326) %>%
  #  st_transform(32736)
  
  #source_coords <- st_coordinates(source_sf)
  #target_coords <- st_coordinates(target_sf)
  
  # st_transform preserves row order, so source_coords rows align with source_pts rows
  #source_df <- data.frame(X=source_coords[,1], Y=source_coords[,2],
                          #value=source_pts$value)
  #target_df <- data.frame(X=target_coords[,1], Y=target_coords[,2])
  
  source_df <- source_pts |> 
    rename(X = x_coord, Y = y_coord) |> 
    select(!Cell_ID)
  
  target_df = target_pts |> 
    rename(X = x_coord, Y = y_coord) |> 
    select(!Cell_ID)
  
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

site_grids <- indices_spatial %>%
  distinct(site, grid_number) %>%
  arrange(site, grid_number)

CELL_LEVEL_VARS = as.vector(colnames(indices_spatial)[c(3)])

all_interpolated <- map_dfr(CELL_LEVEL_VARS, function(var) {
  cat("Interpolating:", var, "\n")
  map2_dfr(site_grids$site, site_grids$grid_number, function(sc, gn) {
    interpolate_grid_idw(sc, gn, indices_spatial, meter_centroids,
                         rock_lookup, var)
  })
})

#tomst_interpolated <- all_interpolated %>%
#  group_by(Cell_ID) %>%
#  summarise(across(all_of(CELL_LEVEL_VARS),
#                   ~first(.[!is.na(.)])),
#            .groups="drop")

cat("\nInterpolated surface:", nrow(all_interpolated), "cells\n")

# High-rock cells will be absent from tomst_interpolated — add them back as NA
# so the output has all cells and rock cells simply have NA for TOMST variables
all_cell_ids <- meter_centroids %>%
  distinct(Cell_ID)

tomst_interpolated <- all_cell_ids %>%
  left_join(all_interpolated, by="Cell_ID")

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
    filter(site == !!site_code,
           grid_number      == grid_num,
           (is.na(high_rock) | !high_rock)) %>%
    select(Cell_ID, x_coord, y_coord, value=all_of(var_name)) %>%
    filter(!is.na(value), !is.na(x_coord), !is.na(y_coord))
  
  if (nrow(grid_loggers) < 3) return(NULL)
  
  map_dfr(seq_len(nrow(grid_loggers)), function(i) {
    held_out <- grid_loggers[i, ]
    training <- grid_loggers[-i, ]
    
    #train_sf <- st_as_sf(training, coords=c("lon","lat"), crs=4326) %>%
    #  st_transform(32736)
    #held_sf  <- st_as_sf(held_out, coords=c("lon","lat"), crs=4326) %>%
    #  st_transform(32736)
    
    # st_transform preserves row order — coordinates align with original data rows
    #train_df <- data.frame(X=st_coordinates(train_sf)[,1],
    #                       Y=st_coordinates(train_sf)[,2],
    #                       value=training$value)
    #held_df  <- data.frame(X=st_coordinates(held_sf)[,1],
    #                       Y=st_coordinates(held_sf)[,2])
    
    train_df <- training |> 
      rename(X = x_coord, Y = y_coord)
    held_df <- held_out |> 
      rename(X = x_coord, Y = y_coord)
    
    obs_min <- min(training$value)
    obs_max <- max(training$value)
    
    idw_mod <- gstat(formula=value~1, locations=~X+Y,
                     data=train_df, nmax=nrow(train_df), set=list(idp=IDP))
    pred    <- predict(idw_mod, newdata=held_df)
    pred_val <- pmax(pmin(pred$var1.pred, obs_max), obs_min)
    
    tibble(site=site_code, grid_number=grid_num, variable=var_name,
           Cell_ID=held_out$Cell_ID, observed=held_out$value,
           predicted=pred_val, error=pred_val - held_out$value)
  })
}

validation_results <- map_dfr(CELL_LEVEL_VARS, function(var) {
  cat("Validating:", var, "\n")
  map2_dfr(site_grids$site, site_grids$grid, function(sc, gn) {
    run_loloo_validation(sc, gn, indices_spatial, var)
  })
})

# =============================================================================
# SECTION 7 Validation summary
# =============================================================================

validation_grid <- validation_results %>%
  group_by(site, grid_number, variable) %>%
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