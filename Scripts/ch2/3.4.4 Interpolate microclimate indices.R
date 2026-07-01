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
