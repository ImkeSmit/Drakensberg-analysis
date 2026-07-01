#Script for interpolating TOMSt indices####



imke_indices <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices.csv", row.names = 1)
env_data   <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv")
centroids  <-read.xlsx("All_data/clean_data/Environmental data/all_cell_coords.xlsx") |> 
  mutate(Cell_ID = paste0(Site, "_G", GRID, "_", COLUMN, ROW)) 


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
  left_join(centroids %>% select(Cell_ID, CellX, CellY), by="Cell_ID") %>%
  # Flag loggers sitting on high-rock cells
  left_join(rock_lookup %>% select(Cell_ID, high_rock), by="Cell_ID") 

n_rock_loggers <- sum(indices_spatial$high_rock, na.rm=TRUE)
cat("Loggers on high-rock cells (excluded from interpolation):",
    n_rock_loggers, "\n\n")
