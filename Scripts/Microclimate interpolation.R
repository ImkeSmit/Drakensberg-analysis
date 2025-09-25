#Script to interpolate microclimate data

library(sp)
library(gstat)
library(raster)
library(tidyverse)
library(tidylog)


#Get microclimate data ready
#We will start with only grid 1 from WH

surface_temp <- read.delim("All_data/clean_data/Tomst_data/Witsieshoek_tomst_data.txt") |> 
  filter(grid == 1) |> 
  group_by(cellref) |> 
  mutate(mean_surface_temp = mean(surface_temp)) |> 
  ungroup() |> 
  distinct(mean_surface_temp, grid, cell, cellref) |> 
  mutate(column = str_match(cell, "^([A-H])([0-9]+)$")[,2], 
         row = str_match(cell, "^([A-H])([0-9]+)$")[,3], 
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(row = as.numeric(row))

#now convert to a spatial object
coordinates(surface_temp) <- ~row +ncolumn #x + y
proj4string(surface_temp) <- CRS("+proj=utm +zone=35 +south +datum=WGS84")

# Define full grid (1x1 m resolution)
x_range <- 1:20
y_range <- 1:8

grid <- expand.grid(x = x_range, y = y_range)
coordinates(grid) <- ~x + y
gridded(grid) <- TRUE
proj4string(grid) <- CRS(proj4string(surface_temp))

# Create IDW model and interpolate
idw_result <- idw(mean_surface_temp ~ 1, locations = surface_temp, newdata = grid, idp = 2.0)  # idp = power (2 is common)

# Convert to raster for plotting
raster_surface_temp <- raster(idw_result)

# Plot
plot(raster_surface_temp, main = "Interpolated Temperature (IDW)")
points(surface_temp, pch = 16, col = "black")


####Let's try kriging
#define variogram model
varmod <- vgm()

kriging_result <- krige(formula = mean_surface_temp ~ 1, 
                        locations = surface_temp, newdata = grid, 
                        model = "Exp") #exponential model fitted to variogram
  
