##Map ses values in each grid
library(tidyverse)
library(sp)
library(gstat)
library(raster)

#import SES at cell scale, computed from unscaled RaoQ
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  mutate(grid = substr(cellref, 1,3),
         column = substr(cellref, 4,4), 
         row = substr(cellref, 5,5),
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(row = as.numeric(row)) |> 
  filter(grid == "BK1", trait == "Height_cm")
cell_ses$elevation <- as.factor(cell_ses$elevation)

# Define full grid (1x1 m resolution)
x_range <- 1:20
y_range <- 1:8

grid <- expand.grid(x = x_range, y = y_range)
coordinates(grid) <- ~x + y
gridded(grid) <- TRUE
proj4string(grid) <- CRS(proj4string(surface_temp))

ses_raster <- raster(grid)
ses_raster <- setValues(ses_raster, cell_ses$SES)
plot(ses_raster, main = "SES of height, BK1")
