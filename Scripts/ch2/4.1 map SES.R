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
         row = substr(cellref, 5,6),
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(row = as.numeric(row)) 
cell_ses$elevation <- as.factor(cell_ses$elevation)


###Create maps of SES of all traits in all grids####

#function to make the ratser of each grid
make_ses_raster <- function(dat) {
  
  # full 20 × 8 grid
  x_range <- 1:20
  y_range <- 1:8
  
  grid <- expand.grid(x = x_range, y = y_range)
  coordinates(grid) <- ~x + y
  gridded(grid) <- TRUE
  
  r <- raster(grid)
  r <- setValues(r, dat$SES)
  
  plot(r, main = paste0("SES - ", dat$trait[1], " - ", dat$grid[1]))
}

##Loop over traits and grids
traits <- unique(cell_ses$trait)
grids  <- unique(cell_ses$grid)

for (tr in traits) {
  
  message("Processing trait: ", tr)
  
  # filter data for this trait
  trait_dat <- cell_ses %>% filter(trait == tr)
  
  # get grids that contain this trait
  trait_grids <- unique(trait_dat$grid)
  
  # open PDF
  pdf(paste0("Figures/SES_maps_trait_", tr, ".pdf"), width = 12, height = 10)
  
  # define plotting layout
  n <- length(trait_grids)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  par(mfrow = c(nr, nc))
  
  # loop over grids
  for (g in trait_grids) {
    
    g_dat <- trait_dat %>% filter(grid == g)
    
    # need values sorted to match the raster cell order:
    g_dat <- g_dat[order(g_dat$ncolumn, g_dat$row), ]
    
    # Create spatial grod object
    x_range <- 1:20
    y_range <- 1:8
    
    grid_obj <- expand.grid(x = x_range, y = y_range)
    coordinates(grid_obj) <- ~x + y
    gridded(grid_obj) <- TRUE
    r <- raster(grid_obj) 
    
    #Insert here something to add NA values to any cells that may be missing
    
    # add SES values to raster
    r <- setValues(r, g_dat$SES)
    
    plot(r, main = paste0("SES - ", tr, " - ", g))
  }
  
  dev.off()
}


