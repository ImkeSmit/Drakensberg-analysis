##Map ses values in each grid
library(tidyverse)
library(sp)
library(gstat)
library(raster)

#import SES at cell scale
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  mutate(site = substr(cellref, 1,2),
        grid = substr(cellref, 1,3),
         column = substr(cellref, 4,4), 
         row = substr(cellref, 5,6),
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(row = as.numeric(row)) 
cell_ses$elevation <- as.factor(cell_ses$elevation)


###Create maps of SES of all traits in all grids####

##Loop over traits and grids
traits <- unique(cell_ses$trait)
grids  <- unique(cell_ses$grid)
sitelist <- c("BK", "GG", "WH")   # site prefixes

for (tr in traits) {
  
  message("Processing trait: ", tr)
  
  # filter data for this trait
  trait_dat <- cell_ses %>% filter(trait == tr)
  
  # get grids that contain this trait
  trait_grids <- unique(trait_dat$grid)
  
  # Compute global min and max SES
  zmin <- min(trait_dat$SES, na.rm = TRUE)
  zmax <- max(trait_dat$SES, na.rm = TRUE)
  
  # open PDF
  pdf(paste0("Figures/SES_maps_trait_", tr, ".pdf"), width = 12, height = 10)
  
  # ---- LOOP OVER SITES: BK → GG → WH ----
  for (s in sitelist) {
    
    site_dat <- trait_dat %>% filter(site == s)
    site_grids <- unique(site_dat$grid)
  
  # define plotting layout
  n <- length(site_grids)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  par(mfrow = c(nr, nc))
  
  # loop over grids
  for (g in site_grids) {
    
    g_dat <- trait_dat %>% filter(grid == g)
  
    #remove extra columns
    g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", "SES"))]
    
    # Create spatial grid object
    x_range <- 1:20
    y_range <- 1:8
    
    grid_obj <- expand.grid(x = x_range, y = y_range)
    coordinates(grid_obj) <- ~x + y
    gridded(grid_obj) <- TRUE
    r <- raster(grid_obj) 
    
    #Check if there are missing cells in data
    possible_cells <- paste(grid_obj$x, grid_obj$y, sep = "_")
    observed_cells <- paste(g_dat$row, g_dat$ncolumn, sep = "_")
    missing <- which(is.na(match(possible_cells, observed_cells)))
  
    if(length(missing) > 0) { #create a filler
    
    filler <- data.frame(row = grid_obj$x[missing], ncolumn = grid_obj$y[missing], SES = NA)
    
    g_dat <- rbind(g_dat, filler) }
    
    # need values sorted to match the raster cell order:
    g_dat <- g_dat[order(g_dat$ncolumn, g_dat$row), ]
    
    # add SES values to raster
    r <- setValues(r, g_dat$SES)
    
    plot(r, main = paste0("SES - ", tr, " - ", g), 
         zlim = c(zmin, zmax))
  }}
  
  dev.off()
}


