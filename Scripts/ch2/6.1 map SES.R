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
  separate_wider_delim(cellref, delim = "_", names = c("site", "grid", "cell"), cols_remove = F) |> 
  mutate(grid = paste0(site, grid), 
        column = substr(cell, 1,1), 
        row = substr(cell, 2,3),
        ncolumn = match(column, LETTERS[1:8]))
cell_ses$elevation <- as.factor(cell_ses$elevation)


###Create maps of SES of all traits in all grids####

##Loop over traits and grids
traits <- unique(cell_ses$trait)
grids  <- unique(cell_ses$grid)
sitelist <- c("BK", "GG", "WH")   # site prefixes

map_grid_variation <- function(data, variable, traits, grids, sitelist) {

for (tr in traits) {
  
  message("Processing trait: ", tr)
  
  # filter data for this trait
  trait_dat <- data %>% filter(trait == tr)
  
  # get grids that contain this trait
  trait_grids <- unique(trait_dat$grid)
  
  # Compute global min and max SES
  zmin <- min(trait_dat[, which(colnames(data) == variable)], na.rm = TRUE)
  zmax <- max(trait_dat[, which(colnames(data) == variable)], na.rm = TRUE)
  
  # open PDF
  pdf(paste0("Figures/", variable , "_maps_trait_", tr, ".pdf"), width = 12, height = 10)
  
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
    g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", variable))]
    
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
    
    filler <- data.frame(row = grid_obj$x[missing], ncolumn = grid_obj$y[missing])
    filler[[variable]] <- NA #the variable is set to NA
    
    g_dat <- rbind(g_dat, filler) }
    
    # need values sorted to match the raster cell order:
    g_dat <- g_dat[order(g_dat$ncolumn, g_dat$row), ]
    
    # add SES values to raster
    r <- setValues(r, g_dat[, which(colnames(g_dat) == variable)])
    
    plot(r, main = paste0(variable, " - ", tr, " - ", g), 
         zlim = c(zmin, zmax))
  }}
  
  dev.off()
}
#end function
  }

#map observed RaoQ
map_grid_variation(data = cell_ses, variable = "RaoQ", traits = unique(cell_ses$trait), 
                  grids = unique(cell_ses$grid), sitelist = unique(cell_ses$site))

#map SES
map_grid_variation(data = cell_ses, variable = "SES", traits = unique(cell_ses$trait), 
                   grids = unique(cell_ses$grid), sitelist = unique(cell_ses$site))

#map mean RaoQ of null communities
map_grid_variation(data = cell_ses, variable = "mean_null", traits = unique(cell_ses$trait), 
                   grids = unique(cell_ses$grid), sitelist = unique(cell_ses$site))



####Create maps of one grid####

map_one_grid_variation <- function(data, variable, g) {
  
    
    # filter data for this grid
    g_dat <- data %>% filter(grid == g)
    
    # Compute global min and max SES
    zmin <- min(g_dat[, which(colnames(g_dat) == variable)], na.rm = TRUE)
    zmax <- max(g_dat[, which(colnames(g_dat) == variable)], na.rm = TRUE)
    
    # open PDF
    png(paste0("Figures/", variable ,"_", g, ".png"))
      
      # define plotting layout
      #n <- length(site_grids)
      #nr <- ceiling(sqrt(n))
      #nc <- ceiling(n / nr)
      #par(mfrow = c(nr, nc))
      
        
        #remove extra columns
        g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", variable))]
        
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
          
          filler <- data.frame(row = grid_obj$x[missing], ncolumn = grid_obj$y[missing])
          filler[[variable]] <- NA #the variable is set to NA
          
          g_dat <- rbind(g_dat, filler) }
        
        # need values sorted to match the raster cell order:
        g_dat <- g_dat[order(g_dat$ncolumn, g_dat$row), ]
        g_dat <- as.data.frame(g_dat)
        
        # add SES values to raster
        r <- setValues(r, g_dat[, which(colnames(g_dat) == variable)])
        
        plot(r, main = paste0(variable, " - ", g), 
             zlim = c(zmin, zmax))
    
    dev.off()
  #end function
}


####Create raster of one grid and one variable####

one_grid_raster <- function(data, variable, g) {

  # filter data for this grid
  g_dat <- data %>% filter(grid == g)
  
  
  #remove extra columns
  g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", variable))]
  
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
    
    filler <- data.frame(row = grid_obj$x[missing], ncolumn = grid_obj$y[missing])
    filler[[variable]] <- NA #the variable is set to NA
    
    g_dat <- rbind(g_dat, filler) }
  
  # need values sorted to match the raster cell order:
  g_dat <- g_dat[order(g_dat$ncolumn, g_dat$row), ]
  g_dat <- as.data.frame(g_dat)
  
  return(g_dat)
  #end function
}



map_data <- cell_ses |> 
  filter(trait == "Height_cm")
map_one_grid_variation(data = map_data, variable = "SES", g = "WHG7")


env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(elevation = case_when(grepl("BK", Cell_ID) == T ~ "3000", #add elevation variable
                               grepl("WH", Cell_ID) == T ~ "2500",
                               grepl("GG", Cell_ID) == T ~ "2000",.default = NA)) |> 
  dplyr::select(-c("site", "grid")) |> 
  separate_wider_delim(Cell_ID, delim = "_", names = c("site", "grid", "cell"), cols_remove = F) |> 
  mutate(grid = paste0(site, grid), 
         column = substr(cell, 1,1), 
         row = substr(cell, 2,3),
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2)) |> 
  dplyr::select(!c(soil_temperature_adj_campaign2, soil_moisture_adj_campaign1, soil_temperature_adj_campaign1,
                   aeolian_process, fluvial_process, slope_process,
                   geology1, geology2, geology3,
                   geology4, geology5, mesotopo, aspect, veg_max_height))


WHG7_rock_r <- one_grid_raster(data = env, variable = "rock_cover", g = "WHG7")

grid_rock_WHG7 <- ggplot(WHG7_rock_r, aes(row, ncolumn, fill = rock_cover))+
geom_raster()+
scale_fill_gradient(low = "white", high = "darkorange2", na.value = "pink")+
theme_void()+
labs(x = " ", y = " ", fill = "Rock cover (%)")  +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20))
ggsave("Figures//rock_cover_WHG7.png", grid_rock_WHG7, device = png, width = 2000, units = "px")


WHG7_SES_r <- one_grid_raster(data = cell_ses, variable = "SES", g = "WHG7")

grid_SES_WHG7<- ggplot(WHG7_SES_r, aes(row, ncolumn, fill = SES))+
  geom_raster()+
  scale_fill_gradient2(low = "blue3", mid = "white", high = "darkorange2", na.value = "darkgray")+
  theme_void()+
  labs(x = " ", y = " ", fill = "SES") +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20))
ggsave("Figures//SES_WHG7.png", grid_SES_WHG7, device = png)

