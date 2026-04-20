#test for spatial autocorrelation of SES
library(tidyverse)
library(tidylog)
library(gstat)
library(raster)
library(spdep)

cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA), 
                                grid = str_sub(cellref, 1, 3), 
                               column = str_sub(cellref, 4,4), 
                               row = as.numeric(str_sub(cellref, 5,6)), 
                               ncolumn = match(column, LETTERS[1:8])) 
cell_ses$elevation <- as.factor(cell_ses$elevation)


#Moran test for every trait in every grid
traits <- unique(cell_ses$trait)
grids  <- unique(cell_ses$grid)
sitelist <- c("BK", "GG", "WH")   # site prefixes

Moran_result <- data.frame(grid = NA, trait = NA, MoranI = NA, p_value = NA)

l = 1
for(g in grids) {
  for(tr in traits) {
   g_dat <- cell_ses %>% filter(grid == "BK1", trait == "Height_cm")
   #remove extra columns
    g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", "SES"))] 
    
    ##Now get Moran's I for g_dat
    coords <- cbind(g_dat$row, g_dat$ncolumn)
    knn <- knearneigh(coords, k = 8)   # k = number of neighbors
    nb     <- knn2nb(knn)
    lw     <- nb2listw(nb, style = "W")   # row-standardized weights
    
    MI <- moran.test(g_dat$SES, lw)
    
    Moran_result[l, 1] <- g
    Moran_result[l, 2] <- tr
    Moran_result[l, 3] <- MI$estimate[1]
    Moran_result[l, 4] <- MI$p.value
    
    l = l+1
  }
}





##Now get Moran's I for g_dat
coords <- cbind(g_dat$row, g_dat$ncolumn)
knn <- knearneigh(coords, k = 8)   # k = number of neighbors
nb     <- knn2nb(knn)
lw     <- nb2listw(nb, style = "W")   # row-standardized weights

MI <- moran.test(g_dat$SES, lw)
