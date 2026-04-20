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




#make the grids a spatial object
g_dat <- cell_ses %>% filter(grid == "BK1", trait == "Height_cm")

#remove extra columns
g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", "SES"))]


##Now get Moran's I for r
coords <- cbind(g_dat$row, g_dat$ncolumn)
knn <- knearneigh(coords, k = 8)   # k = number of neighbors
nb     <- knn2nb(knn)
lw     <- nb2listw(nb, style = "W")   # row-standardized weights

MI <- moran.test(g_dat$SES, lw)
