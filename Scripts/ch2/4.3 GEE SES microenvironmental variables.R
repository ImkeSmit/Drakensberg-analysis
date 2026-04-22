###Script to analyse relationship between SES and microenvironmental variables
library(spdep)
library(spatialreg)
library(spind)
library(remotes)
library(dataDownloader)
library(osfr)
library(tidyverse)
library(tidylog)

#imort SES data
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA), 
         site = case_when(grepl("BK", cellref) == T ~ "BK", #add elevation variable
                          grepl("WH", cellref) == T ~ "WH",
                          grepl("GG", cellref) == T ~ "GG",.default = NA),
         grid = str_sub(cellref, 1, 3), 
         column = str_sub(cellref, 4,4), 
         row = as.numeric(str_sub(cellref, 5,6)), 
         ncolumn = match(column, LETTERS[1:8])) |> 
          mutate(Cell_ID = paste0(site, "_G", str_sub(cellref, 3, 3), "_", column, row), 
                 elevation = as.numeric(elevation)) |> 
          select(-cellref)

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth))


##Combine SES and environmental data
comb <- env |> 
  select(-c(site, grid, column,row)) |> 
  inner_join(cell_ses, by = "Cell_ID") |> 
##Now we need to change the coordinates to reflect the spatial structure of the whole dataset
##Coordinates within a site differ by 80
##Coordinates between sites differ by 1000
  mutate(y_new = case_when(site == "GG" ~ row, 
                           site == "WH" ~ row+160, 
                           site == "BK" ~ row+140, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                           grepl("2", grid) ~ y_new+20, 
                           grepl("3", grid) ~ y_new+20*2,
                           grepl("4", grid) ~ y_new+20*3, 
                           grepl("5", grid) ~ y_new+20*4, 
                           grepl("6", grid) ~ y_new+20*5, 
                           grepl("7", grid) ~ y_new+20*6, 
                           grepl("8", grid) ~ y_new+20*7, .default = NA)) |> 
  mutate(x_coord = ncolumn, 
         rock_cover = as.numeric(rock_cover), 
         mean_soil_depth = as.numeric(mean_soil_depth)) |> 
  ungroup() 



###model for SES of height
#isolate trait
Hdat <- comb |> filter(trait == "Height_cm") |> 
                 drop_na(SES, rock_cover, mean_soil_depth) |> #remove rows with na's
                select(Cell_ID,grid, x_coord, y_coord, SES, rock_cover, mean_soil_depth) 

###We have to create the missing coordinates (between the girds and sites)
y_sequence <- unique(Hdat$y_coord)[order(unique(Hdat$y_coord))]
complete_y <- c(y_sequence[1]:y_sequence[length(y_sequence)])
missing_ycoord <- setdiff(complete_y, y_sequence)
missing_ycoord <- rep(missing_ycoord, each = 8) #repeat each ycoord 8 times, because there are 8 possible x coords
missing_xcoord <- rep(c(1:8), length(unique(missing_ycoord)))
#Create table to bind to Hdat
addcoords <- tibble(Cell_ID = NA, x_coord = missing_xcoord, y_coord = missing_ycoord, 
                    SES = NA, rock_cover = NA, mean_soil_depth = NA)
#Hdat <- bind_rows(Hdat, addcoords) |> #add missing coordinates
#  arrange(y_coord) 

Hdat2 <- Hdat |> filter(!is.na(SES)) |> 
  filter(grid %in% c("BK1", "BK2")) |>  #lets first only look at these grids because they are all the same size
  arrange(grid, y_coord) |> 
  mutate(grid = as.factor(grid))
  
#Select 3 complete grids to model, and the spaces between them
#startrow <- which(Hdat$grid == "BK1")[1]
#endrow <- which(Hdat$grid == "BK2")[length(which(Hdat$grid == "BK2"))]

#Hdat2 <- Hdat[c(startrow:endrow), ]  


coords <- cbind(Hdat2$x_coord, Hdat2$y_coord)

# Compute pairwise distances WITHIN each grid
# First, get one representative grid to build R from
# (assuming all grids have the same 8x20 layout)
Hdat |> group_by(grid) |> 
  summarise(length(unique(Cell_ID)))

grid1_dat <- Hdat2 |> filter(grid == unique(grid)[1])
grid_coords <- cbind(grid1_dat$x_coord, grid1_dat$y_coord)

dist_within <- as.matrix(dist(grid_coords))

cutoff <- 5  # your autocorrelation range from variogram

# Option A: Binary cutoff
R <- ifelse(dist_within <= cutoff, 1, 0)
diag(R) <- 1  # gee::gee expects 1s on diagonal
#This identifies which cells belong to the cluster to compute spatial weights within

testlm <- lm(SES ~ rock_cover + mean_soil_depth,
   family = "gaussian", data = Hdat2)

gee2 <- gee::gee(SES ~ rock_cover + mean_soil_depth,
            family = gaussian, data = Hdat2,
            id = grid,
            corstr = "fixed",
            R= R, 
            scale.fix = FALSE)

summary(gee2)


gee1 <- spind::GEE(SES ~ mean_soil_depth, 
            family = gaussian, data = Hdat,
            coord = coords, 
            corstr = "exchangeable", #can be fixed, exchangeable or quadratic
            #fixed means the same autocorrelation structure is applied to the whole dataset
            cluster = 3,
            scale.fix = FALSE)
summary(gee1, printAutoCorPars = TRUE)

plot(gee1)
