###Script to analyse relationship between SES and microenvironmental variables
library(spdep)
library(spatialreg)
library(spind)
library(remotes)
library(dataDownloader)
library(osfr)
library(tidyverse)
library(tidylog)

#import SES data
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
          select(-c(cellref, site, grid, column, row)) |> 
  filter(trait == "Height_cm")

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2)) |> 
  select(!c(soil_temperature_adj_campaign2, soil_moisture_adj_campaign1, soil_temperature_adj_campaign1,
            aeolian_process, fluvial_process, slope_process,
            geology1, geology2, geology3,
            geology4, geology5, mesotopo, aspect, veg_max_height))


##Combine SES and environmental data
comb <- env |> 
  #select(-c(site, grid, column,row)) |> 
  full_join(cell_ses, by = "Cell_ID") |> 
##Now we need to change the coordinates to reflect the spatial structure of the whole dataset
#make the grids contiguous, differing by 20m along the y axis
  mutate(y_new = case_when(site == "GG" ~ row, 
                           site == "WH" ~ row+160, 
                           site == "BK" ~ row+160+140, .default = NA)) |> 
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


###Fill in the x and y coordinates of cells that do not have SES
comb2 <- comb |> 
  mutate(ncolumn = match(column, LETTERS[1:8]),
         elevation = case_when(grepl("BK", Cell_ID) == T ~ "3000", #add elevation variable
                               grepl("WH", Cell_ID) == T ~ "2500",
                               grepl("GG", Cell_ID) == T ~ "2000", .default = NA),
         y_new = case_when(site == "GG" ~ row, 
                           site == "WH" ~ row+160, 
                           site == "BK" ~ row+160 +140, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                             grepl("2", grid) ~ y_new+20, 
                             grepl("3", grid) ~ y_new+20*2,
                             grepl("4", grid) ~ y_new+20*3, 
                             grepl("5", grid) ~ y_new+20*4, 
                             grepl("6", grid) ~ y_new+20*5, 
                             grepl("7", grid) ~ y_new+20*6, 
                             grepl("8", grid) ~ y_new+20*7, .default = NA)) |> 
  mutate(x_coord = ncolumn) |> 
  ungroup()


#All grids must have 160 cells, check this
comb2 |> group_by(trait, site, grid) |> 
  summarise(n = n()) 


###model for SES of height
#isolate trait
Hdat <- comb |> filter(trait == "Height_cm") |> 
                 drop_na(SES, rock_cover, mean_soil_depth) |> #remove rows with na's
                select(Cell_ID,grid, x_coord, y_coord, SES, rock_cover, mean_soil_depth, elevation) 

##We may have to impute values for missing cells to be able to run the gee model

#isolate two grids to model for now
Hdat2 <- Hdat |> #filter(!is.na(SES)) |> 
  #filter(grid %in% c("BK1", "BK2")) |>  #lets first only look at these grids because they are all the same size
  arrange(grid, y_coord) |> 
  mutate(grid = as.factor(grid))

#Coordinates to compute spatial relationships from
coords <- cbind(Hdat2$x_coord, Hdat2$y_coord)

# Compute pairwise distances WITHIN each grid
# First, get one representative grid to build R from
# (assuming all grids have the same 8x20 layout)
grid1_dat <- Hdat2 |> filter(grid == unique(grid)[1])
grid_coords <- cbind(grid1_dat$x_coord, grid1_dat$y_coord)
dist_within <- as.matrix(dist(grid_coords))


###Draw a variogram
# Load libraries
library(gstat)
library(sp)

coordinates(Hdat2) = ~x_coord + y_coord

# 2. Compute the empirical variogram
# The formula 'log(zinc) ~ 1' assumes a constant mean (no spatial trend)
v_empirical <- variogram(SES ~ 1, data = Hdat2)

# 3. Plot the empirical variogram
plot(v_empirical)

# 4. Define an initial model (e.g., Spherical with sill=1, range=900, nugget=0.1)
initial_model <- vgm(psill = 1, model = "Sph", range = 10, nugget = 0.1)

# 5. Fit the model to empirical data
v_fitted <- fit.variogram(v_empirical, initial_model)

# 6. Plot empirical points and the fitted model together
plot(v_empirical, model = v_fitted)




#Draw a correlogram to see at which distance the spatial autocorrelation in model residuals decreases
testlm <- lm(SES ~ rock_cover + mean_soil_depth, data = Hdat2)
resid <- resid(testlm)

neibs <- knn2nb(knearneigh(coords, k = 4))
spcor <- sp.correlogram(neibs, resid, method = "I", order = 20)
plot(spcor)

cutoff <- 9  # your autocorrelation range from correlogram

# Option A: Binary cutoff
R <- ifelse(dist_within <= cutoff, 1, 0)
diag(R) <- 1  # gee::gee expects 1s on diagonal
#This identifies which cells belong to the cluster to compute spatial weights within


##Option 2 build in distance decay based on correlogram
R2 <- exp(-dist_within / cutoff)
diag(R2) <- 1


gee2 <- gee::gee(SES ~ rock_cover + mean_soil_depth,
            family = gaussian, data = Hdat2,
            id = grid,
            corstr = "fixed",
            R= R2, 
            scale.fix = FALSE)

summary(gee2)

#Get p values
coefs <- summary(gee2)$coefficients
p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
cbind(coefs, p_value = round(p_values, 4))


