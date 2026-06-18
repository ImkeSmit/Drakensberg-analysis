####MODEL MICROCLIMATE INDICES OFR THE WHOLE GRID####
###USING ALL ENVIRONMENTAL VARIABLES####
library(tidyverse)
library(tidylog)
library(corrplot)

###Import microclimate indices
ind <- read.csv("all_data/clean_data/Environmental data/Imke_microclimate_indices.csv", row.names = 1) |> 
  mutate(grid = str_split_i(Cell_ID, "_", 2), 
         column = str_sub(Cell_ID, 7,7), 
         row = as.numeric(str_sub(Cell_ID, 8,9)), 
         ncolumn = match(column, LETTERS[1:8]))


###Import environmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2), 
         soil_depth_CV = as.numeric(soil_depth_CV), 
         mesotopo = as.factor(mesotopo)) |> 
  #variables we are interested in
  dplyr::select(c(Cell_ID:row, vascular_cover, rock_cover, soil_cover, northness, eastness, 
                  mesotopo, veg_max_height, mean_soil_depth, slope_height)) 


###Correlation of env variables
cordf <- env |> 
  select(!mesotopo) |> 
  drop_na()
cormat <- cor(cordf[, c(6:13)])
corrplot(cormat, method = "number", type = "lower", 
         tl.cex = 0.8,    # shrink variable name text
         tl.srt = 45)
##Only vascular cover and rock cover are highly correlated

