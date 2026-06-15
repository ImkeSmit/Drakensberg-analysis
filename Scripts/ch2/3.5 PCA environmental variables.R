####PCA of all microenvironmental variables####
library(tidyverse)
library(tidylog)

env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2), 
         soil_depth_CV = as.numeric(soil_depth_CV)) |> 
  select(c(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, soil_depth_CV, mean_soil_depth, slope_height))
  
