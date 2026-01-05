###Cluster cells according to environmental similarity###
library(tidyverse)
library(tidylog)

#import environmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth))

sdepth <- env |> 
  dplyr::select(Cell_ID, mean_soil_depth) |> 
  filter(!is.na(mean_soil_depth)) |> 
  column_to_rownames("Cell_ID")

sd_dist <- dist(sdepth, method = "euclidean")

sd_hclust <- hclust(sd_dist, method = "ward.D")

#cut the tree to have 5 clusters
cut <- cutree(sd_hclust, k = 5)
plot(sd_hclust)
rect.hclust(sd_hclust, k = 5, col = "red")