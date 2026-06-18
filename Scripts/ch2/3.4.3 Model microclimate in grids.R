####MODEL MICROCLIMATE INDICES OFR THE WHOLE GRID####
###USING ALL ENVIRONMENTAL VARIABLES####
library(tidyverse)
library(tidylog)

###Import microclimate indices
ind <- read.csv("all_data/clean_data/Environmental data/Imke_microclimate_indices.csv", row.names = 1) |> 
  mutate(grid = str_split_i(Cell_ID, "_", 2), 
         column = str_sub(Cell_ID, 7,7), 
         row = as.numeric(str_sub(Cell_ID, 8,9)), 
         ncolumn = match(column, LETTERS[1:8]))
