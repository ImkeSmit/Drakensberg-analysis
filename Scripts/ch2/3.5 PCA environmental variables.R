####PCA of all microenvironmental variables####
library(tidyverse)
library(tidylog)

###Microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2), 
         soil_depth_CV = as.numeric(soil_depth_CV)) |> 
  #variables we are interested in
  select(c(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, 
           soil_depth_CV, mean_soil_depth, slope_height)) 

###Microclimate indices for the summer
ind <- read.csv("C:\\Users\\imke6\\Documents\\PhD 2025\\Ch2 niche modelling\\All_data\\Clean_data\\microclimate\\Pekka_cleaned_data\\TOMST data\\transfer_402058_files_56a24364\\microclimate_indices.csv") |> 
  mutate(area = case_when(area == "BOK" ~ "BK", 
                          area == "GOL" ~ "GG", 
                          area == "WIT" ~ "WH"), 
         grid = str_split_i(grid, "_", 2), 
         cell = toupper(str_split_i(site, "_", 3)), 
         Cell_ID = paste0(area, "_G", grid, "_", cell)) |> 
  select(!c(site, impprop)) |> 
  rename(site = "area") |> 
  #pick variables we are interested in
  filter(stat %in% c("avg", "fdh", "gdh5"), sensor %in% c("T1", "moist")) |>  #T1 is the soil temperature, this is the most relaible temperature as it isn't excesively heated by solar radiation
  mutate(variable = paste(sensor, stat, sep = "_")) |> 
  select(!c(sensor, stat)) |> 
  pivot_wider(names_from = variable, values_from = value)

#look at ranges
hist(ind$value_T1_gdh5) 
hist(ind$value_T1_fdh)
hist(ind$value_T1_avg)
hist(ind$value_moist_avg)


####Remote sensing derived variables
rms <- read.csv("All_data/clean_data/Environmental data/Zonal_stats_all.csv") |> 
  rename(Cell_ID = CELL_ID)


###Now we merge them all together
all_env <- env |> 
  full_join(ind, by = "Cell_ID") |> 
  full_join(rms, by = "Cell_ID")
