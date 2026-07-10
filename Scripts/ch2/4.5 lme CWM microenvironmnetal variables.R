###Models of CWm ~ microenvironmental variables with lme####
###And predictions of CWm traits over environmental gradients####
library(tidyverse)
library(tidylog)
library(nlme)
library(MuMIn)
library(DHARMa)
library(tictoc)
library(ggridges)
library(emmeans)
library(multcomp)
library(conflicted)
library(performance)
library(see)
conflict_prefer_all("tidylog", quiet = TRUE)

#import trait and abundance data
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1) |> 
  mutate(log_Height = log(Height_cm), 
         log_LA = log(Leaf_area_mm2), 
         log_LDMC = log(LDMC))


#compute CWM of each trait for each cell
library(FD)
cwm <- functcomp(x = mean_traits, a = as.matrix(abun_matrix))
cwm2 <- cwm |>
  rownames_to_column(var = "cellref") |> 
  rename(Cell_ID = cellref) |> 
  pivot_longer(cols = !Cell_ID, names_to = "trait", values_to = "cwm_value")

#import microenvironmental data
env2 <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, 
         soil_depth_CV, mean_soil_depth, slope_height) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA))

#import remote sensing derived variables
rms <- read.csv("All_data/clean_data/Environmental data/Zonal_stats_all.csv") |> 
  select(CELL_ID, STD) |> 
  rename(Cell_ID = CELL_ID)

#import interpolated microclimate indices
micro_idw <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices_idw_interpolated.csv", row.names = 1)


##Combine CWM and environmental data
cwm_comb <- env2 |> 
  #join to microclimate indices |> 
  full_join(micro_idw, by = "Cell_ID") |> 
  #join to remote sensing data |> 
  full_join(rms, by = "Cell_ID") |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cwm2, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)


####model CWM ~ microenvironmental variables####
traitlist <- c("Height_cm","Leaf_area_mm2", "log_LDMC","SLA") #the traits we chose to log here are different to the traits logged in the elevation model.
sitelist <- c("GG", "WH", "BK")
variables <- c("mean_T1_growing_season" , "mean_moist_growing_season" , "rock_cover", "mean_soil_depth" , "STD")



for (t in 1:length(traitlist)) {
  for(s in 1:length(sitelist)) {
  modeldat <-  cwm_comb |> 
    filter(trait == traitlist[t], 
           site == sitelist[s]) |> 
    mutate(elevation = as.factor(elevation), 
           grid = as.factor(grid)) |> 
    drop_na()
  
  
  model<- lme(cwm_value ~ mean_T1_growing_season + mean_moist_growing_season + rock_cover+ mean_soil_depth + STD ,
              random = ~1|grid, 
              correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
              data = modeldat) #only gaussian family possible
  
  for(v in 1:length(variables)) {
    pred_var <- variables[v]
    mean_var <- setdiff(variables, variables[v])
    
    pred_dat <- modeldat |> 
      mutate(accross(mean_var), ~ mean(.x, na.rm = TRUE))
    
    
    
    
  }
  
  
}}