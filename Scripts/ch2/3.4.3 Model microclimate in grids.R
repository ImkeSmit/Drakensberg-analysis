####MODEL MICROCLIMATE INDICES OFR THE WHOLE GRID####
###USING ALL ENVIRONMENTAL VARIABLES####
library(tidyverse)
library(tidylog)
library(corrplot)
library(MASS)

###Import microclimate indices
ind <- read.csv("all_data/clean_data/Environmental data/Imke_microclimate_indices.csv", row.names = 1)# |> 
  #mutate(grid = str_split_i(Cell_ID, "_", 2), 
  #       column = str_sub(Cell_ID, 7,7), 
  #       row = as.numeric(str_sub(Cell_ID, 8,9)), 
  #       ncolumn = match(column, LETTERS[1:8]))


###Import environmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2), 
         soil_depth_CV = as.numeric(soil_depth_CV), 
         mesotopo = as.factor(mesotopo)) |> 
  #variables we are interested in
  dplyr::select(c(Cell_ID:row, vascular_cover, rock_cover, soil_cover, northness, eastness, 
                  mesotopo, veg_max_height, mean_soil_depth, slope_height, MAX)) 


###Correlation of env variables
cordf <- env |> 
  dplyr::select(!mesotopo) |> 
  drop_na()
cormat <- cor(cordf[, c(6:14)])
corrplot(cormat, method = "number", type = "lower", 
         tl.cex = 0.8,    # shrink variable name text
         tl.srt = 45)
##Only vascular cover and rock cover are highly correlated


###Join microclimate and env variables
micro_env <- ind |> 
  left_join(env, by = "Cell_ID") |>
  rename(site = site.x) |> 
  dplyr::select(!site.y) |> 
  mutate(grid = paste0(site, grid))

###Golden gate data
GG_micro_env <- micro_env |>  filter(site.x == "GG") 
WH_micro_env <- micro_env |>  filter(site.x == "WH")
BK_micro_env <- micro_env |>  filter(site.x == "BK")


###PREDICT MICROCLIMATE INDICES FROM ENVIRONMENTAL VARIABLES###
#Do this per site. 
#Too few datapoints to do it per grid

#list of sites
gridlist <- c(unique(micro_env$site))

#start loop
for(g in 1:length(gridlist)) {
  one_grid <- micro_env |>  filter(site == gridlist[g])
  
  #full models
  grid_fullmod_T1 <- lm(mean_T1_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
                      veg_max_height+ mean_soil_depth+ slope_height +MAX, data = one_grid)
  #does not contain mesotopo because the model tries to predict to novel mesotopo codes
  
  grid_fullmod_moist <- lm(mean_moist_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
                         veg_max_height+ mean_soil_depth+ slope_height + MAX, data = one_grid)
  
  #backward selection
  T1_bestmod <- stepAIC(grid_fullmod_T1, direction = "backward")
  
  # Store all models visited during stepping, because moistmod sometimes selcts the null model
  visited_models <- list()
  
  moist_bestmod <- stepAIC(
    grid_fullmod_moist,
    direction = "backward",
    trace = TRUE,
    keep = function(model, aic) {
      visited_models[[length(visited_models) + 1]] <<- model  # <<- assigns to parent env
      list(formula = formula(model), aic = aic)
    }
  )
  
  # Compare formula as string to check if best model is the null model
  if (deparse(formula(moist_bestmod)) == deparse(as.formula(mean_moist_growing_season ~ 1))) {
    
    # Extract AICs from visited models
    all_aics     <- sapply(visited_models, AIC)
    
    # Sort and get second lowest
    sorted_idx       <- order(all_aics)
    second_best_form <- formula(visited_models[[sorted_idx[2]]])
    
    # Refit
    moist_bestmod <- update(moist_bestmod, formula = second_best_form)
  }
  
  
  #predictions
  #data for prediction
  pred_dat <- env |> 
    filter(site == gridlist[g]) |> 
    left_join(ind, by = "Cell_ID") |>  #add measured microclimate indices
    dplyr::select(!c(site.y, start_growing_season, end_growing_season)) |> 
    rename(site = site.x)
  
  #Mean T1
  predicted_mean_T1_growing_season <- predict(T1_bestmod, pred_dat, type = "response")
  
  #mean moist
  predicted_mean_moist_growing_season <- predict(moist_bestmod, pred_dat, type = "response")
  
  pred_dat_save <- pred_dat |> 
    mutate(predicted_mean_T1_growing_season = predicted_mean_T1_growing_season, 
           predicted_mean_moist_growing_season = predicted_mean_moist_growing_season,
           T1_bestmod_formula = as.character(summary(T1_bestmod)$call)[2], 
           moist_bestmod_formula = as.character(summary(moist_bestmod)$call)[2], 
           T1_bestmod_rsq = summary(T1_bestmod)$adj.r.squared, 
           moist_bestmod_rsq = summary(moist_bestmod)$adj.r.squared) 
  
  if(g == 1) {
    result <- pred_dat_save
  } else {
    temp_result <- pred_dat_save
    result <- bind_rows(result, temp_result)
  }
}










  




