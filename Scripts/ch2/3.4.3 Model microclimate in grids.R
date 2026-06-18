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
                  mesotopo, veg_max_height, mean_soil_depth, slope_height)) 


###Correlation of env variables
cordf <- env |> 
  dplyr::select(!mesotopo) |> 
  drop_na()
cormat <- cor(cordf[, c(6:13)])
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
                      veg_max_height+ mean_soil_depth+ slope_height, data = one_grid)
  #does not contain mesotopo because the model tries to predict to novel mesotopo codes
  
  grid_fullmod_moist <- lm(mean_moist_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
                         veg_max_height+ mean_soil_depth+ slope_height, data = one_grid)
  
  #backward selection
  T1_bestmod <- stepAIC(grid_fullmod_T1, direction = "backward")
  
  moist_bestmod <- stepAIC(grid_fullmod_moist, direction = "backward")
  
  
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
           T1_bestmod_formula = paste0(as.character(formula(T1_bestmod))), 
           moist_bestmod_formula = as.character(formula (moist_bestmod))) 
  
  if(g == 1) {
    result <- pred_dat_save
  } else {
    temp_result <- pred_dat_save
    result <- bind_rows(result, temp_result)
  }
    
}







###===MODEL SELECTION####
###GOLDEN GATE####
##mean T1
fullmod_GG <- lm(mean_T1_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
              mesotopo+ veg_max_height+ mean_soil_depth+ slope_height, data = GG_micro_env)
par(mfrow = c(2,2))
plot(fullmod_GG) ##Looks ok

GG_bestmod <- stepAIC(fullmod_GG, direction = "backward")

##mean moist
fullmod_GG_moist <- lm(mean_moist_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
                   mesotopo+ veg_max_height+ mean_soil_depth+ slope_height, data = GG_micro_env)
par(mfrow = c(2,2))
plot(fullmod_GG_moist) ##Looks ok

GG_bestmod_moist <- stepAIC(fullmod_GG_moist, direction = "backward")


###WITSIESHOEK####
#mean T1
fullmod_WH <- lm(mean_T1_growing_season ~ rock_cover+soil_cover+northness+eastness+ 
                 mesotopo+veg_max_height+mean_soil_depth+slope_height, data = WH_micro_env)

par(mfrow = c(2,2))
plot(fullmod_WH)##A little heavy tailed

bc <- boxcox(fullmod_WH, plotit = T)
# Extract the best lambda
best_lambda <- bc$x[which.max(bc$y)]
#apply boxcox transformation
WH_micro_env$y_transformed <- (WH_micro_env$mean_T1_growing_season^best_lambda - 1) / best_lambda

fullmod_WH2 <- lm(y_transformed ~ rock_cover+soil_cover+northness+eastness+ 
                   mesotopo+veg_max_height+mean_soil_depth+slope_height, data = WH_micro_env)

par(mfrow = c(2,2))
plot(fullmod_WH2)##Not really better

#backward selection
WH_bestmod <- stepAIC(fullmod_WH, direction = "backward")
par(mfrow = c(2,2))
plot(WH_bestmod)



##mean moist
fullmod_WH_moist <- lm(mean_moist_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
                         mesotopo+ veg_max_height+ mean_soil_depth+ slope_height, data = WH_micro_env)
par(mfrow = c(2,2))
plot(fullmod_WH_moist) ##Looks ok

WH_bestmod_moist <- stepAIC(fullmod_WH_moist, direction = "backward")
#null model is best??


##BOKONG####
#mean T1
fullmod_BK <- lm(mean_T1_growing_season ~ rock_cover+soil_cover+northness+eastness+
                 mesotopo+veg_max_height+mean_soil_depth+slope_height, data = BK_micro_env)
par(mfrow = c(2,2))
plot(fullmod_BK)


BK_bestmod <- stepAIC(fullmod_BK, direction = "backward")
par(mfrow = c(2,2))
plot(BK_bestmod)


#mean moist
fullmod_BK_moist <- lm(mean_moist_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
                         mesotopo+ veg_max_height+ mean_soil_depth+ slope_height, data = BK_micro_env)
par(mfrow = c(2,2))
plot(fullmod_BK_moist) 

BK_bestmod_moist <- stepAIC(fullmod_BK_moist, direction = "backward")


####===PREDICTION===####
###GOLDEN GATE####

formula(GG_bestmod)

#get data to predict over
GG_pred_dat <- env |> 
  filter(site == "GG") |> 
  left_join(ind, by = "Cell_ID") |>  #add measured microclimate indices
  dplyr::select(!c(site.y, start_growing_season, end_growing_season)) |> 
  rename(site = site.x)

#check that grids are complete
GG_pred_dat |>  group_by(grid) |> 
summarise(ncells = n())

#Mean T1
predicted_mean_T1_growing_season <- predict(GG_bestmod, GG_pred_dat, type = "response")

#mean moist
predicted_mean_moist_growing_season <- predict(GG_bestmod_moist, GG_pred_dat, type = "response")

GG_pred_dat <- GG_pred_dat |> 
  mutate(predicted_mean_T1_growing_season = predicted_mean_T1_growing_season, 
         predicted_mean_moist_growing_season = predicted_mean_moist_growing_season)

###Validation
#How to do??



###WITSIESHOEK####
#get data to predict over
WH_pred_dat <- env |> 
  filter(site == "WH") |> 
  left_join(ind, by = "Cell_ID") |>  #add measured microclimate indices
  dplyr::select(!c(site.y, start_growing_season, end_growing_season)) |> 
  rename(site = site.x)

#check that grids are complete
WH_pred_dat |>  group_by(grid) |> 
  summarise(ncells = n())

#Mean T1
WH_predicted_mean_T1_growing_season <- predict(WH_bestmod, WH_pred_dat, type = "response")

#mean moist
WH_predicted_mean_moist_growing_season <- predict(lm(mean_moist_growing_season ~ slope_height, data = WH_micro_env),
                                                  #put next best model in here because the best model is ~1
                                                  WH_pred_dat, type = "response")

WH_pred_dat <- WH_pred_dat |> 
  mutate(predicted_mean_T1_growing_season = WH_predicted_mean_T1_growing_season, 
         predicted_mean_moist_growing_season = WH_predicted_mean_moist_growing_season)

###Validation
#How to do??



###BOKONG####
#get data to predict over
BK_pred_dat <- env |> 
  filter(site == "BK") |> 
  left_join(ind, by = "Cell_ID") |>  #add measured microclimate indices
  dplyr::select(!c(site.y, start_growing_season, end_growing_season)) |> 
  rename(site = site.x)

#check that grids are complete
BK_pred_dat |>  group_by(grid) |> 
  summarise(ncells = n())

#Mean T1
BK_predicted_mean_T1_growing_season <- predict(BK_bestmod, BK_pred_dat, type = "response")

#mean moist
BK_predicted_mean_moist_growing_season <- predict(BK_bestmod_moist,
                                                  BK_pred_dat, type = "response")

WH_pred_dat <- WH_pred_dat |> 
  mutate(predicted_mean_T1_growing_season = WH_predicted_mean_T1_growing_season, 
         predicted_mean_moist_growing_season = WH_predicted_mean_moist_growing_season)

###Validation
#How to do??


  




