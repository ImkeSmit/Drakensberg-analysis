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
  left_join(env, by = "Cell_ID")

###Golden gate data
GG_micro_env <- micro_env |>  filter(site.x == "GG") 
WH_micro_env <- micro_env |>  filter(site.x == "WH")
BK_micro_env <- micro_env |>  filter(site.x == "BK")




###===MODEL SELECTION===####
###GOLDEN GATE####
fullmod_GG <- lm(mean_T1_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
              mesotopo+ veg_max_height+ mean_soil_depth+ slope_height, data = GG_micro_env)
par(mfrow = c(2,2))
plot(fullmod_GG) ##Looks ok

GG_bestmod <- stepAIC(fullmod_GG, direction = "backward")


###WITSIESHOEK####
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



##BOKONG####
fullmod_BK <- lm(mean_T1_growing_season ~ rock_cover+soil_cover+northness+eastness+
                 mesotopo+veg_max_height+mean_soil_depth+slope_height, data = BK_micro_env)
par(mfrow = c(2,2))
plot(fullmod_BK)


BK_bestmod <- stepAIC(fullmod_BK, direction = "backward")
par(mfrow = c(2,2))
plot(BK_bestmod)


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

predicted_mean_T1_growing_season <- predict(GG_bestmod, GG_pred_dat, type = "response")
GG_pred_dat <- GG_pred_dat |> 
  mutate(predicted_mean_T1_growing_season = predicted_mean_T1_growing_season)

###Validation
#How to do??


  




