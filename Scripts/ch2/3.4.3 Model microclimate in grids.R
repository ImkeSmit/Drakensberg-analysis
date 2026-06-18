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
  select(!mesotopo) |> 
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




###Build model for predicting microclimate per site
###GOLDEN GATE####
fullmod_GG <- lm(mean_T1_growing_season ~ rock_cover+ soil_cover+ northness+ eastness+ 
              mesotopo+ veg_max_height+ mean_soil_depth+ slope_height, data = GG_micro_env)
plot(fullmod_GG) ##Looks ok

GG_bestmod <- stepAIC(fullmod_GG, direction = "backward")


###WITSIESHOEK####
fullmod_WH <- glm(mean_T1_growing_season ~ rock_cover+soil_cover+northness+eastness+ 
                 mesotopo+veg_max_height+mean_soil_depth+slope_height, data = WH_micro_env, family = "Gamma")
plot(fullmod_WH)##A little heavy tailed

WH_bestmod <- stepAIC(fullmod_WH, direction = "backward")
plot(WH_bestmod)



##BOKONG####
fullmod_BK <- glm(mean_T1_growing_season ~ rock_cover+soil_cover+northness+eastness+
                 mesotopo+veg_max_height+mean_soil_depth+slope_height, data = BK_micro_env, family = "Gamma")
plot(fullmod_BK)


BK_bestmod <- stepAIC(fullmod_BK, direction = "backward")
plot(WH_bestmod)
