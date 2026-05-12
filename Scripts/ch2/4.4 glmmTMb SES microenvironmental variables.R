####MOdelling with glmmTMB####
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(tictoc)
#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  rename(Cell_ID = cellref)

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2)) |> 
  dplyr::select(!c(soil_temperature_adj_campaign2, soil_moisture_adj_campaign1, soil_temperature_adj_campaign1,
                   aeolian_process, fluvial_process, slope_process,
                   geology1, geology2, geology3,
                   geology4, geology5, mesotopo, aspect, veg_max_height))


##Combine SES and environmental data
comb <- env |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cell_ses, by = "Cell_ID", relationship = "one-to-many") |>
  ##Now we need to change the coordinates to reflect the spatial structure of the whole dataset
  #make the grids contiguous, differing by 20m along the y axis
  mutate(elevation = as.factor(case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000, .default = NA))) |> 
  mutate(y_new = case_when(site == "GG" ~ row, 
                           site == "WH" ~ row+160, 
                           site == "BK" ~ row+160+140, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                             grepl("2", grid) ~ y_new+20, 
                             grepl("3", grid) ~ y_new+20*2,
                             grepl("4", grid) ~ y_new+20*3, 
                             grepl("5", grid) ~ y_new+20*4, 
                             grepl("6", grid) ~ y_new+20*5, 
                             grepl("7", grid) ~ y_new+20*6, 
                             grepl("8", grid) ~ y_new+20*7, .default = NA)) |> 
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         x_coord = ncolumn, 
         rock_cover = as.numeric(rock_cover), 
         mean_soil_depth = as.numeric(mean_soil_depth)) |> 
  ungroup() 


###Fill in the x and y coordinates of cells that do not have SES
comb2 <- comb |> 
  mutate(ncolumn = match(column, LETTERS[1:8]),
         elevation = as.factor(case_when(grepl("BK", Cell_ID) == T ~ "3000", #add elevation variable
                               grepl("WH", Cell_ID) == T ~ "2500",
                               grepl("GG", Cell_ID) == T ~ "2000", .default = NA)),
         y_new = case_when(site == "GG" ~ row, 
                           site == "WH" ~ row+160, 
                           site == "BK" ~ row+160 +140, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                             grepl("2", grid) ~ y_new+20, 
                             grepl("3", grid) ~ y_new+20*2,
                             grepl("4", grid) ~ y_new+20*3, 
                             grepl("5", grid) ~ y_new+20*4, 
                             grepl("6", grid) ~ y_new+20*5, 
                             grepl("7", grid) ~ y_new+20*6, 
                             grepl("8", grid) ~ y_new+20*7, .default = NA)) |> 
  mutate(x_coord = ncolumn) |> 
  ungroup()


Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid))) 


Hdat$pos<- numFactor(Hdat$x_coord, Hdat$y_coord)
Hdat$grid <- as.factor(Hdat$grid)
Hdat$elevation <- as.factor(Hdat$elevation)

Hdat2<- Hdat[-which(is.na(Hdat$SES)), ]
Hdat2<- Hdat2[which(Hdat2$site =="GG"), ]

##Gaussian distribution
test1<- glmmTMB(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height+
                + exp(pos +0|grid), 
                family = gaussian, data = Hdat2)
test1
summary(test1)
test1_res<- simulateResiduals(test1)
plot(test1_res)


##Gamma distribution
Hdat2$SESplus <- Hdat2$SES+ abs(min(Hdat2$SES, na.rm = T)) + 1

test2<- glmmTMB(SESplus ~ rock_cover +
                  (1|grid) + exp(pos +0|grid), 
                family = Gamma(link = "identity"), data = Hdat2) #takes 4 hours
test2
summary(test2)
test2_res <- simulateResiduals(test2)
plot(test2_res) ##doesn't want to simulate these residuals


####t distribution for all grids in GG
test5<- glmmTMB(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height+
                + exp(pos +0|grid), #only add the spatial random effect since variation between grids are likely due to random spatial factors only
                family = t_family(link = "identity"), data = Hdat2) #start 14:39, finish some time before 15:50
test5 #false conversion
summary(test5)
test5_res <- simulateResiduals(test5)
plot(test5_res)
diagnose(test5, check_hessian = F)
#Is it necessary to add the random effect and the spatial random effect?


####skewnormal for all grids in GG

#get starting parameters from the gaussian model
pl <- test1$obj$env$parList()
pl <- pl[lengths(pl) > 0] # Filter out empty components

tic()
test6<- glmmTMB(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height+
                  + exp(pos+0|grid), #only add the spatial random effect since variation between grids are likely due to random spatial factors only
                family = skewnormal(link = "identity"), data = Hdat2, start = pl) 
toc() #564.14 sec elapsed, 9 min
test6
summary(test6)
test6_res <- simulateResiduals(test6)
plot(test6_res)
diagnose(test6, check_hessian = T)#
#Is it necessary to add the random effect and the spatial random effect?


ggplot(Hdat2, aes(x = rock_cover, y = SES)) +
  geom_point()+
  theme_classic()


ggplot(Hdat2, aes(x = northness, y = SES)) +
  geom_point()+
  theme_classic()

ggplot(Hdat2, aes(x = soil_moisture_adj_campaign2, y = SES)) +
  geom_point()+
  theme_classic()

ggplot(Hdat2, aes(x = mean_soil_depth, y = SES)) +
  geom_point()+
  theme_classic()

ggplot(Hdat2, aes(x = slope_height, y = SES)) +
  geom_point()+
  theme_classic()

####t distribution
#only for one grid
Hdat3 <- Hdat2[Hdat2$grid == "GG1", ]
test3<- glmmTMB(SES ~ rock_cover +
                  (1|grid) + exp(pos +0|grid), 
                family = t_family(link = "identity"), data = Hdat3) 
summary(test3)
test3_res <- simulateResiduals(test3)
plot(test3_res) #looks the best so far


###tweedie distribution
test4<- glmmTMB(SESplus ~ rock_cover +
                  (1|grid) + exp(pos+0|grid), 
                family = tweedie(link = "log"), data = Hdat3) 
summary(test4) #doesn't converge!
test4_res <- simulateResiduals(test4)
plot(test4_res)
