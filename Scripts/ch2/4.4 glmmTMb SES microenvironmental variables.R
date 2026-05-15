####MOdelling with glmmTMB####
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(tictoc)
library(ggplot2)
library(ggridges)
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
                           site == "WH" ~ row+160+1000, 
                           site == "BK" ~ row+160+140+2000, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                             grepl("2", grid) ~ y_new+20+100, 
                             grepl("3", grid) ~ y_new+20*2+200,
                             grepl("4", grid) ~ y_new+20*3+300, 
                             grepl("5", grid) ~ y_new+20*4+400, 
                             grepl("6", grid) ~ y_new+20*5+500, 
                             grepl("7", grid) ~ y_new+20*6+600, 
                             grepl("8", grid) ~ y_new+20*7+700, .default = NA)) |> 
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
                           site == "WH" ~ row+160+1000, 
                           site == "BK" ~ row+160 +140+2000, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                             grepl("2", grid) ~ y_new+20+100, 
                             grepl("3", grid) ~ y_new+20*2+200,
                             grepl("4", grid) ~ y_new+20*3+300, 
                             grepl("5", grid) ~ y_new+20*4+400, 
                             grepl("6", grid) ~ y_new+20*5+500, 
                             grepl("7", grid) ~ y_new+20*6+600, 
                             grepl("8", grid) ~ y_new+20*7+700, .default = NA)) |> 
  mutate(x_coord = ncolumn, 
         zrock_cover = scale(rock_cover), #standardise variables
         znorthness = scale(northness), 
         zsoil_moist = scale(soil_moisture_adj_campaign2), 
         zsoil_depth = scale(mean_soil_depth), 
         zslope_height = scale(slope_height)) |> 
  ungroup()


####===========================####
#=======POSTER ANALYSIS===========#

####SES HEIGHT####
#isolate SES of height
#leave heavy tail
Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA), #also select cells which have no SES measurement. This is necessary to make the grid complete
         !is.na(SES)) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))


tic()
tmod1<- glmmTMB(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                      zslope_height + (1|grid), 
                    family = t_family(link = "identity"), data = Hdat)
toc()
summary(tmod1)
tmod1_res <- simulateResiduals(tmod1)
plot(tmod1_res) #looks ok...

write.csv(summary(tmod1)$coefficients$cond, "All_data/comm_assembly_results/SES_height_env_model_results.csv")

em_tmod1 <- emmeans(tmod1, specs = "elevation", type = "response")
cld(em_tmod1, Letters = letters, adjust = "Tukey")
#2500 elevation has lower SES than other two

####SES SLA####
#isolate SES of SLA
#leave heavy tail
SLAdat <- comb2 |> 
  filter(trait %in% c("SLA", NA), #also select cells which have no SES measurement. This is necessary to make the grid complete
         !is.na(SES)) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))


tic()
tmod2<- glmmTMB(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height + (1|grid), 
                family = t_family(link = "identity"), data = SLAdat)
toc()
summary(tmod2)
tmod2_res <- simulateResiduals(tmod2)
plot(tmod2_res) #looks good
write.csv(summary(tmod2)$coefficients$cond, "All_data/comm_assembly_results/SES_SLA_env_model_results.csv")


em_tmod2 <- emmeans(tmod2, specs = "elevation", type = "response")
cld(em_tmod2, Letters = letters, adjust = "Tukey")
#2500 elevation has higher SES than other two


###Figures###
ses_ridges <- comb |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()




#=================================#
#OLD CODE BELOW- TRYING OUT MODELS#
Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA), #also select cells which have no SES measurement. This is necessary to make the grid complete
         !is.na(SES), 
         SES<=2) |> ##remove very large SES
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))


###sample half the data from each grid####
positions <- Hdat |> 
  dplyr::select(x_coord, y_coord, grid) |> 
  mutate(pos = numFactor(x_coord, y_coord))  

  for(g in c(unique(Hdat$grid))) {
    one_grid <- positions |> filter(grid == g)
    
    if(nrow(one_grid)>80) {
    keep <- sample(one_grid$pos, 80)
    }else{keep <- one_grid$pos}
    
    if(g == unique(Hdat$grid[1])) {
      keep_vector <- keep
    }else {
      temp_vector <- keep
      keep_vector <- c(keep_vector, temp_vector)
    }
  }

Hdat2 <- Hdat |> 
  filter(pos %in% keep_vector)


##First run Gaussian, without spatial decay
#WORKS#
tic()
gausmod<- glmmTMB(SES ~ elevation + rock_cover+ northness + soil_moisture_adj_campaign2 + mean_soil_depth + 
                    slope_height+ (1|grid), 
                family = gaussian, data = Hdat)
toc()
summary(gausmod)
gausmod_res <- simulateResiduals(gausmod)
plot(gausmod_res)
#test for spatial autocorrelation in residuals
used_rows <- as.integer(rownames(model.frame(gausmod))) #get data actually used in model
dat_used  <- Hdat[used_rows, ]
testSpatialAutocorrelation(gausmod_res, x = dat_used$x_coord, y = dat_used$y_coord)
#normality and HOV assumption looks ok
#spatial autocorrelation is significant

#get starting parameters from the gaussian model
pl <- gausmod$obj$env$parList()
pl <- pl[lengths(pl) > 0] # Filter out empty components



##Run Gaussian, WITH spatial decay
#DOES NOT WORK#
Hdat2$SESplus = Hdat2$SES+ abs(min(Hdat2$SES))
tic()
gausmod_spat<- glmmTMB(SESplus ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                    zslope_height + exp(pos+0|grid), 
                  family = gaussian, data = Hdat2)
toc()
summary(gausmod_spat)
gausmod_spat_res <- simulateResiduals(gausmod_spat)
plot(gausmod_spat_res)
#test for spatial autocorrelation in residuals
used_rows <- as.integer(rownames(model.frame(gausmod_spat))) #get data actually used in model
dat_used  <- Hdat2[used_rows, ]
testSpatialAutocorrelation(gausmod_res, x = dat_used$x_coord, y = dat_used$y_coord)
#Error in fitTMB(TMBStruc) : 
#negative log-likelihood is NaN at starting parameter values



###Now run skewnormal
##DOES NOT WORK###
tic()
snmod<- glmmTMB(SES ~ elevation + zrock_cover+ znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height + (1|grid)+ exp(pos+0|grid), 
                  family = skewnormal(link = "identity"), data = Hdat2)
toc()
summary(snmod) #doesnt converge
sn_res <- simulateResiduals(snmod)
plot(sn_res)
#Error in fitTMB(TMBStruc) : 
#negative log-likelihood is NaN at starting parameter values
#does not converge with any combo of random effects or link function


##Run t distribution
tic()
tmod_spat<- glmmTMB(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                         zslope_height + (1|grid), 
                       family = t_family(link = "identity"), data = Hdat2)
toc()
summary(tmod_spat)
tmod_spat_res <- simulateResiduals(tmod_spat)
plot(tmod_spat_res)
##THIS MODEL HAS THE BEST DIAGNOSIC PLOTS SO FAR
#spatial model with T distribution does not work





##Gaussian distribution
test1<- glmmTMB(SES ~ northness + soil_moisture_adj_campaign2 + mean_soil_depth + 
                 exp(pos +0|grid), 
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
test6<- glmmTMB(SES ~ northness + soil_moisture_adj_campaign2 + mean_soil_depth 
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
