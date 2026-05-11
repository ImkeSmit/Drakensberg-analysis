####MOdelling with glmmTMB####
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(DHArma)
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
         elevation = case_when(grepl("BK", Cell_ID) == T ~ "3000", #add elevation variable
                               grepl("WH", Cell_ID) == T ~ "2500",
                               grepl("GG", Cell_ID) == T ~ "2000", .default = NA),
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

test1<- glmmTMB(SES ~ rock_cover +
                  (1|grid) + exp(pos +0|grid), 
                family = gaussian, data = Hdat2)
test1
summary(test1)
res<- data.frame(residuals = residuals(test1))
ggplot(res, aes(sample = residuals)) +
  stat_qq(
    color = "#2C7BB6",
    alpha = 0.7,
    size  = 1.8) +
  stat_qq_line(
    color    = "#D7191C",
    linewidth = 1,
    linetype = "dashed") +
  labs( title = "SESplus, Gamma(link = log)",
        x        = "Theoretical Quantiles",
        y        = " Residuals")+
  theme_bw(base_size = 13) 



Hdat2$SESplus <- Hdat2$SES+ abs(min(Hdat2$SES, na.rm = T)) + 1

test2<- glmmTMB(SESplus ~ rock_cover +
                  (1|grid) + exp(pos +0|grid), 
                family = Gamma(link = "identity"), data = Hdat2) #takes 4 hours
test2
summary(test2)
res2<- data.frame(residuals = residuals(test2))

ggplot(res2, aes(sample = residuals)) +
  stat_qq(
    color = "#2C7BB6",
    alpha = 0.7,
    size  = 1.8) +
  stat_qq_line(
    color    = "#D7191C",
    linewidth = 1,
    linetype = "dashed") +
  labs( title = "SESplus, Gamma(link = identity)",
        x        = "Theoretical Quantiles",
        y        = " Residuals")+
  theme_bw(base_size = 13) 
