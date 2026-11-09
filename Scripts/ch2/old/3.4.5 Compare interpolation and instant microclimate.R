###script to investigate hoe well instantaneous and interpolated soil moisture correspond

#import microenvironmental data containing microenv variables
env3 <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row, soil_moisture_adj_campaign2, soil_moisture_adj_campaign1, 
         soil_temperature_adj_campaign1, soil_temperature_adj_campaign2) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA)) |> 
  mutate(ncolumn = match(column, LETTERS[1:8]),  #also add x and y coordinates
         grid = paste0(site, grid), 
         elevation = as.factor(elevation), 
         grid = as.factor(grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row) 

#import interpolated microclimate indices
micro_idw <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices_idw_interpolated.csv", row.names = 1) 

###Correlation between soil moisture campaigns and interpolated soil moisture of the growing season
moisture_dat <- env3 |> 
  inner_join(micro_idw, by = "Cell_ID") |> 
  drop_na()

#campaign 1
cor(moisture_dat$soil_moisture_adj_campaign1, moisture_dat$mean_moist_growing_season)
#0.357 not good!

#campaign 2
cor(moisture_dat$soil_moisture_adj_campaign2, moisture_dat$mean_moist_growing_season)
#0.6347639 a bit better

#campaign 1 and 2
cor(moisture_dat$soil_moisture_adj_campaign2, moisture_dat$soil_moisture_adj_campaign1)
#0.8113965 


###Correlation between soil moisture TEMPERATURE and interpolated soil temperature of the growing season

#campaign 1
cor(moisture_dat$soil_temperature_adj_campaign1, moisture_dat$mean_T1_growing_season)
#0.276079 not good!

#campaign 2
cor(moisture_dat$soil_temperature_adj_campaign2, moisture_dat$mean_T1_growing_season)
#0.6850881 a bit better

#campaign 1 and 2
cor(moisture_dat$soil_temperature_adj_campaign2, moisture_dat$soil_temperature_adj_campaign1)
#0.6857863


##Graph
#temperature
moisture_dat |> 
  pivot_longer(cols = contains(c("campaign", "growing")), names_to = "variable", values_to = "value") |> 
  filter(variable %in% c("mean_T1_growing_season", "soil_temperature_adj_campaign1", "soil_temperature_adj_campaign2")) |> 
ggplot(aes(x = Cell_ID, y = value, color = variable)) +
  geom_point()

#moisture
moisture_dat |> 
  pivot_longer(cols = contains(c("campaign", "growing")), names_to = "variable", values_to = "value") |> 
  filter(variable %in% c("mean_moist_growing_season", "soil_moisture_adj_campaign1", "soil_moisture_adj_campaign2")) |> 
  ggplot(aes(x = Cell_ID, y = value, color = variable)) +
  geom_point() +
  ylim(0, 100)





