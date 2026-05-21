###Modelling with glmmTMB on subsampled data####
###Data is subsampled to eliminate spatial autocorrelation####
library(tidyverse)
library(tidylog)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(tictoc)
library(ggplot2)
library(ggridges)
library(emmeans)
library(multcomp)
library(MuMIn)

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
         zrock_cover = c(scale(rock_cover)), #standardise variables
         znorthness = c(scale(northness)), 
         zsoil_moist = c(scale(soil_moisture_adj_campaign2)), 
         zsoil_depth = c(scale(mean_soil_depth)), 
         zslope_height = c(scale(slope_height))) |> 
  ungroup()



####===SES HEIGHT====####
#isolate SES of height
Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))
  
  
###===Sample uncorrelated cells===####  
#Run Function_select_independent_cells.R
Hdat_subs<- select_independent_cells(Hdat, grid_var = "grid", x = "x_coord", y = "y_coord", value_col = "SES",
                                max_search_radius = 3, lag_threshold = 4)
#between 10 and 5 cells per grid
#lets look at the ones with few cells

##GG6
GG6_subs <- Hdat_subs |> 
  filter(grid == "GG6") |> 
  dplyr::select(Cell_ID)

GG6_full <- Hdat |> 
  filter(grid == "GG6") |> 
  mutate(chosen = case_when(Cell_ID %in% c(GG6_subs$Cell_ID) ~ "yes",
                            is.na(SES) ~ "NA",
                            .default = "no"))

ggplot() +
  geom_tile(data = GG6_full, aes(row, ncolumn, fill = chosen), 
            colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "no"             = "orange",
    "NA"               = "grey",
    "yes"              = "green"))


##BK4
BK4_subs <- Hdat_subs |> 
  filter(grid == "BK4") |> 
  dplyr::select(Cell_ID)

BK4_full <- Hdat |> 
  filter(grid == "BK4") |> 
  mutate(chosen = case_when(Cell_ID %in% c(BK4_subs$Cell_ID) ~ "yes",
                            is.na(SES) ~ "NA",
                            .default = "no"))

ggplot() +
  geom_tile(data = BK4_full, aes(row, ncolumn, fill = chosen), 
            colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "no"             = "orange",
    "NA"               = "grey",
    "yes"              = "green"))


###===Model===####
tmod1<- glmmTMB(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height + (1|grid), 
                family = t_family(link = "identity"), data = Hdat_subs)
toc()
summary(tmod1)
tmod1_res <- simulateResiduals(tmod1)
plot(tmod1_res) #looks ok...

r.squaredGLMM(tmod1)

write.csv(summary(tmod1)$coefficients$cond, "All_data/comm_assembly_results/SES_height_env_model_results.csv")

em_tmod1 <- emmeans(tmod1, specs = "elevation", type = "response")
cld(em_tmod1, Letters = letters, adjust = "Tukey")
#2500 elevation has lower SES than other two


##Variable importance:##
R2full<- r.squaredGLMM(tmod1)[[1]]

predictors <- c("elevation", "zrock_cover", "znorthness","zsoil_moist","zsoil_depth" ,"zslope_height" )

importance <- sapply(predictors, function(var) {
  # Refit without this variable
  f <- as.formula(paste("SES ~", paste(setdiff(predictors, var), collapse = " + "), "+ (1|grid)"))
  m_drop <- glmmTMB(f, data = Hdat, family = t_family(link = "identity"), REML = FALSE)
  
  r2_drop <- r.squaredGLMM(m_drop)[,"R2m"]
  
  R2full - r2_drop  # importance = R² lost by removing this variable
  
})
sort(importance, decreasing = T)




