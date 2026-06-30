###Modelling with nlme, no subsampling####
library(tidyverse)
library(tidylog)
library(nlme)
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
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2), 
         soil_depth_CV = as.numeric(soil_depth_CV)) |> 
  #variables we are interested in
  dplyr::select(c(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, 
           soil_depth_CV, mean_soil_depth, slope_height)) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000, .default = NA))


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


###Check collinearity#### 
library(corrplot)
cordf <- comb2 |> 
  filter(trait == "Height_cm") |> #look at just one set of env data, it repeats for every trait
  dplyr::select(c(zrock_cover, znorthness, zsoil_moist, zsoil_depth, zslope_height)) |> 
  drop_na()
cormat<- cor(cordf)
#cormat[cormat > 0.7]
#cormat[cormat <-0.7]
#none are highly correlated
corrplot(cormat, type = "lower", method = "number")



####SES HEIGHT####
#isolate SES of height
#leave heavy tail
Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA), #also select cells which have no SES measurement. This is necessary to make the grid complete
         !is.na(SES)) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid))) |> 
  drop_na()

#descriptive stats
#how many cells
nrow(Hdat) #2880
Hdat |> group_by(site) |> 
  summarise(n = n())


tic()
tmod1<- lme(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height,
            random = ~1|grid, 
            correlation = corExp(form = ~ x_coord + y_coord | grid, nugget = TRUE), #exponential correlation structure
            data = Hdat) #only gaussian family possible
toc()
summary(tmod1)
anova(tmod1)



# Check the estimated range and nugget effect
tmod1$modelStruct$corStruct
intervals(tmod1, which = "var-cov")


#Compare against a model without spatial structure to see if it improves fit
tic()
tmod1_nonspat<- lme(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
              zslope_height,
            random = ~1|grid, 
            data = Hdat) #only gaussian family possible
toc()

anova(tmod1_nonspat, tmod1)
#spatial model has lower AIC and is a significantly better fit than the nonspatial model


# Residual diagnostics
plot(tmod1)
qqnorm(tmod1, ~ resid(., type = "normalized"))

# Optional: variogram of normalized residuals to visually check
# whether spatial autocorrelation has been adequately captured
plot(Variogram(tmod1, resType = "normalized"))
plot(Variogram(tmod1_nonspat, resType = "normalized"))




####SES SLA####
#isolate SES of SLA
#leave heavy tail
SLAdat <- comb2 |> 
  filter(trait %in% c("SLA", NA), #also select cells which have no SES measurement. This is necessary to make the grid complete
         !is.na(SES)) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = as.factor(paste0(site, grid))) |> 
  drop_na()

#descriptive stats
#how many cells
nrow(SLAdat) #2880
SLAdat |> group_by(site) |> 
  summarise(n = n())


tic()
tmod2<- lme(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height,
            random = ~1|grid, 
            correlation = corExp(form = ~ x_coord + y_coord | grid, nugget = TRUE),
            data = SLAdat)
toc()
summary(tmod2)
anova(tmod2)


#Compare against a model without spatial structure to see if it improves fit
tic()
tmod2_nonspat<- lme(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                      zslope_height,
                    random = ~1|grid, 
                    data = SLAdat) #only gaussian family possible
toc()

anova(tmod2_nonspat, tmod2)
#spatial model has lower AIC and is a significantly better fit than the nonspatial model


# Residual diagnostics
plot(tmod2) #looks pretty good
qqnorm(tmod2, ~ resid(., type = "normalized")) #pretty ok

# Optional: variogram of normalized residuals to visually check
# whether spatial autocorrelation has been adequately captured
plot(Variogram(tmod2, resType = "normalized"))
plot(Variogram(tmod2_nonspat, resType = "normalized"))




