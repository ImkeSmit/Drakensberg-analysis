###Modelling with nlme, no subsampling####
library(tidyverse)
library(tidylog)
library(nlme)
library(MuMIn)
library(DHARMa)
library(tictoc)
library(ggridges)
library(emmeans)
library(multcomp)
library(conflicted)
library(performance)
library(see)
conflict_prefer_all("tidylog", quiet = TRUE)


#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  rename(Cell_ID = cellref)

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, 
           soil_depth_CV, mean_soil_depth, slope_height) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA))

#import remote sensing derived variables
rms <- read.csv("All_data/clean_data/Environmental data/Zonal_stats_all.csv") |> 
  select(CELL_ID, STD) |> 
  rename(Cell_ID = CELL_ID)

#import interpolated microclimate indices
micro_idw <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices_idw_interpolated.csv", row.names = 1)


##Combine SES and environmental data
comb <- env |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cell_ses, by = "Cell_ID", relationship = "one-to-many") |>
  #join to microclimate indices |> 
  full_join(micro_idw, by = "Cell_ID") |> 
  #join to remote sensing data |> 
  full_join(rms, by = "Cell_ID") |> 
  mutate(ncolumn = match(column, LETTERS[1:8])) |> 
  rename(x_coord = ncolumn, 
         y_coord = row)


###Check collinearity#### 
library(corrplot)
cordf <- comb |> 
  filter(trait == "Height_cm") |> #look at just one set of env data, it repeats for every trait
  select(mean_T1_growing_season, mean_moist_growing_season, STD, rock_cover, northness, mean_soil_depth, slope_height) |> 
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
  filter(trait == "Height_cm") |> 
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
            correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #exponential correlation structure
            data = Hdat) #only gaussian family possible
toc()
summary(tmod1)
anova(tmod1)



# Check the estimated range and nugget effect
tmod1$modelStruct$corStruct


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

tmod1_res <- simulateResiduals(tmod1)
plot(tmod1_res)

check_model(tmod1)



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

check_model(tmod2)


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




