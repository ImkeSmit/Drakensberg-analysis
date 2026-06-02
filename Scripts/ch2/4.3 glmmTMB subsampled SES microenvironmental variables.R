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
                                max_search_radius = 2)

#Save the cells that we use
Height_incl_cells <- data.frame(trait = "Height_cm", included_cells = c(unique(Hdat_subs$Cell_ID)))
write.csv(Height_incl_cells, "All_data/comm_assembly_results/included_cells_Height.csv")


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


##BK7
BK7_subs <- Hdat_subs |> 
  filter(grid == "BK7") |> 
  dplyr::select(Cell_ID)

BK7_full <- Hdat |> 
  filter(grid == "BK7") |> 
  mutate(chosen = case_when(Cell_ID %in% c(BK7_subs$Cell_ID) ~ "yes",
                            is.na(SES) ~ "NA",
                            .default = "no"))

ggplot() +
  geom_tile(data = BK7_full, aes(row, ncolumn, fill = chosen), 
            colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "no"             = "orange",
    "NA"               = "grey",
    "yes"              = "green"))


####Graph SES~env variables
Hdat_subs |> 
  dplyr::select(c(SES, rock_cover , northness , soil_moisture_adj_campaign2 , mean_soil_depth , slope_height)) |> 
  pivot_longer(cols = !SES, names_to = "variable", values_to = "value") |> 
  ggplot(aes(y = SES, x = value)) +
  geom_point()+
  facet_wrap(~variable, scales = "free_x")+
  theme_bw()


###===Descriptive stats SES Height===####
#ncells
length(unique(Hdat_subs$Cell_ID)) #196

#nspecies
included_cells <- (unique(Hdat_subs$Cell_ID))

abun_matrix <- read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)
mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1)

subs_abun_matrix<- abun_matrix[rownames(abun_matrix) %in% c(included_cells), ]
subs_abun_matrix2 <- subs_abun_matrix[, colSums(subs_abun_matrix) > 0, drop = F]

nsp <- ncol(subs_abun_matrix2) #138

#how many sp from each site? 
incl_sp <- c(colnames(subs_abun_matrix2))



###===Model===####
tmod1<- glmmTMB(SES ~ elevation + zrock_cover +  znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height + (1|grid), 
                family = t_family(link = "identity"), 
                data = Hdat_subs) 
#performance::check_singularity(tmod1) #singular fit uh oh

#update the modle with gamme priors to see if that helps
#prior <- data.frame(
#  prior = "gamma(100, 2.5)",  # mean can be 1, but even 1e8
#  class = "ranef" )          # for random effects

#tmod1_update <- update(tmod1, priors = prior) 
#check_singularity(tmod1_update) #not singular anymore


summary(tmod1)
Anova(tmod1)
tmod1_res <- simulateResiduals(tmod1)
plot(tmod1_res) #looks ok...
#sample size large enough with t family if lag threshold = 4 and search radius = 2

r.squaredGLMM(tmod1) #0.09365263 0.1094613


em_tmod1 <- emmeans(tmod1, specs = "elevation", type = "response")
cld(em_tmod1, Letters = letters, adjust = "Tukey")

#test for spatial autocorrelation
used_rows <- as.integer(rownames(model.frame(tmod1)))
dat_used  <- Hdat_subs[used_rows, ]
testSpatialAutocorrelation(tmod1_res, x = dat_used$x_coord, y = dat_used$y_coord)
#Still have autocorrelation... Check if we are simulating residuals correctly

##Variable importance:##
R2full<- r.squaredGLMM(tmod1)[[1]]

predictors <- c("elevation", "zrock_cover", "znorthness","zsoil_moist","zsoil_depth" ,"zslope_height" )
importance <- c(rep(NA, length(predictors)))
names(importance) <- predictors

for(var in predictors) {
  f <- as.formula(paste("SES ~", paste(setdiff(predictors, var), collapse = " + "), "+ (1|grid)"))
  
  m_drop <- glmmTMB(f, data = Hdat_subs, family = t_family(link = "identity"), REML = FALSE)
  
  r2_drop <- r.squaredGLMM(m_drop)[,"R2m"]
  
  imp <- R2full - r2_drop  # importance = R² lost by removing this variable
  
  importance[which(names(importance) == var)] <- imp
}
sort(importance, decreasing = T)
##Some models not converging, how to fix???

#save results
tmod1_results <-as.data.frame(summary(tmod1)$coefficients$cond)
tmod1_results$variable_importance <- c(0,0,importance)
tmod1_results$emmean <- c(-0.238,-0.437, -0.301, 0,0,0,0,0)
tmod1_results$em_std_er <- c(0.0527,0.0290,0.0279,0,0,0,0,0)
model <- paste(as.character(formula(tmod1))[[2]], as.character(formula(tmod1))[[1]],as.character(formula(tmod1))[[3]],
               "; ", "family = " , 
               as.character(family(tmod1)[[1]]), "link = ", 
               as.character(family(tmod1)[[3]]))
tmod1_results$model <- model
write.csv(tmod1_results, "All_data/comm_assembly_results/glmmTMB_subsampled_SES_height_env_model_results.csv")



####===SES of SLA====####
#isolate SES of SLA
#leave heavy tail
SLAdat <- comb2 |> 
  filter(trait %in% c("SLA", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))


#=====sample uncorrelated cells===####
#Run Function_select_independent_cells.R
SLAdat_subs<- select_independent_cells(SLAdat, grid_var = "grid", x = "x_coord", y = "y_coord", value_col = "SES",
                                     max_search_radius = 2)

#Save the cells that we use
SLA_incl_cells <- data.frame(trait = "SLA", included_cells = c(unique(SLAdat_subs$Cell_ID)))
write.csv(SLA_incl_cells, "All_data/comm_assembly_results/included_cells_SLA.csv")

#check sampling:
data_filled <- impute_cells(df = SLAdat, 
                            cols_to_impute = colnames(SLAdat)[c(25, 31:35)], 
                            neighbours = 4)

#WH5 has 64 cells, lets look at the correlation structure
grid_correlation_structure(grid_vector = c(unique(data_filled$grid)), data_filled, 
                           formula = "SES ~ zrock_cover + znorthness + zsoil_moist + zsoil_depth + zslope_height", 
                          k_specified = 4)
#Jup, WH5 already has nonsigificant lag at 1

##WH5
WH5_subs <- SLAdat_subs |> 
  filter(grid == "WH5") |> 
  dplyr::select(Cell_ID)

WH5_full <- SLAdat |> 
  filter(grid == "WH5") |> 
  mutate(chosen = case_when(Cell_ID %in% c(WH5_subs$Cell_ID) ~ "yes",
                            is.na(SES) ~ "NA",
                            .default = "no"))

ggplot() +
  geom_tile(data = WH5_full, aes(row, ncolumn, fill = chosen), 
            colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "no"             = "orange",
    "NA"               = "grey",
    "yes"              = "green"))


###===Descriptive stats SES SLA===####
#ncells
length(unique(SLAdat_subs$Cell_ID)) #345

#nspecies
included_cells <- (unique(SLAdat_subs$Cell_ID))

abun_matrix <- read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)
mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1)

subs_abun_matrix<- abun_matrix[rownames(abun_matrix) %in% c(included_cells), ]
subs_abun_matrix2 <- subs_abun_matrix[, colSums(subs_abun_matrix) > 0, drop = F]

nsp <- ncol(subs_abun_matrix2) #144



#====Model====#####
tmod2<- glmmTMB(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height + (1|grid), 
                family = t_family(link = "identity"), data = SLAdat_subs)
summary(tmod2)
tmod2_res <- simulateResiduals(tmod2)
plot(tmod2_res) #looks ok!

#test for spatial autocorrelation
used_rows <- as.integer(rownames(model.frame(tmod2)))
dat_used  <- SLAdat_subs[used_rows, ]
testSpatialAutocorrelation(tmod2_res, x = dat_used$x_coord, y = dat_used$y_coord)
#still has significant spatial autocorrelation, but may need to simulate residuals differently??

em_tmod2 <- emmeans(tmod2, specs = "elevation", type = "response")
cld(em_tmod2, Letters = letters, adjust = "Tukey")
#2500 elevation has higher SES than other two

r.squaredGLMM(tmod2) #[1,] 0.221175 0.2699643

##Variable importance:##
R2full<- r.squaredGLMM(tmod2)[[1]]

predictors <- c("elevation", "zrock_cover", "znorthness","zsoil_moist","zsoil_depth" ,"zslope_height" )
importance_SLA <- c(rep(NA, length(predictors)))
names(importance_SLA) <- predictors

for(var in predictors) {
  f <- as.formula(paste("SES ~", paste(setdiff(predictors, var), collapse = " + "), "+ (1|grid)"))
  
  m_drop <- glmmTMB(f, data = SLAdat_subs, family = t_family(link = "identity"), REML = FALSE)
  
  r2_drop <- r.squaredGLMM(m_drop)[,"R2m"]
  
  imp <- R2full - r2_drop  # importance = R² lost by removing this variable
  
  importance_SLA[which(names(importance_SLA) == var)] <- imp
}
sort(importance_SLA, decreasing = T)
#save results
tmod2_results <-as.data.frame(summary(tmod2)$coefficients$cond)
tmod2_results$variable_importance <- c(0,0,importance_SLA)
write.csv(tmod2_results, "All_data/comm_assembly_results/glmmTMB_subsampled_SES_SLA_env_model_results.csv")




