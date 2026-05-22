###Script to analyse relationship between SES and microenvironmental variables
library(spdep)
library(tidyverse)
library(tidylog)
library(reshape2)
library(ggplot2)
library(corrplot)
library(Matrix)
library(gee)

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
  mutate(x_coord = ncolumn, 
         zrock_cover = c(scale(rock_cover)), #standardise variables
         znorthness = c(scale(northness)), 
         zsoil_moist = c(scale(soil_moisture_adj_campaign2)), 
         zsoil_depth = c(scale(mean_soil_depth)), 
         zslope_height = c(scale(slope_height))) |> 
  ungroup()

#=====================================
#####GEE FOR SES OF HEIGHT####
#=====================================
Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = paste0(site, grid))  #grid variable must be unique accross sites for the impute function
  
#All grids must have 160 cells, check this
check <- Hdat |> group_by(site, grid) |> 
  summarise(n = n()) #all ok

#===============#
###Imputation####
#===============#
#now we need to impute the missing SES or predictor variables because the Gee won't work if there are NA's
#for now, we fill fill the NA cells with the mean of it's 8 nearest neighbours
#run Function_impute_cells.R
Hdat_filled <- impute_cells(df = Hdat, 
                            cols_to_impute = colnames(Hdat)[c(25, 31:35)], 
                            neighbours = 4)
Hdat_filled <- Hdat_filled |> 
  arrange(y_coord) |> 
  mutate(grid = as.factor(grid), #cluster ID must be a factor
         elevation = as.factor(elevation))

#Check what was filled
cols <- c(colnames(Hdat)[c(25, 31:35)])

cat("=== NA summary before and after imputation ===\n\n")
for (col in cols) {
  if (!col %in% names(Hdat)) next
  n_before <- sum(is.na(Hdat[[col]]))
  n_after  <- sum(is.na(Hdat_filled[[col]]))
  cat(sprintf(
    "%-35s  NAs before: %3d  |  NAs after: %3d  |  Filled: %3d\n",
    col, n_before, n_after, n_before - n_after
  ))
}


###Check collinearity#### 
library(corrplot)
cordf <- Hdat_filled |> dplyr::select(c(zrock_cover, znorthness, zsoil_moist, zsoil_depth, zslope_height)) |> 
  drop_na()
cormat<- cor(cordf)
#cormat[cormat > 0.7]
#cormat[cormat <-0.7]
#none are highly correlated
corrplot(cormat, type = "lower", method = "number")

#=====================================================#
####Spatial autocorrelation structure for each grid####
#=====================================================#
#Run Function_grid_correlation_structure.R
##Cannot put elevation in this formula, becaus eit looks at correlation per grid
decay_df <- grid_correlation_structure(grid_vector = c(unique(Hdat_filled$grid)), 
                                   data = Hdat_filled, 
                                   formula = "SES ~ zrock_cover + znorthness + zsoil_moist + zsoil_depth + zslope_height", 
                                   k_specified = 4)
##GG7 doesn't print, possibly because it has too many NA's


###Build 3 correlation matrices####
#with lowest b, mean b and highest b
#function for correlation matrices of one grid
make_grid_corr <- function(b, first_nonsig_lag, coords) {
  
  # Euclidean distance matrix in grid-step units
  dmat <- as.matrix(dist(coords[, c("x_coord", "row")])) #row is in steps of 1-20
  
  # Exponential correlation, zeroed beyond first_nonsig_lag steps
  corr <- exp(-b * dmat)
  corr[dmat > first_nonsig_lag] <- 0 #replace with range_dist to make it more conservative
  diag(corr) <- 1
  corr[corr < 0] <- 0 #replace all negative values with zero, we do not model negative spatial autocorrelation here
  
  corr
}#end make grid corr

#get lowest, highest and mean distance decay
lowest_b <- decay_df |> #slowest distance decay
  filter(b == min(b, na.rm = T))

highest_b <- decay_df |> #fastest distance decay
  filter(b == max(b, na.rm = T))

mean_b <- decay_df |> 
  summarise(b = mean(b, na.rm = T), 
            first_nonsig_lag = ceiling(mean(first_nonsig_lag, na.rm = T)))

coordinates_one_grid <- Hdat_filled |> #get the coordinates of one grid 
  filter(grid == "GG1") |> 
  dplyr::select(x_coord, row)

#lowest distance decay
lowest_R <- make_grid_corr(b = lowest_b$b, 
                           first_nonsig_lag = lowest_b$first_nonsig_lag, 
                           coords = coordinates_one_grid)

highest_R <- make_grid_corr(b = highest_b$b, 
                            first_nonsig_lag = highest_b$first_nonsig_lag, 
                            coords = coordinates_one_grid)

mean_R <- make_grid_corr(b = mean_b$b, 
                         first_nonsig_lag = mean_b$first_nonsig_lag, 
                         coords = coordinates_one_grid)

#map R as a sanity check
df <- melt(highest_R)
colnames(df) <- c("row", "col", "correlation")

pdf("Figures/R_plot.pdf")
ggplot(df, aes(x = col, y = row, fill = correlation)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_flip() +
  theme_minimal() 
dev.off()

#=============#
####Run GEE####
#=============#
#MEAN DECAY RATE
Hdat_filled$f_elevation<- as.factor(Hdat_filled$elevation)
gee_mean <- gee::gee(SES ~ f_elevation+ zrock_cover + znorthness + zsoil_moist + zsoil_depth + zslope_height,
                        family = gaussian, 
                        data = Hdat_filled,
                        id = grid,
                        corstr = "fixed",
                        R = mean_R, #needs to be the same dimension as one group
                        scale.fix = T, scale.value = 1, #this is what Pete used in his code 
                        silent = F) 

summary(gee_mean)
#Get p values
coefs <- summary(gee_mean)$coefficients
p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
gee_mean_results <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
gee_mean_results$variable <- row.names(gee_mean_results)
row.names(gee_mean_results) <- NULL
gee_mean_results <- gee_mean_results[, c(7,1,2,3,4,5,6)] #reorder columns
write.csv(gee_mean_results, "All_data/comm_assembly_results//GEE_SES_height_env_model_results.csv")

resid_df<- data.frame(residuals = gee_mean$residuals)
hist(resid_df$residuals)

ggplot(resid_df, aes(sample = residuals)) +
  stat_qq(
    color = "#2C7BB6",
    alpha = 0.7,
    size  = 1.8) +
  stat_qq_line(
    color    = "#D7191C",
    linewidth = 1,
    linetype = "dashed") +
  labs( title = "SES, gaussian",
    x        = "Theoretical Quantiles",
    y        = " Residuals")+
  theme_bw(base_size = 13) 


#===============================#
######PERMUTATION TEST###########
#===============================#

#Build 999 randomised grids
#Run Function_randomise_grids

random_list <- randomise_grids(data = Hdat_filled, 
                               var = c("rock_cover", "northness", "soil_moisture_adj_campaign2", 
                                       "mean_soil_depth", "slope_height"),
                               iterations = 999)
#We randomise the env variables becuase our null hypothesis is that env conditions do not affect the SES
saveRDS(random_list, file = "All_data//comm_assembly_results//randomised_env_grids.rds")

#Loop through random_list, performing gee and extracting p value for each one
random_list <- readRDS("All_data//comm_assembly_results//randomised_env_grids.rds")
for (l in 1:length(random_list)) {
  data <- random_list[[l]]
  
  one_gee <- gee::gee(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
                       family = gaussian, data = data,
                       id = grid,
                       corstr = "independence", #because spatial autocorrelation was removed during randomisation
                       scale.fix = T, scale.value = 1, #this is what Pete used in his code 
                       silent = F) 
  
  coefs <- summary(one_gee)$coefficients
  p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
  
  if(l == 1) {
  results_random <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
  results_random$l = l
  results_random$variable <- row.names(results_random)
  row.names(results_random) <- NULL
  } else {
  
  results_temp <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
  results_temp$l = l
  results_temp$variable <- row.names(results_temp)
  row.names(results_temp) <- NULL
  
  results_random <- rbind(results_random, results_temp)
  }
}


###Compute p values from permutation test
vars <- c(unique(results_random$variable))
for (v in vars) {
  one_var <- results_random |> 
    filter(variable == v)
  
  obs <- gee_mean_results[which(gee_mean_results$variable == v), which(colnames(gee_mean_results) == "Estimate")]
  se <- gee_mean_results[which(gee_mean_results$variable == v), which(colnames(gee_mean_results) == "Robust S.E.")]
  
  #What is the probability that the estimate of the randomised models different from the observed model?
  #mean gets the proportion of values that are greater or equal to the observed model estimate
  #0.12/12% were as extreme or more extreme than the observed estimate
  p_val <- length(which(abs(one_var$Estimate) >= abs(obs)) == T)/(nrow(one_var)+1) #add one so that p value cannot be zero
  perc_2.5 <- quantile(one_var$Estimate, probs = 0.025)
  perc_97.5 <- quantile(one_var$Estimate, probs = 0.975)
  significance <- obs < perc_2.5 | obs > perc_97.5
  
  
  if(v == "(Intercept)") {
    p_table <- data.frame(variable = v, observed_estimate = obs,
                          observed_se = se, p_value = p_val, 
                          perc_2.5 = perc_2.5, perc_97.5 = perc_97.5, 
                          significant = significance)
  }else {
    p_temp <- data.frame(variable = v, observed_estimate = obs,
                         observed_se = se, p_value = p_val, 
                         perc_2.5 = perc_2.5, perc_97.5 = perc_97.5, 
                         significant = significance)
    p_table <- rbind(p_table, p_temp)
  }
}#end loop through vars
#save results
write.csv(p_table, "All_data//comm_assembly_results//SES_height_env_results.csv")


###Graph the null distributions and the observed values
facet_names <- c("Intercept", "Rock cover", "Northness", "Soil moisture", "soil depth", "slope height")
names(facet_names) <- c(p_table$variable)

perm_test <- ggplot(results_random, aes(x = Estimate)) +
  geom_histogram() +
  geom_vline(data = p_table, aes(xintercept = observed_estimate), color = "red")+
  facet_wrap(~variable, scales = "free", labeller = as_labeller(facet_names))

ggsave("Figures//permutation_test.png" ,perm_test)




#==============================#
#####GEE FOR SES OF SLA#########
#==============================#
SLAdat <- comb2 |> 
  filter(trait %in% c("SLA", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = paste0(site, grid))  #grid variable must be unique accross sites for the impute function

#All grids must have 160 cells, check this
check <- SLAdat |> group_by(site, grid) |> 
  summarise(n = n()) #all ok

#===============#
###Imputation####
#===============#
#Everything ran with the SES NA's, maybe imputation is not necessary?
#now we need to impute the missing SES or predictor variables because the Gee won't work if there are NA's
#for now, we fill fill the NA cells with the mean of it's 8 nearest neighbours
#run Function_impute_cells.R
SLAdat_filled <- impute_cells(df = SLAdat, 
                            cols_to_impute = colnames(SLAdat)[c(25, 31:35)])
SLAdat_filled <- SLAdat_filled |> 
  arrange(y_coord) |> 
  mutate(grid = as.factor(grid), #cluster ID must be a factor
         elevation = as.factor(elevation))

#Check what was filled
cols <- c(colnames(SLAdat)[c(25, 31:35)])

cat("=== NA summary before and after imputation ===\n\n")
for (col in cols) {
  if (!col %in% names(SLAdat)) next
  n_before <- sum(is.na(SLAdat[[col]]))
  n_after  <- sum(is.na(SLAdat_filled[[col]]))
  cat(sprintf(
    "%-35s  NAs before: %3d  |  NAs after: %3d  |  Filled: %3d\n",
    col, n_before, n_after, n_before - n_after
  ))
} ###3 NA values remain because all of their neighbours are NA


#=====================================================#
####Spatial autocorrelation structure for each grid####
#=====================================================#
#Run Function_grid_correlation_structure.R
decay_df <- grid_correlation_structure(grid_vector = c(unique(SLAdat_filled$grid)), 
                                       data = SLAdat_filled, 
                                       formula = "SES ~ zrock_cover + znorthness + zsoil_moist + zsoil_depth + zslope_height", 
                                       k_specified = 4)
#again GG7 not printed, probably too many Na's

###Build correlation matrix####
#function for correlation matrices of one grid
make_grid_corr <- function(b, first_nonsig_lag, coords) {
  
  # Euclidean distance matrix in grid-step units
  dmat <- as.matrix(dist(coords[, c("x_coord", "row")])) #row is in steps of 1-20
  
  # Exponential correlation, zeroed beyond first_nonsig_lag steps
  corr <- exp(-b * dmat)
  corr[dmat > first_nonsig_lag] <- 0 #replace with range_dist to make it more conservative
  diag(corr) <- 1
  corr[corr < 0] <- 0 #replace all negative values with zero, we do not model negative spatial autocorrelation here
  
  corr
}#end make grid corr

#mean distance decay
mean_b <- decay_df |> 
  summarise(b = mean(b, na.rm = T), 
            first_nonsig_lag = ceiling(mean(first_nonsig_lag, na.rm = T)))

coordinates_one_grid <- SLAdat_filled |> #get the coordinates of one grid 
  filter(grid == "GG1") |> 
  dplyr::select(x_coord, row)


mean_R <- make_grid_corr(b = mean_b$b, 
                         first_nonsig_lag = mean_b$first_nonsig_lag, 
                         coords = coordinates_one_grid)

#map R as a sanity check
df <- melt(mean_R)
colnames(df) <- c("row", "col", "correlation")

pdf("Figures/R_plot.pdf")
ggplot(df, aes(x = col, y = row, fill = correlation)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_flip() +
  theme_minimal() 
dev.off()

#=============#
####Run GEE####
#=============#
gee_SLA <- gee::gee(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + zslope_height,
                       family = gaussian, data = SLAdat_filled,
                       id = grid,
                       corstr = "fixed",
                       R = mean_R, #needs to be the same dimension as one group
                       scale.fix = T, scale.value = 1, #this is what Pete used in his code 
                       silent = F) 

#Get p values
coefs <- summary(gee_SLA)$coefficients
p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
gee_SLA_results <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
gee_SLA_results$variable <- row.names(gee_SLA_results)
row.names(gee_SLA_results) <- NULL
gee_SLA_results <- gee_SLA_results[, c(7,1,2,3,4,5,6)] #reorder columns
write.csv(gee_SLA_results, "All_data/comm_assembly_results//GEE_SES_SLA_env_model_results.csv")

resid_df<- data.frame(residuals = gee_SLA$residuals)
hist(resid_df$residuals)

ggplot(resid_df, aes(sample = residuals)) +
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


#==================================#
#        Permutation test          #
#==================================#

#Build 999 randomised grids
#Run Function_randomise_grids

random_list <- randomise_grids(data = SLAdat_filled, 
                               var = c("rock_cover", "northness", "soil_moisture_adj_campaign2", 
                                       "mean_soil_depth", "slope_height"),
                               iterations = 999)
#We randomise the env variables becuase our null hypothesis is that env conditions do not affect the SES
saveRDS(random_list, file = "All_data//comm_assembly_results//randomised_env_grids_SLA.rds")

#Loop through random_list, performing gee and extracting p value for each one
random_list <- readRDS("All_data//comm_assembly_results//randomised_env_grids_SLA.rds")
for (l in 1:length(random_list)) {
  data <- random_list[[l]]
  
  one_gee <- gee::gee(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
                      family = gaussian, data = data,
                      id = grid,
                      corstr = "independence", #because spatial autocorrelation was removed during randomisation
                      scale.fix = T, scale.value = 1, #this is what Pete used in his code 
                      silent = F) 
  
  coefs <- summary(one_gee)$coefficients
  p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
  
  if(l == 1) {
    results_random <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
    results_random$l = l
    results_random$variable <- row.names(results_random)
    row.names(results_random) <- NULL
  } else {
    
    results_temp <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
    results_temp$l = l
    results_temp$variable <- row.names(results_temp)
    row.names(results_temp) <- NULL
    
    results_random <- rbind(results_random, results_temp)
  }
}


###Compute p values from permutation test
vars <- c(unique(results_random$variable))
for (v in vars) {
  one_var <- results_random |> 
    filter(variable == v)
  
  obs <- gee_mean_results[which(gee_mean_results$variable == v), which(colnames(gee_mean_results) == "Estimate")]
  se <- gee_mean_results[which(gee_mean_results$variable == v), which(colnames(gee_mean_results) == "Robust S.E.")]
  
  #What is the probability that the estimate of the randomised models different from the observed model?
  #mean gets the proportion of values that are greater or equal to the observed model estimate
  #0.12/12% were as extreme or more extreme than the observed estimate
  p_val <- length(which(abs(one_var$Estimate) >= abs(obs)) == T)/(nrow(one_var)+1) #add one so that p value cannot be zero
  perc_2.5 <- quantile(one_var$Estimate, probs = 0.025)
  perc_97.5 <- quantile(one_var$Estimate, probs = 0.975)
  significance <- obs < perc_2.5 | obs > perc_97.5
  
  
  if(v == "(Intercept)") {
    p_table <- data.frame(variable = v, observed_estimate = obs,
                          observed_se = se, p_value = p_val, 
                          perc_2.5 = perc_2.5, perc_97.5 = perc_97.5, 
                          significant = significance)
  }else {
    p_temp <- data.frame(variable = v, observed_estimate = obs,
                         observed_se = se, p_value = p_val, 
                         perc_2.5 = perc_2.5, perc_97.5 = perc_97.5, 
                         significant = significance)
    p_table <- rbind(p_table, p_temp)
  }
}#end loop through vars

#save results
write.csv(p_table, "All_data//comm_assembly_results//SES_SLA_env_results.csv")


###Graph the null distributions and the observed values
facet_names <- c("Intercept", "Rock cover", "Northness", "Soil moisture", "soil depth", "slope height")
names(facet_names) <- c(p_table$variable)

perm_test <- ggplot(results_random, aes(x = Estimate)) +
  geom_histogram() +
  geom_vline(data = p_table, aes(xintercept = observed_estimate), color = "red")+
  facet_wrap(~variable, scales = "free", labeller = as_labeller(facet_names))

ggsave("Figures//permutation_test.png" ,perm_test)


  

