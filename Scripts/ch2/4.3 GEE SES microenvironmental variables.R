###Script to analyse relationship between SES and microenvironmental variables
library(spdep)
#library(spatialreg)
#library(spind)
library(tidyverse)
library(tidylog)
library(reshape2)
library(ggplot2)
library(corrplot)
library(Matrix)
library(gee)

#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA), 
         site = case_when(grepl("BK", cellref) == T ~ "BK", #add elevation variable
                          grepl("WH", cellref) == T ~ "WH",
                          grepl("GG", cellref) == T ~ "GG",.default = NA),
         grid = str_sub(cellref, 1, 3), 
         column = str_sub(cellref, 4,4), 
         row = as.numeric(str_sub(cellref, 5,6)), 
         ncolumn = match(column, LETTERS[1:8])) |> 
          mutate(Cell_ID = paste0(site, "_G", str_sub(cellref, 3, 3), "_", column, row), 
                 elevation = as.numeric(elevation)) |> 
          select(-c(cellref, site, grid, column, row)) 

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2)) |> 
  select(!c(soil_temperature_adj_campaign2, soil_moisture_adj_campaign1, soil_temperature_adj_campaign1,
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
  mutate(x_coord = ncolumn, 
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
                            cols_to_impute = colnames(Hdat)[c(6:20, 25)])
Hdat_filled <- Hdat_filled |> 
  arrange(y_coord) |> 
  mutate(grid = as.factor(grid), #cluster ID must be a factor
         elevation = as.numeric(elevation))

#Check what was filled
cols <- c(colnames(Hdat)[c(6:20, 25)])

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
cordf <- Hdat_filled |> select(c(colnames(Hdat)[c(6:13, 15:20, 25)]))
cormat<- cor(cordf)
#cormat[cormat >0.7]
#cormat[cormat <-0.7]
#none are highly correlated
corrplot(cormat, type = "lower", method = "number")

#=====================================================#
####Spatial autocorrelation structure for each grid####
#=====================================================#
#Run Function_grid_correlation_structure.R
decay_df <- grid_correlation_structure(grid_vector = c(unique(Hdat_filled$grid)), 
                                   data = Hdat_filled, 
                                   formula = "SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height", 
                                   k_specified = 4)



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
select(x_coord, row)

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
#LOWEST DECAY RATE
gee_lowest <- gee::gee(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
            family = gaussian, data = Hdat_filled,
            id = grid,
            corstr = "fixed",
            R = lowest_R, #needs to be the same dimension as one group
            scale.fix = T, scale.value = 1, #this is what Pete used in his code 
            silent = F) 

summary(gee_lowest)

#Get p values
coefs <- summary(gee_lowest)$coefficients
p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
cbind(coefs, p_value = round(p_values, 4))


#HIGHEST DECAY RATE
gee_highest <- gee::gee(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
                       family = gaussian, data = Hdat_filled,
                       id = grid,
                       corstr = "fixed",
                       R = highest_R, #needs to be the same dimension as one group
                       scale.fix = T, scale.value = 1, #this is what Pete used in his code 
                       silent = F) 

summary(gee_highest)

#Get p values
coefs <- summary(gee_highest)$coefficients
p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
cbind(coefs, p_value = round(p_values, 4))


#MEAN DECAY RATE
gee_mean <- gee::gee(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
                        family = gaussian, data = Hdat_filled,
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


#===============================#
######PERMUTATION TEST###########
#===============================#

#Build 999 randomised grids
#for now, we will randomise only SES
#Run Function_randomise_grids

random_list <- randomise_grids(data = Hdat_filled, 
                               var = c("rock_cover", "northness", "soil_moisture_adj_campaign2", 
                                       "mean_soil_depth", "slope_height"),
                               iterations = 999)
##Perhaps we should rather be randomising the env variables
#Becuase our null hypothesis is that env conditions do not affect the SES


#Loop through random_list, performing gee and extracting p value for each one
for (l in 1:length(random_list)) {
  data <- random_list[[l]]
  
  one_gee <- gee::gee(randomised_SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
                       family = gaussian, data = data,
                       id = grid,
                       corstr = "fixed",
                       R = mean_R, #needs to be the same dimension as one group
                       scale.fix = T, scale.value = 1, #this is what Pete used in his code 
                       silent = F) 
  
  coefs <- summary(one_gee)$coefficients
  p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
  
  if(l == 1) {
  results_table <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
  results_table$l = l
  results_table$variable <- row.names(results_table)
  row.names(results_table) <- NULL
  } else {
  
  results_temp <- as.data.frame(cbind(coefs, p_value = round(p_values, 4)))
  results_temp$l = l
  results_temp$variable <- row.names(results_temp)
  row.names(results_temp) <- NULL
  
  results_table <- rbind(results_table, results_temp)
  }
}


#Let's look at the p value for rock_cover
p_rock <- results_table |> 
  filter(variable == "rock_cover")

library(ggplot2)
ggplot(p_rock, aes(x = Estimate)) +
  geom_histogram() +
  geom_vline(xintercept = 0.004)

#What is the probability that the estimate of the randomised models different from the observed model?
mean(abs(p_rock$Estimate) >= abs(gee_mean_results$Estimate[2]))
#0.45/45% were as extreme or more extreme than the observed estimate

#Thus, the probability of obtaining this estimate when the null hypothesis is true is large
#too large to reject H0
