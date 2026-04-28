###Script to analyse relationship between SES and microenvironmental variables
library(spdep)
library(spatialreg)
library(spind)
library(tidyverse)
library(tidylog)
library(corrplot)
library(Matrix)

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

###Imputation####
#now we need to impute the missing SES or predictor variables because the Gee won't work if there are NA's
#for now, we fill fill the NA cells with the mean of it's 8 nearest neighbours
#run Function_impute_cells.R
Hdat_filled <- impute_cells(df = Hdat, 
                            cols_to_impute = colnames(Hdat)[c(6:20, 25)])

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
cormat[cormat >0.7]
cormat[cormat <-0.7]
#none are highly correlated
corrplot(as.matrix(cormat), type = "lower", method = "number")


####Spatial autocorrelation structure for each grid####
#Run Function_grid_correlation_structure.R
decay_df <- grid_correlation_structure(grid_vector = c(unique(Hdat_filled$grid)), 
                                   data = Hdat_filled, 
                                   formula = "SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height", 
                                   k_specified = 4)




####Build correlation matrix, R####
grid_params <- decay_df
#For some grids, the negative exponential curve could not be fitted succesfully
#Give these grids the average decay constant from all the grids
mean_b   <- mean(grid_params$b,na.rm = TRUE)
mean_lag <- round(mean(grid_params$first_nonsig_lag, na.rm = TRUE))
mean_range_dist <- round(mean(grid_params$range_dist, na.rm = TRUE))

grid_params <- grid_params %>%
  mutate(b = ifelse(is.na(b),mean_b,b),
         first_nonsig_lag = ifelse(is.na(first_nonsig_lag), mean_lag, first_nonsig_lag), 
         range_dist = ifelse(is.na(range_dist), mean_range_dist, range_dist))

Rtest <- build_R(data = Hdat_filled, 
                grid_params = grid_params, 
                cutoff = "range_dist")

#map R as a sanity check
df <- melt(Rtest)
colnames(df) <- c("row", "col", "correlation")

pdf("Figures/R_plot.pdf")
ggplot(df, aes(x = col, y = row, fill = correlation)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() 
dev.off()



####Run GEE####
Hdat_filled$grid <- as.factor(Hdat_filled$grid)
gee2 <- gee::gee(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
            family = gaussian, data = Hdat_filled,
            id = grid,
            corstr = "fixed",
            R= Rtest[1:160, 1:160], #needs to be the same dimension as one group
            scale.fix = FALSE) #start 12:15

summary(gee2)

#Get p values
coefs <- summary(gee2)$coefficients
p_values <- 2 * pnorm(abs(coefs[, "Robust z"]), lower.tail = FALSE)
cbind(coefs, p_value = round(p_values, 4))

Hdat_final <- Hdat_filled |> 
  arrange(y_coord) |> 
  group_by(grid) |> 
  mutate(order_in_grid = c(1:160))
  

zcor <- fixed2Zcor(cor.fixed = Rtest, 
                   id = Hdat_final$grid, 
                   waves = Hdat_final$order_in_grid)

test2<- geepack::geeglm(SES ~ rock_cover + northness + soil_moisture_adj_campaign2 + mean_soil_depth + slope_height,
                        family = "gaussian", 
                        data = Hdat_final, 
                        id = c(Hdat_final$grid), 
                        zcor = zcor, 
                        corstr = "userdefined")
