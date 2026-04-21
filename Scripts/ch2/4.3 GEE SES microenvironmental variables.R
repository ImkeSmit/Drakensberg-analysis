###Script to analyse relationship between SES and microenvironmental variables
library(spdep)
library(spatialreg)
library(spind)
library(remotes)
library(dataDownloader)
library(osfr)
library(tidyverse)
library(tidylog)

#imort SES data
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
          mutate(Cell_ID = paste0(site, "_G", str_sub(cellref, 3, 3), "_", column, row)) |> 
          select(-cellref)
cell_ses$elevation <- as.factor(cell_ses$elevation) 

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth))


##Combine SES and environmental data
comb <- env |> 
  select(-c(site, grid, column,row)) |> 
  inner_join(cell_ses, by = "Cell_ID") |> 
##Now we need to change the coordinates to reflect the spatial structure of the whole dataset
##Coordinates within a site differ by 80
##Coordinates between sites differ by 1000
  mutate(y_new = case_when(site == "GG" ~ ncolumn, 
                           site == "WH" ~ ncolumn+1000, 
                           site == "BK" ~ ncolumn+2000, .default = NA)) |> 
  group_by(site) |> 
  mutate(y_coord = case_when(grepl("1", grid) ~ y_new, 
                           grepl("2", grid) ~ y_new+80, 
                           grepl("3", grid) ~ y_new+160,
                           grepl("4", grid) ~ y_new+240, 
                           grepl("5", grid) ~ y_new+320, 
                           grepl("6", grid) ~ y_new+400, 
                           grepl("7", grid) ~ y_new+480, 
                           grepl("8", grid) ~ y_new+560, .default = NA)) |> 
  mutate(x_coord = row)




comb_H <- cell_ses |> 
  filter(trait == "Height_cm", grid == "WH1") |> #right now we can only model one grid at a time
  inner_join(env, by = "Cell_ID")

coords <- cbind(comb_H$row.x, comb_H$ncolumn)

gee1 <- GEE(SES ~ rock_cover + mean_soil_depth, 
            family = "gaussian", data = comb_H,
            coord = coords, corstr = "fixed", scale.fix = FALSE)
summary(gee1, printAutoCorPars = TRUE)

plot(gee1)
