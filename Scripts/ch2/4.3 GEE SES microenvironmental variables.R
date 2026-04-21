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
  mutate(x_coord = row, 
         rock_cover = as.numeric(rock_cover), 
         mean_soil_depth = as.numeric(mean_soil_depth)) |> 
  ungroup() 



###model for SES of height
#isolate trait
Hdat <- comb |> filter(trait == "Height_cm") |> 
                 #drop_na(SES, rock_cover, mean_soil_depth) |> #remove rows with na's
                select(Cell_ID, x_coord, y_coord, SES, rock_cover, mean_soil_depth)

###We have to create the missing coordinates (between the girds and sites)
y_sequence <- unique(Hdat$y_coord)[order(unique(Hdat$y_coord))]
complete_y <- c(y_sequence[1]:y_sequence[length(y_sequence)])
missing_ycoord <- setdiff(complete_y, y_sequence)
missing_ycoord <- rep(missing_ycoord, each = 8) #repeat each ycoord 8 times, because there are 8 possible x coords
missing_xcoord <- rep(c(1:8), length(unique(missing_ycoord)))
#Create table to bind to Hdat
addcoords <- tibble(Cell_ID = NA, x_coord = missing_xcoord, y_coord = missing_ycoord, 
                    SES = NA, rock_cover = NA, mean_soil_depth = NA)
Hdat <- bind_rows(Hdat, addcoords) |> #add missing coordinates
  arrange(y_coord) 

Hdat <- Hdat |> filter(!is.na(SES))

coords <- cbind(Hdat$x_coord, Hdat$y_coord)


#Let's specify a cutoff distance above which sptaial autocorrelation is not a problem anymore
dist_matrix <- as.matrix(dist(coords))

# Set your cutoff distance (in same units as your coordinates)
cutoff <- 80

# Binary: 1 if within cutoff, 0 if beyond
weight_matrix <- ifelse(dist_matrix <= cutoff, 1, 0)
diag(weight_matrix) <- 0  # no self-correlation

gee2 <- GEE(SES ~ rock_cover + mean_soil_depth,
            family = "gaussian", data = Hdat,
            coord = weight_matrix,
            corstr = "fixed",
            scale.fix = FALSE)


gee1 <- spind::GEE(SES ~ mean_soil_depth, 
            family = "gaussian", data = Hdat,
            coord = coords, 
            corstr = "exchangeable", #can be fixed, exchangeable or quadratic
            #fixed means the same autocorrelation structure is applied to the whole dataset
            cluster = 3,
            scale.fix = FALSE)
summary(gee1, printAutoCorPars = TRUE)

plot(gee1)
