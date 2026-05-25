####Script to generate figures for posters and publications####
library(ggplot2)
library(ggridges)
library(tidyverse)
library(ggh4x)
library(ggpubr)
library(scales)

####Hypothesis figure SES~elevation####
#generate fake data
SESlow <- dnorm(200, mean = 1, sd = 0.2)
SESmid <- dnorm(200, mean = 0, sd = 0.2)
SEShigh <- dnorm(200, mean = -1, sd = 0.2)
SES_all <- c(SESlow, SESmid, SEShigh)

hp_dat <- tibble(elevation = c(rep("2000", 200), rep("2500", 200), rep("3000", 200)), 
                 SES = SES_all)

hypothesis_ses_ridges <- hp_dat |>
  ggplot(aes(x = SES, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  labs(x = "SES", y = "Elevation (m a.s.l.)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 20), 
        axis.text = element_text(size = 18), 
        panel.grid = element_blank()) 
ggsave(hypothesis_ses_ridges, filename = "hypothesis_SES_elevation_poster.png", path = "Figures", 
       width = 1100, height = 1100, units = "px")



####Hypothesis grid maps####
##Function to retrieve data of one grid
one_grid_raster <- function(data, variable, g) {
  
  library(sp)
  library(raster)
  
  # filter data for this grid
  g_dat <- data %>% filter(grid == g)
  
  
  #remove extra columns
  g_dat <- g_dat[, which(colnames(g_dat) %in% c("row", "ncolumn", variable))]
  
  # Create spatial grid object
  x_range <- 1:20
  y_range <- 1:8
  
  grid_obj <- expand.grid(x = x_range, y = y_range)
  coordinates(grid_obj) <- ~x + y
  gridded(grid_obj) <- TRUE
  r <- raster(grid_obj) 
  
  #Check if there are missing cells in data
  possible_cells <- paste(grid_obj$x, grid_obj$y, sep = "_")
  observed_cells <- paste(g_dat$row, g_dat$ncolumn, sep = "_")
  missing <- which(is.na(match(possible_cells, observed_cells)))
  
  if(length(missing) > 0) { #create a filler
    
    filler <- data.frame(row = grid_obj$x[missing], ncolumn = grid_obj$y[missing])
    filler[[variable]] <- NA #the variable is set to NA
    
    g_dat <- rbind(g_dat, filler) }
  
  # need values sorted to match the raster cell order:
  g_dat <- g_dat[order(g_dat$ncolumn, g_dat$row), ]
  g_dat <- as.data.frame(g_dat)
  
  return(g_dat)
  #end function
}


###Environmental grid####
##import environmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(elevation = case_when(grepl("BK", Cell_ID) == T ~ "3000", #add elevation variable
                               grepl("WH", Cell_ID) == T ~ "2500",
                               grepl("GG", Cell_ID) == T ~ "2000",.default = NA)) |> 
  dplyr::select(-c("site", "grid")) |> 
  separate_wider_delim(Cell_ID, delim = "_", names = c("site", "grid", "cell"), cols_remove = F) |> 
  mutate(grid = paste0(site, grid), 
         column = substr(cell, 1,1), 
         row = substr(cell, 2,3),
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth), 
         soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2)) |> 
  dplyr::select(!c(soil_temperature_adj_campaign2, soil_moisture_adj_campaign1, soil_temperature_adj_campaign1,
                   aeolian_process, fluvial_process, slope_process,
                   geology1, geology2, geology3,
                   geology4, geology5, mesotopo, aspect, veg_max_height))


WHG7_rock_r <- one_grid_raster(data = env, variable = "rock_cover", g = "WHG7")

grid_rock_WHG7 <- ggplot(WHG7_rock_r, aes(y = row, x = ncolumn, fill = rock_cover))+
  geom_tile(color = "black", linewidth = 0.5)+
  scale_fill_gradient(low = "white", high = "darkorange2", na.value = "pink")+
  coord_fixed(ratio = 1)+
  theme_void()+
  labs(x = " ", y = " ", fill = " ", title = "Rock cover (%)")  +
  theme(legend.text = element_text(size = 17), legend.title = element_text(size = 20), 
        title = element_text(size = 18))
ggsave("Figures//rock_cover_WHG7.png", grid_rock_WHG7, device = png, height = 900, units = "px")


####SES grid####
cell_ses_height <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  separate_wider_delim(cellref, delim = "_", names = c("site", "grid", "cell"), cols_remove = F) |> 
  mutate(grid = paste0(site, grid), 
         column = substr(cell, 1,1), 
         row = substr(cell, 2,3),
         ncolumn = match(column, LETTERS[1:8])) |> 
  filter(trait == "Height_cm") #map SES of height


WHG7_SES_r <- one_grid_raster(data = cell_ses_height, variable = "SES", g = "WHG7")

grid_SES_WHG7 <- ggplot(WHG7_SES_r, aes(y = row, x = ncolumn, fill = SES))+
  geom_tile(color = "black", linewidth = 0.5)+
  scale_fill_gradient2(low = "blue3", mid = "white", high = "darkorange2", na.value = "black")+
  coord_fixed(ratio = 1)+
  theme_void()+
  labs(x = " ", y = " ",fill = " ", title = "SES")  +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20), 
        title = element_text(size = 18))
ggsave("Figures//SES_WHG7.png", grid_SES_WHG7, device = png, height = 900, units = "px")



####Variable importance####
#import results from subsampled glmmTMB models
height_mod <- read.csv("All_data/comm_assembly_results/glmmTMB_subsampled_SES_height_env_model_results.csv")
colnames(height_mod) = c("variable", "estimate", "std_eror","z_value", "p_value", "var_imp")
SLA_mod <- read.csv("All_data/comm_assembly_results/glmmTMB_subsampled_SES_SLA_env_model_results.csv")
colnames(SLA_mod) = c("variable", "estimate", "std_eror","z_value", "p_value", "var_imp")

imp_H <- height_mod$var_imp[3:8]
imp_SLA <- SLA_mod$var_imp[3:8]

#make the figure by combining two figures:
imp_SLA_only <- data.frame(var = c("elevation", "rock cover", "northness","soil moisture","soil depth" ,"slope height" ), 
                         imp_SLA = imp_SLA, row.names = NULL) |> 
  ggplot(aes(x = var, y = imp_SLA)) +
  geom_bar(stat = "identity")+
  labs(x = "", y = "Variable importance", title = "SES of SLA")+
  scale_y_continuous(labels = label_number(accuracy = 0.01))+
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 18), 
        axis.text.y = element_blank(), axis.text.x = element_text(size = 16),
        title = element_text(size = 20), 
        panel.grid = element_blank())

var_imp_panes <- ggarrange(imp_H_only, imp_SLA_only, ncol = 2, nrow = 1, widths = c(1.5,1))
ggsave(var_imp_panes, filename = "variable_importance_poster.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")




#####SES ~ elevation ridges####
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

l1 <- c("Height_cm" = "SES~of~plant~height", "SLA" = "SES~of~SLA")
ridges_letters <- data.frame(trait = c(rep("Height_cm", 3), rep("SLA", 3)), 
                             elevation = as.factor(c(rep(c("2000", "2500", "3000"), 2))),
                             letters = c("a", "b", "a", "a", "b", "a"), 
                             x_pos = c(7.5,7.5,7.5,4.5,4.5,4.5))

ses_ridges <- comb2 |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = SES, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), scale = "free_x", 
             strip.position = "bottom")+
  labs(x = " ", y = "Elevation (m a.s.l.)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = ridges_letters, aes(x = x_pos, y = elevation, label = letters, size = 16))+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 18), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(ses_ridges, filename = "SES_elevation_poster.png", path = "Figures")


####SES ridges on subsampled data####
##subsample cells
Hdat <- comb2 |> 
  filter(trait %in% c("Height_cm", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))

Hdat_subs<- select_independent_cells(Hdat, grid_var = "grid", x = "x_coord", y = "y_coord", value_col = "SES",
                                     max_search_radius = 2)

#also for SLA
SLAdat <- comb2 |> 
  filter(trait %in% c("SLA", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord))

SLAdat_subs<- select_independent_cells(SLAdat, grid_var = "grid", x = "x_coord", y = "y_coord", value_col = "SES",
                                       max_search_radius = 2)

##Bind together
all_subs <- bind_rows(Hdat_subs, SLAdat_subs)

#labels
l1 <- c("Height_cm" = "SES~of~plant~height", "SLA" = "SES~of~SLA")
ridges_letters <- data.frame(trait = c(rep("Height_cm", 3), rep("SLA", 3)), 
                             elevation = as.factor(c(rep(c("2000", "2500", "3000"), 2))),
                             letters = c("a", "b", "a", "a", "b", "a"), 
                             x_pos = c(6.5,6.5,6.5,4,4,4))


ses_ridges_subsampled <- all_subs |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = SES, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), scale = "free_x", 
             strip.position = "bottom")+
  labs(x = " ", y = "Elevation (m a.s.l.)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = ridges_letters, aes(x = x_pos, y = elevation, label = letters, size = 18))+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 20), 
        axis.text = element_text(size = 18), strip.text = element_text(size = 20), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(ses_ridges_subsampled, filename = "SES_elevation_subsampled.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")



####Variable importance bar graph####
l3 <- c("imp_H" = "SES of plant height", "imp_SLA" = "SES of SLA")

var_imp <- data.frame(var = c("elevation", "rock cover", "northness","soil moisture","soil depth" ,"slope height" ), 
                      imp_H = importance, imp_SLA = importance_SLA, row.names = NULL) |> 
  arrange(imp_H) |> 
  pivot_longer(!var, names_to = "trait", values_to = "var_imp") |>
  arrange(trait, var_imp) |> 
  ggplot(aes(x = var, y = var_imp)) +
  geom_bar(stat = "identity")+
  facet_wrap(~trait, strip.position = "top", labeller = as_labeller(l3), scales = "free")+
  labs(x = "", y = "Variable importance")+
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 18), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank())
ggsave(var_imp, filename = "variable_importance_poster.png", path = "Figures", width = 1800, units = "px")



####CWm ridges####
#import trait and abundance data
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1)


#compute CWM of each trait for each cell
library(FD)
cwm <- functcomp(x = mean_traits, a = as.matrix(abun_matrix))
cwm <- cwm |>
  rownames_to_column(var = "cellref") |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  pivot_longer(cols = c("Height_cm", "LDMC", "Leaf_area_mm2", "SLA", "Thickness_mm"), names_to = "trait", values_to = "cwm_value")
cwm$elevation <- as.factor(cwm$elevation) 

##Graph
l1 <- c("Height_cm" = "Plant~height~(cm)", "SLA" = "SLA~(mm^2/mg)")

trait_ridges <- cwm |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = cwm_value, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), scale = "free_x", strip.position = "bottom")+
  labs(x = "", y = "Elevation (m a.s.l.)") +
  #geom_text(data = ridges_letters, aes(x = x_pos, y = elevation, label = letters, size = 16))+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 18), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(trait_ridges, filename = "trait_elevation_poster.png", path = "Figures")


###Subset for the cells used in analysis
included <- all_subs |> 
  group_by(trait) |> 
  distinct(Cell_ID) |> 
  rename(cellref = Cell_ID)

cwm_subs <- cwm |> 
  inner_join(included, by = c("trait", "cellref"))


##Graph
l1 <- c("Height_cm" = "Plant~height~(cm)", "SLA" = "SLA~(mm^2/mg)")

trait_ridges_subsampled <- cwm_subs |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = cwm_value, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), scale = "free_x", strip.position = "bottom")+
  labs(x = "", y = "Elevation (m a.s.l.)") +
  #geom_text(data = ridges_letters, aes(x = x_pos, y = elevation, label = letters, size = 16))+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 20), 
        axis.text = element_text(size = 18), strip.text = element_text(size = 20), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(trait_ridges_subsampled, filename = "trait_elevation_poster_subsampled.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")
