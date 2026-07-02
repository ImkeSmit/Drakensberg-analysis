####Script to generate figures for posters and publications####
library(ggplot2)
library(ggridges)
library(tidyverse)
library(ggh4x)
library(ggpubr)
library(scales)
library(glmmTMB)

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
imp_H_only <- data.frame(var = c("elevation", "rock cover", "northness","soil moisture","soil depth" ,"slope height" ), 
                           imp_H = imp_H, row.names = NULL) |> 
                          mutate(imp_H = case_when(imp_H<0 ~ 0, .default = imp_H)) |> 
  ggplot(aes(x = var, y = imp_H)) +
  geom_bar(stat = "identity")+
  labs(x = "", y = "Variable importance", title = "SES of plant height")+
  scale_y_continuous(labels = label_number(accuracy = 0.01), n.breaks = 3)+
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 16),
        title = element_text(size = 16), 
        panel.grid = element_blank())


imp_SLA_only <- data.frame(var = c("elevation", "rock cover", "northness","soil moisture","soil depth" ,"slope height" ), 
                         imp_SLA = imp_SLA, row.names = NULL) |> 
                          mutate(imp_SLA = case_when(imp_SLA<0 ~ 0, .default = imp_SLA)) |>
  ggplot(aes(x = var, y = imp_SLA)) +
  geom_bar(stat = "identity")+
  labs(x = "", y = "Variable importance", title = "SES of SLA")+
  scale_y_continuous(labels = label_number(accuracy = 0.01))+
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 18), 
        axis.text.y = element_blank(), axis.text.x = element_text(size = 16),
        title = element_text(size = 16), 
        panel.grid = element_blank())

var_imp_panes <- ggarrange(imp_H_only, imp_SLA_only, ncol = 2, nrow = 1, widths = c(1,0.75))
ggsave(var_imp_panes, filename = "variable_importance_poster.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")




#####SES ~ elevation ridges####
#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/SES_RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
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
  #join to microclimate indices |> 
  full_join(micro_idw, by = "Cell_ID") |> 
  #join to remote sensing data |> 
  full_join(rms, by = "Cell_ID") |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cell_ses, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid), 
         elevation = as.factor(elevation)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)

l1 <- c("log_Height" = "SES~of~plant~height", "log_SLA" = "SES~of~SLA", "log_LDMC" = "SES~of~LDMC", "log_LA" = "SES~of~LA")
ridges_letters <- data.frame(trait = c(rep("log_Height", 3), rep("log_SLA", 3), rep("log_LDMC", 3),rep("log_LA", 3)), 
                             elevation = as.factor(c(rep(c("2000", "2500", "3000"), 4))),
                             letters = c(SES_ele_cld$log_Height$.group, 
                                         SES_ele_cld$log_SLA$.group, 
                                         SES_ele_cld$log_LDMC$.group, 
                                         SES_ele_cld$log_LA$.group), 
                             x_pos = c(4.5, 4.5, 4.5,
                                       6,6,6, 
                                       9,9,9,
                                       8.1,8.1,8.1), 
                             emmean = c(SES_ele_cld$log_Height$emmean,
                                        SES_ele_cld$log_SLA$emmean,
                                        SES_ele_cld$log_LDMC$emmean,
                                        SES_ele_cld$log_LA$emmean), 
                             emmean_SE = c(SES_ele_cld$log_Height$SE,
                                           SES_ele_cld$log_SLA$SE,
                                           SES_ele_cld$log_LDMC$SE,
                                           SES_ele_cld$log_LA$SE))

ses_ridges <- comb |>
  filter(trait %in% c("log_Height", "log_SLA", "log_LDMC", "log_LA")) |> 
  ggplot(aes(x = SES, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), scale = "free_x", 
             strip.position = "bottom")+
  labs(x = " ", y = "Elevation (m a.s.l.)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = ridges_letters, aes(x = x_pos, y = elevation, label = letters, size = 16))+
  geom_point(
    data = ridges_letters,
    aes(x = emmean, y = as.numeric(elevation)),
    size = 2, inherit.aes = FALSE) +
  geom_errorbarh(
    data = ridges_letters,
    aes(xmin = emmean - emmean_SE,
        xmax = emmean + emmean_SE,
        y = as.numeric(elevation)),
    height = 0.08,
    linewidth = 1, inherit.aes = FALSE)+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 18), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(ses_ridges, filename = "SES_elevation_poster.png", path = "Figures")


####SES ridges on subsampled data####
##subsample cells
Height_incl_cells <- read.csv("All_data/comm_assembly_results/included_cells_Height.csv", row.names = 1) |> 
  rename(Cell_ID = included_cells)

Hdat_subs <- comb2 |> 
  filter(trait %in% c("Height_cm", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "Height_cm",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord)) |> 
  inner_join(Height_incl_cells, by = c("trait", "Cell_ID"))


#also for SLA
SLA_incl_cells <- read.csv("All_data/comm_assembly_results/included_cells_SLA.csv", row.names = 1) |> 
  rename(Cell_ID = included_cells)

SLAdat_subs <- comb2 |> 
  filter(trait %in% c("SLA", NA)) |>  #also select cells which have no SES measurement. This is necessary to make the grid complete
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = as.factor(paste0(site, grid)), 
         pos = numFactor(x_coord, y_coord)) |> 
  inner_join(SLA_incl_cells, by = c("trait", "Cell_ID"))


##Bind together
all_subs <- bind_rows(Hdat_subs, SLAdat_subs)

#labels
l1 <- c("Height_cm" = "SES~of~plant~height", "SLA" = "SES~of~SLA")
ridges_letters <- data.frame(trait = c(rep("Height_cm", 3), rep("SLA", 3)), 
                             elevation = as.factor(c(rep(c("2000", "2500", "3000"), 2))),
                             letters = c("a", "b", "a", "a", "b", "a"), 
                             x_pos = c(6.5,6.5,6.5,4,4,4), 
                             emmean = c(-0.238, -0.437, -0.301, -0.8908, -0.0923, -0.7040), #type over from results table
                             emmean_SE = c(0.0527, 0.0290, 0.0279, 0.114, 0.125, 0.106))

# Interpolate density height at each group's mean
segments_df <- all_subs %>%
  filter(trait %in% c("Height_cm", "SLA")) %>%
  group_by(trait, elevation) %>%
  summarise(
    mean_val = mean(SES),
    dens_x = list(density(SES)$x),
    dens_y = list(density(SES)$y),
    .groups = "drop"
  ) %>%
  mutate(
    # Interpolate height at the mean
    mean_y_raw = mapply(function(dx, dy, mx) approx(dx, dy, xout = mx)$y,
                        dens_x, dens_y, mean_val),
    # Get the max density per ridge — same normalisation ggridges uses
    max_y_raw  = sapply(dens_y, max),
    # Scale so that max density = scale parameter (set scale = 1 in geom_density_ridges)
    mean_y_scaled = (mean_y_raw / max_y_raw) * 1,  # 0.9 gives a small safety margin
    emmean = c(-0.238, -0.437, -0.301, -0.8908, -0.0923, -0.7040),
    emmean_SE = c(0.0527, 0.0290, 0.0279, 0.114, 0.125, 0.106)
  )


ses_ridges_subsampled <- all_subs |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = SES, y = elevation)) +
  geom_density_ridges(alpha = 0.5, scale = 1) +
  geom_point(
    data = segments_df,
    aes(x = emmean, y = as.numeric(elevation)),
    size = 2, inherit.aes = FALSE) +
  geom_errorbarh(
    data = segments_df,
    aes(xmin = emmean - emmean_SE,
      xmax = emmean + emmean_SE,
      y = as.numeric(elevation)),
    height = 0.08,
    linewidth = 1, inherit.aes = FALSE) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), 
             scale = "free_x", strip.position = "bottom") +
  labs(x = " ", y = "Elevation (m a.s.l.)") +
  geom_vline(xintercept = 0, colour = "red") +
  geom_text(data = ridges_letters, 
            aes(x = x_pos, y = elevation, label = letters, size = 18)) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title = element_text(size = 20), 
        axis.text = element_text(size = 18), 
        strip.text = element_text(size = 20), 
        strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid = element_blank())
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
  scale_x_continuous(breaks = c(0,25,50))+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14), strip.text = element_text(size = 18), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(trait_ridges, filename = "trait_elevation_poster.png", path = "Figures")


###Subset for the cells used in analysis
SLA_incl_cells <- read.csv("All_data/comm_assembly_results/included_cells_SLA.csv", row.names = 1) |> 
  rename(Cell_ID = included_cells)
H_incl_cells <- read.csv("All_data/comm_assembly_results/included_cells_Height.csv", row.names = 1) |> 
  rename(Cell_ID = included_cells)

incl_cells <- bind_rows(SLA_incl_cells, H_incl_cells)

cwm_subs <- cwm |> 
  rename(Cell_ID = cellref) |> 
  inner_join (incl_cells, by = c("trait", "Cell_ID")) #subset


##Graph
l1 <- c("Height_cm" = "CWM~Plant~height~(cm)", "SLA" = "CWM~SLA~(mm^2/mg)")
ridges_letters <- data.frame(trait = c(rep("Height_cm", 3), rep("SLA", 3)), 
                             elevation = as.factor(c(rep(c("2000", "2500", "3000"), 2))),
                             letters = c("a", "b", "c", "a", "b", "c"), 
                             x_pos = c(60,60,60,60,60,60), 
                             emmean = c(19.9, 11.1, 14.4, 18.2, 38.4, 15.5), #type over from results table
                             emmean_SE = c(0.567, 0.619, 0.602, 0.849, 1.060, 0.908))
                            

trait_ridges_subsampled <- cwm_subs |>
  filter(trait %in% c("Height_cm", "SLA")) |> 
  ggplot(aes(x = cwm_value, y = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  geom_point(
    data = ridges_letters,
    aes(x = emmean, y = as.numeric(elevation)),
    size = 2, inherit.aes = FALSE) +
  geom_errorbarh(
    data = ridges_letters,
    aes(xmin = emmean - emmean_SE,
        xmax = emmean + emmean_SE,
        y = as.numeric(elevation)),
    height = 0.08,
    linewidth = 1, inherit.aes = FALSE) +
  facet_wrap(~trait, labeller = as_labeller(l1, default = label_parsed), scale = "free_x", strip.position = "bottom")+
  geom_text(data = ridges_letters, aes(x = x_pos, y = elevation, label = letters, size = 18))+
  labs(x = "", y = "Elevation (m a.s.l.)") +
  scale_x_continuous(breaks = c(0,25,50))+
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 20), 
        axis.text = element_text(size = 18), strip.text = element_text(size = 20), 
        strip.background = element_blank(),
        strip.placement = "outside", panel.grid = element_blank()) 
ggsave(trait_ridges_subsampled, filename = "trait_elevation_poster_subsampled.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")



####Scatterplots of SES Height ~microenvironmental variables####
SLA_incl_cells <- read.csv("All_data/comm_assembly_results/included_cells_SLA.csv", row.names = 1) |> 
  rename(Cell_ID = included_cells)
H_incl_cells <- read.csv("All_data/comm_assembly_results/included_cells_Height.csv", row.names = 1) |> 
  rename(Cell_ID = included_cells)

incl_cells <- bind_rows(SLA_incl_cells, H_incl_cells)

plotdat <- comb2 |> 
  inner_join(incl_cells, by = c("trait", "Cell_ID"))


###Generate predictions for lines of best fit###
Hdat_subs <- plotdat |> filter(trait == "Height_cm")
tmod1<- glmmTMB(SES ~ elevation + zrock_cover +  znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height, 
                family = t_family(link = "identity"), 
                data = Hdat_subs) 

rock_dat <- Hdat_subs |> 
  dplyr::select(elevation , zrock_cover, znorthness, zsoil_moist, zsoil_depth, zslope_height) |> 
  mutate(elevation = "2000", 
         znorthness = mean(znorthness), 
         zsoil_moist = mean(zsoil_moist, na.rm = T), 
         zsoil_depth = mean(zsoil_depth), 
         zslope_height = mean(zslope_height))
rock_predictions <- predict(tmod1, rock_dat, type = "response")
rock_pred_dat <- data.frame(microenv_var = "zrock_cover", predictions = rock_predictions, value = rock_dat$zrock_cover)


moist_dat <- Hdat_subs |> 
  dplyr::select(elevation , zrock_cover, znorthness, zsoil_moist, zsoil_depth, zslope_height) |> 
  drop_na() |> 
  mutate(elevation = "2000", 
         znorthness = mean(znorthness), 
         zrock_cover = mean(zrock_cover), 
         zsoil_depth = mean(zsoil_depth), 
         zslope_height = mean(zslope_height))
moist_predictions <- predict(tmod1, moist_dat, type = "response")
moist_pred_dat <- data.frame(microenv_var = "zsoil_moist", predictions = moist_predictions, value = moist_dat$zsoil_moist)

prediction_lines <- rbind(rock_pred_dat, moist_pred_dat)

l1 <- c("zrock_cover" = "Rock~~cover~(standardised)", "zsoil_moist" = "Soil~moisture~(standardised)")

Hplots <- plotdat |> 
  filter(trait == "Height_cm") |> 
  dplyr::select(Cell_ID, SES, zrock_cover, zsoil_moist) |> 
  pivot_longer(cols = c(zrock_cover, zsoil_moist), names_to = "microenv_var", values_to = "value") |> 
  ggplot(aes(x = value, y = SES)) +
  geom_point()+
  geom_line(data = prediction_lines, aes(x = value, y = predictions), color = "red", linewidth = 1) +
  facet_wrap(~ microenv_var, scale = "free_x", labeller = as_labeller(l1, default = label_parsed), 
             strip.position = "bottom") +
  labs(x = " ", y = "SES of Height") +
  theme_bw() +
  theme( 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14), 
        strip.text = element_text(size = 16), 
        strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid = element_blank())
ggsave(Hplots, filename = "SES_Height_microenv_subsampled.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")




####Scatterplots of SES SLA ~microenvironmental variables####

###Generate predictions for lines of best fit###
SLAdat_subs <- plotdat |> filter(trait == "SLA")
tmod2<- glmmTMB(SES ~ elevation + zrock_cover +  znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height, 
                family = t_family(link = "identity"), 
                data = SLAdat_subs) 

north_dat <- SLAdat_subs |> 
  dplyr::select(elevation , zrock_cover, znorthness, zsoil_moist, zsoil_depth, zslope_height) |> 
  drop_na() |> 
  mutate(elevation = "2000", 
         zrock_cover = mean(zrock_cover), 
         zsoil_moist = mean(zsoil_moist), 
         zsoil_depth = mean(zsoil_depth), 
         zslope_height = mean(zslope_height))
north_predictions <- predict(tmod2, north_dat, type = "response")
north_pred_dat <- data.frame(microenv_var = "znorthness", predictions = north_predictions, value = north_dat$znorthness)


moist_dat <- SLAdat_subs |> 
  dplyr::select(elevation , zrock_cover, znorthness, zsoil_moist, zsoil_depth, zslope_height) |> 
  drop_na() |> 
  mutate(elevation = "2000", 
         znorthness = mean(znorthness), 
         zrock_cover = mean(zrock_cover), 
         zsoil_depth = mean(zsoil_depth), 
         zslope_height = mean(zslope_height))
moist_predictions <- predict(tmod2, moist_dat, type = "response")
moist_pred_dat <- data.frame(microenv_var = "zsoil_moist", predictions = moist_predictions, value = moist_dat$zsoil_moist)

prediction_lines <- rbind(north_pred_dat, moist_pred_dat)

l1 <- c("znorthness" = "Northness~(standardised)", "zsoil_moist" = "Soil~moisture~(standardised)")

SLAplots <- plotdat |> 
  filter(trait == "SLA") |> 
  dplyr::select(Cell_ID, SES, znorthness, zsoil_moist) |> 
  pivot_longer(cols = c(znorthness, zsoil_moist), names_to = "microenv_var", values_to = "value") |> 
  ggplot(aes(x = value, y = SES)) +
  geom_point()+
  geom_line(data = prediction_lines, aes(x = value, y = predictions), color = "red", linewidth = 1) +
  facet_wrap(~ microenv_var, scale = "free_x", labeller = as_labeller(l1, default = label_parsed), 
             strip.position = "bottom") +
  labs(x = " ", y = "SES of SLA") +
  theme_bw() +
  theme( 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    strip.text = element_text(size = 16), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    panel.grid = element_blank())
ggsave(SLAplots, filename = "SES_SLA_microenv_subsampled.png", path = "Figures", 
       width = 2300, height = 1400, units = "px")
