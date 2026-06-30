####PCA of all microenvironmental variables####
library(tidyverse)
library(tidylog)
library(ggbiplot)

###Microenvironmental data
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

###Microclimate indices for the summer
ind <- read.csv("C:\\Users\\imke6\\Documents\\PhD 2025\\Ch2 niche modelling\\All_data\\Clean_data\\microclimate\\Pekka_cleaned_data\\TOMST data\\transfer_402058_files_56a24364\\microclimate_indices.csv") |> 
  mutate(area = case_when(area == "BOK" ~ "BK", 
                          area == "GOL" ~ "GG", 
                          area == "WIT" ~ "WH"), 
         grid = str_split_i(grid, "_", 2), 
         cell = toupper(str_split_i(site, "_", 3)), 
         Cell_ID = paste0(area, "_G", grid, "_", cell)) |> 
  dplyr::select(!c(site, impprop)) |> 
  rename(site = "area") |> 
  #pick variables we are interested in
  filter(stat %in% c("avg", "fdh", "gdh5"), sensor %in% c("T1", "moist"), period %in% c("annual", "summer")) |>  #T1 is the soil temperature, this is the most relaible temperature as it isn't excesively heated by solar radiation
  mutate(variable = paste(sensor, stat, period, sep = "_")) |> 
  dplyr::select(!c(sensor, stat, period)) |> 
  pivot_wider(names_from = variable, values_from = value) |> 
  arrange(site, grid)

ind %>% 
  group_by(Cell_ID) %>% 
  filter(n() > 1) %>% 
  ungroup() #no duplicates here


####Remote sensing derived variables
rms <- read.csv("All_data/clean_data/Environmental data/Zonal_stats_all.csv") |> 
  rename(Cell_ID = CELL_ID)

dups <- rms %>% 
  group_by(Cell_ID) %>% 
  filter(n() > 1) %>% 
  ungroup() #There are some duplicates. For now we will just remove them

rms <- rms |> 
  filter(!(Cell_ID %in% c(dups$Cell_ID)))


###Now we merge them all together
all_env <- env |> 
  full_join(ind, by = "Cell_ID") |> 
  full_join(rms, by = "Cell_ID")  |> 
  dplyr::select(c(Cell_ID, site.x, elevation, rock_cover:slope_height, T1_avg_annual:moist_avg_summer, STD)) |> 
  column_to_rownames(var = "Cell_ID") |> 
  mutate(across(c(rock_cover: STD), ~ as.numeric(scale(.x))))
  
####Do the principal component analysis####
all_env_subs <- all_env |> 
  drop_na()


all_env_pca <- princomp(all_env_subs[, c(3:ncol(all_env_subs))], scores = T) #not including elevation
summary(all_env_pca) #proportion of variance is the variance explained by the PC
all_env_pca$scores #
all_env_pca$loadings #How much each var contributed to building the PC
all_env_pca$scale #scaling applied to each variable
all_env_pca$center #means

#make biplot
biplot(all_env_pca, choices = c("Comp.1", "Comp.2"))


###GGplot biplot
env_pca <- ggbiplot(all_env_pca, choices = c(1,2), 
                      varname.size = 4, varname.color = "black", 
                      groups = c(all_env_subs$site.x)) +
  geom_point(aes(color = all_env_subs$site.x), alpha = 0.8)+
  scale_color_manual(values = c("blue", "red", "green"))+
  theme_classic() 
env_pca$layers <- c(env_pca$layers, env_pca$layers[[2]], env_pca$layers[[3]]) #move the arrows to plot in the foreground
#ggsave(plot = trait_pca, filename = "trait.pca.png", path = "Figures") 


###correlation
cormat <- cor(all_env_subs[, c(3:ncol(all_env_subs))])
corrplot(cormat, method = "number", type = "lower", 
         tl.cex = 0.8,    # shrink variable name text
         tl.srt = 45)




