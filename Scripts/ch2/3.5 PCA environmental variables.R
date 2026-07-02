####PCA of all microenvironmental variables####
library(tidyverse)
library(tidylog)
library(ggbiplot)

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
  mutate(ncolumn = match(column, LETTERS[1:8])) |> 
  rename(x_coord = ncolumn, 
         y_coord = row)


#check for duplicates
dups <- comb %>% 
  group_by(Cell_ID) %>% 
  filter(n() > 1) %>% 
  ungroup() #no duplicates


  
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




