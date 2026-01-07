###Cluster cells according to environmental similarity###
library(tidyverse)
library(tidylog)
library(ggplot2)
library(ggbiplot)
library(multcomp)
library(multcompView)
library(rcompanion)

#import environmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  mutate(mean_soil_depth = as.numeric(mean_soil_depth))

sdepth <- env |> 
  dplyr::select(Cell_ID, mean_soil_depth, northness, vascular_cover) |> 
  filter(!is.na(mean_soil_depth)) |> 
  column_to_rownames("Cell_ID")
#should add elevation in here

sd_dist <- dist(sdepth, method = "euclidean")

sd_hclust <- hclust(sd_dist, method = "ward.D")

#cut the tree to have 5 clusters
cut <- cutree(sd_hclust, k = 5)
plot(sd_hclust)
rect.hclust(sd_hclust, k = 5, border  = "red")

#we could also divide the cells into categories of deviation from the mean
#clustering only really makes sense if we use multiple variables

#make a dataframe of the cell, the sdepth cluster and the SES

cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA)) |> 
  mutate(site = substr(cellref, 1,2),
         grid = substr(cellref, 3,3),
         column = substr(cellref, 4,4), 
         row = substr(cellref, 5,6),
         ncolumn = match(column, LETTERS[1:8])) |> 
  mutate(row = as.numeric(row), 
         Cell_ID = paste0(site, "_G", grid, "_", column, row))

clusters <- tibble(Cell_ID = names(cut), 
                   cluster = as.factor(cut)) |> 
  left_join(cell_ses, by = "Cell_ID") |> 
  left_join(env, by = "Cell_ID")

ggplot(clusters, aes(x = cluster, y = mean_soil_depth)) +
  geom_boxplot()

ggplot(clusters, aes(x = cluster, y = northness)) +
  geom_boxplot()

ggplot(clusters, aes(x = cluster, y = vascular_cover)) +
  geom_boxplot()

ggplot(clusters, aes(x = cluster, y = SES)) +
  geom_boxplot()

#Are these clusters MEANINGFUL? 
#May be better to create groups based on deviation from the mean for each variable separately. That will be easier to interpret than these multivariate clusters


###PCA of environmental variables
env_subset <- env |> 
  select(Cell_ID, vascular_cover, rock_cover, northness, soil_moisture_adj_campaign2, 
         soil_temperature_adj_campaign1, veg_median_height, mean_soil_depth, slope_height) |> 
  mutate(soil_moisture_adj_campaign2 = as.numeric(soil_moisture_adj_campaign2), 
         soil_temperature_adj_campaign1 = as.numeric(soil_temperature_adj_campaign1)) |> 
  mutate(elevation = case_when(grepl("BK", Cell_ID) == T ~ 3000, #add elevation variable
                               grepl("WH", Cell_ID) == T ~ 2500,
                               grepl("GG", Cell_ID) == T ~ 2000,.default = NA)) |> 
  column_to_rownames(var = "Cell_ID") |> 
  filter(vascular_cover < 110) |> 
  drop_na() 

env_pca <- prcomp(env_subset, scale = T)
biplot(env_pca)
env_pca$loadings

#perform kmeans clustering
env_scaled <- scale(env_subset)
env_kmeans <- kmeans(env_scaled, centers = 4, iter.max = 20, nstart = 1)

#visualise clusters
library(factoextra)
fviz_cluster(env_kmeans, data = env_scaled, show_labels = T)

kmeans_clusters <- env_subset |> 
  rownames_to_column(var = "Cell_ID") |> 
  inner_join(tibble(k_cluster = as.factor(env_kmeans$cluster), 
                          Cell_ID = names(env_kmeans$cluster)), by = "Cell_ID") |> 
  left_join(cell_ses, by = "Cell_ID")

##another way of visualising
env_pca2 <- prcomp(kmeans_clusters[, 2:10], scale = T)

ggbiplot(env_pca2,
         groups = kmeans_clusters$k_cluster,
         labels = kmeans_clusters$Cell_ID,
         labels.size = 2,
         var.factor = 1.4,
         ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1,
         circle = F,
         varname.size = 4,
         varname.color = "black") +
  labs(fill = "Region", color = "Region") +
  theme(legend.direction = 'horizontal', legend.position = 'top')



#visualise variation in SES of height accross clusters
kmeans_clusters |> 
  filter(trait == "Height_cm") |> 
ggplot(aes(x = k_cluster, y = SES)) +
  geom_boxplot()

#test for differences in SES between clusters
#HEIGHT
#test if median ses of height differ between elevation groups
modeldat <- kmeans_clusters[which(kmeans_clusters$trait == "Height_cm"), ]
kr_height <- kruskal.test(SES ~ k_cluster, data = modeldat) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_height <- pairwise.wilcox.test(modeldat$SES,  modeldat$k_cluster, p.adjust.method = "bonferroni")
height_cld <- multcompLetters(fullPTable(wx_height$p.value))
#only cluster 1 and 4 do not differ


#SLA
kmeans_clusters |> 
  filter(trait == "SLA") |> 
  ggplot(aes(x = k_cluster, y = SES)) +
  geom_boxplot()

#test if median ses of SLA differ between elevation groups
modeldat <- kmeans_clusters[which(kmeans_clusters$trait == "SLA"), ]
kr_SLA <- kruskal.test(SES ~ k_cluster, data = modeldat) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_SLA <- pairwise.wilcox.test(modeldat$SES,  modeldat$k_cluster, p.adjust.method = "bonferroni")
SLA_cld <- multcompLetters(fullPTable(wx_SLA$p.value))
#cluster 1 is lower than the other clusters
