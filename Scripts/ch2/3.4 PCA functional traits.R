###PCA OF ALL TRAITS####

####Import community and trait data####
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1) |> 
  mutate(zheight = c(scale(Height_cm)), 
        zLDMC = c(scale(LDMC)), 
         zLA = c(scale(Leaf_area_mm2)), 
         zSLA = c(scale(SLA)), 
         zthickness = c(scale(Thickness_mm))) |> 
  dplyr::select(contains(c("z"))) |> 
  drop_na()



####Do the principal component analysis####
all_FT_pca <- princomp(mean_traits, scores = T)
summary(all_FT_pca) #proportion of variance is the variance explained by the PC
all_FT_pca$scores #
all_FT_pca$loadings #How much each var contributed to building the PC
all_FT_pca$scale #scaling applied to each variable
all_FT_pca$center #means

#make biplot
biplot(all_FT_pca, choices = c("Comp.1", "Comp.2"))
plot(all_FT_pca$scores[, 1], all_FT_pca$scores[, 2])



###GGplot biplot
trait_pca <- ggbiplot(all_FT_pca, choices = c(1,2), 
                      varname.size = 4, varname.color = "black") +
  geom_point(colour = "azure3", alpha = 0.8)+
  theme_classic() 
trait_pca$layers <- c(trait_pca$layers, trait_pca$layers[[2]], trait_pca$layers[[3]]) #move the arrows to plot in the foreground
ggsave(plot = trait_pca, filename = "trait.pca.png", path = "Figures")


###correlation
cormat <- cor(mean_traits)
corrplot(cormat, method = "number", type = "lower")
