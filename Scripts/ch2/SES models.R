#Does community assembly mechanisms vary with scale and with elevation?#
library(tidyverse)
library(tidylog)
library(openxlsx)
library(ggplot2)
library(ggridges)
library(gamlss)
library(rcompanion)
library(multcompView)

###Import original community data####
#occurrence data
drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))

#trait data
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) |> 
  pivot_longer(cols = c(Wet_mass_mg, Dry_mass_mg, Chlorophyll_mg_per_m2, Ft, Height_cm, 
                        Thickness_mm, Leaf_area_mm2, SLA, LDMC), names_to = "trait", values_to = "value")

#combine occurrence and trait data
FT_join <- drak |> 
  inner_join(FT, by = c("site", "grid", "cell", "taxon")) |> #inner join to only work with taxa that have trait data
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

#create abundance matrix
abun_matrix <- FT_join |> 
  select(cellref, taxon, cover) |> 
  distinct(cellref, taxon, .keep_all = T) |> 
  ungroup() |> 
  arrange(taxon) |> 
  mutate(cover = ceiling(cover)) |> #change cover values to integer to use in null models
  pivot_wider(names_from = taxon, values_from = cover) 

abun_matrix <- as.data.frame(abun_matrix)
row.names(abun_matrix) <- abun_matrix$cellref
abun_matrix <- abun_matrix[, -1]

#replace NA values with 0
for(r in 1:nrow(abun_matrix)) {
  for(c in 1:ncol(abun_matrix)) {
    
    if(is.na(abun_matrix[r,c])) {
      abun_matrix[r,c] <- 0
    }
  }
}

###TEST FOR EF####
###Cell Scale, C5, pool = entire####
#import SES at cell scale, computed from unscaled RaoQ
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA))
cell_ses$elevation <- as.factor(cell_ses$elevation)  

####Descriptive figures####

RQ_ridges <- cell_ses |> 
  ggplot(aes(x = RaoQ, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()
#RaoQ distributions also very heavy tailed. Thus very few cells with high Fdiv. Is this normal?
#Could using another Fdiv measure help? probably not...
ggsave(path = "Figures",plot = RQ_ridges, filename = "RaoQ_elevation_by_traits.png")

ses_ridges <- cell_ses |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()


###Descriptive statistics####
#How many cells are significantly less or more diverse than expected under null model?
sign_positive_ses <- cell_ses |> 
  group_by(trait) |> 
  filter(SES > 2) |> 
  summarise(n = n())

sign_negative_ses <- cell_ses |> 
  group_by(trait) |> 
  filter(SES < -2) |> 
  summarise(n = n()) #hmm none have ses lower than -2, so none are significantly less diverse than null?


###Models of SES ~ elevation####
###Height####
modeldat <- cell_ses[which(cell_ses$trait == "Height_cm"), ]

#USe skew SHASHo distribution as recommended by chatgpt.Can handle skewed and heavy tailed data with nonzero values
#other distribution to try is skew normal SN()
test2 <- gamlss(SES ~ elevation, data = modeldat, family = SHASHo()) #skew t does not converge
#only one categorical predictor, thus it is equivalent to an anova under the specified distribution.
#chatgpt says this is a totally sensible approach due to the right skewed nature of my response variable
#warning message, algorthm RS has not yet converged
summary(test2)
plot(test2) #looks much better

#Kruskal-wallis test
kr_height <- kruskal.test(SES ~ elevation, data = modeldat) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_height <- pairwise.wilcox.test(modeldat$SES,  modeldat$elevation, p.adjust.method = "bonferroni")
height_cld <- multcompLetters(fullPTable(wx_height$p.value))

#get medians
median_height <- cell_ses |> 
  filter(trait == "Height_cm") |> 
  group_by(elevation) |> 
  summarise(median = median(SES,  na.rm = T))

#2500 lower than 2000, 3000 lower than 2000, 3000 not lower than 2500
#functional convergence increases with elevation (SES decreases with elevation)


####LDMC####
modeldat_LDMC <- cell_ses[which(cell_ses$trait == "LDMC"), ]

test_LDMC <- gamlss(SES ~ elevation, data = modeldat_LDMC, family = SEP2()) #sshasho does not converge
#only one categorical predictor, thus it is equivalent to an anova under the specified distribution.
#chatgpt says this is a totally sensible approach due to the right skewed nature of my response variable
#warning message, algorthm RS has not yet converged
summary(test_LDMC)
plot(test_LDMC) #looks much better

#Kruskal-wallis test
kr_LDMC <- kruskal.test(SES ~ elevation, data = modeldat_LDMC) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_LDMC <- pairwise.wilcox.test(modeldat_LDMC$SES,  modeldat_LDMC$elevation, p.adjust.method = "bonferroni")
LDMC_cld <- multcompLetters(fullPTable(wx_LDMC$p.value))


#get medians
median_LDMC <- cell_ses |> 
  filter(trait == "LDMC") |> 
  group_by(elevation) |> 
  summarise(median = median(SES,  na.rm = T))

#2500 not different from 2000, 3000 higher than 2000, 3000 higher than 2500
#Functional convergence decreases with elevation (median SES switches from negative to positive)


###Leaf area####
modeldat_LA <- cell_ses[which(cell_ses$trait == "Leaf_area_mm2"), ]

test_LA <- gamlss(SES ~ elevation, data = modeldat_LA, family = ST3()) #sshasho does not converge
summary(test_LA)
plot(test_LA) 

#Kruskal-wallis test
kr_LA <- kruskal.test(SES ~ elevation, data = modeldat_LA) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_LA <- pairwise.wilcox.test(modeldat_LA$SES,  modeldat_LA$elevation, p.adjust.method = "bonferroni")
#all medians differ significantly
LA_cld <- multcompLetters(fullPTable(wx_LA$p.value))


#get medians
median_LA <- cell_ses |> 
  filter(trait == "Leaf_area_mm2") |> 
  group_by(elevation) |> 
  summarise(median = median(SES,  na.rm = T))

#2500 higher than 2000, 3000 lower than 2000, 3000 lower than 2500
#Functional convergence decreases between 2000 and 2500, but increases at 3000.
#EF is strongest at 3000, and weakest at 2500


####SLA####
modeldat_SLA <- cell_ses[which(cell_ses$trait == "SLA"), ]

test_SLA <- gamlss(SES ~ elevation, data = modeldat_SLA, family = ST4()) #sshasho does not converge

summary(test_SLA)
plot(test_SLA) #pretty good

#Kruskal-wallis test
kr_SLA <- kruskal.test(SES ~ elevation, data = modeldat_SLA) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_SLA <- pairwise.wilcox.test(modeldat_SLA$SES,  modeldat_SLA$elevation, p.adjust.method = "bonferroni")
SLA_cld <- multcompLetters(fullPTable(wx_SLA$p.value))

#all are significantly different

#get medians
median_SLA <- cell_ses |> 
  filter(trait == "SLA") |> 
  group_by(elevation) |> 
  summarise(median = median(SES,  na.rm = T))

#2500 higher than 2000, 3000 lower than 2000, 3000 lower than 2500
#Functional convergence decreases between 200 and 2500, then increases at 3000
#strongets EF at 3000, weakest at 2500

####SES~elevation summary figure####
ses_ridges <- cell_ses |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic() 

letters_df <- tibble(trait = c(rep("Height", 3), rep("LDMC", 3), rep("Leaf_area_mm2", 3), rep("SLA", 3)), 
                     elevation = rep(c("2000", "2500", "3000"), 4), 
                     letters = c(multcompLetters(fullPTable(wx_height$p.value)), multcompLetters(fullPTable(wx_LDMC$p.value)), 
                                 multcompLetters(fullPTable(wx_LA$p.value)), multcompLetters(fullPTable(wx_SLA$p.value))))
  



###What is going on in the cells with high SES values?
high <- cell_ses |> 
  slice_max(SES, n = 20)

high_plots <- cell_ses |> 
  slice_max(SES, n = 20) |> 
  distinct(cellref)

high_mat <- abun_matrix[which(row.names(abun_matrix) %in% c(high_plots$cellref)) , ]
specnumber(high_mat) #all have between 2 and 6 species

all_abun <- as.data.frame(specnumber(abun_matrix))
all_abun$cellref <- row.names(all_abun)
row.names(all_abun) <- NULL
colnames(all_abun) <- c( "sprichness", "cellref")
#NB!This is the number of sp in a plot that have trait values!!

cell_ses2 <- cell_ses |> 
  left_join(all_abun, by = "cellref") 

#Look at distribution of SES when we remove plots with 2 sp
ses2_ridges <- cell_ses2 |> 
  filter(sprichness > 2) |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()
#does not really improve heavy tails
#also removes 42% of plots


###Cell Scale, C2, pool = site####
#import SES at cell scale, computed from unscaled RaoQ
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_cells_C2_site.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA))
cell_ses$elevation <- as.factor(cell_ses$elevation)  



ses_ridges <- cell_ses |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()
#only GG sites retained! need to force each sp to have ant least one occurrence


###Grid scale####
grid_ses <- read.csv("All_data/comm_assembly_results/RQ_grids_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ 3000, #add elevation variable
                               grepl("WH", cellref) == T ~ 2500,
                               grepl("GG", cellref) == T ~ 2000,.default = NA))


grid_ses_ridges <- grid_ses |> 
  mutate(elevation_char = as.character(elevation)) |> 
  ggplot(aes(x = SES, y = elevation_char, fill = elevation_char)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()
