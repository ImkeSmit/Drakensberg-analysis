#Does community assembly mechanisms vary with scale and with elevation?#
#ALL SPECIES###
library(tidyverse)
library(tidylog)
library(openxlsx)
library(ggplot2)
library(ggridges)
library(gamlss)
library(rcompanion)
library(multcompView)

###Import original community data####
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1)

###TEST FOR EF####
###All sp - Cell Scale, C5, pool = entire####
#import SES at cell scale, computed from unscaled RaoQ
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA))
cell_ses$elevation <- as.factor(cell_ses$elevation)  

####Descriptive figures####

RQ_ridges <- cell_ses |> 
  ggplot(aes(x = RaoQ, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~trait, scale = "free_x", nrow = 3, ncol = 2) +
  theme_classic()
#RaoQ distributions also very heavy tailed. Thus very few cells with high Fdiv. Is this normal?
#Could using another Fdiv measure help? probably not...
ggsave(path = "Figures",plot = RQ_ridges, filename = "RaoQ_elevation_by_traits.png",
       width = 1200, height = 1400, units = "px")

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
  summarise(n = n()) #very few have significant negative ses


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
#test if median ses differ between elevation groups
kr_height <- kruskal.test(SES ~ elevation, data = modeldat) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_height <- pairwise.wilcox.test(modeldat$SES,  modeldat$elevation, p.adjust.method = "bonferroni")
height_cld <- multcompLetters(fullPTable(wx_height$p.value))

#get medians
median_height <- cell_ses |> 
  filter(trait == "Height_cm") |> 
  group_by(elevation) |> 
  summarise(median = median(SES,  na.rm = T))

#2500 lower than 2000, 3000 lower than 2000, 3000 higher than 2500
#functional convergence increases with elevation (SES decreases with elevation)

#test if medians differ from zero - wilcoxon signed rank test
height_2000 <- wilcox.test(modeldat[which(modeldat$elevation == "2000"), ]$SES, mu = 0, alternative = "two.sided")
height_2500 <- wilcox.test(modeldat[which(modeldat$elevation == "2500"), ]$SES, mu = 0, alternative = "two.sided")
height_3000 <- wilcox.test(modeldat[which(modeldat$elevation == "3000"), ]$SES, mu = 0, alternative = "two.sided")
height_stars <- c(" ", "*", "*")

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
#Functional convergence decreases with elevation (median SES switches from negative to positive)

#test if medians differ from zero - wilcoxon signed rank test
LDMC_2000 <- wilcox.test(modeldat_LDMC[which(modeldat_LDMC$elevation == "2000"), ]$SES, mu = 0, alternative = "two.sided")
LDMC_2500 <- wilcox.test(modeldat_LDMC[which(modeldat_LDMC$elevation == "2500"), ]$SES, mu = 0, alternative = "two.sided")
LDMC_3000 <- wilcox.test(modeldat_LDMC[which(modeldat_LDMC$elevation == "3000"), ]$SES, mu = 0, alternative = "two.sided")
LDMC_stars <- c("*", "*", "*")



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

#test if medians differ from zero - wilcoxon signed rank test
LA_2000 <- wilcox.test(modeldat_LA[which(modeldat_LA$elevation == "2000"), ]$SES, mu = 0, alternative = "two.sided")
LA_2500 <- wilcox.test(modeldat_LA[which(modeldat_LA$elevation == "2500"), ]$SES, mu = 0, alternative = "two.sided")
LA_3000 <- wilcox.test(modeldat_LA[which(modeldat_LA$elevation == "3000"), ]$SES, mu = 0, alternative = "two.sided")
LA_stars <- c("*", "*", "*")



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

#2500 higher than 2000, 3000 higher than 2000, 3000 lower than 2500
#Functional convergence decreases between 200 and 2500, then increases at 3000
#strongets EF at 3000, weakest at 2500

#test if medians differ from zero - wilcoxon signed rank test
SLA_2000 <- wilcox.test(modeldat_SLA[which(modeldat_SLA$elevation == "2000"), ]$SES, mu = 0, alternative = "two.sided")
SLA_2500 <- wilcox.test(modeldat_SLA[which(modeldat_SLA$elevation == "2500"), ]$SES, mu = 0, alternative = "two.sided")
SLA_3000 <- wilcox.test(modeldat_SLA[which(modeldat_SLA$elevation == "3000"), ]$SES, mu = 0, alternative = "two.sided")
SLA_stars <- c("*", "*", "*")


####Leaf thickness####
modeldat_LT <- cell_ses[which(cell_ses$trait == "Thickness_mm"), ]

#Kruskal-wallis test
kr_LT <- kruskal.test(SES ~ elevation, data = modeldat_LT) #medians of at least two groups differ
# Conduct pairwise comparisons with Wilcoxon rank-sum test
wx_LT <- pairwise.wilcox.test(modeldat_LT$SES,  modeldat_LT$elevation, p.adjust.method = "bonferroni")
LT_cld <- multcompLetters(fullPTable(wx_LT$p.value))

#all are significantly different

#get medians
median_LT <- cell_ses |> 
  filter(trait == "Thickness_mm") |> 
  group_by(elevation) |> 
  summarise(median = median(SES,  na.rm = T))

#2500 higher than 2000, 3000 lower than 2000, 3000 lower than 2500
#Functional divergence at 2500, convergence at low and high elevation. 

#test if medians differ from zero - wilcoxon signed rank test
LT_2000 <- wilcox.test(modeldat_LT[which(modeldat_LT$elevation == "2000"), ]$SES, mu = 0, alternative = "two.sided")
LT_2500 <- wilcox.test(modeldat_LT[which(modeldat_LT$elevation == "2500"), ]$SES, mu = 0, alternative = "two.sided")
LT_3000 <- wilcox.test(modeldat_LT[which(modeldat_LT$elevation == "3000"), ]$SES, mu = 0, alternative = "two.sided")
LT_stars <- c("*", "*", "*")


####SES~elevation summary figure####
ses_ridges <- cell_ses |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~trait, nrow = 3, ncol = 2) +
  theme_classic() 

#significance letters
letters_df <- tibble(trait = c(rep("Height_cm", 3), rep("LDMC", 3), rep("Leaf_area_mm2", 3), 
                               rep("SLA", 3), rep("Thickness_mm", 3)),
                     elevation = rep(c("2000", "2500", "3000"), 5), 
                     letters = c(unlist(height_cld)[1:3], unlist(LDMC_cld)[1:3], 
                                 unlist(LA_cld)[1:3], unlist(SLA_cld)[1:3], unlist(LT_cld)[1:3]))
#medians
med_df <- cell_ses %>%
  group_by(trait, elevation) %>%
  summarize(med = median(SES, na.rm = TRUE))

#stars
stars_df <- tibble(trait = c(rep("Height_cm", 3), rep("LDMC", 3), rep("Leaf_area_mm2", 3), 
                             rep("SLA", 3), rep("Thickness_mm", 3)), 
                     elevation = rep(c("2000", "2500", "3000"), 5), 
                     stars = c(height_stars, LDMC_stars, 
                                 LA_stars, SLA_stars, LT_stars))

#improved figure
ses_ridges2 <- ses_ridges+
  geom_segment(data = med_df,
    aes(x = med, xend = med,
        y = as.numeric(elevation) - 0.01,
        yend = as.numeric(elevation) + 0.1),
    linetype = "solid", size = 0.6) +
  geom_text(data = letters_df,
    aes(x = 7, y = elevation, label = letters),
    color = "black",size = 4) +
  geom_text(data = stars_df,
            aes(x = 6, y = elevation, label = stars),
            color = "black",size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  theme(legend.position = "bottom")

ses_ridges3 <- ses_ridges2 +
  xlim(min(cell_ses$SES),3) #remove heavy tails for now

ggsave(filename = "C5_SES_elevation_by_traits.png", plot = ses_ridges3, path= "Figures", 
       width = 1200, height = 1500, units = "px")



###What is going on in the cells with high SES values?
high <- cell_ses |> 
  slice_max(SES, n = 20)

high_plots <- cell_ses |> 
  slice_max(SES, n = 20) |> 
  distinct(cellref)

high_mat <- abun_matrix[which(row.names(abun_matrix) %in% c(high_plots$cellref)) , ]
specnumber(high_mat) #do not necessarily have few species

all_abun <- as.data.frame(specnumber(abun_matrix))
all_abun$cellref <- row.names(all_abun)
row.names(all_abun) <- NULL
colnames(all_abun) <- c( "sprichness", "cellref")
#NB!This is the number of sp in a plot that have trait values
#lowest is 2, highest is 30

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
#but now only remove 16 plots


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
