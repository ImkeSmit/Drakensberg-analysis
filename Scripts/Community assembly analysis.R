##Commmunity assembly analysis####
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(sads)
library(vegan)

drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx")

ab <- drak |> 
  group_by(taxon) |> 
  summarise(abundance = n(), 
            total_cover = sum(cover)) |> 
  ungroup()
ab$abundance <- as.numeric(ab$abundance)
ab <- as.data.frame(ab)
ab <- ab[order(ab$abundance, decreasing = TRUE), ]

#SAD
SAD <- ggplot(ab, aes(x = abundance)) +
  geom_histogram() +
  ylab("number of species") +
  xlab("Number of cells occupied")

SAD_cover <- ggplot(ab, aes(x = total_cover)) +
  geom_histogram() +
  ylab("Number of species") +
  xlab("Total cover")

#RAD
drak_rad <- rad(ab$abundance)
plot(drak_rad)


#fit lognormal distribution
lognorm <- fitsad(ab$abundance, sad = "lnorm")
lognorm
summary(lognorm)
AIC(lognorm)


#fit fisher logseries
logseries <- fitsad(ab$abundance, sad = "ls")
summary(logseries)
AIC(logseries)

#fit zsm
#ML method of alonso and McKane
neutral_ML <- fitsad(ab$abundance, sad = "mzsm")
summary(neutral_ML)
AIC(neutral_ML)

#ANALYTICAL method of Volkov
neutral_ana <- fitsad(ab$abundance, sad = "volkov") #this will take a long time
#start 13:09


#visually inspect model fit
plot(lognorm)
plot(logseries)
plot(neutral)

#plot fitted models over observed abundance octaves
drak_octaves <- octav(ab$abundance)
plot(drak_octaves)

#get predictions from each model over the octaves
lognorm_pred <- octavpred(lognorm)
logseries_pred <- octavpred(logseries)
neutral_pred <- octavpred(neutral)

#plot predictions
plot(drak_octaves)
lines(lognorm_pred, col = "red")
lines(logseries_pred, col = "blue")
lines(neutral_pred, col = "green")

#compare AIC values
AICtab(lognorm, logseries, neutral, base = T)
#neutral and logseries are equivalent

#likelihood ratio test
test <- anova(neutral, lognorm) #can we do a likelihood ratio test this way? the models are not nested

#extract loglikelihood. Higher means better fit
logLik(neutral) #-2158.167 (df=1)
logLik(lognorm) #-2188.739 (df=2)


#####DO TRAITS PREDICT ABUNDANCES####
#Import trait data and get mean trait for each sp
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  filter(Taxon %in% c(unique(drak$taxon))) |>  #remove sp that aren't in abundance data
  group_by(Taxon) |> 
  summarise(mean_SLA = mean(SLA, na.rm = T), #calculate mean trait per species
            mean_LDMC = mean(LDMC, na.rm = T), #remove NA vals before calculating mean
            mean_Height = mean(Height_cm, na.rm = T), 
            mean_Thickness = mean(Thickness_mm, na.rm = T), 
            mean_LA = mean(Leaf_area_mm2), na.rm = T) |> 
  mutate(z_SLA = as.vector(scale(mean_SLA)), #standardise traits to mean = 0 and sd = 1
         z_LDMC = as.vector(scale(mean_LDMC)), 
         z_Height = as.vector(scale(mean_Height)), 
         z_Thickness =  as.vector(scale(mean_Thickness)), 
         z_LA =  as.vector(scale(mean_LA))) |> 
  select(Taxon, z_LDMC, z_Height, z_Thickness, z_LA, z_SLA) |>
  drop_na() #remove rows that have NA in any of the columns
  

pca <- princomp(FT[, c(2:6)])
summary(pca) #proportion of variance is the variance explained by the PC
loadings(pca) #How much each var contributed to building the PC
plot(pca)
biplot(pca)
pca$scores
pca$scale #scaling applied to each var. Should be 1 because I standardised all variables?

#extract pca scores
pca_scores <- data.frame(taxon = FT$Taxon, 
                         PC1 = pca$scores[, 1], 
                         PC2 = pca$scores[, 2], 
                         PC3 = pca$scores[, 3])



#Now we need to get the relative abundances of each species
ab <- ab |> 
  mutate(rel_abundance = abundance/sum(abundance), 
         rel_cover = total_cover/sum(total_cover))


#join PC scores to relative abundances
modeldat <- ab |> 
  inner_join(pca_scores, by = "taxon")
#in this dataset there is only one singleton, probably because we removed the sp that were not sampled for traits
#many singletons in the occurrence data were likely not found again for trait sampling
#problem???

###model relative abundance ~ trait score
hist(modeldat$rel_abundance)
mod1 <- glmmTMB(rel_abundance ~ PC1 +PC2 + PC3, data = modeldat, family = lognormal(link = "log"))
summary(mod1)

#can also use standard effect size approach and not model
