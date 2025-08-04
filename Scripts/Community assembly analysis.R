##Commmunity assembly analysis####
library(openxlsx)
library(tidyverse)
library(tidylog)
library(ggplot2)
library(sads)

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
#package can fit nuetral models from volkov or from Alonso and McKane. 
#read papers to determine which is most appropriate
neutral <- fitsad(ab$abundance, sad = "mzsm")
summary(neutral)
AIC(neutral)


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

#likelihood ratio test
test <- anova(neutral, lognorm) #can we do a likelihood ratio test this way? the models are not nested

#extract loglikelihood. Higher means better fit
logLik(neutral) #-2158.167 (df=1)
logLik(lognorm) #-2188.739 (df=2)
