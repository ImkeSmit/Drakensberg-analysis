##Commmunity assembly analysis####
library(openxlsx)
library(tidyverse)
library(tidylog)
library(untb)
library(ggplot2)
library(vegan)
library(VGAM)
library(elliptic)
Sys.setenv(gp_binary = "C:/Program Files/Pari64-2-17-2/gp.exe")
gp_binary = "C:/Program Files/Pari64-2-17-2/gp.exe"

data("saunders")
summary(saunders.tot)

#we have to calculate logkda with method R (not in optimal.params bcause Pari won't work)
k <- logkda(saunders.tot, method = "R")
zsm_params <- optimal.params(D = saunders.tot, log.kda = k)
#theta          m 
#31.5741581  0.9844249

#now we can fit the zsm
zsm_fit <- zsm(J = sum(saunders.tot), #J = size of local comm, the sum of abundances
               P = 0.001, #Abundance in the metacommunity is correlated with theta
               m = zsm_params[2]) #probability of immigration from metacomm


prob_density <- data.frame(prob = zsm_fit, abundance = seq(from = 1, to = sum(saunders.tot) +1, by = 1))

ggplot(prob_density, aes(x = abundance, y = prob)) +
  geom_line() +
  xlim(1, 100)

exp_ab <- expected.abundance(J = 10, theta = zsm_params[1])
sum(prob_density) #this should be one



###Plot SAD for drakensberg data

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

ggplot(ab, aes(x = log(total_cover))) +
  geom_histogram() +
  ylab("Number of species") +
  xlab("log Total cover")


#Fit lognormal distribution to abundance data
lognormal_fit <- rad.lognormal(c(ab$abundance)) #AIC = 9448.980570
plot(lognormal_fit)

preston <- prestonfit(c(ab$abundance))
plot(preston) #does the weird binning


#Fit fischer logseries
f_alpha <- fisher.alpha(c(ab$abundance))
#alpha = 61.60134
fisher <- fisherfit(c(ab$abundance))
plot(fisher)


#fit zsm
#transfrom to count object required by untb
abundance_count <- c(ab$abundance)
names(abundance_count) <- c(ab$taxon)
abundance_count <- count(abundance_count)

zsm_params <- optimal.params(abundance_count)

zsm_fit <- zsm(abundance_count, )
et <- ettiene(abundance_count)
theta <- optimal.theta(abundance_count)


