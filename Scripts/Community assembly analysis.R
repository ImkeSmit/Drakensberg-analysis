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
preston(saunders.tot, n= 9)
optimal.theta(saunders.tot)
optimal.params(saunders.tot, gp_binary = "C:/Program Files/Pari64-2-17-2/gp.exe")
etienne(saunders.tot)

logkda(saunders.tot, method = "R")


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


