##Commmunity assembly analysis####
library(openxlsx)
library(tidyverse)
library(tidylog)
library(untb)
library(ggplot2)
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
