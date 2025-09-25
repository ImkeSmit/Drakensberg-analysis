#Community assembly analysis 2.0###
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