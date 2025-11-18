###Script to create Three-D species lists
library(tidyverse)
library(openxlsx)

#import clean data
comm <- read.xlsx("All_data/clean_data/threed/community_2025.xlsx")

plotsp <- comm |> 
  group_by(destSiteID, origBlockID, turfID) |> 
  distinct(Species) |> 
  mutate(destplotID = as.numeric(str_split_i(turfID, "_", 3)), 
         survey_date = "March 2025") |> 
  arrange(destplotID) |> 
  select(survey_date, destSiteID, origBlockID, destplotID, turfID, Species)

write.xlsx(plotsp, "All_data/clean_data/threed/plot_composition_2025.xlsx")

##Ceate list of plots that need clipping
metaturfid <- read.xlsx("All_data/clean_data/threed/metaTurfID.xlsx") 

toclip <- metaturfid |> 
  filter(grazing %in% c("M", "I")) |> 
  arrange(destPlotID)

write.xlsx(toclip, "All_data/clean_data/threed/plots_to_clip.xlsx")

