#import and clean community data from Three_D experiment
library(tidyverse)
library(readxl)
library(openxlsx)
library(stringi)

####CREATE METADATA####
#Create metadata file detailing the treatments associated with each turfID
# high and mid
origSiteID <-  c("high", "mid")
origBlockID <-  c(1:10)
origPlotID <- tibble(origPlotID = 1:160)
warming <-  c("A", "W")
grazing <-  c("C", "M", "I", "N")
# Nitrogen level needs to be in a certain order
nitrogen <- tibble(Nlevel = rep(rep(c(1,6,5,3,10,7,4,8,9,2), each = 8), 2))
# add corresponding N amount in kg per ha and year
NitrogenDictionary <- tibble(Nlevel = c(1,6,5,3,10,7,4,8,9,2),
                             Namount_kg_ha_y = c(0, 5, 1, 0, 150, 10, 0.5, 50, 100, 0))

# cross site, block warm and grazing treatment
meta <- crossing(origSiteID, origBlockID, warming, grazing) %>% 
  bind_cols(nitrogen)

# Vik (is done separately because it is only destination site)
vik <- tibble(
  origSiteID = factor("low", levels = c("high", "mid", "low")),
  origBlockID = rep(1:10, each = 4),
  origPlotID = 161:200,
  destSiteID = factor(NA, levels = c("high", "mid", "low")),
  Nlevel = rep(c(1,6,5,3,10,7,4,8,9,2), each = 4),
  warming = "W",
  grazing = rep(c("notN", "notN", "notN", "N"), 10),
  fence = if_else(grazing == "N", "out", "in"))

# randomize warming and grazing treatment
set.seed(32) # seed is needed to replicate sample_frac
meta2 <- meta %>% 
  # create variable for grazing treatment inside or outside fence
  mutate(fence = if_else(grazing == "N", "out", "in")) %>% 
  mutate(origSiteID = factor(origSiteID, levels = c("high", "mid", "low"))) %>%
  arrange(origSiteID) %>% # site needs to be arranged, because transplant goes only in one direction
  group_by(origSiteID, origBlockID, Nlevel, fence) %>%
  sample_frac() %>% # randomization
  ungroup() %>% 
  bind_cols(origPlotID) %>% # add plotID
  mutate(destSiteID = case_when(
    origSiteID == "high" & warming == "A" ~ "high",
    origSiteID == "mid" & warming == "W" ~ "low",
    TRUE ~ "mid")) %>%
  mutate(destSiteID = factor(destSiteID, levels = c("high", "mid", "low"))) %>%
  bind_rows(vik) %>% # add Vik
  group_by(origSiteID, origBlockID, warming, fence) %>% 
  mutate(rownr = row_number())


# Join meta2 to warmed plots
metaTurfID <- left_join(
  meta2 %>% filter(origPlotID < 161), # remove plots from vik
  # only warmed plots, remove unused rows
  meta2 %>% filter(warming == "W") %>% select(-grazing, -destSiteID, destPlotID = origPlotID), 
  by = c("destSiteID" = "origSiteID", "origBlockID" = "origBlockID", "rownr" = "rownr", "fence" = "fence", "Nlevel" = "Nlevel", "warming" = "warming"), 
  suffix = c("", "_dest")) %>% 
  mutate(destBlockID = origBlockID,
         destPlotID = ifelse(is.na(destPlotID), origPlotID, destPlotID),
         turfID = paste0(origPlotID, "_", warming, "N", Nlevel, grazing,  "_", destPlotID)) %>% 
  ungroup() %>% 
  select(-fence, -rownr) %>% 
  #CHANGE PLOTID 23-103 TO 23 AMBIENT, AND 24 TO 24-103 WARMING (wrong turf was transplanted!)
  mutate(warming = ifelse(origSiteID == "high" & origPlotID == 23, "A", warming),
         destPlotID = ifelse(origSiteID == "high" & origPlotID == 23, 23, destPlotID),
         turfID = ifelse(origSiteID == "high" & origPlotID == 23, "23_AN5N_23", turfID),
         
         warming = ifelse(origSiteID == "high" & origPlotID == 24, "W", warming),
         destPlotID = ifelse(origSiteID == "high" & origPlotID == 24, 103, destPlotID),
         turfID = ifelse(origSiteID == "high" & origPlotID == 24, "24_WN5N_103", turfID)) %>% 
  mutate(destSiteID = as.character(destSiteID)) %>% 
  mutate(destSiteID = case_when(turfID == "23_AN5N_23" ~ "high",
                                turfID == "24_WN5N_103" ~ "mid",
                                TRUE ~ destSiteID))

write.xlsx(metaTurfID, file = "All_data/clean_data/threed/metaTurfID.xlsx", colNames = TRUE)


####IMPORT VEG SURVEY DATA####
#run import_community script first
#CHeck import community script again before importing 2026 data
#script is custom for 2025 data
metadat <- read.xlsx("All_data/clean_data/threed/metaTurfID.xlsx", colNames = T)

veg2025 <- import_community(metadat, filepath = "All_data/raw_data/raw_threed_data/2025")

veg_only <- veg2025 |> 
  filter(!Species %in% c("Total Cover (%)","Vascular plants","Bryophyes","Lichen", "Litter","Bare soil",
                         "Bare rock","Poop","Height / depth (cm)","Vascular plant layer","Moss layer" ))

#There are a few plots with no total cover measurements. Add them from looking at the pictures. 

write.xlsx(veg_only, "All_data/clean_data/threed/community_2025.xlsx")

abiotic_only <- veg2025 |> 
  filter(Species %in% c("Vascular plants","Bryophyes","Lichen", "Litter","Bare soil",
                         "Bare rock","Poop","Vascular plant layer","Moss layer" )) |>
  select(!Cover) |> 
  rename(Variable = Species) |> 
  #rename variable names
  mutate(Variable = case_when(Variable == "Vascular plants" ~ "Vascular plant cover", 
                              Variable == "Bryophyes" ~ "Bryophyte cover", 
                              Variable == "Lichen" ~ "Lichen cover",
                              Variable == "Litter" ~ "Litter cover",
                              Variable == "Bare soil" ~ "Bare soil cover",
                              Variable == "Bare rock" ~ "Bare rock cover", 
                              Variable == "Poop" ~ "Poop cover",
                              Variable == "Vascular plant layer" ~ "Vascular plant layer height",
                              Variable == "Moss layer" ~ "Moss layer height", 
                              .default = Variable)) |>
  mutate(across(c(`1`:`25`), as.numeric)) #make all numeric

#create table of vascular height
vasc_height <- abiotic_only |> 
  filter(Variable == "Vascular plant layer height") |> 
  select(turfID, Variable, `1`:`4`) |> 
  rename(Vascular_plant_height1 = `1`, 
         Vascular_plant_height2 = `2`,
         Vascular_plant_height3 = `3`,
         Vascular_plant_height4 = `4`) |> 
  select(!Variable)

#create table of moss height
moss_height <- abiotic_only |> 
  filter(Variable == "Moss layer height") |> 
  select(turfID, Variable, `1`:`4`) |> 
  rename(Moss_layer_height1 = `1`, 
         Moss_layer_height2 = `2`,
         Moss_layer_height3 = `3`,
         Moss_layer_height4 = `4`) |> 
  select(!Variable)

#join moss and vascular height to other data
abiotic_only2 <- abiotic_only |> 
  filter(!Variable %in% c("Vascular plant layer height", "Moss layer height")) |>
  left_join(vasc_height, by = "turfID") |> 
  left_join(moss_height, by = "turfID")


#identify problematic cover values
summary(abiotic_only2) #cols1 -25 must have max cover of 100.5
#col 10 has max of 999

prob10 <- abiotic_only2 |> 
  slice_max(`10`, n = 10) #problem value is at 93_AN6M_93
#fix it
abiotic_only2[which(abiotic_only2$turfID == "93_AN6M_93" & abiotic_only2$Variable == "Vascular plant cover"), 
              which(colnames(abiotic_only2) == "10")] <- 99






#replace NA's with 0 where appropriate
replace_cols = colnames(abiotic_only2)[14:38]
for(r in 1:nrow(abiotic_only2)) {
    
    for(i in 1:length(replace_cols)) {
      
      element <- abiotic_only2[r, which(colnames(abiotic_only2) == replace_cols[i])]
      
      if(is.na(element)) {
        abiotic_only2[r, which(colnames(abiotic_only2) == replace_cols[i])] <- 0
      }
    }
  }
  

  
  
  
