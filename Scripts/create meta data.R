###################################
### CREATE TURFID AND RANDOMIZE ###
###################################

#source("R/Load packages.R")

# Create meta data
# Lia and Joa
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
#set.seed(33) # for Chinese experiment
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
         turfID = paste0(origPlotID, " ", warming, "N", Nlevel, grazing,  " ", destPlotID)) %>% 
  ungroup() %>% 
  select(-fence, -rownr) %>% 
  #CHANGE PLOTID 23-103 TO 23 AMBIENT, AND 24 TO 24-103 WARMING (wrong turf was transplanted!)
  mutate(warming = ifelse(origSiteID == "high" & origPlotID == 23, "A", warming),
         destPlotID = ifelse(origSiteID == "high" & origPlotID == 23, 23, destPlotID),
         turfID = ifelse(origSiteID == "high" & origPlotID == 23, "23 AN5N 23", turfID),
         
         warming = ifelse(origSiteID == "high" & origPlotID == 24, "W", warming),
         destPlotID = ifelse(origSiteID == "high" & origPlotID == 24, 103, destPlotID),
         turfID = ifelse(origSiteID == "high" & origPlotID == 24, "24 WN5N 103", turfID)) %>% 
  mutate(destSiteID = as.character(destSiteID)) %>% 
  mutate(destSiteID = case_when(turfID == "23 AN5N 23" ~ "high",
                                turfID == "24 WN5N 103" ~ "mid",
                                TRUE ~ destSiteID))

write.xlsx(metaTurfID, file = "All_data/clean_data/threed/metaTurfID.xlsx", colNames = TRUE)



#########################################################################
#### CHINA ####

# Create meta data
# Middle and Top
origSiteID <-  c("Top", "Middle")
origBlockID <-  c(1:10)
origPlotID <- tibble(origPlotID = 1:160)
warming <-  c("A", "W")
grazing <-  c("C", "M", "I", "N")
# Nitrogen level needs to be in a certain order
#nitrogen <- tibble(Nlevel = rep(rep(c(4,1,2,3,6,5,10,9,7,8), each = 8), 2))
nitrogen <- tibble(Nlevel = rep(rep(c(1,6,5,3,10,7,4,8,9,2), each = 8), 2))
# add corresponding N amount in kg per ha and year
# NitrogenDictionary <- tibble(Nlevel = c(4,1,2,3,6,5,10,9,7,8),
#                              Namount_kg_ha_y = c(0.5,0,0,0,5,1,150,100,10,50))
NitrogenDictionary <- tibble(Nlevel = c(1,6,5,3,10,7,4,8,9,2),
                             Namount_kg_ha_y = c(0,5,1,0,150,10,0.5,50,100,0))




# cross site, block warm and grazing treatment
meta <- crossing(origSiteID, origBlockID, warming, grazing) %>% 
  bind_cols(nitrogen)

# Low (is done separately because it is only destination site)
low <- tibble(
  origSiteID = factor("Low", levels = c("Top", "Middle", "Low")),
  origBlockID = rep(1:10, each = 4),
  origPlotID = 161:200,
  destSiteID = factor(NA, levels = c("Top", "Middle", "Low")),
  #Nlevel = rep(c(4,1,2,3,6,5,10,9,7,8), each = 4),
  Nlevel = rep(c(1,6,5,3,10,7,4,8,9,2), each = 4),
  warming = "W",
  grazing = rep(c("notN", "notN", "notN", "N"), 10),
  fence = if_else(grazing == "N", "out", "in"))

# randomize warming and grazing treatment
set.seed(33) # seed is needed to replicate sample_frac
meta2 <- meta %>% 
  # create variable for grazing treatment inside or outside fence
  mutate(fence = if_else(grazing == "N", "out", "in")) %>% 
  mutate(origSiteID = factor(origSiteID, levels = c("Top", "Middle", "Low"))) %>%
  arrange(origSiteID) %>% # site needs to be arranged, because transplant goes only in one direction
  group_by(origSiteID, origBlockID, Nlevel, fence) %>%
  sample_frac() %>% # randomization
  ungroup() %>% 
  bind_cols(origPlotID) %>% # add plotID
  mutate(destSiteID = case_when(
    origSiteID == "Top" & warming == "A" ~ "Top",
    origSiteID == "Middle" & warming == "W" ~ "Low",
    TRUE ~ "Middle")) %>% 
  mutate(destSiteID = factor(destSiteID, levels = c("Top", "Middle", "Low"))) %>%
  bind_rows(low) %>% # add low
  group_by(origSiteID, origBlockID, warming, fence) %>% 
  mutate(rownr = row_number())


# Join meta2 to warmed plots
metaTurfID_China <- left_join(
  meta2 %>% filter(origPlotID < 161), # remove plots from vik
  # only warmed plots, remove unused rows
  meta2 %>% filter(warming == "W") %>% select(-grazing, -destSiteID, destPlotID = origPlotID), 
  by = c("destSiteID" = "origSiteID", "origBlockID" = "origBlockID", "rownr" = "rownr", "fence" = "fence", "Nlevel" = "Nlevel", "warming" = "warming"), 
  suffix = c("", "_dest")) %>% 
  mutate(destBlockID = origBlockID,
         destPlotID = ifelse(is.na(destPlotID), origPlotID, destPlotID),
         turfID = paste0(origPlotID, " ", warming, "N", Nlevel, grazing,  " ", destPlotID)) %>% 
  ungroup() %>% 
  select(-fence, -rownr)
