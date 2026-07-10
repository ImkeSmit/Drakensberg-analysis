###Models of CWm ~ microenvironmental variables with lme####
library(tidyverse)
library(tidylog)
library(nlme)
library(MuMIn)
library(DHARMa)
library(tictoc)
library(ggridges)
library(emmeans)
library(multcomp)
library(conflicted)
library(performance)
library(see)
conflict_prefer_all("tidylog", quiet = TRUE)

#import trait and abundance data
abun_matrix <-read.csv("All_data/comm_assembly_results/abun_matrix.csv", row.names = 1)

mean_traits <- read.csv("All_data/comm_assembly_results/mean_traits.csv", row.names = 1) |> 
  mutate(log_Height = log(Height_cm), 
         log_LA = log(Leaf_area_mm2), 
         log_LDMC = log(LDMC)) |> #use the same log variables as we do in the models
  select(!c(LDMC, Leaf_area_mm2, Height_cm, Thickness_mm))


#compute CWM of each trait for each cell
library(FD)
cwm <- functcomp(x = mean_traits, a = as.matrix(abun_matrix))
cwm2 <- cwm |>
  rownames_to_column(var = "cellref") |> 
  rename(Cell_ID = cellref) |> 
  mutate(elevation = case_when(grepl("BK", Cell_ID) == T ~ "3000", #add elevation variable
                               grepl("WH", Cell_ID) == T ~ "2500",
                               grepl("GG", Cell_ID) == T ~ "2000",.default = NA)) |> 
  pivot_longer(cols = c("log_Height", "log_LDMC", "log_LA", "SLA"), names_to = "trait", values_to = "cwm_value")

#import env data so that we can join x and y coordinates
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row) 

cwm_comb <- env |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  inner_join(cwm2, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)


####model CWM ~ elevation####
traitlist <- c("log_Height", "log_LDMC", "log_LA", "SLA")
#lists to store results in
CWM_ele_cld<- vector(mode= "list", length = length(traitlist))
names(CWM_ele_cld) = traitlist

for (t in 1:length(traitlist)) {
  modeldat <-  cwm_comb |> 
    filter(trait == traitlist[t]) |> 
    mutate(elevation = as.factor(elevation), 
           grid = as.factor(grid)) |> 
    drop_na()
  
  
  model<- lme(cwm_value ~ elevation ,
              random = ~1|grid, 
              correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
              data = modeldat) #only gaussian family possible
  
  em_model <- emmeans(model, specs = "elevation", type = "response")
  comp_letters <-cld(em_model, Letters = letters, adjust = "Tukey", sort = FALSE)
  
  #Also save results in a list object so we can call it with quarto
  CWM_ele_cld[[t]] <- comp_letters
}