###Modelling with nlme, no subsampling####
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

###=========================================###
###Loop to run SES~elevation for all traits####
###=====C5 NULL MODEL, pool = entire========####
###=========================================###
#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/SES_RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  rename(Cell_ID = cellref) 

#import microenvironmental data containing x and y coordinates
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA))

##Combine SES and environmental data
comb <- env |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  inner_join(cell_ses, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)
  
#These variables have the best diagnostics
#LA and SLA still have TERRIBLE DIAGNOSTICS
traitlist <- c("log_Height", "log_LDMC", "log_LA", "SLA")
#lists to store results in
SES_ele_summary <- vector(mode= "list", length = length(traitlist))
names(SES_ele_summary) = traitlist

SES_ele_Rsq<- vector(mode= "list", length = length(traitlist))
names(SES_ele_Rsq) = traitlist

SES_ele_anova<- vector(mode= "list", length = length(traitlist))
names(SES_ele_anova) = traitlist

SES_ele_cld<- vector(mode= "list", length = length(traitlist))
names(SES_ele_cld) = traitlist

for (t in 1:length(traitlist)) {
  modeldat <-  comb |> 
    filter(trait == traitlist[t]) |> 
    mutate(elevation = as.factor(elevation), 
           grid = as.factor(grid)) |> 
    drop_na()
  
  
  model<- lme(SES ~ elevation ,
                  random = ~1|grid, 
                  correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
                  data = modeldat) #only gaussian family possible
  
  ###Save check_model plot
  plot_file <- paste0("All_data/comm_assembly_results/C5_null_model_results/checkmodel_SES_elevation/checkmodel_lme_SES_", traitlist[t], "_elevation.png")
  png(plot_file, width = 1600, height = 1200, res = 150)
  print(check_model(model))   # print() forces the plot to actually draw to the device
  dev.off()
  
  ###Save model results
  output_file <- paste0("All_data/comm_assembly_results/C5_null_model_results/lme_results_SES_elevation/lme_SES_" ,traitlist[t], "_elevation_results.txt")
  sink(output_file)
  
  # ── 1.Trait ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  TRAIT\n")
  cat("===========================================\n")
  print(traitlist[t])
  cat("\n\n")
  
  
  # ── 1. Model Formula ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  MODEL FORMULA\n")
  cat("===========================================\n")
  print(formula(model))
  cat("\n\n")
  
  # ── 2. Summary Table ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  MODEL SUMMARY\n")
  cat("===========================================\n")
  print(summary(model))
  cat("\n\n")
  
  
  # ── 2. R squared ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  R SQUARED\n")
  cat("===========================================\n")
  print(r.squaredGLMM(model))
  cat("\n\n")
  
  
  # ── 3. ANOVA Table ────────────────────────────────────────────
  cat("===========================================\n")
  cat("  ANOVA TABLE\n")
  cat("===========================================\n")
  print(anova(model))
  cat("\n\n")
  
  # ── 4. EMmeans Table ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  ESTIMATED MARGINAL MEANS (emmeans)\n")
  cat("===========================================\n")
  em_model <- emmeans(model, specs = "elevation", type = "response")
  comp_letters <-cld(em_model, Letters = letters, adjust = "Tukey", sort = FALSE)
  print(comp_letters)
  cat("\n")
  
  # --- Close the sink ---
  sink()
  
  #Also save results in a list object so we can call it with quarto
  SES_ele_summary[[t]] <- summary(model)
  SES_ele_Rsq[[t]] <- r.squaredGLMM(model)
  SES_ele_anova[[t]] <- anova(model)
  SES_ele_cld[[t]] <- comp_letters
  
}

###Using the logged traits are generally better for model diagnostics than the raw traits
##BUT the leaf area diagnostics are still very bad




###===========================================================###
###Loop to run SES~microenvironmental variables for all traits####
###=============C5 NULL MODEL, pool = site====================####
###===========================================================###
#import ses
cell_ses_poolsite <- read.csv("All_data/comm_assembly_results/SES_RQ_weighted_cells_C5_poolsite.csv", row.names = 1) |> 
  rename(Cell_ID = cellref) 

#import microenvironmental data
env2 <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, 
         soil_depth_CV, mean_soil_depth, slope_height) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA))

#import remote sensing derived variables
rms <- read.csv("All_data/clean_data/Environmental data/Zonal_stats_all.csv") |> 
  select(CELL_ID, STD) |> 
  rename(Cell_ID = CELL_ID)

#import interpolated microclimate indices
micro_idw <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices_idw_interpolated.csv", row.names = 1)


##Combine SES and environmental data
comb2 <- env2 |> 
  #join to microclimate indices |> 
  full_join(micro_idw, by = "Cell_ID") |> 
  #join to remote sensing data |> 
  full_join(rms, by = "Cell_ID") |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cell_ses_poolsite, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)


###Check collinearity#### 
library(corrplot)
cordf <- comb2 |> 
  filter(trait == "Height_cm") |> #look at just one set of env data, it repeats for every trait
  select(mean_T1_growing_season, mean_moist_growing_season, STD, rock_cover, northness, mean_soil_depth, slope_height) |> 
  drop_na()
cormat<- cor(cordf)
corrplot(cormat, type = "lower", method = "number")


###Loop starts here####
#These variables have the best diagnostics after comparing diagnostic plots of log and raw variables
traitlist <- c("Height_cm","Leaf_area_mm2", "log_LDMC","SLA")
sitelist <- c("GG", "WH", "BK")
#lists to store results in
SES_microenv_summary <- vector(mode= "list", length = length(traitlist) * length(sitelist))

SES_microenv_Rsq<- vector(mode= "list", length = length(traitlist) * length(sitelist))

SES_microenv_anova<- vector(mode= "list", length = length(traitlist) * length(sitelist))

l = 1
for (t in 1:length(traitlist)) {
  for(s in 1:length(sitelist)) {
  modeldat <-  comb2 |> 
    filter(trait == traitlist[t], 
           site == sitelist[s]) |> 
    mutate(elevation = as.factor(elevation), 
           grid = as.factor(grid)) |> 
    drop_na()
  
  
  model<- lme(SES ~ mean_T1_growing_season + mean_moist_growing_season + rock_cover+ mean_soil_depth + STD,
              random = ~1|grid, 
              correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
              data = modeldat) #only gaussian family possible
  
  ###Save check_model plot
  plot_file <- paste0("All_data/comm_assembly_results/checkmodel_SES_microenv/checkmodel_lme_SES_", traitlist[t], "_microenv_", sitelist[s], ".png")
  png(plot_file, width = 1600, height = 1200, res = 150)
  print(check_model(model))   # print() forces the plot to actually draw to the device
  dev.off()
  
  ###Save model results
  output_file <- paste0("All_data/comm_assembly_results/lme_results_SES_microenv/lme_SES_" ,traitlist[t], "_microenv_", sitelist[s], ".txt")
  sink(output_file)
  
  # ── 1.Trait ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  TRAIT\n")
  cat("===========================================\n")
  print(traitlist[t])
  cat("\n\n")
  
  
  # ── 1. Model Formula ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  MODEL FORMULA\n")
  cat("===========================================\n")
  print(formula(model))
  cat("\n\n")
  
  # ── 2. Summary Table ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  MODEL SUMMARY\n")
  cat("===========================================\n")
  print(summary(model))
  cat("\n\n")
  
  
  # ── 2. R squared ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  R SQUARED\n")
  cat("===========================================\n")
  print(r.squaredGLMM(model))
  cat("\n\n")
  
  
  # ── 3. ANOVA Table ────────────────────────────────────────────
  cat("===========================================\n")
  cat("  ANOVA TABLE\n")
  cat("===========================================\n")
  print(anova(model))
  cat("\n\n")
  
  # --- Close the sink ---
  sink()
  
  #Also save results in a list object so we can call it with quarto
  names(SES_microenv_summary)[[l]] <- paste(traitlist[t], sitelist[s], sep = "_")
  SES_microenv_summary[[l]] <- summary(model)
  
  names(SES_microenv_Rsq)[[l]] <- paste(traitlist[t], sitelist[s], sep = "_")
  SES_microenv_Rsq[[l]] <- r.squaredGLMM(model)
  
  names(SES_microenv_anova)[[l]] <- paste(traitlist[t], sitelist[s], sep = "_")
  SES_microenv_anova[[l]] <- anova(model)
  
  l = l+1

}}
#plot(checkmodel()) returns the following error: Check it out when you have internet
#Converting missing values (`NA`) into regular values currently not possible for variables of class `NULL`.




###=========================================###
###Loop to run SES~elevation for all traits####
###=====C2 NULL MODEL, pool = site=========####
###=========================================###
#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/SES_RQ_weighted_cells_C2_poolsite.csv", row.names = 1) |> 
  rename(Cell_ID = cellref)

#import microenvironmental data containing x and y coordinates
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA))

##Combine SES and environmental data
comb <- env |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  inner_join(cell_ses, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid), #each grid must have a unique id 
         elevation = as.factor(elevation), 
         grid = as.factor(grid)) |> 
  rename(x_coord = ncolumn, 
         y_coord = row)

#Run the loop for all traits
traitlist <- c("log_Height", "log_LDMC", "log_LA", "log_SLA", "Height_cm", "LDMC", "Leaf_area_mm2", "SLA")
#lists to store results in
C2_SES_ele_summary <- vector(mode= "list", length = length(traitlist))
names(C2_SES_ele_summary) = traitlist

C2_SES_ele_Rsq<- vector(mode= "list", length = length(traitlist))
names(C2_SES_ele_Rsq) = traitlist

C2_SES_ele_anova<- vector(mode= "list", length = length(traitlist))
names(C2_SES_ele_anova) = traitlist

C2_SES_ele_cld<- vector(mode= "list", length = length(traitlist))
names(C2_SES_ele_cld) = traitlist

for (t in 1:length(traitlist)) {
  modeldat <-  comb |> 
    filter(trait == traitlist[t]) |> 
    drop_na()
  
  
  model<- lme(SES ~ elevation ,
              random = ~1|grid, 
              correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
              data = modeldat) #only gaussian family possible
  
  ###Save check_model plot
  plot_file <- paste0("All_data/comm_assembly_results/C2_null_model_results/checkmodel_SES_elevation/checkmodel_lme_SES_", traitlist[t], "_elevation.png")
  png(plot_file, width = 1600, height = 1200, res = 150)
  print(check_model(model))   # print() forces the plot to actually draw to the device
  dev.off()
  
  ###Save model results
  output_file <- paste0("All_data/comm_assembly_results/C2_null_model_results/lme_results_SES_elevation/lme_SES_" ,traitlist[t], "_elevation_results.txt")
  sink(output_file)
  
  # ── 1.Trait ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  TRAIT\n")
  cat("===========================================\n")
  print(traitlist[t])
  cat("\n\n")
  
  
  # ── 1. Model Formula ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  MODEL FORMULA\n")
  cat("===========================================\n")
  print(formula(model))
  cat("\n\n")
  
  # ── 2. Summary Table ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  MODEL SUMMARY\n")
  cat("===========================================\n")
  print(summary(model))
  cat("\n\n")
  
  
  # ── 2. R squared ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  R SQUARED\n")
  cat("===========================================\n")
  print(r.squaredGLMM(model))
  cat("\n\n")
  
  
  # ── 3. ANOVA Table ────────────────────────────────────────────
  cat("===========================================\n")
  cat("  ANOVA TABLE\n")
  cat("===========================================\n")
  print(anova(model))
  cat("\n\n")
  
  # ── 4. EMmeans Table ──────────────────────────────────────────
  cat("===========================================\n")
  cat("  ESTIMATED MARGINAL MEANS (emmeans)\n")
  cat("===========================================\n")
  em_model <- emmeans(model, specs = "elevation", type = "response")
  comp_letters <-cld(em_model, Letters = letters, adjust = "Tukey", sort = FALSE)
  print(comp_letters)
  cat("\n")
  
  # --- Close the sink ---
  sink()
  
  #Also save results in a list object so we can call it with quarto
  C2_SES_ele_summary[[t]] <- summary(model)
  C2_SES_ele_Rsq[[t]] <- r.squaredGLMM(model)
  C2_SES_ele_anova[[t]] <- anova(model)
  C2_SES_ele_cld[[t]] <- comp_letters
  
}



###===========================================================###
###Loop to run SES~microenvironmental variables for all traits####
###=============C2 NULL MODEL, pool = site====================####
###===========================================================###
#import ses
cell_ses_poolsite <- read.csv("All_data/comm_assembly_results/SES_RQ_weighted_cells_C2_poolsite.csv", row.names = 1) |> 
  rename(Cell_ID = cellref) 

#import microenvironmental data
env2 <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
  #variables we are interested in
  select(Cell_ID:row, rock_cover, northness, soil_moisture_adj_campaign2, 
         soil_depth_CV, mean_soil_depth, slope_height) |> 
  #add elevation variables
  mutate(elevation = case_when(site == "GG" ~ 2000, 
                               site == "WH" ~ 2500, 
                               site == "BK" ~ 3000,
                               .default = NA))

#import remote sensing derived variables
rms <- read.csv("All_data/clean_data/Environmental data/Zonal_stats_all.csv") |> 
  select(CELL_ID, STD) |> 
  rename(Cell_ID = CELL_ID)

#import interpolated microclimate indices
micro_idw <- read.csv("All_data/clean_data/Environmental data/Imke_microclimate_indices_idw_interpolated.csv", row.names = 1)


##Combine SES and environmental data
comb2 <- env2 |> 
  #join to microclimate indices |> 
  full_join(micro_idw, by = "Cell_ID") |> 
  #join to remote sensing data |> 
  full_join(rms, by = "Cell_ID") |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cell_ses_poolsite, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)


###Loop starts here####
#Run the loop for all traits
traitlist <- c("log_Height", "log_LDMC", "log_LA", "log_SLA", "Height_cm", "LDMC", "Leaf_area_mm2", "SLA")
sitelist <- c("GG", "WH", "BK")
#lists to store results in
C2_SES_microenv_summary <- vector(mode= "list", length = length(traitlist) * length(sitelist))

C2_SES_microenv_Rsq<- vector(mode= "list", length = length(traitlist) * length(sitelist))

C2_SES_microenv_anova<- vector(mode= "list", length = length(traitlist) * length(sitelist))

l = 1
for (t in 1:length(traitlist)) {
  for(s in 1:length(sitelist)) {
    modeldat <-  comb2 |> 
      filter(trait == traitlist[t], 
             site == sitelist[s]) |> 
      mutate(elevation = as.factor(elevation), 
             grid = as.factor(grid)) |> 
      drop_na()
    
    
    model<- lme(SES ~ mean_T1_growing_season + mean_moist_growing_season + rock_cover+ mean_soil_depth + STD,
                random = ~1|grid, 
                correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
                data = modeldat) #only gaussian family possible
    
    ###Save check_model plot
    plot_file <- paste0("All_data/comm_assembly_results/C2_null_model_results/checkmodel_SES_microenv/checkmodel_lme_SES_", traitlist[t], "_microenv_", sitelist[s], ".png")
    png(plot_file, width = 1600, height = 1200, res = 150)
    print(check_model(model))   # print() forces the plot to actually draw to the device
    dev.off()
    
    ###Save model results
    output_file <- paste0("All_data/comm_assembly_results/C2_null_model_results/lme_results_SES_microenv/lme_SES_" ,traitlist[t], "_microenv_", sitelist[s], ".txt")
    sink(output_file)
    
    # ── 1.Trait ──────────────────────────────────────────
    cat("===========================================\n")
    cat("  TRAIT\n")
    cat("===========================================\n")
    print(traitlist[t])
    cat("\n\n")
    
    
    # ── 1. Model Formula ──────────────────────────────────────────
    cat("===========================================\n")
    cat("  MODEL FORMULA\n")
    cat("===========================================\n")
    print(formula(model))
    cat("\n\n")
    
    # ── 2. Summary Table ──────────────────────────────────────────
    cat("===========================================\n")
    cat("  MODEL SUMMARY\n")
    cat("===========================================\n")
    print(summary(model))
    cat("\n\n")
    
    
    # ── 2. R squared ──────────────────────────────────────────
    cat("===========================================\n")
    cat("  R SQUARED\n")
    cat("===========================================\n")
    print(r.squaredGLMM(model))
    cat("\n\n")
    
    
    # ── 3. ANOVA Table ────────────────────────────────────────────
    cat("===========================================\n")
    cat("  ANOVA TABLE\n")
    cat("===========================================\n")
    print(anova(model))
    cat("\n\n")
    
    # --- Close the sink ---
    sink()
    
    #Also save results in a list object so we can call it with quarto
    names(C2_SES_microenv_summary)[[l]] <- paste(traitlist[t], sitelist[s], sep = "_")
    C2_SES_microenv_summary[[l]] <- summary(model)
    
    names(C2_SES_microenv_Rsq)[[l]] <- paste(traitlist[t], sitelist[s], sep = "_")
    C2_SES_microenv_Rsq[[l]] <- r.squaredGLMM(model)
    
    names(C2_SES_microenv_anova)[[l]] <- paste(traitlist[t], sitelist[s], sep = "_")
    C2_SES_microenv_anova[[l]] <- anova(model)
    
    l = l+1
    
  }}

