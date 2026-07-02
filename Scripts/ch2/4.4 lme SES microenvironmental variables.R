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


#import SES data
cell_ses <- read.csv("All_data/comm_assembly_results/SES_RQ_weighted_cells_C5_entire.csv", row.names = 1) |> 
  rename(Cell_ID = cellref)

#import microenvironmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv") |> 
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
comb <- env |> 
  #join to microclimate indices |> 
  full_join(micro_idw, by = "Cell_ID") |> 
  #join to remote sensing data |> 
  full_join(rms, by = "Cell_ID") |> 
  #join, one row in env matches many rows in cell_ses due to it containing ses of different traits
  full_join(cell_ses, by = "Cell_ID", relationship = "one-to-many") |>
  mutate(ncolumn = match(column, LETTERS[1:8]), 
         grid = paste0(site, grid)) |> #each grid must have a unique id 
  rename(x_coord = ncolumn, 
         y_coord = row)


###Check collinearity#### 
library(corrplot)
cordf <- comb |> 
  filter(trait == "Height_cm") |> #look at just one set of env data, it repeats for every trait
  select(mean_T1_growing_season, mean_moist_growing_season, STD, rock_cover, northness, mean_soil_depth, slope_height) |> 
  drop_na()
cormat<- cor(cordf)
#cormat[cormat > 0.7]
#cormat[cormat <-0.7]
#none are highly correlated
corrplot(cormat, type = "lower", method = "number")


###############################
#####SES ~ ELEVATION MODELS####
##C5 NULLMODEL, POOL = ENTIRE##


####SES HEIGHT####
Hdat <- comb |> 
  filter(trait == "log_Height") |> 
  mutate(elevation = as.factor(elevation), 
         grid = as.factor(grid)) |> 
  drop_na()


tic()
H_ele_mod<- lme(SES ~ elevation ,
            random = ~1|grid, 
            correlation = corSpher(form = ~ x_coord + y_coord|grid, nugget = TRUE), #spherical structure
            data = Hdat) #only gaussian family possible
toc()
summary(H_ele_mod)
anova(H_ele_mod)
emmeans(H_ele_mod, specs = "elevation")

#diagnostics
check_model(H_ele_mod)
#looks better!



###Write a loop for all models
traitlist <- c("log_Height", "log_SLA", "log_LDMC", "log_Thickness")

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
  
  output_file <- paste0("All_data/comm_assembly_results/lme_SES_" ,traitlist[t], "_elevation_results.txt")
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
  comp_letters <-cld(em_model, Letters = letters, adjust = "Tukey")
  print(comp_letters)
  cat("\n")
  
  # --- Close the sink ---
  sink()
  
}



####SES SLA####
#isolate SES of SLA
#leave heavy tail
SLAdat <- comb2 |> 
  filter(trait %in% c("SLA", NA), #also select cells which have no SES measurement. This is necessary to make the grid complete
         !is.na(SES)) |> 
  arrange(y_coord, x_coord) |> 
  mutate(trait = "SLA",  #give all records a trait
         grid = as.factor(paste0(site, grid))) |> 
  drop_na()

#descriptive stats
#how many cells
nrow(SLAdat) #2880
SLAdat |> group_by(site) |> 
  summarise(n = n())


tic()
tmod2<- lme(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                  zslope_height,
            random = ~1|grid, 
            correlation = corExp(form = ~ x_coord + y_coord | grid, nugget = TRUE),
            data = SLAdat)
toc()
summary(tmod2)
anova(tmod2)

check_model(tmod2)


#Compare against a model without spatial structure to see if it improves fit
tic()
tmod2_nonspat<- lme(SES ~ elevation + zrock_cover + znorthness + zsoil_moist + zsoil_depth + 
                      zslope_height,
                    random = ~1|grid, 
                    data = SLAdat) #only gaussian family possible
toc()

anova(tmod2_nonspat, tmod2)
#spatial model has lower AIC and is a significantly better fit than the nonspatial model


# Residual diagnostics
plot(tmod2) #looks pretty good
qqnorm(tmod2, ~ resid(., type = "normalized")) #pretty ok

# Optional: variogram of normalized residuals to visually check
# whether spatial autocorrelation has been adequately captured
plot(Variogram(tmod2, resType = "normalized"))
plot(Variogram(tmod2_nonspat, resType = "normalized"))




