#Does community assembly mechanisms vary with scale and with elevation?#
library(tidyverse)
library(tidylog)
library(openxlsx)
library(ggplot2)
library(ggridges)

###Cell Scale####
#import SES at cell scale
cell_ses <- read.csv("All_data/comm_assembly_results/RQ_cells_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ "3000", #add elevation variable
                               grepl("WH", cellref) == T ~ "2500",
                               grepl("GG", cellref) == T ~ "2000",.default = NA))
cell_ses$elevation <- as.factor(cell_ses$elevation)  



ses_ridges <- cell_ses |> 
  ggplot(aes(x = SES, y = elevation, fill = elevation)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()


#models of ses~ elevation
hist(cell_ses$SES)
modeldat <- cell_ses[which(cell_ses$trait == "Height_cm"), ]
test <- lm(SES ~ elevation, data = modeldat)
summary(test)
anova(test)

#how badly are the assumptions violated?
plot(test)#ja very badly

#Let's try a GAM
#USe skew t distribution as recommended by chatgpt.Can handle skewed and heavy tailed data with nonzero values
test2 <- gamlss(SES ~ elevation, data = modeldat, family = ST1())
summary(test2)
plot(test2) #looks much better


###Grid scale####
grid_ses <- read.csv("All_data/comm_assembly_results/RQ_grids_C5_entire.csv", row.names = 1) |> 
  mutate(elevation = case_when(grepl("BK", cellref) == T ~ 3000, #add elevation variable
                               grepl("WH", cellref) == T ~ 2500,
                               grepl("GG", cellref) == T ~ 2000,.default = NA))


grid_ses_ridges <- grid_ses |> 
  mutate(elevation_char = as.character(elevation)) |> 
  ggplot(aes(x = SES, y = elevation_char, fill = elevation_char)) +
  geom_density_ridges(alpha = 0.5) +
  facet_wrap(~trait) +
  theme_classic()
