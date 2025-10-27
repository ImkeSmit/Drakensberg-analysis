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


###What is going on in the cells with high SES values?
high <- cell_ses |> 
  slice_max(SES, n = 20)

high_plots <- cell_ses |> 
  slice_max(SES, n = 20) |> 
  distinct(cellref)

#occurrence data
drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))

#trait data
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) |> 
  pivot_longer(cols = c(Wet_mass_mg, Dry_mass_mg, Chlorophyll_mg_per_m2, Ft, Height_cm, 
                        Thickness_mm, Leaf_area_mm2, SLA, LDMC), names_to = "trait", values_to = "value")

#combine occurrence and trait data
FT_join <- drak |> 
  inner_join(FT, by = c("site", "grid", "cell", "taxon")) |> #inner join to only work with taxa that have trait data
  mutate(cellref = paste0(site, grid, cell)) |> 
  select(!c(column, row))

#create abundance matrix
abun_matrix <- FT_join |> 
  select(cellref, taxon, cover) |> 
  distinct(cellref, taxon, .keep_all = T) |> 
  ungroup() |> 
  arrange(taxon) |> 
  mutate(cover = ceiling(cover)) |> #change cover values to integer to use in null models
  pivot_wider(names_from = taxon, values_from = cover) 

abun_matrix <- as.data.frame(abun_matrix)
row.names(abun_matrix) <- abun_matrix$cellref
abun_matrix <- abun_matrix[, -1]

#replace NA values with 0
for(r in 1:nrow(abun_matrix)) {
  for(c in 1:ncol(abun_matrix)) {
    
    if(is.na(abun_matrix[r,c])) {
      abun_matrix[r,c] <- 0
    }
  }
}

high_mat <- abun_matrix[which(row.names(abun_matrix) %in% c(high_plots$cellref)) , ]
specnumber(high_mat) #all have between 2 and 6 species

all_abun <- as.data.frame(specnumber(abun_matrix))
all_abun$cellref <- row.names(all_abun)
row.names(all_abun) <- NULL

cell_ses2 <- cell_ses |> 
  left_join(all_abun, by = "cellref")


#GG4F18 has the highest SES, look at it's diversity
s1 <- abun_matrix[which(rownames(abun_matrix) == "GG4F18") ,, drop = F]
s2 <- s1[, -which(colSums(abun_matrix) == 0)]

abun_matrix |> 
  rownames_to_column(var = "cellref") |> 
  filter(cellref == "GG4F18") |> 
  select(where(~ sum(.) != 0))


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
