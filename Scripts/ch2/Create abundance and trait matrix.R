###This is the script used to create the abundance and trait matrix used in further analyses
#we only do this once to avoid inconsistencies in downstream analyses
library(openxlsx)
library(tidyverse)
library(tidylog)
library(vegan)

###Matrices for all life forms####
####Import community and trait data###
#occurrence data
drak <- read.xlsx("All_data/clean_data/micro_climb_occurrence.xlsx") |> 
  mutate(cell = paste0(column, row))

#trait data
FT <- read.xlsx("All_data/clean_data/micro-climb_traits.xlsx") |> 
  rename(taxon = Taxon, 
         site = Site, grid = Grid, cell = Cell) |> 
  pivot_longer(cols = c(Wet_mass_mg, Dry_mass_mg, Chlorophyll_mg_per_m2, Ft, Height_cm, 
                        Thickness_mm, Leaf_area_mm2, SLA, LDMC), names_to = "trait", values_to = "value")


#Get mean traits for species
mean_traits_all <- FT |> 
  filter(trait %in% c("Height_cm", "Leaf_area_mm2","Thickness_mm", "SLA", "LDMC")) |> 
  group_by(taxon, trait) |> 
  summarise(mean_trait = mean(value, na.rm = T)) |> 
  pivot_wider(names_from = trait, values_from = mean_trait) |> 
  ungroup() |> 
  arrange(taxon)

#Do inner join between trait and cover data to get sp that match between the two
FT_join <- drak |> 
  inner_join(mean_traits_all, by = "taxon") |> #inner join to only work with taxa that have trait data
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

#remove sites that have only one species 
no <- specnumber(abun_matrix)
abun_matrix <- abun_matrix[-which(no == 1), ]

#order species alphabetically
abun_matrix <- abun_matrix[, order(colnames(abun_matrix))]
colnames(abun_matrix)[79] <- "helichrysum_albo_brunneum" #this names causes problems in csv format

#Save abundance matrix:
write.csv(abun_matrix, "All_data/comm_assembly_results/abun_matrix.csv")


###Polish trait matrix:
#get row and columns names right
mean_traits <- as.data.frame(mean_traits_all)
row.names(mean_traits) <- mean_traits$taxon
mean_traits <- mean_traits[, -1]

#order rownames alphabetically
mean_traits <- mean_traits[order(row.names(mean_traits)), ]
#fix the helichrysum name to be the same as in abun_matrix
row.names(mean_traits)[which(row.names(mean_traits) == "helichrysum_albo-brunneum")] <- "helichrysum_albo_brunneum"

#remove species that aren't also in the abun_matrix
mean_traits <- mean_traits[which(row.names(mean_traits) %in% colnames(abun_matrix)) , ]

#save trait matrix
write.csv(mean_traits, "All_data/comm_assembly_results/mean_traits.csv")




####Matrices for forbs and graminoids separately####
#create life form classification
lf <- data.frame(taxon = colnames(abun_matrix), 
                 growth_form = NA)

graminoids <- c("andropogon_filiformes", "anthoxanthum_ecklonii","aristida_adscensionis", "aristida_junciformis","brachiaria_serrata","bulbostylis_humilis",
                "carex_zuluensis","cymbopogon_dieterlenii","cymbopogon_prolixus","cyperus_rupestris",
                "cyperus_semitrifidus","digitaria_monodactyla","digitaria_sanguinium","diheteropogon_filifolius",
                "ehrharta_longiglumous","elionurus_muticus","eragrostis_capensis","eragrostis_curvula",
                "eragrostis_plana","eragrostis_racemosa","eulalia_villosa","festuca_caprina","festuca_costata",
                "festuca_scabra",
                "ficinia_cinammomea","ficinia_stolonifera","harpochloa_falx","heteropogon_contortus","melinis_nerviglumis",
                "merxmuellera_drakensbergensis","microchloa_caffra","miscanthus_capensis","pentameris_airoides",
                "pentameris_exserta","pentameris_setifolia","poa_binata","rendlia_altera","sporobolus_centrifugus",
                "stiburus_alopecuroides","tenaxia_disticha","themeda_triandra","tristachya_leucothrix")

#classify each sp as a graminoid or not
for(i in 1:nrow(lf)) {
  sp <- lf[i,1]
  
  if(any(sp == graminoids)) {
    lf[i,2] <- "graminoid"
  }else {lf[i,2] <- "forb" }
}

#Now subset the abundance and trait matrices for forbs and graminoids

###graminoid abundance
abun_graminoids <- abun_matrix[ , graminoids]
#remove sites that have only one species 
no <- specnumber(abun_graminoids)
abun_graminoids <- abun_graminoids[-which(no < 2), ]

###graminoid traits
traits_graminoids <- mean_traits[graminoids, ]


###forb abundance
abun_forbs <- abun_matrix[ , -which(colnames(abun_matrix) %in% graminoids)]
#remove sites that have only one species 
no <- specnumber(abun_forbs)
abun_forbs <- abun_forbs[-which(no < 2), ]

###forb traits
traits_forbs <- mean_traits[-which(rownames(mean_traits) %in% graminoids), ]

#write to csv
write.csv(abun_graminoids, "All_data/comm_assembly_results/abun_graminoids.csv")
write.csv(traits_graminoids, "All_data/comm_assembly_results/traits_graminoids.csv")
write.csv(abun_forbs, "All_data/comm_assembly_results/abun_forbs.csv")
write.csv(traits_forbs, "All_data/comm_assembly_results/traits_forbs.csv")

