####Script to clean Micro-Climb FT data####
library(tidyverse)
library(tidylog)
library(readxl)
library(janitor)
library(ggplot2)
library(openxlsx)

#load name standardising function:
standardise_names <- function(data, #dataframe containing species names that need to be corrected
                              data_species_column, #column in data contain the speciesnames, must be a string
                              #species names must be in the format genus_species.
                              #I.e. no capitals, with an underscore separating genus and specific epithet
                              naming_system,  #dataframe containing the old names that need to be changed, and the names they should be changed to
                              correct_name,  #column in naming_system that has the correct name
                              synonym) #columns in naming_system that has the synonyms
{
  #add change tracker column to keep track of names changed
  data$change_tracker <- NA
  
  #remove names that do not have synonyms
  naming_system <- naming_system |> 
    filter(!is.na(synonym1))
  
  for (i in 1:nrow(data)) {
    old_name <- data[i, which(colnames(data) == data_species_column)]
    new_name <- NA
    
    found <- FALSE
    for (j in 1:nrow(naming_system)) { # looks whether species name should be corrected and replaces it with the new_heli_name_system in case
      #found <- grepl(old_name, naming_system[j, which(colnames(naming_system) %in% synonym)])
      found <- any(old_name == as.character(naming_system[j, synonym]))
      
      # if (is.na(found)){ # only runs if the species is missing
      #    found <- FALSE
      #  }
      
      if (TRUE %in% found){ # only runs if the species is a synonym
        new_name <- naming_system[j, which(colnames(naming_system) %in% correct_name)] # finds the true name of the species and saves it
        break
      }
    }
    
    if (TRUE %in% found) { # replaces the species in the trait database with the saved true name if "found" is "TRUE"
      data[i, which(colnames(data) == data_species_column)] <- new_name
      
      #add column to keep track of which names changed
      data[i, which(colnames(data) == "change_tracker")] <- paste0(old_name, " -> ", new_name)
      
    }
  }#end loop through rows
  return(data)
}#end function

####GOLDEN GATE####
GG <- read.xlsx("All_data/raw_trait_data/GG_dataset_Functional_traits.xlsx") |> 
  select(Ref.number:Notes) |> 
  mutate(Cell = str_split_i(Grid_Cell, "_", 2)) #this creates 2 NA's
#the NA's are created because the Grid_Cell column does not follow the same format as the others
#fix it manually
GG[which(GG$Grid_Cell == "G7B14"), which(colnames(GG) == "Cell")] <- "B14"
GG[which(GG$Grid_Cell == "G4-c16"), which(colnames(GG) == "Cell")] <- "C16"

#clean up column names
names(GG) <- gsub("\\.", "_", names(GG)) #replace full stops in column names
names(GG) <- gsub("[()]", "", names(GG)) #remove brackets

#Rename columns, and drop the Grid_Cell_column
GG <- GG |> 
  select(!c(Grid_Cell, `Dry/Wet`)) |> 
  rename(Taxon = Species, 
         Total_leaf_area_cm2 = Total_Area_cm2, 
         Chlorophyll_mg_per_m2 = 'Chlor_mg/m2', 
         Average_leaf_area_cm2 = Average_Area_cm2, 
         Scan_name = Image_number, 
         Width_at_tear_mm = Width_mm, 
         Sample_ID = Ref_number) |> 
  mutate(Taxon = str_to_lower(Taxon), #change speciesnames to lower case and replace spaces with underscores
         Taxon = str_squish(Taxon),
         Taxon = str_replace_all(Taxon, " ", "_")) |> 
  mutate(Site = "GG", #add site variable
         Sample_ID = paste0(Site, Sample_ID)) 


#Are there NA's where there shouldn't be?
GG[which(is.na(GG$Grid)), ]
GG[which(is.na(GG$Cell)), ]
GG[which(is.na(GG$Taxon)), ]

#Check that all trait values are numeric
summary(GG)

unique(GG$FTT_N) #there are NA values, they will get NA
GG$FTT_N <- as.numeric(GG$FTT_N)

unique(GG$Width_at_tear_mm) #There are values of BA and NA
GG[which(GG$Width_at_tear_mm == "BA"), ] #they will get NA
GG$Width_at_tear_mm <- as.numeric(GG$Width_at_tear_mm)

#standardise names
name_trail <- read.xlsx("All_data/clean_data/micro_climb_ALL_names_editing.xlsx", sheet = "editing")
GG_clean_names <- standardise_names(data = GG, data_species_column = "Taxon", 
                                      naming_system = name_trail, correct_name = "taxon", 
                                    synonym = c("synonym1", "synonym2", "synonym3"))
unique(GG_clean_names$change_tracker) #all in order


###Let's look at the notes
GG_probs <- GG_clean_names |> 
  filter(!is.na(Area_Notes) | !is.na(Notes))

  
####WITSIESHOEK####
WH <- read.xlsx("All_data/raw_trait_data/WH_Functional_traits_dataset.xlsx", sheet = "Data_entry") |> 
  mutate(Cell = paste0(Column, Row), 
         Site = "WH") |> 
  select(!c(Date, Column, Row))

#clean up column names
names(WH) <- gsub("\\.", "_", names(WH)) #replace full stops in column names
names(WH) <- gsub("[()]", "", names(WH)) #remove brackets

WH <- WH |> #try to standardise column names between the datasets
  rename(Taxon = Species, 
         Wet_mass_mg = Wet_Mass, 
         Dry_mass_mg = Dry_Mass, 
         Chlorophyll_mg_per_m2 = Chlorophyll, 
         Total_leaf_area_cm2 = Total_Leaf_Area, 
         Average_leaf_area_cm2  = Average_Leaf_Area, 
         Scan_name = Scan_Name, 
         Sample_ID = Envelope_Number, 
         Width_at_tear_mm = Width_mm, 
         Notes2 = X18) |> 
  mutate(Taxon = str_to_lower(Taxon), #change speciesnames to lower case and replace spaces with underscores
         Taxon = str_squish(Taxon),
         Taxon = str_replace_all(Taxon, " ", "_"), 
         Sample_ID = paste0(Site, Sample_ID)) 

#Are there NA's where there shouldn't be?
WH[which(is.na(WH$Grid)), ]
WH[which(is.na(WH$Cell)), ]
WH[which(is.na(WH$Taxon)), ]

#Check that all trait values are numeric
summary(WH)
unique(WH$Height_cm) #there is a value of "S"
WH[which(WH$Height_cm == "S") , ] #This will get an NA when we change it to numeric
WH$Height_cm <- as.numeric(WH$Height_cm)

unique(WH$Width_at_tear_mm) #there are NA and " NA" values, they will get NA values
WH[which(is.na(WH$Width_at_tear_mm)) , ]
WH[which(WH$Width_at_tear_mm == " NA") , ]
WH$Width_at_tear_mm <- as.numeric(WH$Width_at_tear_mm)

unique(WH$Dry_mass_mg) #there are scan names in here
WH[grep("jpg", WH$Dry_mass_mg), ] #these will get NA, what are the numbers in the notes column??
WH$Dry_mass_mg <- as.numeric(WH$Dry_mass_mg)

#now we can calculate SLA
WH <- WH |> 
  mutate(SLA = Total_leaf_area_cm2/Dry_mass_mg)

#standardise names
name_trail <- read.xlsx("All_data/clean_data/micro_climb_ALL_names_editing.xlsx", sheet = "editing")
WH_clean_names <- standardise_names(data = WH, data_species_column = "Taxon", 
                                    naming_system = name_trail, correct_name = "taxon", 
                                    synonym = c("synonym1", "synonym2", "synonym3"))
unique(WH_clean_names$change_tracker) #all in order


####BOKONG####
BK <- read.xlsx("All_data/raw_trait_data/FT_Bokong_19Nov.xlsx", sheet = "FT measurements") |> 
  select(!Date) |> 
  mutate(Site = "BK")

#clean up column names
names(BK) <- gsub("\\.", "_", names(BK)) #replace full stops in column names
names(BK) <- gsub("[()]", "", names(BK)) #remove brackets

BK <- BK |> 
  rename(Taxon = Species, 
         Chlorophyll_mg_per_m2 = Chlorofil, 
         Scan_name = Image_reference, 
         FTT_N = FTT) |> 
  mutate(Taxon = str_to_lower(Taxon), #change speciesnames to lower case and replace spaces with underscores
         Taxon = str_squish(Taxon),
         Taxon = str_replace_all(Taxon, " ", "_"), 
         Sample_ID = paste0(Site, Sample_ID))

#Check that all trait values are numeric
summary(BK)
unique(BK$Chlorophyll_mg_per_m2) #there is a value of " ", NA, NA-no signal
BK[which(BK$Chlorophyll_mg_per_m2 == " ") , ] #This will get an NA when we change it to numeric
BK[which(BK$Chlorophyll_mg_per_m2 == "NA - no signal") , ]
BK[which(is.na(BK$Chlorophyll_mg_per_m2)) , ]
BK$Chlorophyll_mg_per_m2 <- as.numeric(BK$Chlorophyll_mg_per_m2)

unique(BK$FTT_N)#there are values of "cannot take measurement" These will get NA
BK[which(BK$FTT_N == "Cannot take measurement") , ]
BK$FTT_N <- as.numeric(BK$FTT_N)

unique(BK$Width_at_tear_mm) #again values of cannot take measurement
BK[which(BK$Width_at_tear_mm == "Cannot take measurement") , ]
BK$Width_at_tear_mm <- as.numeric(BK$Width_at_tear_mm)

#standardise names
name_trail <- read.xlsx("All_data/clean_data/micro_climb_ALL_names_editing.xlsx", sheet = "editing")
BK_clean_names <- standardise_names(data = BK, data_species_column = "Taxon", 
                                    naming_system = name_trail, correct_name = "taxon", 
                                    synonym = c("synonym1", "synonym2", "synonym3"))
unique(BK_clean_names$change_tracker) #all in order



####Checking Trait values####
#Combine trait data of all three sites
FT_allsites <- GG_clean_names |> 
  bind_rows(WH_clean_names) |> 
  bind_rows(BK_clean_names)

#make a copy to apply corrections to
FT_checked <- FT_allsites

FT_long <- FT_allsites |> 
  pivot_longer(cols = Wet_mass_mg:SLA, names_to = "Trait", values_to = "Value")

#check the distributions of all variables
ggplot(FT_long) + 
  geom_histogram(aes(x = Value)) +
  facet_wrap(~Trait, scales = "free") #none seem to have outrageous values
#but lets check what species are in the tails of the distribution

# Define tail threshold (e.g., 5%)
tail_threshold <- 0.05

# Function to identify tail records
get_trait_tails <- function(trait_vector, threshold = 0.05) {
  # Remove NAs for quantile calculation
  lower_cutoff <- quantile(trait_vector, probs = threshold, na.rm = TRUE)
  upper_cutoff <- quantile(trait_vector, probs = 1 - threshold, na.rm = TRUE)
  
  tail_status <- sapply(trait_vector, function(x) {
    if (is.na(x)) {
      return(NA)
    } else if (x < lower_cutoff) {
      return("low_tail")
    } else if (x > upper_cutoff) {
      return("high_tail")
    } else {
      return("middle")
    }
  })
  
  return(tail_status)
}


FT_allsites$LA_tail <- get_trait_tails(FT_allsites$Average_leaf_area_cm2, tail_threshold)
LA_outlier <- FT_allsites[which(FT_allsites$LA_tail == "high_tail"), ] #those with very high LA are all large leaved sp

FT_allsites$Chlor_tail <- get_trait_tails(FT_allsites$Chlorophyll_mg_per_m2, tail_threshold)
Chlor_outlier <- FT_allsites[which(FT_allsites$Chlor_tail == "low_tail"), ] #all seem in order, no notes about sp with low chlorophyll readings

FT_allsites$Dry_mass_tail <- get_trait_tails(FT_allsites$Dry_mass, tail_threshold)
Dry_mass_outlier <- FT_allsites[which(FT_allsites$Dry_mass_tail == "high_tail"), ] #heavy leaves make sense for these sp


FT_allsites$FTT_tail <- get_trait_tails(FT_allsites$FTT_N, tail_threshold)
FTT_outlier <- FT_allsites[which(FT_allsites$FTT_tail == "high_tail"), ] 
#two have notes of FTT-split in GG, what does it mean? Otherwise nothing alarming

FT_allsites$H_tail <- get_trait_tails(FT_allsites$Height_cm, tail_threshold)
H_outlier <- FT_allsites[which(FT_allsites$H_tail == "high_tail"), ] 
#two have notes of FTT-split in GG, what does it mean? Otherwise nothing alarming

FT_allsites$SLA_tail <- get_trait_tails(FT_allsites$SLA, tail_threshold)
SLA_outlier <- FT_allsites[which(FT_allsites$SLA_tail == "high_tail"), ] 
#many succulens which make sense, will have to test this another way too

FT_allsites$Thickness_tail <- get_trait_tails(FT_allsites$Thickness_mm, tail_threshold)
Thickness_outlier <- FT_allsites[which(FT_allsites$Thickness_tail == "high_tail"), ] 
Thickness_problems <- Thickness_outlier[order(Thickness_outlier$Thickness_mm), ][c(309,310), which(colnames(Thickness_outlier) == "Sample_ID")]
#The two highest values are a little wack
#Tenaxia disticha with thickness = 26mm
#Gazania krebsiana with thickness = 95mm
#remove these values, correct in FT_checked
FT_checked[which(FT_checked$Sample_ID %in% c(Thickness_problems)), which(colnames(FT_checked) == "Thickness_mm")] <- NA

FT_allsites$TLA_tail <- get_trait_tails(FT_allsites$Total_leaf_area_cm2, tail_threshold)
TLA_outlier <- FT_allsites[which(FT_allsites$TLA_tail == "high_tail"), ] #all looks fine
#the largest areas come from very large leaves

FT_allsites$Wet_mass_tail <- get_trait_tails(FT_allsites$Wet_mass_mg, tail_threshold)
Wet_mass_outlier <- FT_allsites[which(FT_allsites$Wet_mass_tail == "high_tail"), ] #all looks fine
#The leaves with the highest masses make sense

FT_allsites$Width_tail <- get_trait_tails(FT_allsites$Width_at_tear_mm, tail_threshold)
Width_outlier <- FT_allsites[which(FT_allsites$Width_tail == "high_tail"), ] 
#There is an elionurus muticus with a width of 7.86, that is wack
#the 3 highest widths are improbable
Width_outlier[which(Width_outlier$Taxon == "elionurus_muticus") , ]
Width_outlier[which(Width_outlier$Taxon == "pentameris_setifolia") , ] #this one is probably fine

#Fix the elionurus problem
Width_problem1 <- Width_outlier[which(Width_outlier$Taxon == "elionurus_muticus") , which(colnames(Width_outlier) == "Sample_ID")]
FT_checked[which(FT_checked$Sample_ID == Width_problem1), which(colnames(FT_checked) == "Width_at_tear_mm")] <- NA

#fix the 3 highest widths
Width_problem2 <- Width_outlier[order(Width_outlier$Width_at_tear_mm), ][c(131:133), which(colnames(Width_outlier) == "Sample_ID")]
FT_checked[which(FT_checked$Sample_ID %in% Width_problem2), which(colnames(FT_checked) == "Width_at_tear_mm")] <- NA


###When finished, write to file
write.xlsx(GG_clean_names, "All_data/clean_data/Golden_Gate_clean_traits.xlsx")
write.xlsx(WH_clean_names, "All_data/clean_data/Witsieshoek_clean_traits.xlsx")
write.xlsx(BK_clean_names, "All_data/clean_data/Bokong_clean_traits.xlsx")
