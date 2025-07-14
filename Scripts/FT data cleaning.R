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
  select(!Grid_Cell) |> 
  mutate(Site = "GG") |>  #add site variable
  rename(Taxon = Species, 
         Total_leaf_area_cm2 = Total_Area_cm2, 
         Chlorophyll_mg_per_m2 = 'Chlor_mg/m2', 
         Average_leaf_area_cm2 = Average_Area_cm2, 
         Scan_name = Image_number, 
         Width_at_tear_mm = Width_mm) |> 
  mutate(Taxon = str_to_lower(Taxon), #change speciesnames to lower case and replace spaces with underscores
         Taxon = str_squish(Taxon),
         Taxon = str_replace_all(Taxon, " ", "_"))


#Are there NA's where there shouldn't be?
GG[which(is.na(GG$Grid)), ]
GG[which(is.na(GG$Cell)), ]
GG[which(is.na(GG$Taxon)), ]

#Check that all trait values are numeric
summary(GG)

unique(GG$FTT_N) #there are NA values, they will get NA
GG$FTT_N <- as.numeric(GG$FTT_N)

unique(GG$Width_mm) #There are values of BA and NA
GG[which(GG$Width_mm == "BA"), ] #they will get NA
GG$Width_mm <- as.numeric(GG$Width_mm)

#standardise names
name_trail <- read.xlsx("All_data/clean_data/micro_climb_ALL_names_editing.xlsx", sheet = "editing")
GG_clean_names <- standardise_names(data = GG, data_species_column = "Taxon", 
                                      naming_system = name_trail, correct_name = "taxon", 
                                    synonym = c("synonym1", "synonym2", "synonym3"))
unique(GG_clean_names$change_tracker) #all in order

  
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
         Width_at_tear_mm = Width_mm) |> 
  mutate(Taxon = str_to_lower(Taxon), #change speciesnames to lower case and replace spaces with underscores
         Taxon = str_squish(Taxon),
         Taxon = str_replace_all(Taxon, " ", "_"))

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
         Taxon = str_replace_all(Taxon, " ", "_"))

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
