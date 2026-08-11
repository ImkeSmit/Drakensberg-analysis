####Function to standardise names according to a key with synonyms####
#create function to change incorrect names
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
  
  data <- tibble(data)
  
  #remove names that do not have synonyms
  naming_system <- naming_system |> 
    filter(!is.na(synonym1))
  
  for (i in 1:nrow(data)) {
    old_name <- data[i, which(colnames(data) == data_species_column)]
    new_name <- NA
    
    found <- FALSE
    for (j in 1:nrow(naming_system)) { # looks whether species name should be corrected and replaces it with the new name
      #found <- grepl(old_name, naming_system[j, which(colnames(naming_system) %in% synonym)])
      found <- any(old_name$data_species_column == as.character(naming_system[j, synonym]))
      
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