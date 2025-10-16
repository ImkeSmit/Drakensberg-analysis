# import community
library(tidyverse)
library(stringi)
library(readxl)
import_community <- function(metaTurfID, filepath){ #e.g., "All_data/raw_data/raw_threed_data/2025"
  
  #### COMMUNITY DATA ####
  ### Read in files
  files <- dir(path = filepath, pattern = "\\.xlsx$", full.names = TRUE, recursive = TRUE)
  
  #Function to read in meta data
  metaComm_raw <- map_df(set_names(files), function(file) {
    print(file)
    file %>% 
      excel_sheets() %>% 
      set_names() %>% 
      # exclude sheet to check data and taxonomy file
      discard(. %in% c("CHECK", "taxonomy", "empty", "Species list", "veg survey protocol")) %>% 
      map_df(~ read_xlsx(path = file, sheet = .x, n_max = 1, col_types = c("text", rep("text", 29))), .id = "sheet_name")
  }, .id = "file") %>%
    mutate(
      Date = dplyr::case_when(
        suppressWarnings(!is.na(as.numeric(Date))) ~ as.Date(as.numeric(Date), origin = "1899-12-30"),
        TRUE ~ lubridate::ymd(Date)
      )
    ) #there are two dates that are NA, these sheets have no dates on them
  
  # need to break the workflow here, otherwise tedious to find problems
  metaComm <- metaComm_raw %>% 
    select(sheet_name, Date, turfID, destSiteID, destBlockID, destPlotID, Recorder, Scribe)  %>% 
    
    # make date
    mutate(Date = ymd(Date),
           Year = year(Date),
           destBlockID = as.numeric(destBlockID),
           destPlotID = as.numeric(destPlotID)) %>%
    mutate(Year = case_when(is.na(Year)~ 2025, .default = Year)) %>% #add a year for the two sheets that don't have dates
    #There is a turf naming mistake in the data
    #turf 73_WN2I_158 should be 78_WN2I_158
    mutate(turfID = ifelse(turfID == "73_WN2I_158", "78_WN2I_158", turfID), 
           destSiteID = ifelse(turfID == "73_WN2I_158", "mid", destSiteID), 
           destBlockID =  ifelse(turfID == "73_WN2I_158", 10, destBlockID), 
           destPlotID = ifelse(turfID == "73_WN2I_158", 158, destPlotID)) %>%
    
    #join to metaTurfID
    left_join(metaTurfID , by = "turfID") %>% 
    mutate(destSiteID = coalesce(destSiteID.x, destSiteID.y),
           destBlockID = coalesce(destBlockID.x, destBlockID.y),
           destPlotID = coalesce(destPlotID.x, destPlotID.y)) %>%
  select(- c(destSiteID.x, destSiteID.y, destBlockID.x, destBlockID.y, destPlotID.x, destPlotID.y))
  
  
  # validate input
  # rules <- validator(Date = is.Date(Date),
  #                    Rec = is.character(Recorder),
  #                    Scr = is.character(Scribe))
  # out <- confront(metaComm, rules)
  # summary(out)
  
  
  # Function to read in data
  comm  <- map_df(set_names(files), function(file) {
    file %>% 
      excel_sheets() %>% 
      set_names() %>% 
      discard(. %in% c("CHECK", "taxonomy", "empty", "Species list", "veg survey protocol")) %>% 
      map_df(~ read_xlsx(path = file, sheet = .x, skip = 2, n_max = 61, col_types = "text"), .id = "sheet_name")
  }, .id = "file") %>% 
    select(file:Remark) %>% 
    rename("Cover" = `%`) %>% 
    #mutate(Year = 2025)
    mutate(Year = as.numeric(stri_extract_last_regex(file, "\\d{4}")))
  

  community <- metaComm %>% 
    left_join(comm, by = c("sheet_name", "Year")) %>% 
    select(origSiteID:origPlotID, destSiteID:turfID, warming:Nlevel, Date, Year, Species:Cover, Recorder, Scribe, Remark, file)
  
    
}

