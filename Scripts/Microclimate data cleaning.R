##Cleaning microclimate data and getting it ready to use
library(tidyverse)
library(tidylog)
library(lubridate)
library(data.table)
library(patchwork)
library(openxlsx)

###create function to read in files
readdata <- function(i){
  nn <- sum(grepl(i, fi$file))
  
  if(nn > 1){
    
    fi2 <- fi %>% filter(grepl(i, fi$file))
    
    df2 <- data.frame()
    for(ii in fi2$file){
      print(ii)
      d <- fread(ii)
      
      d %>% select(V2,V3,V4,V5,V6,V7) -> d
      
      d %>% filter(!duplicated(.$V2, fromLast = T)) -> d
      
      df2 <- bind_rows(df2, d)
    }
    
    df2 %>% filter(!duplicated(.$V2, fromLast = T)) -> df2
    
    df2$plot <- fi[which(fi$file2 == ii),"plot"]
    
    df2 %>% mutate(across(V4:V6, ~as.numeric(gsub(",",".\\",.)))) -> df2
    
    df2 %>% mutate(V2 = ymd_hm(V2, tz = "UTC")) %>% 
      mutate(V2 = with_tz(V2, tzone = "Etc/GMT-2")) -> df2
    
    
    return(df2)
    
  } else {
    
    print(i)
    d <- fread(fi$file[fi$file == i])
    
    d %>% select(V2,V3,V4,V5,V6,V7) -> d
    
    d %>% filter(!duplicated(.$V2, fromLast = T)) -> d
    
    d %>% mutate(across(V4:V6, ~as.numeric(gsub(",",".\\",.)))) -> d
    
    d$plot <- fi[which(fi$file == i),"plot"]
    
    d %>% mutate(V2 = ymd_hm(V2, tz = "UTC")) %>% 
      mutate(V2 = with_tz(V2, tzone = "Etc/GMT-2")) -> d
    
    return(d)
    
  }
  
}


####WITSIESHOEK####
datadir <- "All_data/raw_microclimate_data/Witsieshoek/Witsies_Tomst_data_final_reading_Feb2024"

# read names x tomst id table
ids <- read_csv(paste0(datadir, "/", "tomst_ids.csv")) %>% 
  rename(date_installed = `Date installed`, 
         date_removed = `Date removed`) %>%
  mutate(plot = paste0(Site, "_",Grid,Cell), 
         date_installed = dmy(date_installed), 
         date_removed = ymd(date_removed))


####BOKONG####
datadir <- "All_data/raw_microclimate_data/Bokong/Tomst_data_Final_reading_February2024"

# read names x tomst id table
ids <- read_excel(paste0(datadir, "/", "Bokong_Logger_notes.xlsx")) %>% 
  rename(date_installed = `Date installed`, 
         date_removed = `Date removed`) %>%
  mutate(plot = paste0(Site, "_",Grid,Cell), 
         date_installed = dmy(date_installed), 
         date_removed = ymd(date_removed))

# List logger data files to read
f <- list.files(datadir, pattern = "data_[0-9]", full.names = T, recursive = T)

fi <- data.frame(file = f)

fi$TMS <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",rev(strsplit(x, "/")[[1]])[1]), "_")[[1]][1])))

fi <- full_join(fi, ids, by = join_by(TMS)) %>% 
  mutate(plot = toupper(plot)) %>% 
  filter(!is.na(file))

fi$file2 <- gsub("_..csv", "", fi$file)

fi <- fi[order(fi$plot),]

# This should be zero
fi %>% filter(is.na(plot)) %>% nrow #there is one file (logger) That doesn't have a plot
#lets remove this row for now
fi <- fi %>% filter(!is.na(plot))

# read the Tomst data in
mylist <- lapply(fi$file, readdata)
bk_microclim <- rbindlist( mylist )

# Rename columns
bk_microclim <- bk_microclim %>% rename(datetime = V2,
              zone = V3,
              soil_temp = V4, #T1
              surface_temp = V5, #T2
              air_temp = V6, #T3
              moist = V7) 

bk_microclim <- bk_microclim %>% arrange(plot, datetime) 

# Remove implausible dates
mind <- ymd("2023-01-09") #the date of insertion is in a text file in the Bokong microclimate folder
maxd <- ymd("2024-02-04")
bk_microclim2 <- bk_microclim |> filter(between(datetime, mind, maxd))

