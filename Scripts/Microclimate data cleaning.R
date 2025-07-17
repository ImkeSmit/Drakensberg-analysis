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

datadir <- ""
datadir <- "All_data/raw_microclimate_data/Bokong/Tomst_data_Final_reading_February2024"

#################################################################################3

# read names x tomst id table

ids <- read_excel(paste0(datadir, "/", "Bokong_Logger_notes.xlsx")) %>% 
  mutate(plot = paste0(Site, "_",Grid,Cell))

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


# read the Tomst data in
mylist <- lapply(fi$file, readdata)
df <- rbindlist( mylist )