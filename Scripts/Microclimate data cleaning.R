##Cleaning microclimate data and getting it ready to use
library(tidyverse)
library(tidylog)
library(lubridate)
library(data.table)
library(patchwork)
library(openxlsx)
library(readxl)

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
    df2$site <- fi[which(fi$file2 == ii),"site"]
    df2$grid <- fi[which(fi$file2 == ii),"grid"]
    df2$cell <- fi[which(fi$file2 == ii),"Cell"]
    
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
    d$site <- fi[which(fi$file == i),"site"]
    d$grid <- fi[which(fi$file == i),"grid"]
    d$cell <- fi[which(fi$file == i),"Cell"]
    
    d %>% mutate(V2 = ymd_hm(V2, tz = "UTC")) %>% 
      mutate(V2 = with_tz(V2, tzone = "Etc/GMT-2")) -> d
    
    return(d)
    
  }
  
}


####GOLDEN GATE####
datadir <- "All_data/raw_data/raw_microclimate_data"

# read names x tomst id table
ids <- read.delim(paste0(datadir, "/", "tomst_ids.csv"), sep = ";") |>  
  filter(site == "goldengate") |> 
  rename(Cell = plot) |> 
  mutate(Cell = toupper(Cell),
         site = toupper(site),
         plot = paste0(site, "_",grid,Cell))

# List logger data files to read
f <- list.files(paste0(datadir, "/", "GoldenGate/Golden_Gate_Tomsp_Final_reading_Feb2024"), pattern = "data_[0-9]", full.names = T, recursive = T)

fi <- data.frame(file = f)

fi$TMS <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",rev(strsplit(x, "/")[[1]])[1]), "_")[[1]][1])))

fi <- full_join(fi, ids, by = join_by(TMS)) %>% 
  mutate(plot = toupper(plot)) %>% 
  filter(!is.na(file))

fi$file2 <- gsub("_..csv", "", fi$file)

fi <- fi[order(fi$plot),]

# This should be zero
fi %>% filter(is.na(plot)) %>% nrow 
#there is one logger file that doesn't have info on which cell it was in. remove this file
fi <- fi |> 
  filter(!is.na(plot))

# read the Tomst data in
mylist <- lapply(fi$file, readdata)
gg_microclim <- rbindlist( mylist )

# Rename columns
gg_microclim <- gg_microclim %>% rename(datetime = V2,
                                        zone = V3,
                                        soil_temp = V4, #T1
                                        surface_temp = V5, #T2
                                        air_temp = V6, #T3
                                        moist = V7) 

gg_microclim <- gg_microclim %>% arrange(plot, datetime) 

# Remove implausible dates
mind <- ymd("2023-01-01") #can't find the installation or extraction date. Assume these dates
maxd <- ymd("2024-02-28")
gg_microclim2 <- gg_microclim |> filter(between(datetime, mind, maxd)) |> 
  distinct(datetime, plot, .keep_all = T) #remove duplicate rows

#look at the duplicate rows for quality control
gg_microclim |> 
  group_by_all() |> 
  filter(n() >1)


####WITSIESHOEK####
datadir <- "All_data/raw_data/raw_microclimate_data"

# read names x tomst id table
ids <- read.delim(paste0(datadir, "/", "tomst_ids.csv"), sep = ";") |>  
  filter(site == "witsieshoek") |> 
  rename(Cell = plot) |> 
  mutate(Cell = toupper(Cell),
         site = toupper(site),
         plot = paste0(site, "_",grid,Cell))

# List logger data files to read
f <- list.files(paste0(datadir, "/", "Witsieshoek/Witsies_Tomst_data_final_reading_Feb2024"), pattern = "data_[0-9]", full.names = T, recursive = T)

fi <- data.frame(file = f)

fi$TMS <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",rev(strsplit(x, "/")[[1]])[1]), "_")[[1]][1])))

fi <- full_join(fi, ids, by = join_by(TMS)) %>% 
  mutate(plot = toupper(plot)) %>% 
  filter(!is.na(file))

fi$file2 <- gsub("_..csv", "", fi$file)

fi <- fi[order(fi$plot),]

# This should be zero
fi %>% filter(is.na(plot)) %>% nrow 

# read the Tomst data in
mylist <- lapply(fi$file, readdata)
wh_microclim <- rbindlist( mylist )

# Rename columns
wh_microclim <- wh_microclim %>% rename(datetime = V2,
                                        zone = V3,
                                        soil_temp = V4, #T1
                                        surface_temp = V5, #T2
                                        air_temp = V6, #T3
                                        moist = V7) 

wh_microclim <- wh_microclim %>% arrange(plot, datetime) 

# Remove implausible dates
mind <- ymd("2023-01-01") #can't find the installation or extraction date. Assume these dates
maxd <- ymd("2024-02-28")
wh_microclim2 <- wh_microclim |> filter(between(datetime, mind, maxd)) |> 
  distinct(datetime, plot, .keep_all = T)


####BOKONG####
datadir <- "All_data/raw_data/raw_microclimate_data"

# read names x tomst id table
ids <- read.delim(paste0(datadir, "/", "tomst_ids.csv"), sep = ";") |> 
  filter(site == "bokong") |> 
  rename(Cell = plot) |> 
  mutate(Cell = toupper(Cell),
         site = toupper(site),
         plot = paste0(site, "_",grid,Cell))

# List logger data files to read
f <- list.files(paste0(datadir, "/", "Bokong/Tomst_data_Final_reading_February2024"), pattern = "data_[0-9]", full.names = T, recursive = T)

fi <- data.frame(file = f)

fi$TMS <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",rev(strsplit(x, "/")[[1]])[1]), "_")[[1]][1])))

fi <- full_join(fi, ids, by = join_by(TMS)) %>% 
  mutate(plot = toupper(plot)) %>% 
  filter(!is.na(file))

fi$file2 <- gsub("_..csv", "", fi$file)

fi <- fi[order(fi$plot),]

# This should be zero
fi %>% filter(is.na(plot)) %>% nrow 

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
maxd <- ymd("2024-02-04") #from bokong logger notes file
bk_microclim2 <- bk_microclim |> filter(between(datetime, mind, maxd)) |> 
  distinct(datetime, plot, .keep_all = T) #remove duplicate entries

#look at the duplicate rows for quality control
bk_microclim |> 
  group_by_all() |> 
  filter(n() >1)

####Write the cleaned microclimate data to text files####
#Excel cna't open a table with this many rows, so we write it to a text file
write.table(bk_microclim2, file = "All_data/clean_data/Tomst_data/Bokong_tomst_data.txt", sep = "\t",
            row.names = TRUE, col.names = c(colnames(bk_microclim2)))

write.table(wh_microclim2, file = "All_data/clean_data/Tomst_data/Witsieshoek_tomst_data.txt", sep = "\t",
            row.names = TRUE, col.names = c(colnames(wh_microclim2)))

write.table(gg_microclim2, file = "All_data/clean_data/Tomst_data/GoldenGate_tomst_data.txt", sep = "\t",
            row.names = TRUE, col.names = c(colnames(gg_microclim2)))


