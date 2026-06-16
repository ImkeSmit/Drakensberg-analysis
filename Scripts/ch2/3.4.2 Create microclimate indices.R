####Identify growing season from TOMSt data####
###And calculate mean temperature for the growing season####
library(tidyverse)
library(tidylog)
library(lubridate)

#T1 = soil temp
#T3 = air temp


#import raw tomst data 
tms <- read.csv("C:\\Users\\imke6\\Documents\\PhD 2025\\Ch2 niche modelling\\All_data\\Clean_data\\microclimate\\Pekka_cleaned_data\\TOMST data\\transfer_402058_files_56a24364\\tomst_data_imputed.csv") |> 
  rename(cellref = site) |> 
  filter(!is.na(cellref)) |> 
  mutate(site = str_split_i(cellref, "_", 1),
         site = case_when(site == "BOK" ~ "BK", 
                          site == "GOL" ~ "GG", 
                          site == "WIT" ~ "WH"),
         grid = str_split_i(cellref, "_", 2),
          cell = toupper(str_split_i(cellref, "_", 3)), 
         Cell_ID = paste0(site, "_G", grid, "_", cell))
  
  
mean_daily_temp <- tms |> 
  mutate(date = as_date(datetime)) %>%
  group_by(site, Cell_ID, date) %>%
  summarise(mean_daily_T1 = mean(T1, na.rm = TRUE), )


