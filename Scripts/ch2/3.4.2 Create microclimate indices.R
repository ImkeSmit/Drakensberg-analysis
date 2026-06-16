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
  

#mean daily soil temp per logger  
mean_daily_temp <- tms |> 
  mutate(date = as_date(datetime)) |> 
  group_by(site, Cell_ID, date)|> 
  summarise(mean_daily_T1 = mean(T1, na.rm = TRUE)) |> 
  ungroup()

#min daily soil temp per logger
min_daily_temp <- tms |> 
  mutate(date = as_date(datetime))|> 
  group_by(site, Cell_ID, date) |> 
  summarise(min_daily_T1 = min(T1, na.rm = TRUE)) |> 
  ungroup()


##plot
ggplot(mean_daily_temp, aes( x = date, y = mean_daily_T1)) +
  geom_line()+
  facet_wrap(~site, nrow = 3, ncol = 1)

##plot
ggplot(min_daily_temp, aes( x = date, y = min_daily_T1)) +
  geom_line()+
  facet_wrap(~site, nrow = 3, ncol = 1)


###USING MEAN TEMP####
###Period where mean annual temp is consistently above 5 degrees
#first date after June 21 where all logers at a site had a mean  temp above 5
start_grw <- mean_daily_temp |>
  filter(date > ymd("2023-06-21")) |>           # only dates after 21 Jun 2023
  group_by(site, date) %>%
  summarise(all_above_5 = all(mean_daily_T1 > 5),    # TRUE only if every Cell_ID > 5
            .groups = "drop") |> 
  filter(all_above_5) |> #only keep dates where all loggers had temps above 5
  group_by(site) |> 
  slice_min(date, n = 1) |>   #earliest date per site
  rename(start_grw_date = date) |> 
  select(!all_above_5)

##latest date before 21 June where all loggers at a site had a mean temp above 5
end_grw <- mean_daily_temp |>
  filter(date < ymd("2023-06-21")) |>           # only dates before 21 Jun 2023
  group_by(site, date) %>%
  summarise(all_above_5 = all(mean_daily_T1 > 5),    # TRUE only if every Cell_ID > 5
            .groups = "drop") |> 
  filter(all_above_5) |> #only keep dates where all loggers had temps above 5
  group_by(site) |> 
  slice_max(date, n = 1) |>   #latest date per site
  rename(end_grw_date = date) |> 
  select(!all_above_5)

growing_season_mean_temps <- start_grw |> 
  left_join(end_grw, by = "site") #do not consider the years, only the days!



###USING MIN TEMP####
###Period where mean annual temp is consistently above 5 degrees
#first date after June 21 where all logers at a site had a mean  temp above 5
start_grw <- min_daily_temp |>
  filter(date > ymd("2023-06-21")) |>           # only dates after 21 Jun 2023
  group_by(site, date) %>%
  summarise(all_above_5 = all(min_daily_T1 > 5),    # TRUE only if every Cell_ID > 5
            .groups = "drop") |> 
  filter(all_above_5) |> #only keep dates where all loggers had temps above 5
  group_by(site) |> 
  slice_min(date, n = 1) |>   #earliest date per site
  rename(start_grw_date = date) |> 
  select(!all_above_5)

##latest date before 21 June where all loggers at a site had a mean temp above 5
end_grw <- min_daily_temp |>
  filter(date < ymd("2023-06-21")) |>           # only dates before 21 Jun 2023
  group_by(site, date) %>%
  summarise(all_above_5 = all(min_daily_T1 > 5),    # TRUE only if every Cell_ID > 5
            .groups = "drop") |> 
  filter(all_above_5) |> #only keep dates where all loggers had temps above 5
  group_by(site) |> 
  slice_max(date, n = 1) |>   #latest date per site
  rename(end_grw_date = date) |> 
  select(!all_above_5)

growing_season_min_temps <- start_grw |> 
  left_join(end_grw, by = "site") #do not consider the years, only the days!

