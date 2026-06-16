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

####===DEFINE GROWING SEASON===####
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


###====Calculate indices for growing season===####
##GG#
GG_start <- c(growing_season_min_temps[which(growing_season_min_temps$site == "GG"), "start_grw_date"])[[1]]
GG_end <- c(growing_season_min_temps[which(growing_season_min_temps$site == "GG"), "end_grw_date"])[[1]]

GG_ind <- tms |> 
  mutate(date = as_date(datetime)) |>
  filter(site == "GG", 
         !between(date, ymd(GG_end), ymd(GG_start)), #exclude non-growing season
          between(date, ymd("2023-01-01"), ymd("2023-12-31"))) |>  #only use from 2023 growing season
  group_by(site, Cell_ID) |> 
  summarise(mean_T1_growing_season = mean(T1, na.rm = T)) |> 
  mutate(start_growing_season = ymd(GG_start), 
         end_growing_season = ymd(GG_end))
  

##WH#
WH_start <- c(growing_season_min_temps[which(growing_season_min_temps$site == "WH"), "start_grw_date"])[[1]]
WH_end <- c(growing_season_min_temps[which(growing_season_min_temps$site == "WH"), "end_grw_date"])[[1]]

WH_ind <- tms |> 
  mutate(date = as_date(datetime)) |>
  filter(site == "WH", 
         !between(date, ymd(WH_end), ymd(WH_start)), #exclude non-growing season
         between(date, ymd("2023-01-01"), ymd("2023-12-31"))) |>  #only use from 2023 growing season
  group_by(site, Cell_ID) |> 
  summarise(mean_T1_growing_season = mean(T1, na.rm = T)) |> 
  mutate(start_growing_season = ymd(WH_start), 
         end_growing_season = ymd(WH_end))


##BK#
BK_start <- c(growing_season_min_temps[which(growing_season_min_temps$site == "BK"), "start_grw_date"])[[1]]
BK_end <- c(growing_season_min_temps[which(growing_season_min_temps$site == "BK"), "end_grw_date"])[[1]]

BK_ind <- tms |> 
  mutate(date = as_date(datetime)) |>
  filter(site == "BK", 
         !between(date, ymd(BK_end), ymd(BK_start)), #exclude non-growing season
         between(date, ymd("2023-01-01"), ymd("2023-12-31"))) |>  #only use from 2023 growing season
  group_by(site, Cell_ID) |> 
  summarise(mean_T1_growing_season = mean(T1, na.rm = T)) |> 
  mutate(start_growing_season = ymd(BK_start), 
         end_growing_season = ymd(BK_end))


ind_mean_T1 <- bind_rows(GG_ind, WH_ind, BK_ind)
write.csv(ind_mean_T1, "all_data/clean_data/Environmental data/Imke_microclimate_indices.csv")

ggplot(ind_mean_T1, aes(x = site, y = mean_T1_growing_season)) +
  geom_boxplot() 
