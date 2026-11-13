###Working with weather station data####
library(tidyverse)
library(tidylog)
library(lubridate)
library(ggplot2)

###Witsieshoek weather station, behind WH mountain lodge####
#Import daily data from multiple years

witsies_longtime <- read.csv("All_data\\weather station data\\UFS2_Witsies_TableDay.dat-upto Oct 2026.csv", sep = ";", header = FALSE, row.names = NULL)
colnames(witsies_longtime) <- c(witsies_longtime[2,])
witsies_longtime <- witsies_longtime[-c(1:4), ]

witsies_daily <- witsies_longtime |> 
  tibble() |> 
  mutate(
    TIMESTAMP = ymd_hm(TIMESTAMP, tz = "Africa/Johannesburg"),
    across(WSpd_Min:Rain_Tot, as.numeric),
    Date = as_date(TIMESTAMP)) |> 
  filter(between(Date, as_date("2023-01-01"), as_date("2025-12-31"))) |> 
  rename(AirTemp_mean = AirTemp_Avg, 
         AirTemp_max = AirTemp_Max, 
         AirTemp_min = AirTemp_Min, 
         Rain_total = Rain_Tot)|>
  filter(AirTemp_min > -15) |> 
  select(Date, AirTemp_min, AirTemp_mean, AirTemp_max, Rain_total)

ggplot(witsies_daily, aes(y = AirTemp_max, x = Date)) +
  geom_point()

ggplot(witsies_daily, aes(y = AirTemp_min, x = Date)) +
  geom_point()

ggplot(witsies_daily, aes(y = AirTemp_mean, x = Date)) +
  geom_point()

ggplot(witsies_daily, aes(y = Rain_Tot, x = Date)) +
  geom_point()

write.csv(witsies_daily, "All_data\\weather station data\\Witsieshoek_daily_2023-2025.csv")


####Sentinel car park####

###Import hourly data from multiple years
scarpark_longtime <- read.csv("All_data\\weather station data\\UFS3_Sentinel_TableHour-2023-2025.csv", sep = ";", header = FALSE, row.names = NULL)
colnames(scarpark_longtime) <- c(scarpark_longtime[2,])
scarpark_longtime <- scarpark_longtime[-c(1:4), ]


#compute daiy averages, maxima, minima
scarpark_daily <- scarpark_longtime |> 
  tibble() |> 
  select(TIMESTAMP, StationID, AirTemp_Avg, Rain_Tot) |> 
  mutate(
    TIMESTAMP = ymd_hm(TIMESTAMP, tz = "Africa/Johannesburg"),
    across(AirTemp_Avg:Rain_Tot, as.numeric),
    Date = as_date(TIMESTAMP)) |> 
  filter(AirTemp_Avg > -15) |>
  group_by(Date) |> 
  summarise(
    AirTemp_mean = mean(AirTemp_Avg, na.rm = TRUE),
    AirTemp_min  = min(AirTemp_Avg, na.rm = TRUE),
    AirTemp_max  = max(AirTemp_Avg, na.rm = TRUE),
    Rain_total   = sum(Rain_Tot, na.rm = TRUE),
    n_obs        = n())


ggplot(scarpark_daily, aes(y = AirTemp_max, x = Date)) +
  geom_point()

ggplot(scarpark_daily, aes(y = AirTemp_min, x = Date)) +
  geom_point()

ggplot(scarpark_daily, aes(y = AirTemp_mean, x = Date)) +
  geom_point()


#save data
write.csv(scarpark_daily, "All_data\\weather station data\\Sentinel_carpark_daily_2023-2025.csv")


####Alpine base weather station####
###Import hourly data from multiple years
alpbase_longtime <- read.csv("All_data\\weather station data\\UFS4 direct CR300_TableHour-2023-2025.csv", sep = ";", header = FALSE, row.names = NULL)
colnames(alpbase_longtime) <- c(alpbase_longtime[2,])
alpbase_longtime <- alpbase_longtime[-c(1:4), which(colnames(alpbase_longtime) %in% c("TIMESTAMP", "StationID", "AirTemp_Avg", "Rain_Tot"))]

#compute daiy averages, maxima, minima
alpbase_daily <- alpbase_longtime |> 
  tibble() |> 
  mutate(
    TIMESTAMP = ymd_hm(TIMESTAMP, tz = "Africa/Johannesburg"),
    across(AirTemp_Avg:Rain_Tot, as.numeric),
    Date = as_date(TIMESTAMP)) |>
  filter(AirTemp_Avg > -15) |>
  filter(!between(Date, as_date("2025-05-01"), as_date("2025-07-10"))) |> #takes no measurements here
  group_by(Date) |> 
  summarise(
    AirTemp_mean = mean(AirTemp_Avg, na.rm = TRUE),
    AirTemp_min  = min(AirTemp_Avg, na.rm = TRUE),
    AirTemp_max  = max(AirTemp_Avg, na.rm = TRUE),
    Rain_total   = sum(Rain_Tot, na.rm = TRUE),
    n_obs        = n())

ggplot(alpbase_daily, aes(y = AirTemp_max, x = Date)) +
  geom_point()

ggplot(alpbase_daily, aes(y = AirTemp_min, x = Date)) +
  geom_point()

ggplot(alpbase_daily, aes(y = AirTemp_mean, x = Date)) +
  geom_point()

#save data
write.csv(alpbase_daily, "All_data\\weather station data\\Alpine_base_daily_2023-2025.csv")
