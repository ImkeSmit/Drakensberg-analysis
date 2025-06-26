###Working with weather station data####
library(tidyverse)
library(tidylog)
library(lubridate)
library(ggplot2)

###Witsieshoek weather station, behind WH mountain lodge####
witsies <- read.delim("All_data\\weather station data\\Witsieshoek_lodge_daily.txt", sep = ",", header = T, row.names = NULL) 

colnames(witsies) <- c("Timestamp","StationID","WSpd_Min","WSpd_TMn",          
                        "WSpd_Max","WSpd_TMx","WSpd_Avg","WSpd_Std","WDir_Avg",          
                        "WDir_Std","AirTemp_Min","AirTemp_TMn","AirTemp_Max","AirTemp_TMx",       
                        "AirTemp_Avg","RH_Min","RH_TMn","RH_Max","RH_TMx",            
                        "DewPointTemp_Min","DewPointTemp_TMn","DewPointTemp_Max","DewPointTemp_TMx","DewPointTemp_Avg",  
                        "SlrW_Max","SlrW_TMx","SlrW_Avg","SlrMJ_Tot","Rain_Tot",          
                        "BPress_Min","BPress_TMn","BPress_Max","BPress_TMx","BPress_Avg",        
                        "LoggerSerialNumber","ProgramName","ProgramSignature","LoggerBattery_Min","LoggerBattery_TMn", 
                        "LoggerBattery_Max","LoggerBattery_TMx","LoggerBattery_Avg","LoggerTemp_Min","LoggerTemp_TMn",    
                        "LoggerTemp_Max","LoggerTemp_TMx","LoggerTemp_Avg","LoggerLithiumBatt","PingTime_Avg",      
                        "ScanCount","ETos","Rso","SunHrs_Tot","WC_C_Min",          
                        "WC_C_TMn", "extra")

witsies2 <- witsies |> #the scancount column was separated into two because there is a comma in the value.
  mutate(ScanCount = paste(ScanCount, ETos)) |> 
  select(!ETos) 
colnames(witsies2)[51:55] <- c("ETos","Rso","SunHrs_Tot","WC_C_Min",          
                              "WC_C_TMn")
#fix the date format
witsies2 <- witsies2 |> 
  mutate(Timestamp = gsub('"', '', Timestamp)) |> 
  mutate(Timestamp = mdy_hms(Timestamp))


#mean daily temperature
mat_witsies <- witsies2 |> 
  filter(!is.na(AirTemp_Avg)) |> 
  summarise(mat = mean(AirTemp_Avg)) #12.2989

#annual rainfall
Rai_witsies <- witsies2 |> 
  filter(!is.na(Rain_Tot)) |> 
  summarise(rai = sum(Rain_Tot)) #3780.1



####Sentinel car park

scarpark <- read.delim("All_data\\weather station data\\Sentinel_car_park_daily.txt", sep = ",", header = T, row.names = NULL)
colnames(scarpark) <- c("Timestamp","StationID","WSpd_Min","WSpd_TMn",          
                        "WSpd_Max","WSpd_TMx","WSpd_Avg","WSpd_Std","WDir_Avg",          
                        "WDir_Std","AirTemp_Min","AirTemp_TMn","AirTemp_Max","AirTemp_TMx",       
                        "AirTemp_Avg","RH_Min","RH_TMn","RH_Max","RH_TMx",            
                        "DewPointTemp_Min","DewPointTemp_TMn","DewPointTemp_Max","DewPointTemp_TMx","DewPointTemp_Avg",  
                        "SlrW_Max","SlrW_TMx","SlrW_Avg","SlrMJ_Tot","Rain_Tot",          
                        "BPress_Min","BPress_TMn","BPress_Max","BPress_TMx","BPress_Avg",        
                        "LoggerSerialNumber","ProgramName","ProgramSignature","LoggerBattery_Min","LoggerBattery_TMn", 
                        "LoggerBattery_Max","LoggerBattery_TMx","LoggerBattery_Avg","LoggerTemp_Min","LoggerTemp_TMn",    
                        "LoggerTemp_Max","LoggerTemp_TMx","LoggerTemp_Avg","LoggerLithiumBatt","PingTime_Avg",      
                        "ScanCount","ETos","Rso","SunHrs_Tot","WC_C_Min",          
                        "WC_C_TMn" , "extra")

scarpark2 <- scarpark |> #the scancount column was separated into two because there is a comma in the value.
  mutate(ScanCount = paste(ScanCount, ETos)) |> 
  select(!ETos) 
colnames(scarpark2)[51:55] <- c("ETos","Rso","SunHrs_Tot","WC_C_Min",          
                               "WC_C_TMn")
#fix the date format
scarpark2 <- scarpark2 |> 
  mutate(Timestamp = dmy(Timestamp))

#date range
range(scarpark2$Timestamp) #"2024-07-09" "2025-06-25"

#mean daily temperature
mat_witsies <- scarpark2 |> 
  filter(!is.na(AirTemp_Avg)) |> 
  summarise(mat = mean(AirTemp_Avg)) #10.52712

#annual rainfall
Rai_witsies <- scarpark2 |> 
  filter(!is.na(Rain_Tot)) |> 
  summarise(rai = sum(Rain_Tot)) #4071.2


