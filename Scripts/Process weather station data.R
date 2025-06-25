###Working with weather station data####

library(tidyverse)
library(lubridate)
library(ggplot2)

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
                        "WC_C_TMn" )
scarpark$Timestamp <- dmy(scarpark$Timestamp) #change to date format

ggplot(scarpark) +
  geom_bar(aes(x = Timestamp, y = Rain_Tot), stat = "identity")

ggplot(scarpark) +
  geom_line(aes(x = Timestamp, y = AirTemp_Avg))
