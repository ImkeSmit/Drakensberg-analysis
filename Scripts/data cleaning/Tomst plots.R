###Create timeseries plots of each Tomst logger

library(ggplot)
library(ggpubr)

#import data
GG <- read.table("All_data/clean_data/Tomst_data/GoldenGate_tomst_data.txt", sep = "\t")

TMS<- unique(GG$TMS)
one_logger <- GG |> 
  filter(TMS == TMS[1])

soil_temp<-   ggplot(one_logger, aes(x = as.Date(datetime), y = soil_temp)) +
  geom_line(col = "black") +
  theme_classic()

surface_temp <- ggplot(one_logger, aes(x = as.Date(datetime), y = surface_temp)) +
  geom_line(col = "darkgreen") +
  theme_classic()

air_temp <- ggplot(one_logger, aes(x = as.Date(datetime), y = air_temp)) +
  geom_line( col = "blue") +
  theme_classic()

soil_moist <- ggplot(one_logger, aes(x = as.Date(datetime), y = raw_moisture)) +
  geom_line( col = "red") +
  theme_classic()



