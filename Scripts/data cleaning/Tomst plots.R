###Create timeseries plots of each Tomst logger

library(ggplot)
library(ggpubr)

#import data
GG <- read.table("All_data/clean_data/Tomst_data/GoldenGate_tomst_data.txt", sep = "\t")

TMS_codes<- unique(GG$TMS)


pdf("Figures//goldengate_loggers.pdf")
for(t in 1:length(TMS_codes)) {
  
  one_logger <- GG[which(GG$TMS == TMS_codes[t]), ]
  
  soil_temp<- ggplot(one_logger, aes(x = as.Date(datetime), y = soil_temp)) +
    geom_line(col = "black") +
    theme_classic()+
    labs(title = as.character(TMS_codes[t]))
  
  surface_temp <- ggplot(one_logger, aes(x = as.Date(datetime), y = surface_temp)) +
    geom_line(col = "darkgreen") +
    theme_classic()+
    labs(title = as.character(TMS_codes[t]))
  
  air_temp <- ggplot(one_logger, aes(x = as.Date(datetime), y = air_temp)) +
    geom_line( col = "blue") +
    theme_classic()+
    labs(title = as.character(TMS_codes[t]))
  
  soil_moist <- ggplot(one_logger, aes(x = as.Date(datetime), y = raw_moisture)) +
    geom_line( col = "red") +
    theme_classic()+
    labs(title = as.character(TMS_codes[t]))
  
  print(soil_temp/surface_temp)
  print(air_temp/soil_moist)
  
  #gridExtra::grid.arrange(soil_temp, surface_temp, air_temp, soil_moist, ncol = 1, nrow = 2)
  
}

dev.off()


prob <- GG[which(GG$TMS == 95223917), ]
max(prob$datetime)
min(prob$datetime)
