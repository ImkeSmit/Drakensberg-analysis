###########################################################################3
# Read & plot Tomst microclimate data

library(tidyverse)
library(tidylog)
library(lubridate)
library(data.table)
library(patchwork)
library(readxl)
library(ggplot2)

# Set date limits to remove implausible dates
mind <- as.Date("2023-11-10") #get this installment dates from the BOkong logger notes file
maxd <- as.Date("2024-02-04")

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

# Rename columns
df %>% rename(datetime = V2,
              zone = V3,
              soil_temp = V4, #T1
              surface_temp = V5, #T2
              air_temp = V6, #T3
              moist = V7) -> df

df %>% arrange(plot, datetime) -> df

# Remove implausible dates
df2 <- df |> filter(datetime >= mind,
              datetime <= maxd)


############################################################################
# PLOTTINGS
############################################################################

# Select the plot names to be plotted
plots <- unique(df$plot)

# Plot temperatures
pdf("tomst.pdf", 12, 10)
for(i in plots){
  #i <- plots[3]
  print(i)
  
  g1 <- df2 |>filter(plot == i) |> 
    ggplot() +
    geom_line(aes(x = datetime, y = soil_temp, col = 'soil_temp')) +
    geom_line(aes(x = datetime, y = surface_temp, col = 'surface_temp'))+
    geom_line(aes(x = datetime, y = air_temp, col = 'air_temp')) +
    scale_color_manual(name=' ',
                       breaks=c('soil_temp', 'surface_temp', 'air_temp'),
                       values=c('soil_temp'='darkgoldenrod', 
                                'surface_temp'='brown1', 
                                'air_temp'='cornflowerblue'))+
    theme_minimal() +
    ylab("Temperature") + xlab("Date")+
    scale_y_continuous(limits = c(-5, 35))+
    ggtitle(i) +
    theme(legend.position = "bottom")
  
  g2 <- df2 |> filter(plot == i) |> 
    ggplot(aes(x = datetime, y = moist), col = "cornflowerblue") +
    geom_line(col = "cornflowerblue") +
    theme_minimal() +
    ylab("Moisture count") + xlab("Date")+
    scale_y_continuous(limits = c(0, 3900))
  
  
  print(g1 / g2)
  
}
dev.off()


#let's close in on a week in BOKONG_1C14
testplot <- df2 |> filter(plot == "BOKONG_1C14", 
              datetime >= as.Date("2023-11-15"), datetime <= as.Date("2023-11-22")) |> 
  ggplot() +
  geom_line(aes(x = datetime, y = soil_temp, col = 'soil_temp')) +
  geom_line(aes(x = datetime, y = surface_temp, col = 'surface_temp'))+
  geom_line(aes(x = datetime, y = air_temp, col = 'air_temp')) +
  scale_color_manual(name=' ',
                     breaks=c('soil_temp', 'surface_temp', 'air_temp'),
                     values=c('soil_temp'='darkgoldenrod', 
                              'surface_temp'='brown1', 
                              'air_temp'='cornflowerblue'))+
  theme_minimal() +
  ylab("Temperature") + xlab("Date")+
  scale_y_continuous(limits = c(-5, 50))+
  ggtitle("BOKONG_1C14") +
  theme(legend.position = "bottom")
  
