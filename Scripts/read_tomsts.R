###########################################################################3
# Read & plot Tomst microclimate data

library(tidyverse)
library(lubridate)
library(data.table)
library(patchwork)

# Set date limits to remove implausible dates
mind <- "2023-01-12"
maxd <- "2023-01-21"

datadir <- ""
datadir <- "All_data/raw_microclimate_data/Bokong/Tomst_data_Final_reading_February2024"

#################################################################################3

# read names x tomst id table

ids <- read_excel(paste0(datadir, "/", "Bokong_Logger_notes.xlsx")) %>% 
  mutate(plot = paste0(Site, "_",Grid,Cell))

# List logger data files to read
f <- list.files(datadir, pattern = "data_[0-9]", full.names = T, recursive = T)

fi <- data.frame(file = f)

fi$tomst_id <- unlist(lapply(fi$file, function(x) as.numeric(strsplit(gsub("data_","",rev(strsplit(x, "/")[[1]])[1]), "_")[[1]][1])))

fi <- full_join(fi, ids) %>% 
  mutate(plot = toupper(plot)) %>% 
  filter(!is.na(file))

fi$file2 <- gsub("_..csv", "", fi$file)

fi <- fi[order(fi$plot),]

# This should be zero
fi %>% filter(is.na(plot)) %>% nrow

readdata <- function(i){
  nn <- sum(grepl(i, fi$file2))
  
  if(nn > 1){
    
    fi2 <- fi %>% filter(grepl(i, fi$file2))
    
    df2 <- data.frame()
    for(ii in fi2$file2){
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
    d <- fread(fi$file[fi$file2 == i])
    
    d %>% select(V2,V3,V4,V5,V6,V7) -> d
    
    d %>% filter(!duplicated(.$V2, fromLast = T)) -> d
    
    d %>% mutate(across(V4:V6, ~as.numeric(gsub(",",".\\",.)))) -> d
    
    d$plot <- fi[which(fi$file2 == i),"plot"]
    
    d %>% mutate(V2 = ymd_hm(V2, tz = "UTC")) %>% 
      mutate(V2 = with_tz(V2, tzone = "Etc/GMT-2")) -> d
    
    return(d)
    
  }
  
}


# read the Tomst data in
mylist <- lapply(fi$file2, readdata)
df <- rbindlist( mylist )

# Rename columns
df %>% rename(datetime = V2,
              zone = V3,
              T1 = V4,
              T2 = V5,
              T3 = V6,
              moist = V7) -> df

df %>% arrange(plot, datetime) -> df

# Remove implausible dates
df %>% filter(datetime >= mind,
              datetime <= maxd) -> df

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
  df %>% filter(plot == i) %>% 
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "T3"), col = "cornflowerblue") +
    geom_line(aes_string(y = "T2"), col = "brown1") +
    geom_line(aes_string(y = "T1"), col = "darkgoldenrod") +
    theme_minimal() +
    ylab("Temperature") + xlab("Date")+
    scale_y_continuous(limits = c(-5, 35))+
    ggtitle(i) -> g1
  
  df %>% filter(plot == i) %>% 
    ggplot(aes_string(x="datetime")) +
    geom_line(aes_string(y = "moist"), col = "cornflowerblue") +
    theme_minimal() +
    ylab("Moisture count") + xlab("Date")+
    scale_y_continuous(limits = c(0, 3900)) -> g2
  
  
  print(g1 / g2)
  
}
dev.off()
