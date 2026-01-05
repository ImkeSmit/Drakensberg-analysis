###Cluster cells according to environmental similarity###
library(tidyverse)
library(tidylog)

#import environmental data
env <- read.csv("All_data/clean_data/Environmental data/All_Sites_Environmental_Data.csv")
