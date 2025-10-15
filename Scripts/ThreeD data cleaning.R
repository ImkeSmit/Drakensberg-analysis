#import and clean community data from Three_D experiment
library(tidyverse)
library(readxl)
library(openxlsx)


#run import_community script first

high2025 <- import_community()
