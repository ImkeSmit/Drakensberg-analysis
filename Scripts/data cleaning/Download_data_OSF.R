#Example script to download data from OSF
install.packages("remotes")
library(remotes)
remotes::install_github("Between-the-Fjords/dataDownloader") #install the dataDownloader package from github
library(dataDownloader)
install.packages("osfr")
library(osfr)

#Generate a personal access token on OSF and authenticate your PAT
#here is a vignette on how to do that: https://cran.r-project.org/web/packages/osfr/vignettes/auth.html
#You can put your token in a .Renviron file to automatically autheticate when you open your Rproject:
#install.packages("usethis)
#library(usethis)
#usethis::::edit_r_environ("project")
#Then edit the .Renviron in notepad or R with OSF_PAT=yourpersonalaccestoken
#Restart R for it to take effect

get_file(node = "374bk", #node of the Micro-climb OSF project
         file = "Differential_GPS_August_2023.xlsx", #name of file you want to download
         path = "All_data/raw_data", #where you want to save this dataset
         remote_path = "Clean_data/Geospatial_data") #Path on OSF to where this file is stored

#Now you can import the file as normal