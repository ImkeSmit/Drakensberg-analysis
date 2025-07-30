##Commmunity assembly analysis####
install.packages("untb")
library(untb)
Sys.setenv(gp_binary = "C:/Program Files/Pari64-2-17-2/gp.exe")
gp_binary = "C:/Program Files/Pari64-2-17-2/gp.exe"

data("saunders")
summary(saunders.tot)
preston(saunders.tot, n= 9)
optimal.theta(saunders.tot)
optimal.params(saunders.tot, gp_binary = "C:/Program Files/Pari64-2-17-2/gp.exe")
etienne(saunders.tot)

logkda(saunders.tot, method = "R")
