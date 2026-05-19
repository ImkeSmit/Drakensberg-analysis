######################################
# GEE and GLM models for Saana data  #
######################################

setwd("D:\\UP_research\\Bridgette_McMillan\\GEE")
source("Code_geefunctions_modified_v12.R")
source("Source_built_all_subset_models_v6.r")

OutputFileName <- "GEE_test"
ResponseVariable <- "Vasc_cov"                                 # ***
PredictorsColumns <- c(26, 15, 16, 23, 17, 25, 22)               # ***
family1 <- "binomial"                      # "gaussian" or "poisson" or "binomial"
ModelProportion.value <- "ResponseVariable"      # set to NA if not modelling a proportion (i.e. cover); else set to "ResponseVariable" [with quotation marks]

PredictorCombinations2 <- list()

###########################################################
# If only species models requried, list each individually #
###########################################################

# PredictorCombinations2[[1]] <- "richness_1 ~ sm_mean + st_mean + mesotopo + distall + ph"

###########################
# Copy and paste all code #
###########################

 
###############
# Import data #
###############

library(spdep)
library(MASS)
library(gee)
library(utils)
require(geepack)
require(fields)
library(ncf)
library(car)

winDialog(type = "ok", "Please select the input file (.txt)")
FileName <- file.choose()
setwd(dirname(FileName))
data1 <- read.delim(FileName)

increm = 1  

#** convert all values to numeric     ## probably not necessary to use these three lines
all.columns <- 1 : length(data1[1, ])
all.numeric.cols <- all.columns[-which(colnames(data1) == "column")]
for (i in all.numeric.cols) {data1[, i] <- as.numeric(as.character(data1[, i]))}


#** Check for NAs in variables, and subset dataset to drop NAs
data1_all <- data1
number_of_incomplete_rows <- length(which(complete.cases(data1) == FALSE))
number_of_incomplete_rows
data1 <- data1[complete.cases(data1),]

##############################
# Start looping through data #
##############################                                           

ResponseVariable.column <- which(colnames(data1) == ResponseVariable)
ResponseVariable.name <- ResponseVariable

if (!is.na(ModelProportion.value)) {
    ResponseVariable <- data1[, ResponseVariable.column]/100
} else {
    ResponseVariable <- data1[, ResponseVariable.column]
    if (family1 == "poisson" & any(is.wholenumber(ResponseVariable) != TRUE)) {
        print("Data multiplied by 10 and converted to integers")
        ResponseVariable <- round(ResponseVariable * 10, 0)}}

data1$zero.group <- factor(rep("a", nrow(data1)))
coords <- round(data.frame(data1$longitude, data1$latitude))

attach(data1)
nn <- nrow(data1)
id1 <- rep(1, nn)
id2 <- 1 : nrow(data1)
nb <- nb2listw(dnearneigh(cbind(data1$longitude, data1$latitude), 0, 23, longlat = FALSE), glist = NULL, style = "W", zero.policy = TRUE)

# check there are the right number of groups (i.e. grids)
plot(nb, coords)

########################################
# Create all combinations of variables #
########################################

PredictorCombinations2 <- AllSubsets(ResponseVariableColumn = ResponseVariable.column,                                        PredictorsColumns = PredictorsColumns, 
                             data.source = data1, 
                             ModelProportion = ModelProportion.value, 
                             Add.PolynomialTerms = TRUE, 
                             Polynom.exclude = 26, 
                             Polynom.order = 2)
PredictorCombinations2.mat <- as.matrix(PredictorCombinations2)
colnames(PredictorCombinations2.mat) <- c("Model")
write.table(PredictorCombinations2.mat, paste(OutputFileName, "_AllModels.txt", sep = ""), quote = FALSE)


########################
# run null model first #
########################

ResponseVariable.original <- ResponseVariable

model.counter = 0
    Mod_0 <- glm(ResponseVariable ~ 1, family = family1)
    ExtractRegRes(Mod_0, max.terms = length(PredictorsColumns))
    Mod_0_nonGEE <- GEE(model.name = "Mod0", model = Mod_0, family1 = family1, data2 = data1, coord = coords, corstr = "independence", plot.ac = FALSE, graph = FALSE, output = TRUE)
    Mod_0_GEE <- GEE(model.name = "Mod0", model = ResponseVariable ~ 1, family1 = family1, data2 = data1, coord = coords, corstr = "fixed", plot.ac = FALSE, graph = FALSE, output = TRUE)

null.quasidev <- (-2) * QIC(Mod0_IndepMod, Mod0_IndepMod, family = family1)[2]
QIC(Mod0_IndepMod, Mod0_FixedMod, family = family1)
ExtractGEE(Mod0_FixedMod)
temp.objects <- as.character(0)
temp.objects <- c(temp.objects, (paste("Mod",model.counter,"_FixedMod",sep="")), as.character(paste("Mod",model.counter,"_IndepMod",sep="")), as.character(paste("Mod_",model.counter,sep="")), as.character(paste("Mod_", model.counter, "_nonGEE", sep="")), as.character(paste("Mod_", model.counter, "_GEE", sep="")), as.character(paste("Mod", model.counter, "_autocorrelation", sep="")))
remove(list = temp.objects[-1])

##############################
# run all GLM and GEE models #
##############################

for (model.counter in 1 : length(PredictorCombinations2)) {
	assign(paste("Mod_", model.counter, sep=""), glm(eval(parse(text = PredictorCombinations2[[model.counter]])), family = family1, data = data1))

	VIFoutput <- c(1, model.counter, length(data1[, 1]), AIC(eval(parse(text=paste("Mod_", model.counter, sep = "")))))   # in previous version "1" here SubNumber ??
	if (eval(parse(text = paste("length(attr(Mod_", model.counter ,"$terms, 'term.labels'))", sep = ""))) > 1) {
		temp.vif <- vif(eval(parse(text = paste("Mod_", model.counter, sep = ""))))
		if (is.vector(temp.vif)) {VIFoutput <- t(c(VIFoutput, temp.vif))} else {VIFoutput <- t(c(VIFoutput, temp.vif[, 3]))}
		} else {VIFoutput <- t(c(VIFoutput, NA))}
	if (file.exists(paste(OutputFileName, "_VIF.txt", sep = ""))) {
			write.table(VIFoutput, file = paste(OutputFileName,"_VIF.txt", sep = ""), append = TRUE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = FALSE, row.names = FALSE)
			} else {write.table(VIFoutput, file = paste(OutputFileName, "_VIF.txt", sep = ""), append = FALSE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = TRUE, row.names = FALSE)}

	ExtractRegRes(eval(parse(text=paste("Mod_", model.counter, sep = ""))))
	assign(paste("Mod_", model.counter, "_nonGEE", sep=""), GEE(eval(paste("Mod", model.counter, sep = "")), eval(parse(text = paste("Mod_", model.counter, sep = ""))), family1 = family1, data2 = data1, coord = coords, corstr = "independence", plot.ac = FALSE, graph = FALSE, output = TRUE))
	assign(paste("Mod_", model.counter, "_GEE", sep = ""), GEE(eval(paste("Mod", model.counter, sep = "")), eval(parse(text = paste("Mod_", model.counter, sep = ""))), family1 = family1, data2 = data1, coord = coords, corstr = "fixed", plot.ac = FALSE, graph = FALSE, output = TRUE))
	if (separation.error == FALSE) {
		assign(paste("Mod_", model.counter, "_QIC", sep=""), QIC(eval(parse(text=paste("Mod", model.counter, "_IndepMod", sep=""))), eval(parse(text = paste("Mod", model.counter, "_FixedMod", sep = ""))), family = family1))}
	ExtractGEE(eval(parse(text=paste("Mod", model.counter, "_FixedMod", sep=""))))

	temp.objects<-as.character(0)
	temp.objects<-c(temp.objects,as.character(paste("Mod",model.counter,"_FixedMod",sep="")), as.character(paste("Mod",model.counter,"_IndepMod",sep="")), as.character(paste("Mod_",model.counter,sep="")), as.character(paste("Mod_", model.counter, "_nonGEE", sep="")), as.character(paste("Mod_", model.counter, "_GEE", sep="")), as.character(paste("Mod_", model.counter, "_QIC", sep="")), as.character(paste("Mod", model.counter, "_autocorrelation", sep="")))  # ,as.character(paste("Mod",model.counter,"_FixedMod_wcs",sep=""))
	remove(list=temp.objects[-1])
	
	print("-----------")
	print(paste("Completed model #", model.counter, " at ", Sys.time(), sep = ""))
	print("-----------")
	
	save.image(file = paste(OutputFileName, ".RData", sep = ""))
	# counter.start = counter.start + 1
}

detach(data1)

