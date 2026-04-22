# increase nsim and reps from 9 to 19 or 99!
# alter "increment" to determine distance bins


#################################################
# Create function to extract regression results #
#################################################

ExtractRegRes <- function (ModelOutput, OutputFile = OutputFileName) {
	ModelOutputSummary <- summary(ModelOutput)
	if (dim(ModelOutputSummary$coeff)[1] == 1) {temp.VIF <- 999} else {temp.VIF <- NA} 	# VIF set to NA due to use of GVIF in other section of code
	TempResult <- data.frame("Result" = NA)
		TempResult[1] <- "empty"
		TempResult[2] <- noquote(as.character(eval(parse(text=ModelOutput$call[2]))[3]))
		TempResult[3] <- length(data1[,1])
		TempResult[4] <- ModelOutput$df.null
		TempResult[5] <- ModelOutput$df.residual
		TempResult[6] <- ModelOutput$null.deviance - ModelOutput$deviance	#model stat - a.k.a. likelihood ratio statistic
		TempResult[7] <- 1 - pchisq(ModelOutput$null.deviance - ModelOutput$deviance, ModelOutput$df.null - ModelOutput$df.residual)
		TempResult[8] <- ModelOutput$aic
		TempResult[9] <- ModelOutput$deviance				# G2 / Deviance = measure of goodness of fit
		TempResult[10] <- 1 - pchisq(ModelOutput$deviance, ModelOutput$df.residual)
		TempResult[11] <- sum(residuals(ModelOutput, type="pearson")^2)		# Pearson Chi2 = measure of goodness of fit
		TempResult[12] <- 1 - pchisq(sum(residuals(ModelOutput, type="pearson")^2), ModelOutput$df.residual)
  	TempResult[13] <- 999
		TempResult[14] <- length(ModelOutput$coeff)
		TempResult[15] <- ModelOutput$null.deviance
		TempResult[16] <- ModelOutput$deviance
		TempResult[17] <- (ModelOutput$null.deviance-ModelOutput$deviance)/ModelOutput$null.deviance
		TempResult[18] <- 1-(((length(data1[,1])-1)/(length(data1[,1])-length(ModelOutput$coeff)))*(1-((ModelOutput$null.deviance-ModelOutput$deviance)/ModelOutput$null.deviance)))
		TempResult[19] <- noquote(as.character(ModelOutput$call[1]))
		TempResult[20] <- noquote(as.character(ModelOutput$call[3]))
		for (ParameterNumber in 1 : length(ModelOutput$coeff)) {
			if (dim(ModelOutputSummary$coeff)[1] < ParameterNumber) {break}
			TempResult[21 + ((ParameterNumber - 1) * 6)] <- noquote(names(ModelOutput$coeff)[ParameterNumber])
			TempResult[22 + ((ParameterNumber - 1) * 6)] <- ModelOutputSummary$coeff[ParameterNumber,1]
			TempResult[23 + ((ParameterNumber - 1) * 6)] <- ModelOutputSummary$coeff[ParameterNumber,2]
			TempResult[24 + ((ParameterNumber - 1) * 6)] <- ModelOutputSummary$coeff[ParameterNumber,3]
			TempResult[25 + ((ParameterNumber - 1) * 6)] <- ModelOutputSummary$coeff[ParameterNumber,4]
			if(ParameterNumber==1) {TempResult[26 + ((ParameterNumber - 1) * 6)] <- NA} else {TempResult[26 + ((ParameterNumber - 1) * 6)] <- temp.VIF[ParameterNumber - 1]}}

	if (file.exists(paste(OutputFile,"_GLM.txt", sep=""))) {
		write.table(TempResult, file = paste(OutputFile, "_GLM.txt", sep = ""), append = TRUE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = FALSE, row.names = FALSE)
		} else {
		Results.lab1 <- c("Species", "Model", "n", "null d.f.", "residual d.f.", "Model stat.", "Model p", "AIC", "G2 g.o.f.", "g2 p", "Pearson Chi g.o.f.", "Pearson Chi p", "Empty", "Number parameters", "Null dev.", "Residual dev.", "%dev explained", "adj. D2", "Method")
		Results.lab3 <- "Distribution"
		Results.lab2 <- c("Variable", "Coefficient", "SE", "Stat.", "p-value", "VIF")
    counter2 <- length(ModelOutput$coeff)
			for (ParameterNumber in 1:counter2) {
			for (counter in 1 : 6) {
			Results.lab3 <- c(Results.lab3, paste(Results.lab2[counter], "_", ParameterNumber, sep = ""))}}
		colnames(TempResult) <- c(Results.lab1, Results.lab3)
		write.table(TempResult, file = paste(OutputFile, "_GLM.txt", sep = ""), append = FALSE, sep = "\t", dec = ".", quote = FALSE, eol = "\n", na = "NA", col.names = TRUE, row.names = FALSE)}

} # ends function definition




################################################
# Function to extract GEE stats, including QIC #
################################################

ExtractGEE <- function (ModelOutput, QIC.value = QIC.values, OutputFile = OutputFileName) {
	if(separation.error == FALSE){
		TempResult <- data.frame("Result" = NA)
		TempResult[1] <- NA
		TempResult[2] <- noquote(as.character(paste("GEE + ", ModelOutput[[1]], sep="")))
		TempResult[3] <- model.counter
		TempResult[4] <- ModelOutput[[8]]$n
		TempResult[5] <- ModelOutput[[11]][1]
		TempResult[6] <- ModelOutput[[11]][2]
		TempResult[7] <- UsingRobust
		TempResult[8] <- round(ModelOutput[[10]][1] * (exp(-5 * ModelOutput[[10]][2])), 3)
		TempResult[9] <- ModelOutput[[8]]$model$link
		TempResult[10] <- ModelOutput[[8]]$model$varfun
		TempResult[11] <- ModelOutput[[8]]$model$corstr
		TempResult[12] <- ModelOutput[[8]]$scale
		TempResult[13] <- null.quasidev						#null deviance; i.e. without correlation structure or explanatory variables
		TempResult[14] <- (-2) * QIC.values[2]					#residual deviance ???!!! check this value - does not seem correct for poisson disrtibution
		
		TempResult[15] <- (ModelOutput[[8]]$n) - 1					#null d.f.
		TempResult[16] <- (ModelOutput[[8]]$n) - length(ModelOutput[[2]]$coeff)	#residual d.f.
 		TempResult[17] <- (null.quasidev - ((-2) * QIC.values[2])) / null.quasidev	#deviance explained
		TempResult[18] <- 1 - (sum((data1[, ResponseVariable.column] - ModelOutput[[7]])^2) / sum((data1[, ResponseVariable.column] - mean(data1[, ResponseVariable.column])) ^2))
		TempResult[19] <- QIC.values[1]
		TempResult[20] <- QIC.values[2]
		TempResult[21] <- QIC.values[3]
		TempResult[22] <- QIC.values[4]
		TempResult[23] <- length(ModelOutput[[2]]$coeff)
		for (ParameterNumber in 1:length(ModelOutput[[2]]$coeff)) {
			TempResult[24 + ((ParameterNumber - 1) * 8)] <- noquote(rownames(ModelOutput[[9]])[ParameterNumber])
			TempResult[25 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[9]][ParameterNumber,1]
			TempResult[26 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[9]][ParameterNumber,2]
			TempResult[27 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[9]][ParameterNumber,3]
			TempResult[28 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[9]][ParameterNumber,4]
			TempResult[29 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[10]][ParameterNumber,1]
			TempResult[30 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[10]][ParameterNumber,2]
			TempResult[31 + ((ParameterNumber - 1) * 8)] <- ModelOutput[[10]][ParameterNumber,3]}
	} else {
		TempResult <- data.frame("Result"=NA)
		TempResult[1] <- NA
		TempResult[2] <- "GEE"
		TempResult[3] <- model.counter
		TempResult[4] <- "Data separation error"
		TempResult[5] <- ModelOutput[[10]][1]
		TempResult[6] <- ModelOutput[[10]][2]
		TempResult[7] <- ModelOutput[[11]]
		TempResult[8] <- round(ModelOutput[[10]][1] * (exp(-30 * ModelOutput[[10]][2])), 3)
	}

	if (file.exists(paste(OutputFile,"_GEE.txt", sep=""))) {
		write.table(TempResult, file = paste(OutputFile,"_GEE.txt", sep=""), append = TRUE, sep="\t", dec=".", quote = FALSE, eol="\n", na = "NA", col.names=FALSE, row.names=FALSE)
		} else {
		Results.lab1<-c("Species", "Model", "Model number", "n","Decay constant1","Decay constant2","Negative diagonal?","Autocorrel. at 5m", "Link function", "Family", "Correlation structure", "Scale",
			"Null deviance", "Residual deviance", "Null d.f.", "Residual d.f.", "%Dev. explained",
			"Adj. R2", "QIC", "Quasi-likelihood", "Trace", "Naive trace")
		Results.lab3<-"Number of parameters"
		Results.lab2<-c("Variable", "Coefficient", "SE", "Z", "p-value", "Robust SE", "Robust Z", "Robust p")  #check if naive or robust SE and Z used!
		counter2<-length(ModelOutput[[2]]$coeff)
			for (ParameterNumber in 1:counter2) {
			for (counter in 1:8) {
			Results.lab3<-c(Results.lab3,paste(Results.lab2[counter],"_", ParameterNumber,sep=""))}}
		colnames(TempResult)<-c(Results.lab1, Results.lab3)
		write.table(TempResult, file = paste(OutputFile,"_GEE.txt", sep=""), append = FALSE, sep="\t", dec=".", quote = FALSE, eol="\n", na = "NA", col.names=TRUE, row.names=FALSE)}
}



#######################################################################
# R function to create custom correlation matrices for Saana datasets #
#######################################################################

create.Saana.correlation.matrix <- function(data3, formula2, family1, Do.plot = TRUE, Save.plot = TRUE) {
    
    # set parameters
    NumberSites <- 8
    NumberRows <- 20
    NumberColumns <- 8
        
    TotalCells <- length(data3[, 1])
    if(any(class(formula2) == "formula")) {formula1 <- formula2} else {formula1 <- formula2$formula}
    
    # load libraries
    require(ncf)
    
    # create empty matrix
    R <- matrix(data = 0, nrow = TotalCells, ncol = TotalCells)
    DecayValues <- matrix(data = 0, nrow = NumberSites, ncol = 2)
    
    # if plotting...
    if (Save.plot == TRUE) {tiff(paste("Mod_", model.counter, "_exponential_fit.tif", sep = ""), units = "cm", height = 16, width = 24, res = 60)}
    if (Do.plot == TRUE) {par(mfrow=c(2, 4))}
    
    # loop through each site
    for (ii in 1 : NumberSites) {
    
        # extract data subsets for each site
        data4 <- data1[data1$grid == ii, ]
        ResponseVariable <<- ResponseVariable.original[data1$grid == ii]
        
        # calculate GLM model (i.e. null model, assuming no spatial autocorrelation)
        m0.temp <- glm(formula = formula1, family = family1, data = data4)
        res0.temp <- resid(m0.temp)
        
        # calculate autocorrelation from GLM model residuals
        if (all(res0.temp == res0.temp[1])) {
          print("No variation - check response variable"); next
        } else {
          corglm <- correlog(data4$longitude, data4$latitude, res0.temp, na.rm = TRUE, increment = increm, resamp = 9, quiet = TRUE)
        }
        
        #fitting negative exponential curve to data
        if (any(corglm$correlation < 0 | corglm$p > 0.025)) {first.ns.lag <- min(which(corglm$correlation < 0 | corglm$p > 0.025))} else {first.ns.lag <- length(corglm$correlation)}
        if (first.ns.lag < 8) {first.ns.lag <- 8}
        
        nl.data <- cbind(corglm$mean.of.class[1 : length(corglm$correlation)], corglm$correlation[1 : length(corglm$correlation)])
        nl.data[(first.ns.lag : length(corglm$correlation)), 2] <- 0       
           
        exp.formula <- nl.data[,2] ~ beta1 * (exp(1)^(-1 * nl.data[,1] * beta2))
    	  exp.model <- try(nls(exp.formula, start = list(beta1 = 0.75, beta2 = 0.001), trace = FALSE, control = list(maxiter = 500, tol = 1e-04, warnOnly = TRUE)))
        if(inherits(exp.model, "try-error")) {
            nl.data[nl.data[, 2] < 0, 2] <- 0
            exp.model <- try(nls(exp.formula, start = list(beta1 = 0.75, beta2 = 0.001), trace = FALSE, control = list(maxiter = 500, tol = 1e-04, warnOnly = TRUE)))
            if(inherits(exp.model, "try-error")) {print("Adjusted correlation data to ensure exponential curve fit")}
        }
        if (coef(exp.model)[1] >= 1) {
            exp.formula <- nl.data[,2] ~ 1 * (exp(1)^(-1 * nl.data[,1] * beta2))           
            exp.model <- try(nls(exp.formula, start = list(beta2 = 0.001), trace = FALSE, control = list(maxiter = 500, tol = 1e-04, warnOnly = TRUE)))
            if(inherits(exp.model, "try-error")) {
               nl.data[nl.data[, 2] < 0, 2] <- 0
               exp.model <- try(nls(exp.formula, start = list(beta1 = 0.75, beta2 = 0.001), trace = FALSE, control = list(maxiter = 500, tol = 1e-04, warnOnly = TRUE)))
               if(inherits(exp.model, "try-error")) {print("Adjusted correlation data to ensure exponential curve fit")}}
        }
      	DecayConstant <- coef(exp.model)
      	if (length(DecayConstant) == 1) {DecayConstant <- c(1, DecayConstant)}
      	exp.model.error <- exp.model$convInfo$stopMessage

    	  if(DecayConstant[1]<=0) {DecayConstant[1]=0; DecayConstant[2]=0}
    	  if(DecayConstant[2]<=0) {DecayConstant[1]=0; DecayConstant[2]=0}
    	  
    	  print(paste("Site = ", ii, "Decay constants = ", round(DecayConstant[1], 4), " and ", round(DecayConstant[2], 4), sep = ""))
        if (Do.plot == TRUE) { 
           	plot(corglm$mean.of.class[1 : length(corglm$mean.of.class)/2], corglm$correlation[1 : length(corglm$mean.of.class)/2], type = "l", xlab = paste("Distance (", 1 / increm, " m units)", sep = ""), ylab = "Autocorrelation", ylim=c(-1,1), main=ii)
        		points(corglm$mean.of.class[1 : length(corglm$mean.of.class)/2], DecayConstant[1]*(exp(1)^(-1 * corglm$mean.of.class[1 : length(corglm$mean.of.class)/2] * DecayConstant[2])), col="red", type="l")
        		legend("bottomleft", col=c("black", "red"), legend=c("GLM residuals", "Negative exponential fit"), lty=1)
        		text(x = mean(corglm$mean.of.class[1 : length(corglm$mean.of.class)/2]), y = 0.8 ,paste("Decay params = ",round(DecayConstant[1], 2), " & " ,round(DecayConstant[2], 3), "; 1st ns lag = ", first.ns.lag, sep = ""), cex = 0.75)
        		abline(h = 0, lty = 2)
        }
    
        # calculate distance matrix
        D.mat <- (as.matrix(dist(cbind(data4$longitude, data4$latitude))))
        
        # calculate R for data subset
        R.temp <- DecayConstant[1]*(exp(-1 * D.mat * DecayConstant[2]))
      	R.temp[D.mat == 0] <- 1
        if(max(R.temp[D.mat != 0]) > 0) {print(paste("Max. correlation between obs in working correlation structure = ", round(max(R.temp[D.mat != 0]), 3), sep = ""))}
        
        # add R for data subset to R for whole dataset
        R[((1 + length(which(data1$grid <= (ii - 1)))) : ((length(which(data1$grid <= ii))))), ((1 + length(which(data1$grid <= (ii - 1)))) : ((length(which(data1$grid <= ii)))))] <- as.numeric(R.temp)
        
        # save decay constant values to matrix
        DecayValues[ii, ] <- DecayConstant
        
    } # end loop through sites
  if (Save.plot == TRUE) {dev.off()}
  R
}

 
###############
# GEEcode.txt #
###############
# Carl & Kühn: Analyzing Spatial Autocorrelation in Species Distributions using Gaussian and Logit Models

GEE <- function(model.name, model, family1 = "poisson", data2, coord, corstr = "independence", cluster = 3, plot.ac = FALSE, graph = TRUE, output = TRUE, OutputFile = OutputFileName) {

reps = 9
separation.error <<- FALSE
ipl.fail <- FALSE

Rescale.Morans <- "y"
Plot.fit <- "y"

# fix scale if family allows over/under-dispersion
if (family1 == "poisson") {FixedScale <- "n"} else {FixedScale <- "y"}

if(any(class(model) == "formula")) {formula1 <- as.formula(model)} else {formula1 <- model$formula}
# at <- intersect(names(data2), all.vars(formula1))
# if(length(at) == 0) stop("formula: specified notation is missing")
nn <- nrow(data2)
x <- coord[,1]
y <- coord[,2]
if(length(x) != nn) stop("error in dimension")
logic1 <- identical(as.numeric(x), round(x,0))
logic2 <- identical(as.numeric(y), round(y,0))
if(!logic1 | !logic2) stop("coordinates not integer")

m0 <- glm(formula1, family1, data2)
#res0<-resid(m0,type="pearson")						#using default setting for type of residuals (=deviance, not pearson)
res0 <- resid(m0)
corglm <- correlog(x, y, res0, na.rm = T, increment = increm, resamp = 1, quiet = TRUE)	#reps = 1 since only used to generate residuals, not to test significance

if(corstr=="independence"){
	# independence model == normal GLM model
		ashort <- 0
		A <- 0
		fitted <- fitted(m0)
		#resids<-resid(m0,type="pearson")
		resids <- resid(m0)			#using default setting for type of residuals (=deviance, not pearson)
		b <- summary(m0)$coefficients[,1]
		s.e. <- summary(m0)$coefficients[,2]
		z <- summary(m0)$coefficients[,3]
		p <- summary(m0)$coefficients[,4]
		scale <- summary(m0)$dispersion
		beta <- cbind(b,s.e.,z,p)
		if(family1 == "gaussian") {colnames(beta) <- c("Estimate", "Std.Err", "t", "p")}
		if(family1 == "binomial" | family1 == "poisson") {colnames(beta) <- c("Estimate", "Std.Err", "z", "p")}
		cat("---","\n","Coefficients:","\n")
		printCoefmat(beta, has.Pvalue = TRUE, P.values = TRUE)
		cat("---","\n")
		cat("Estimated scale parameter: ",scale,"\n")
		cat("Working correlation parameter(s): ", ashort,"\n" )
	# using geeglm function here to model without any spatial data
		mgee <- geeglm(formula = formula1, family = family1, data = data2, id = id2, corstr = "independence")
    mgee.summary <- summary(mgee)
		print(mgee$geese$vbeta.naiv)
		ind.var <- mgee$geese$vbeta.naiv
		Xmat <- mgee$geese$X
		response <- mgee$y
		ind.beta <- mgee$geese$beta
		DecayConstant <- c(0,0)
		assign(paste(model.name,"IndepMod",sep="_"), list(model.name, mgee, mgee$call, NA, ind.var, ind.beta, mgee$fitted.values, mgee.summary, beta, DecayConstant, "none", Xmat, response), envir = .GlobalEnv)
}


if(corstr == "fixed"){
  R <- create.Saana.correlation.matrix(data3 = data2, formula2 = model, family1 = family1)
  ResponseVariable <<- ResponseVariable.original
  data2 <- data.frame(data2, id1)
	
	# error handler here to account for separation of data
	if (class(model)[1] %in% c("glm", "lm")) {model <- model$formula}
  temp.formula1 <- as.character(model)[-1]
  temp.formula2 <- paste(temp.formula1[1], " ~ ", paste(temp.formula1[-1], collapse= " + "), sep = "")
  model <- as.formula(temp.formula2)
  
  if (FixedScale == "y") {mgee <- try(gee(formula = model, family = family1, data = data2, id = id1, R = R, corstr = "fixed", scale.fix = TRUE, scale.value = 1), silent = TRUE)} else {mgee <- try(gee(formula = model, family = family1, data = data2, id = id1, R = R, corstr = "fixed", scale.fix = FALSE), silent = TRUE)}
	if(inherits(mgee, "try-error")) {
		assign(paste(model.name, "FixedMod", sep="_"), list(model.name, "Error", mgee[1], NA, NA, NA, NA, NA, NA, NA, NA), envir = .GlobalEnv)
		print("ERROR IN GEE")
		cat(mgee[1])
		separation.error<<-TRUE
		UsingRobust <<- FALSE
	} else {
		mgee.summary <- summary(mgee)
		mgee.coefficients <- summary(mgee)$coefficients
		b <- mgee$coeff
		fitted <- mgee$fitted
	# this was my initial code to calculate residuals, but I've reverted to res.gee() function
		#if (family1=="binomial") {resids<-(OccurCode-fitted)/sqrt(fitted*(1-fitted))}
		#if (family1=="poisson") {resids<-(OccurCode-fitted)/sqrt(fitted)}
		#following three lines not necessary when not using groups
  		
		res <- res.gee(model, family1, data2, nn, b = mgee$coeff, R = R)
		#fitted<-res$fitted
		resids <- res$resid
		
		# in case of negative variance values in the diagonal....
		if (any(is.nan(mgee.coefficients[, 2 : 3]))) {
		    NaN.matrix <- which(is.nan(mgee.coefficients[, 1 : 3]), arr.ind = TRUE)
		    for (counter6 in 1 : dim(NaN.matrix)[1]) {
		        mgee.coefficients[NaN.matrix[counter6, 1], NaN.matrix[counter6, 2]] <- mgee.coefficients[NaN.matrix[counter6, 1], 2 + NaN.matrix[counter6, 2]]}
        UsingRobust <<- TRUE} else {UsingRobust <<- FALSE}
		
	# alter here in case of overdispersion
		s.e. = mgee.coefficients[,2]	#change to [,4] for robust SE
		z <- mgee.coefficients[,3]	#change to [,5] for robust Z
		p <- rep(NA,nrow(mgee.coefficients))
		for(ii in 1:nrow(mgee.coefficients)){
			if(z[ii]>=0) {p[ii]<-2*(1-pnorm(z[ii]))} else {p[ii]<-2*(pnorm(z[ii]))}}

	# also calculate robust SE, Z and p
		s.e.robust=mgee.coefficients[,4]
		z.robust<-mgee.coefficients[,5]
		p.robust<-rep(NA,nrow(mgee.coefficients))
		for(ii in 1:nrow(mgee.coefficients)){
		  if(is.nan(z.robust[ii]) | is.na(z.robust[ii])) {p.robust[ii] <- 99} else {
    			if(z.robust[ii]>=0) {p.robust[ii]<-2*(1-pnorm(z.robust[ii]))} else {p.robust[ii]<-2*(pnorm(z.robust[ii]))}}}

		scale<-summary(mgee)[[9]]
		beta<-cbind(b,s.e.,z,p)
		beta.robust<-cbind(s.e.robust,z.robust,p.robust)
		colnames(beta) <- c("Estimate", "Std.Err", "z", "p")
		cat("---","\n","Coefficients:","\n")
		printCoefmat(beta,has.Pvalue=TRUE, P.values = TRUE)
		cat("---","\n")
		cat("Estimated scale parameter: ",scale,"\n")
			
		nonind.var<-mgee$robust.variance
		nonind.var.naive<-mgee$naive.variance
		nonind.beta<-mgee$coeff
		Xmat<-mgee$geese$X
		response<-mgee$y
		assign(paste(model.name, "FixedMod", sep="_"), list(model.name, mgee, mgee$call, nonind.var, nonind.var.naive, nonind.beta, mgee$fitted.values, mgee.summary, beta, beta.robust, NA, NA), envir = .GlobalEnv)
		# assign(paste(model.name, "FixedMod_wcs",sep="_"), as.matrix(mgee$working.correlation), envir = .GlobalEnv)
		}	# closes error handling
	}	# closes calculation of "fixed" correlation structure model

  if(!plot.ac & !graph) {ac0 <- NA; ac <- NA}
  
	if(plot.ac | graph | output) {
		x<-coord[,1]
		y<-coord[,2]
		ac.null.mod <- correlog(x, y, data1[, ResponseVariable.column], na.rm = T, increment = increm, resamp = reps, quiet = TRUE)
		ac.null <- ac.null.mod$correlation
		ac.null.p <- ac.null.mod$p
		ac0 <- corglm$correlation
		ac0.p <- corglm$p
		
		if (separation.error == FALSE) {
			ac.mod <- correlog(x, y, resids, na.rm = T, increment = increm, resamp = reps, quiet = TRUE)
			ac <- ac.mod$correlation
			ac.p <- ac.mod$p}}
	if(plot.ac){
		cat("---","\n")
		cat("Autocorrelation for raw data","\n")
		print(cbind(ac.null, ac.null.p))
		cat("Autocorrelation for glm.residuals","\n")
		print(cbind(ac0, ac0.p))
		if (separation.error==FALSE) {
			cat("Autocorrelation for gee.residuals","\n")
			print(cbind(ac, ac.p))}}

		#calculate global autocorrelation statistic	
		# raw data
		raw.temp.global.moran.gr <- moran.mc(x = data1[, ResponseVariable.column], listw = nb, nsim = reps, zero.policy = FALSE, alternative = "greater", na.action = na.fail, spChk = NULL)
		raw.temp.global.moran.ls <- moran.mc(x = data1[, ResponseVariable.column], listw = nb, nsim = reps, zero.policy = FALSE, alternative = "less", na.action = na.fail, spChk = NULL)
		paste("Moran's I = ", round(raw.temp.global.moran.gr$statistic,2), ", p = ", round(min(raw.temp.global.moran.gr$p.value, raw.temp.global.moran.ls$p.value), 3), sep = "")
		# GLM residuals
		glm.temp.global.moran.gr <- moran.mc(x = res0, listw = nb, nsim = reps, zero.policy = FALSE, alternative = "greater", na.action = na.fail, spChk = NULL)
		glm.temp.global.moran.ls <- moran.mc(x = res0, listw = nb, nsim = reps, zero.policy = FALSE, alternative = "less", na.action = na.fail, spChk = NULL)
		paste("Moran's I = ", round(glm.temp.global.moran.gr$statistic,2), ", p = ", round(min(glm.temp.global.moran.gr$p.value,glm.temp.global.moran.ls$p.value),3), sep="")
		# GEE residuals
		if (separation.error == FALSE) {
			gee.temp.global.moran.gr <- moran.mc(x=resids, listw=nb, nsim=reps, zero.policy=FALSE, alternative="greater", na.action=na.fail, spChk=NULL)
			gee.temp.global.moran.ls <- moran.mc(x=resids, listw=nb, nsim=reps, zero.policy=FALSE, alternative="less", na.action=na.fail, spChk=NULL)
			paste("Moran's I = ", round(gee.temp.global.moran.gr$statistic,2), ", p = ", round(min(gee.temp.global.moran.gr$p.value,gee.temp.global.moran.ls$p.value),3), sep="")
			assign(paste(model.name,"autocorrelation",sep = "_"), data.frame("Species" = rep(1, length(ac.null)), "Model" = rep(noquote(as.character(formula1[3])), length(ac.null)), "Correlation structure" = rep(corstr, length(ac.null)), "Mean distance" = round(ac.null.mod$mean,0), "RawData.AC" = ac.null, "RawData.p" = ac.null.p, "Raw.Morans.I.Stat" = rep(raw.temp.global.moran.gr$statistic, length(ac.null)), "RawData.Morans.I.p" = rep(min(raw.temp.global.moran.gr$p.value, raw.temp.global.moran.ls$p.value), length(ac.null)), "GLM.residuals.AC"=ac0, "GLM.residuals.p"=ac0.p, "GLM.residuals.Morans.I. stat"=rep(glm.temp.global.moran.gr$statistic, length(ac.null)), "GLM.residuals.Morans.I.p" = rep(min(glm.temp.global.moran.gr$p.value, glm.temp.global.moran.ls$p.value), length(ac.null)), "GEE.residuals.AC" = ac,
				 "GEE.residuals.p" = ac.p, "GEE.residuals.Moran.I.stat" = rep(gee.temp.global.moran.gr$statistic,length(ac.null)), "GEE.residuals.Moran.I.p" = rep(min(gee.temp.global.moran.gr$p.value, gee.temp.global.moran.ls$p.value), length(ac.null))), envir = .GlobalEnv)
		} else {
			assign(paste(model.name,"autocorrelation",sep="_"), data.frame("Species"=rep(1, length(ac.null)), "Model number"=model.counter, "Model"=rep(noquote(as.character(formula1[3])), length(ac.null)), "Correlation structure"=rep(corstr,length(ac.null)), "Mean distance"=round(ac.null.mod$mean,0), "RawData.AC"=ac.null, "RawData.p"=ac.null.p, "Raw.Morans.I.Stat"=rep(raw.temp.global.moran.gr$statistic, length(ac.null)), "RawData.Morans.I.p"= rep(min(raw.temp.global.moran.gr$p.value, raw.temp.global.moran.ls$p.value), length(ac.null)), "GLM.residuals.AC"=ac0, "GLM.residuals.p"=ac0.p, "GLM.residuals.Morans.I. stat"=rep(glm.temp.global.moran.gr$statistic, length(ac.null)), "GLM.residuals.Morans.I.p" = rep(min(glm.temp.global.moran.gr$p.value, glm.temp.global.moran.ls$p.value), length(ac.null))), envir = .GlobalEnv)
		}

	if(output){
		if (file.exists(paste(OutputFile,"_autocorrelation.txt", sep=""))) {
			write.table(eval(parse(text=paste(model.name,"autocorrelation",sep="_"))), file = paste(OutputFile,"_autocorrelation.txt", sep=""), append = TRUE, sep="\t", dec=".", quote = FALSE, eol="\n", na = "NA", col.names=FALSE, row.names=FALSE)
			} else {
			write.table(eval(parse(text=paste(model.name,"autocorrelation",sep="_"))), file = paste(OutputFile,"_autocorrelation.txt", sep=""), append = FALSE, sep="\t", dec=".", quote = FALSE, eol="\n", na = "NA", col.names=TRUE, row.names=FALSE)}
		}

	if(graph){
		tiff(paste(OutputFileName,1 , eval(expression(paste("Mod",model.counter, sep=""))), "residual_autocorr.tif",sep="_"), units="cm", height=16, width=24, res=60)
		par(mfrow=c(1,2))
		Xmax <- 300
		plot(ac.null.mod$mean.of.class, ac.null.mod$correlation, typ="b", ylim=c(-1,1), xlim=c(0, Xmax), ylab="Autocorrelation of residuals", xlab=paste("Lag distance (bins of ", 1 / increm," m)",sep=""), pch=22, col="dark green")
		points(corglm$mean.of.class,corglm$correlation, pch=1,type="b", col="black")
		if (separation.error==FALSE) {
			points(ac.mod$mean.of.class, ac.mod$correlation,pch=2,type="b", col="red")}
		abline(h=0)
		leg<-c("Raw data", "GLM residuals","GEE residuals")
		legend("topright",leg, pch=c(22,1,2), bty="n", col=c("dark green", "black", "red"))

		if (separation.error == FALSE) {
			y1<-min(min(ac0),min(ac))-.1; y2<-max(max(ac0),max(ac))+.1
			} else{
			y1<-min(ac0)-.1; y2<-max(ac0)+.1}
		plot(corglm$mean.of.class,corglm$correlation,type="b",ylim=c(y1,y2), xlim=c(0, Xmax), ylab="Autocorrelation of residuals", xlab=paste("Lag distance (bins of ", 1 / increm," m)",sep=""), main=formula1)
		points(corglm$mean.of.class[ac0.p<=0.05], ac0[ac0.p<=0.05], pch=16, type="p", col="black")
		if (separation.error==FALSE) {
			points(ac.mod$mean.of.class, ac.mod$correlation,pch=2,type="b", col="red")
			points(ac.mod$mean.of.class[ac.mod$p<=0.05], ac[ac.mod$p<=0.05], pch=17, type="p", col="red")}
		abline(h=0)
		leg<-c("GLM residuals","GEE residuals")
		legend("topright",leg, pch=c(1,2), bty="n", col=c("black", "red"))
		dev.off()
		}

call<-match.call()
#list(call=call,b=b,s.e.=s.e.,z=z,p=p,
#scale=scale,fitted=fitted,resid=resid,w.ac=ashort,W.ac=A,
#ac.glm=ac0,ac=ac)
}

#########################################################################
res.gee <- function(formula1, family1 = gaussian, data2 , n, clusz = NA, zcor = NA, a = NA, b ,R = NA) {
#########################################################################
# A function to calculate fitted values and residuals
# for Generalized Estimating Equation Models
# for gaussian or binary data (with logit link) or Poisson data (log link)
# Arguments
# formula a formula expression
# family "gaussian", "binomial", "poisson" are allowed
# "binomial" = binary
# data a data frame
# n for maximal cluster size n*n
# clusz an object of class "clus.sz"
# zcor an object of class "genZcor"
# a a vector of correlation parameters
# for clusters only
# as an object of class "a.gee"
# b a vector of regression parameters beta
# R a square matrix of correlation parameters
# for full dimension (=number of observations) only
#########################################################################

  l <- dim(data2)[2]
  ieo <- data2[,l-1]
  if(n != dim(data2)[1]) {
    # irrelevant code deleted here, since using fixed correlation structure
  }
  if(n == dim(data2)[1]) v <- R
  ww<-solve(v)
  s.geese<-svd(ww)
  d.geese<-diag(sqrt(s.geese$d))
  w<-s.geese$u%*%d.geese%*%t(s.geese$u)
  x.matrix<-model.matrix(formula1,data2)
  fitted<-x.matrix%*%b
  fitted<-fitted[1:length(ieo)]
  if(family1 == "poisson") fitted <- exp(fitted)
  if(family1 == "binomial") fitted <- exp(fitted)/(1 + exp(fitted))
  if(family1 == "gaussian") rgeese <- model.frame(formula1, data2)[[1]] - fitted
  if(family1 == "poisson") rgeese <- (model.frame(formula1, data2)[[1]] - fitted)/sqrt(fitted)
  if(family1 == "binomial") rgeese <-(model.frame(formula1, data2)[[1]] - fitted)/sqrt(fitted*(1 - fitted))
  rsgeese <- w%*%rgeese
  resgeeseo <- rsgeese[1 : length(ieo)]
  list(fitted = fitted, resid = resgeeseo)
}


############################################
###QIC FUNCTION FOR MODEL-SELECTION OF SPATIAL MODELS ##########
# gee.ind        GEE with an independence structure
# gee.nonind     GEE fitted (FIXED or EXCHANGEABLE structure)

QIC <- function(gee.ind, gee.nonind, family1 = "binomial") {
	require(MASS)
	response <- gee.ind[[13]]
	Xmat <- gee.ind[[12]]
	if(identical(gee.ind, gee.nonind)) {nonind.var <- gee.nonind[[5]]} else {nonind.var <- gee.nonind[[4]]}
	nonind.var.naive <- gee.nonind[[5]]
	ind.var <- gee.ind[[5]]
	nonind.beta <- gee.nonind[[6]]
	ind.beta <- gee.ind[[6]]

	# Calculate the QIC: Make difference between gee & geepack
	trace.user <- sum(diag(ginv(ind.var) %*% nonind.var))
	trace.user.naive <- sum(diag(ginv(ind.var) %*% nonind.var.naive))

	# Get the fitted values
	mu.ind <- exp(Xmat %*% ind.beta)
  # mu.user <- exp(Xmat %*% nonind.beta)      # original code
  mu.user <- gee.nonind[[7]]		# my modification


	# quasi-likelihood (Poisson, otherwise Binomial)
	if (family1 == "gaussian") {quasi.ind <- sum(((response - mu.user)^2) / -2)}
	if (family1 == "poisson") {quasi.ind <- sum((response * log(mu.user)) - mu.user - (response * (log(response + 0.00001) - 1)))}
	if (family1 == "binomial") {quasi.ind <- sum((response * log(mu.user/(1 - mu.user))) + log(1 - mu.user))}
	if (family1 != "poisson" & family1 != "binomial" & family1 != "gaussian") {print("The quasilikelihood for this family needs to be specified in the QIC function")}

	# Calculate QIC
	qic.user <- 2 * (trace.user - quasi.ind)
	# output values outside of function
	QIC.values <- NULL
	assign("QIC.values", c(QIC = qic.user, quasi.lik = quasi.ind, trace = trace.user, trace.naive = trace.user.naive), envir = .GlobalEnv)
}