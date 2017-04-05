## To convert AIC values from a log response for the purpose of direct comparison with AIC values from a linear response
	## 10 August 2016 by J.D. Muehlbauer
	## Modification of code from Jack Weiss found here: http://www.unc.edu/courses/2008fall/ecol/563/001/docs/lectures/lecture11.htm

##### Input variables #####
	## Model: Saved model output from lm, glm, lme4, or nlme, and possibly others (no guarantees! But if it doesn't work, you'll get an error!)
	## IC: The information criterion to calculate. Accepts "AIC" or "BIC". "AIC is the default.
trueAIC<-function(Model, IC = 'AIC'){
	modclass <- class(Model)
	if('lme' %in% modclass){
		log.y <- getResponse(Model)
		K <- sum(Model$dims$ngrps * Model$dims$ncol)
	} else{
		if('lmerMod' %in% modclass){
			log.y <- attributes(Model)$resp$y
			K <- summary(Model)$ngrps + dim(summary(Model)$coefficients)[1] + 1
		} else{
			log.y <- Model$model[,1]
			K <- dim(summary(Model)$coefficients)[1] + 1
			}
		}
	N <- length(log.y)
	sigma2 <- (sum(residuals(Model)^2)) / N
		## Squared residuals divided by sample size is the MLE of the variance
	loglike <- sum(log(dnorm(log.y, mean = predict(Model), sd = sqrt(sigma2)) * 1 / exp(log.y)))
	if(IC == 'AIC'){
		theAIC <- (-2 * loglike + 2 * K)
		theAICc <- theAIC + ((2 * K * (K + 1))/(N - K - 1))
		outAIC <-cbind(theAIC, theAICc, K, N)
			colnames(outAIC) <- c('AIC', 'AICc', 'k', 'n')
	}
	if(IC == 'BIC'){
		theBIC <- (-2 * loglike + K * log(N))
		outAIC <- cbind(theBIC, K, N)
			colnames(outAIC) <- c('BIC', 'k', 'n')
	}
	rownames(outAIC) <- NULL
	outAIC
}
