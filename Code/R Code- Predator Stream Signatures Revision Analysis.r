##### New analyses based on Ecological Monographs Reviews #####

##### Set working directory, read in data #####

## Working directory
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis')

## Load packages and functions
source('C:/Users/jmuehlbauer/Documents/Misc/Trainings & Info/Stats & Data Analysis/R/AICforLogResponse.r')


## Data
mega8 <- read.csv('R Code & Input Data/PredatorMega8.csv')
breaks <- read.csv('R Code & Input Data/Breaks.csv')

## Remove, rename, and round some columns
neo0 <- subset(mega8, select = c(Code, SiteBank, Group, Trophic, Region, Width, OrderClass, Geomorph, Banks, VegFld, preForest, Forest, Wall, DistShift, Effect, logDist, logEffect))
neo1 <- neo0
colnames(neo1)[match(c('OrderClass', 'VegFld', 'preForest', 'DistShift', 'logEffect', 'logDist'), colnames(neo0))] <- c('Order', 'Veg', 'PreF', 'Dist', 'logE', 'logD')
neo1[, c('Effect', 'logE', 'logD')] <- round(neo1[, c('Effect', 'logE', 'logD')], 4)
 
 
##### Add variables relative to breaks #####

## Dummy variables pre/post breaks
neo1$AnteP <- ifelse(neo1$Dist < neo1$PreF | neo1$PreF %in% NA, 1, 0)
neo1$PostP <- ifelse(neo1$AnteP == 0, 1, 0)
neo1$AnteF <- ifelse(neo1$Dist < neo1$Forest | neo1$Forest %in% NA, 1, 0)
neo1$PostF <- ifelse(neo1$AnteF == 0, 1, 0)
neo1$AnteW <- ifelse(neo1$Dist < neo1$Wall | neo1$Wall %in% NA, 1, 0)
neo1$PostW <- ifelse(neo1$AnteW == 0, 1, 0)

## Distance post-breaks
neo1$DistP <- ifelse(neo1$PostP == 1, neo1$Dist - neo1$PreF, 0)
neo1$logP <- round(ifelse(neo1$PostP == 1, log(neo1$DistP + .05), 0), 4)
neo1$DistF <- ifelse(neo1$PostF == 1, neo1$Dist - neo1$Forest, 0)
neo1$logF <- round(ifelse(neo1$PostF == 1, log(neo1$DistF + .05), 0), 4)
neo1$DistW <- ifelse(neo1$PostW == 1, neo1$Dist - neo1$Wall, 0)
neo1$logW <- round(ifelse(neo1$PostW == 1, log(neo1$DistW + .05), 0), 4)


##### Run simple, fixed effects regressions for all SiteBanks #####

## Put data in lists by SiteBank
list.neo <- list()
list.neo <- lapply(levels(neo1$SiteBank), function(x) neo1[neo1$SiteBank == x, ])
	names(list.neo) <- levels(neo1$SiteBank)

## Create function for running regressions with pre-forest breaks
Pfx <- function(X) {
	sb <- X
	## Compute lm, exp, pow, and break equations for each SiteBank
	lin0 <- lm(Effect ~ Dist, data = sb)
		lin1 <- lm(Effect ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = sb)
		lin2 <- lm(Effect ~ PostP + Dist + PostP:DistP, data = sb)
	exp0 <- lm(logE ~ Dist, data = sb)
		exp1 <- lm(logE ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = sb)
		exp2 <- lm(logE ~ PostP + Dist + PostP:DistP, data = sb)
	pow0 <- lm(logE ~ logD, data = sb)
		pow1 <- lm(logE ~ 0 + AnteP + PostP + AnteP:logD + PostP:logP, data = sb)
		pow2 <- lm(logE ~ PostP + logD + PostP:logP, data = sb)
	## Get AIC values for each regression
	aics <- as.data.frame(c('Model0', 'Model1', 'Model2'))
		colnames(aics) <- 'Model'
		aics$Lin <- sapply(list(lin0, lin1, lin2), AIC)
		aics$Exp <- trueaics <- sapply(list(exp0, exp1, exp2), function(x) trueAIC(x)[1])
		aics$Pow <- trueaics <- sapply(list(pow0, pow1, pow2), function(x) trueAIC(x)[1])
			aics[,-1] <- round(aics[,-1], 2)
	## Make sure line actually decays
	slopes <- aics
		slopes[1, 2:4] <- sapply(list(lin0, exp0, pow0), function(x) ifelse(x$coefficients[2] < 0, 'Neg' , 'Pos'))
		slopes[2, 2:4] <- sapply(list(lin1, exp1, pow1), function(x) ifelse(x$coefficients[3] < 0, 'Neg' , 'Pos'))
		slopes[3, 2:4] <- sapply(list(lin2, exp2, pow2), function(x) ifelse(x$coefficients[3] < 0, 'Neg' , 'Pos'))
		
	## Get best AIC model
	best <- names(unlist(aics))[which(aics == min(aics[, -1]))]
	lst <- list(lin0, lin1, lin2, exp0, exp1, exp2, pow0, pow1, pow2, aics, slopes, best)
	names(lst) <- c('lin0', 'lin1', 'lin2', 'exp0', 'exp1', 'exp2', 'pow0', 'pow1', 'pow2', 'aics', 'slopes', 'best')
	lst
	}
	
## Forest breaks
Ffx <- function(X) {
	sb <- X
	## Compute lm, exp, pow, and break equations for each SiteBank
	lin0 <- lm(Effect ~ Dist, data = sb)
		lin1 <- lm(Effect ~ 0 + AnteF + PostF + AnteF:Dist + PostF:DistF, data = sb)
		lin2 <- lm(Effect ~ PostF + Dist + PostF:DistF, data = sb)
	exp0 <- lm(logE ~ Dist, data = sb)
		exp1 <- lm(logE ~ 0 + AnteF + PostF + AnteF:Dist + PostF:DistF, data = sb)
		exp2 <- lm(logE ~ PostF + Dist + PostF:DistF, data = sb)
	pow0 <- lm(logE ~ logD, data = sb)
		pow1 <- lm(logE ~ 0 + AnteF + PostF + AnteF:logD + PostF:logF, data = sb)
		pow2 <- lm(logE ~ PostF + logD + PostF:logF, data = sb)
	## Get AIC values for each regression
	aics <- as.data.frame(c('Model0', 'Model1', 'Model2'))
		colnames(aics) <- 'Model'
		aics$Lin <- sapply(list(lin0, lin1, lin2), AIC)
		aics$Exp <- trueaics <- sapply(list(exp0, exp1, exp2), function(x) trueAIC(x)[1])
		aics$Pow <- trueaics <- sapply(list(pow0, pow1, pow2), function(x) trueAIC(x)[1])
			aics[,-1] <- round(aics[,-1], 2)
	## Make sure line actually decays
	slopes <- aics
		slopes[1, 2:4] <- sapply(list(lin0, exp0, pow0), function(x) ifelse(x$coefficients[2] < 0, 'Neg' , 'Pos'))
		slopes[2, 2:4] <- sapply(list(lin1, exp1, pow1), function(x) ifelse(x$coefficients[3] < 0, 'Neg' , 'Pos'))
		slopes[3, 2:4] <- sapply(list(lin2, exp2, pow2), function(x) ifelse(x$coefficients[3] < 0, 'Neg' , 'Pos'))
		
	## Get best AIC model
	best <- names(unlist(aics))[which(aics == min(aics[, -1]))]
	lst <- list(lin0, lin1, lin2, exp0, exp1, exp2, pow0, pow1, pow2, aics, slopes, best)
	names(lst) <- c('lin0', 'lin1', 'lin2', 'exp0', 'exp1', 'exp2', 'pow0', 'pow1', 'pow2', 'aics', 'slopes', 'best')
	lst
	}

## Wall breaks	
	Wfx <- function(X) {
	sb <- X
	## Compute lm, exp, pow, and break equations for each SiteBank
	lin0 <- lm(Effect ~ Dist, data = sb)
		lin1 <- lm(Effect ~ 0 + AnteW + PostW + AnteW:Dist + PostW:DistW, data = sb)
		lin2 <- lm(Effect ~ PostW + Dist + PostW:DistW, data = sb)
	exp0 <- lm(logE ~ Dist, data = sb)
		exp1 <- lm(logE ~ 0 + AnteW + PostW + AnteW:Dist + PostW:DistW, data = sb)
		exp2 <- lm(logE ~ PostW + Dist + PostW:DistW, data = sb)
	pow0 <- lm(logE ~ logD, data = sb)
		pow1 <- lm(logE ~ 0 + AnteW + PostW + AnteW:logD + PostW:logW, data = sb)
		pow2 <- lm(logE ~ PostW + logD + PostW:logW, data = sb)
	## Get AIC values for each regression
	aics <- as.data.frame(c('Model0', 'Model1', 'Model2'))
		colnames(aics) <- 'Model'
		aics$Lin <- sapply(list(lin0, lin1, lin2), AIC)
		aics$Exp <- trueaics <- sapply(list(exp0, exp1, exp2), function(x) trueAIC(x)[1])
		aics$Pow <- trueaics <- sapply(list(pow0, pow1, pow2), function(x) trueAIC(x)[1])
			aics[,-1] <- round(aics[,-1], 2)
	## Make sure line actually decays
	slopes <- aics
		slopes[1, 2:4] <- sapply(list(lin0, exp0, pow0), function(x) ifelse(x$coefficients[2] < 0, 'Neg' , 'Pos'))
		slopes[2, 2:4] <- sapply(list(lin1, exp1, pow1), function(x) ifelse(x$coefficients[3] < 0, 'Neg' , 'Pos'))
		slopes[3, 2:4] <- sapply(list(lin2, exp2, pow2), function(x) ifelse(x$coefficients[3] < 0, 'Neg' , 'Pos'))
		
	## Get best AIC model
	best <- names(unlist(aics))[which(aics == min(aics[, -1]))]
	lst <- list(lin0, lin1, lin2, exp0, exp1, exp2, pow0, pow1, pow2, aics, slopes, best)
	names(lst) <- c('lin0', 'lin1', 'lin2', 'exp0', 'exp1', 'exp2', 'pow0', 'pow1', 'pow2', 'aics', 'slopes', 'best')
	lst
	}

## Run functions
lm.P <- lapply(list.neo, function(x) Pfx(x))
lm.F <- lapply(list.neo, function(x) Ffx(x))
lm.W <- lapply(list.neo, function(x) Wfx(x))

## See which functions are the best
BestP <- table(unlist(lapply(lm.P, function(x) x$best)))
BestF <- table(unlist(lapply(lm.F, function(x) x$best)))
BestW <- table(unlist(lapply(lm.W, function(x) x$best)))
	## Linear always sucks, power generally better than exponential, but adding more model complexity tends not to help anything.
	
## See which sites might need to be tossed
Slopes <- lapply(lm.P, function(x) x$slope)
Slope <- names(which(unlist(lapply(lm.P, function(x) x$slope)) == 'Positive'))


##### Try simple models again, using only those SiteBanks with actual decays #####

## Get rid of non-decaying sites (based on simple lm: lin0)
list1.neo <- list.neo[-c(which(names(list.neo) %in% Slope))]
lm1.neo <- lapply(list1.neo, function(x) regfx(x))
Best1 <- table(unlist(lapply(lm1.neo, function(x) x$best)))


##### THIS ISN'T WORKING!!! #####

##### Try fitting regressions just on a great site #####

## Subset data for BOLZ1RB
bolz<-neo1[neo1$SiteBank=='BOLZ1RB',]

## Run function
lm.B <- Pfx(bolz)
plot(bolz$Dist,exp(predict(lm.B$pow1)))
plot(bolz$Dist,exp(predict(lm.B$pow2)))
plot(bolz$Dist,exp(predict(lm.B$pow1)))

## Try quantile regression and logistic quantile regession
library(quantreg)
lmLB <- lm(Effect ~ Dist, data = bolz)
lmLQ <- rq(Effect ~ Dist, data = bolz)
bolz$Logit <- round(log((bolz$Effect-(min(bolz$Effect-.05)))/((max(bolz$Effect)+.05)-bolz$Effect)), 4)
lmLL <- rq(Logit ~ Dist, data = bolz)
lmPB <- lm(logE ~ logD, data = bolz)
lmPQ <- rq(logE ~ logD, data = bolz)
bolz$LBpred <- round(predict(lmLB), 4)
bolz$LQpred <- round(predict(lmLQ), 4)
bolz$LLpred <- round((exp(predict(lmLL))*(max(bolz$Effect)+.05)+min(bolz$Effect)-.05)/(1+exp(predict(lmLL))), 4)
bolz$PBpred <- round(exp(predict(lmPB)), 4)
bolz$PQpred <- round(exp(predict(lmPQ)), 4)
plot(bolz$Dist, bolz$Effect)
lines(bolz$Dist, bolz$LBpred, col = 1, type = 'l')
lines(bolz$Dist, bolz$LQpred, col = 2, type = 'l')
lines(bolz$Dist, bolz$LLpred, col = 3, type = 'l')
lines(bolz$Dist, bolz$PBpred, col = 4, type = 'l')
lines(bolz$Dist, bolz$PQpred, col = 5, type = 'l')

bolz2 <- bolz[bolz$Dist < 100, ]
lm2LB <- lm(Effect ~ Dist, data = bolz2)
lm2LQ <- rq(Effect ~ Dist, data = bolz2)
bolz2$Logit <- round(log((bolz2$Effect-(min(bolz2$Effect-.05)))/((max(bolz2$Effect)+.05)-bolz2$Effect)), 4)
lm2EB <- lm(logE ~ Dist, data = bolz2)
lm2PB <- lm(logE ~ logD, data = bolz2)
lm2EQ <- rq(logE ~ Dist, data = bolz2, tau = .90)
lm2PQ <- rq(logE ~ logD, data = bolz2, tau = .90)
lm2EL <- rq(Logit ~ Dist, data = bolz2, tau = .90)
lm2PL <- rq(Logit ~ logD, data = bolz2, tau = .90)
	explm2EL <- round((exp(predict(lm2EL, newdata = data.frame(Dist = seq(0, 120, .5))))*(max(bolz2$Effect)+.05)+(min(bolz2$Effect)-.05))/(1+exp(predict(lm2EL, newdata = data.frame(Dist = seq(0, 120, .5))))), 4)
	explm2PL <- round((exp(predict(lm2PL, newdata = data.frame(logD = log(seq(.05, 120, .5)))))*(max(bolz2$Effect)+.05)+(min(bolz2$Effect)-.05))/(1+exp(predict(lm2PL, newdata = data.frame(logD = log(seq(.05, 120, .5)))))), 4)
plot(bolz$Dist, bolz$Effect)
lines(seq(.05, 120, .5), exp(predict(lm2PB, newdata = data.frame(logD = log(seq(.05, 120, .5))))), col = 1, type = 'l')
lines(seq(0, 120, .5), exp(predict(lm2EB, newdata = data.frame(Dist = seq(0, 120, .5)))), col = 2, type = 'l')
lines(seq(.05, 120, .5), exp(predict(lm2PQ, newdata = data.frame(logD = log(seq(.05, 120, .5))))), col = 3, type = 'l')
lines(seq(0, 120, .5), exp(predict(lm2EQ, newdata = data.frame(Dist = seq(0, 120, .5)))), col = 4, type = 'l')
lines(seq(0, 120, .5), explm2EL, col = 5, type = 'l')
lines(seq(.05, 120, .5), explm2PL, col = 6, type = 'l')
legend('top', legend = c('lm Power', 'lm Exponential', '90th quantile Power', '90th quantile Exponential', '90th quantile Logistic dist', '90th quantile Logistic logdist'), col = 1:6, lty = 1, bty = 'n')


##### Try fitting beta regressions #####

## Load library
library(betareg)

## Convert any 0 or 1 effects to 0.01 and 0.99
bolz2$Effect2 <- ifelse(bolz2$Effect == 1, .99, ifelse(bolz2$Effect == 0, .01, bolz2$Effect))

## Run beta regression
lm2Beta <- betareg(Effect2 ~ Dist, data = bolz2, link = 'log')
ddat <- data.frame(Dist = seq(.01, 120, .01))
	ddat$logD <- seq(.01, 120, .01)
lines(ddat$Dist, predict(lm2Beta, ddat))
lm2BetaPow <- betareg(Effect2 ~ logD, data = bolz2, link = 'log')
	lines(ddat$logD, predict(lm2BetaPow, ddat))
	

##### Try fitting various functions for pre-break data #####

## Subset data to before break
bolz2 <- bolz[bolz$Dist < 100, ]

## Linear regressions
lin2L <- glm(Effect2 ~ Dist, data = bolz2)
lin2E <- glm(logE2 ~ Dist, data = bolz2)
lin2P <- glm(logE2 ~ logD, data = bolz2)

## Beta regessions
beta2E <- betareg(Effect2 ~ Dist, data = bolz2, link = 'log')
beta2P <- betareg(Effect2 ~ logD, data = bolz2, link = 'log')

## Compare AICs
AICs <- data.frame(Model = c('Linear0', 'Exponential0', 'Power0', 'Beta_Dist0', 'Beta_logDist0', 'Linear1', 'Exponential1', 'Power1', 'Beta_Dist1', 'Beta_logDist1', 'Linear2', 'Exponential2', 'Power2', 'Beta_Dist2', 'Beta_logDist2'))
	AICs$AICpreBreak <- round(c(AIC(lin2L), trueAIC(lin2E)[1], trueAIC(lin2P)[2], AIC(beta2E), AIC(beta2P)))
	## Looks like beta better than the basic linear form, but not the expon or power?
	## Do I trust the AIC transform?

## Plot the various curves
ddat <- data.frame(Dist = seq(.01, 120, .01))
	ddat$logD <- seq(.01, 120, .01)
cols <- c('red', 'cyan', 'orange', 'blue', 'green')
plot(bolz$Dist, bolz$Effect, ylab = 'Portion aquatic (relativized abundance*biomass effect)', xlab = 'Distance from river (m)')
	lines(ddat$Dist, predict(lin2L, ddat, type = 'response'), col = cols[1])
	lines(ddat$logD, exp(predict(lin2E, ddat, type = 'response')), col = cols[2])
	lines(ddat$logD, exp(predict(lin2P, ddat, type = 'response')), col = cols[3])
	lines(ddat$Dist, predict(beta2E, ddat, type = 'response'), col = cols[4])
	lines(ddat$logD, predict(beta2P, ddat, type = 'response'), col = cols[5])
	legend('top', legend = AICs$Model[1:5], title = 'Model', col = cols, lty = 1, bty = 'n')
	## Beta regression with log(Distance) is really the only one that seems decent
	
	
##### Try fitting various functions for all data #####

## Convert any 0 or 1 effects to 0.01 and 0.99
bolz$Effect2 <- ifelse(bolz$Effect == 1, .99, ifelse(bolz$Effect == 0, .01, bolz$Effect))
bolz$logE2 <- log(bolz$Effect2)

## Basic Linear regressions
linL0 <- glm(Effect2 ~ Dist, data = bolz)
linE0 <- glm(logE2 ~ Dist, data = bolz)
linP0 <- glm(logE2 ~ logD, data = bolz)

## Basic Beta regessions
betaE0 <- betareg(Effect2 ~ Dist, data = bolz, link = 'log')
betaP0 <- betareg(Effect2 ~ logD, data = bolz, link = 'log')

## Regressions with Charles' equation 1
linL1 <- glm(Effect2 ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = bolz)
linE1 <- glm(logE2 ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = bolz)
linP1 <- glm(logE2 ~ 0 + AnteP + PostP + AnteP:logD + PostP:logP, data = bolz)
betaE1 <- betareg(Effect2 ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = bolz, link = 'log')
betaP1 <- betareg(Effect2 ~ 0 + AnteP + PostP + AnteP:logD + PostP:logP, data = bolz, link = 'log')

## Regressions with Charles' equation 2
linL2 <- glm(Effect2 ~ PostP + Dist + PostP:DistP, data = bolz)
linE2 <- glm(logE2 ~ PostP + Dist + PostP:DistP, data = bolz)
linP2 <- glm(logE2 ~ PostP + logD + PostP:logP, data = bolz)
betaE2 <- betareg(Effect2 ~ PostP + Dist + PostP:DistP, data = bolz, link = 'log')
betaP2 <- betareg(Effect2 ~ PostP + logD + PostP:logP, data = bolz, link = 'log')

## Compare AICs
AICs$AICoverall <- round(c(AIC(linL0), trueAIC(linE0)[1], trueAIC(linP0)[2], AIC(betaE0), AIC(betaP0), AIC(linL1), trueAIC(linE1)[1], trueAIC(linP1)[2], AIC(betaE1), AIC(betaP1), AIC(linL1), trueAIC(linE1)[1], trueAIC(linP1)[2], AIC(betaE1), AIC(betaP1)))
	## Looks like beta better than the basic linear form, but not the expon or power?
	## Do I trust the AIC transform?
	## Charles' equations aren't helping.

## Plot the various curves
ddat <- data.frame(Dist = seq(.01, 175, .01))
	ddat$logD <- seq(.01, 175, .01)
plot(bolz$Dist, bolz$Effect, ylab = 'Portion aquatic (relativized abundance*biomass effect)', xlab = 'Distance from river (m)')
	lines(ddat$Dist, predict(linL0, ddat, type = 'response'), col = cols[1])
	lines(ddat$logD, exp(predict(linE0, ddat, type = 'response')), col = cols[2])
	lines(ddat$logD, exp(predict(linP0, ddat, type = 'response')), col = cols[3])
	lines(ddat$Dist, predict(betaE0, ddat, type = 'response'), col = cols[4])
	lines(ddat$logD, predict(betaP0, ddat, type = 'response'), col = cols[5])
	legend('top', legend = AICs$Model[1:5], title = 'Model', col = cols, lty = 1, bty = 'n')
	## Beta regression with log(Distance) is really the only one that seems decent
	## STOPPED HERE. Need to repeat code in lines 317-322 for Charles' model 1 and 2 types (probably same color, different line type). Also, WTF with the power curve?