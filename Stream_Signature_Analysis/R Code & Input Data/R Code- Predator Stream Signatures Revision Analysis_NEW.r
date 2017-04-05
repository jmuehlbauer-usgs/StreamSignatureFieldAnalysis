##### New analyses based on Ecological Monographs Reviews #####


##### Set working directory, read in data #####

## Working directory
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis')

## Load packages and functions
library(betareg)
library(quantreg)
source('C:/Users/jmuehlbauer/Documents/Misc/Trainings & Info/Stats & Data Analysis/R/AICforLogResponse.r')

## Data
mega8.0 <- read.csv('R Code & Input Data/PredatorMega8.csv')

## Convert break locations based on DistShift, not Dist
mega8 <- mega8.0
mega8[,c('preForest', 'Forest', 'Wall')] <- mega8.0[,c('preForest', 'Forest', 'Wall')] - (mega8.0$Dist - mega8.0$DistShift)

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


##### Subset data for a characteristic site for ease of viewing #####

## Subset only BOLZ1RB site
bolz<-neo1[neo1$SiteBank=='BOLZ1RB',]

## Convert any 0 or 1 effects to 0.01 and 0.99 for Beta regressions
bolz$Effect2 <- ifelse(bolz$Effect == 1, .99, ifelse(bolz$Effect == 0, .01, bolz$Effect))
bolz$logE2 <- log(bolz$Effect2)

## Calculate empirical logit for logistic quantile regressions
bolz$Logit <- round(log((bolz$Effect-(min(bolz$Effect-.05)))/((max(bolz$Effect)+.05)-bolz$Effect)), 4)

## Subset data to before break
bolz2 <- bolz[bolz$Dist < 128, ]


##### Try fitting various functions for all data #####

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
AICs <- data.frame(Model = c('Linear0', 'Exponential0', 'Power0', 'Beta_Dist0', 'Beta_logDist0', 'Linear1', 'Exponential1', 'Power1', 'Beta_Dist1', 'Beta_logDist1', 'Linear2', 'Exponential2', 'Power2', 'Beta_Dist2', 'Beta_logDist2'))
	AICs$AICoverall <- round(c(AIC(linL0), trueAIC(linE0)[1], trueAIC(linP0)[2], AIC(betaE0), AIC(betaP0), AIC(linL1), trueAIC(linE1)[1], trueAIC(linP1)[2], AIC(betaE1), AIC(betaP1), AIC(linL1), trueAIC(linE1)[1], trueAIC(linP1)[2], AIC(betaE1), AIC(betaP1)))
	## Looks like beta better than the basic linear form, but not the expon or power?
	## Do I trust the AIC transform?
	## Charles' equations aren't helping.


##### Try fitting various functions for pre-break data #####

## Linear regressions
lin2L <- glm(Effect2 ~ Dist, data = bolz2)
lin2E <- glm(logE2 ~ Dist, data = bolz2)
lin2P <- glm(logE2 ~ logD, data = bolz2)

## Beta regessions
beta2E <- betareg(Effect2 ~ Dist, data = bolz2, link = 'log')
beta2P <- betareg(Effect2 ~ logD, data = bolz2, link = 'log')

## Logistic quantile regressions (90% quantile)
quant2E <- rq(Logit ~ Dist, data = bolz2, tau = .90)
quant2P <- rq(Logit ~ logD, data = bolz2, tau = .90)

## Compare AICs
AICs$AICpreBreak <- c(round(c(AIC(lin2L), trueAIC(lin2E)[1], trueAIC(lin2P)[2], AIC(beta2E), AIC(beta2P))), rep(NA, 10))
	## Looks like beta better than the basic linear form, but not the exponential or power?
	## Do I trust the AIC transform?	
	
	
##### Plot the various curves #####

## Prepare predicted values
ddat <- data.frame(Dist = seq(.01, 175, .01))
	ddat$logD <- log(ddat$Dist)
	ddat$AnteP <- ifelse(ddat$Dist < 128, 1, 0)
	ddat$PostP <- ifelse(ddat$AnteP == 0, 1, 0)
	ddat$DistP <- ifelse(ddat$Dist >=  128, ddat$Dist - 128, 0)
	ddat$logP <- ifelse(ddat$Dist >= 128, log(ddat$Dist), 0)
ddat2 <- ddat[ddat$Dist < 120, ]
	
## Set some plotting parameters
par(mfrow = c(2, 2), mar = c(4, 5, .5, .5))
cols <- c(2, 5, 'orange', 4, 3)

## Plot all data, basic models
plot(bolz$Dist, bolz$Effect, ylab = 'Portion aquatic\n(relativized abundance*biomass effect)', xlab = 'Distance from river (m)')
	lines(ddat$Dist, predict(linL0, ddat, type = 'response'), col = cols[1])
	lines(ddat$Dist, exp(predict(linE0, ddat, type = 'response')), col = cols[2])
	lines(ddat$Dist, exp(predict(linP0, ddat, type = 'response')), col = cols[3])
	lines(ddat$Dist, predict(betaE0, ddat, type = 'response'), col = cols[4])
	lines(ddat$Dist, predict(betaP0, ddat, type = 'response'), col = cols[5])
	legend(10, 1, legend = paste(AICs$Model[1:5], ': AIC = ', AICs$AICoverall[1:5], sep = ''), title = 'Model 0', col = cols, lty = 1, bty = 'n')
	## For basic models, beta regression with log(Distance) is really the only one that actually plots a decline, but all are terrible.

## Plot all data, Charles' model 1 formulation
plot(bolz$Dist, bolz$Effect, ylab = 'Portion aquatic\n(relativized abundance*biomass effect)', xlab = 'Distance from river (m)')
	lines(ddat$Dist, predict(linL1, ddat, type = 'response'), col = cols[1], lty = 2)
	lines(ddat$Dist, exp(predict(linE1, ddat, type = 'response')), col = cols[2], lty = 2)
	lines(ddat$Dist, exp(predict(linP1, ddat, type = 'response')), col = cols[3], lty = 2)
	lines(ddat$Dist, predict(betaE1, ddat, type = 'response'), col = cols[4], lty = 2)
	lines(ddat$Dist, predict(betaP1, ddat, type = 'response'), col = cols[5], lty = 2)
	legend(10, 1, legend = paste(AICs$Model[6:10], ': AIC = ', AICs$AICoverall[6:10], sep = ''), title = 'Model 1', col = cols, lty = 2, bty = 'n')
	## Model 1 fits are better, but still vastly overestimate low values and underestimate high ones.

## Plot all data, Charles' model 2 formulation
plot(bolz$Dist, bolz$Effect, ylab = 'Portion aquatic\n(relativized abundance*biomass effect)', xlab = 'Distance from river (m)')
	lines(ddat$Dist, predict(linL2, ddat, type = 'response'), col = cols[1], lty = 3)
	lines(ddat$Dist, exp(predict(linE2, ddat, type = 'response')), col = cols[2], lty = 3)
	lines(ddat$Dist, exp(predict(linP2, ddat, type = 'response')), col = cols[3], lty = 3)
	lines(ddat$Dist, predict(betaE2, ddat, type = 'response'), col = cols[4], lty = 3)
	lines(ddat$Dist, predict(betaP2, ddat, type = 'response'), col = cols[5], lty = 3)
	legend(10, 1, legend = paste(AICs$Model[11:15], ': AIC = ', AICs$AICoverall[11:15], sep = ''), title = 'Model 2', col = cols, lty = 3, bty = 'n')
	## Model 2 fits are essentially equivalent to model 1 fits.

## Plot pre-break data, simple models
par(mar = c(4, 5, .5, .5))
plot(bolz$Dist, bolz$Effect, ylab = 'Portion aquatic\n(relativized abundance*biomass effect)', xlab = 'Distance from river (m)', pch = 19)
	lines(ddat2$Dist, predict(lin2L, ddat2, type = 'response'), col = cols[1])
	lines(ddat2$Dist, exp(predict(lin2E, ddat2, type = 'response')), col = cols[2])
	lines(ddat2$Dist, exp(predict(lin2P, ddat2, type = 'response')), col = cols[3])
	lines(ddat2$Dist, predict(beta2E, ddat2, type = 'response'), col = cols[4])
	lines(ddat2$Dist, predict(beta2P, ddat2, type = 'response'), col = cols[5])
	legend(10, 1, legend = paste(AICs$Model[1:5], ': AIC = ', AICs$AICoverall[1:5], sep = ''), title = 'pre-Forest only', col = cols, lty = 1, bty = 'n')
	## Essentially identical to model 1 and 2 fits before the break. So at lease those models are working before the break, even if they are still terrible.

## Plot a line with un-log-transformed distances
	## Erroenous model usage, but it looks nice.
	## How do I replicate this in a model?
ddat3 <- ddat2
ddat3$logD <- ddat3$Dist
	lines(ddat3$Dist, predict(beta2P, ddat3, type = 'response'), col = 8, lty = 2)
	legend(10, .7, legend = 'Gamed values', col = 8, lty = 2, bty = 'n')
	
## Plot lines for logistic quantile regressions
predE <- round((exp(predict(quant2E, ddat2))*(max(bolz2$Effect)+.05)+(min(bolz2$Effect)-.05))/(1+exp(predict(quant2E, ddat2))), 4)
predP <- round((exp(predict(quant2P, ddat2))*(max(bolz2$Effect)+.05)+(min(bolz2$Effect)-.05))/(1+exp(predict(quant2P, ddat2))), 4)
	lines(ddat2$Dist, predE, col = 'purple', lty = 3)
	lines(ddat2$Dist, predP, col = 6, lty = 3)
	legend(10, .6, legend = c('Logistic 90th quantile, Dist', 'Logistic 90th quantile, logDist'), col = c('purple', 6), lty = 3, bty = 'n')
	