##### New analyses based on Ecological Monographs Reviews #####


##### Set working directory, read in data #####

## Working directory
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream_Signature/Stream_Signature_Analysis')

## Load packages
library(lme4, quietly = TRUE)

## Load functions
source('R Code & Input Data/AICforLogResponse.r')
source('R Code & Input Data/PlottingFunction_wCurves.r')

## Data
mega8.0 <- read.csv('R Code & Input Data/PredatorMega8.csv')
breaks <- read.csv('R Code & Input Data/Breaks.csv')

## Convert break locations based on DistShift, not Dist
mega8 <- mega8.0
mega8[,c('preForest', 'Forest', 'Wall')] <- mega8.0[,c('preForest', 'Forest', 'Wall')] - (mega8.0$Dist - mega8.0$DistShift)

## Remove, rename, add, and round some columns
neo0 <- subset(mega8, select = c(Code, SiteBank, Group, Trophic, Region, Width, OrderClass, Geomorph, Banks, VegFld, preForest, Forest, Wall, DistShift, AbBiomTopCutNoWeirdHighDrop, Effect, logDist, logEffect))
neo1 <- neo0
colnames(neo1)[match(c('OrderClass', 'VegFld', 'preForest', 'DistShift', 'AbBiomTopCutNoWeirdHighDrop', 'logEffect', 'logDist'), colnames(neo0))] <- c('Order', 'Veg', 'PreF', 'Dist', 'AbBiom', 'logE', 'logD')
neo1[, c('Effect', 'logE', 'logD')] <- round(neo1[, c('Effect', 'logE', 'logD')], 4)
neo1$Effect2 <- ifelse(neo1$Effect == 1, .99, ifelse(neo1$Effect == 0, .01, neo1$Effect))
neo1$logE2 <- round(log(neo1$Effect2), 4)
neo1 <- neo1[with(neo1, order(SiteBank, Dist, -Effect)),] 
	rownames(neo1) <- 1:dim(neo1)[1]
	
	
##### Add variables relative to breaks #####

## Dummy variables pre/post breaks
neo1$AnteP <- ifelse(neo1$Dist < neo1$PreF | neo1$PreF %in% NA, 1, 0)
neo1$PostP <- ifelse(neo1$AnteP == 0, 1, 0)
neo1$AnteF <- ifelse(neo1$Dist < neo1$Forest | neo1$Forest %in% NA, 1, 0)
neo1$PostF <- ifelse(neo1$AnteF == 0, 1, 0)
neo1$AnteW <- ifelse(neo1$Dist < neo1$Wall | neo1$Wall %in% NA, 1, 0)
neo1$PostW <- ifelse(neo1$AnteW == 0, 1, 0)
neo1$AnteB <- ifelse(neo1$AnteW == 0 & neo1$AnteP == 0, 0, 1)
neo1$PostB <- ifelse(neo1$AnteB == 0, 1, 0)
neo1$Zone <- ifelse(neo1$PostP == 0 & neo1$PostW == 0, 'A', ifelse(neo1$PostP == 1 & neo1$PostW == 0, 'B', ifelse(neo1$PostP == 0 & neo1$PostW == 1, 'C', 'D'))) 

## Distance post-breaks
neo1$DistP <- ifelse(neo1$PostP == 1, neo1$Dist - neo1$PreF, 0)
neo1$logP <- round(ifelse(neo1$PostP == 1, ifelse(neo1$DistP == 0, -3, log(neo1$DistP)), 0), 4)
neo1$DistF <- ifelse(neo1$PostF == 1, neo1$Dist - neo1$Forest, 0)
neo1$logF <- round(ifelse(neo1$PostF == 1, ifelse(neo1$DistF == 0, -3, log(neo1$DistF)), 0), 4)
neo1$DistW <- ifelse(neo1$PostW == 1, neo1$Dist - neo1$Wall, 0)
neo1$logW <- round(ifelse(neo1$PostW == 1, ifelse(neo1$DistW == 0, -3, log(neo1$DistW)), 0), 4)
neo1$logA <- with(neo1, ifelse(logP != 0 & logW != 0, ifelse(logP > logW, logW, logP), ifelse(logP != 0 & logW == 0, logP, ifelse(logP == 0 & logW != 0, logW, 0))))

## Different potential forms of logD
neo1$logD2 <- ifelse(neo1$PostP == 1, 0, neo1$logD) ## 0 after pre-Forest
neo1$logD3 <- ifelse(neo1$PostW == 1, 0, neo1$logD) ## 0 after Wall
neo1$logD4 <- ifelse(neo1$PostP == 1 | neo1$PostW == 1, 0, neo1$logD) ## 0 after either break
neo1$logD5 <- ifelse(neo1$PostP == 1, neo1$logP, neo1$logD) ## resets after pre-Forest
neo1$logD6 <- ifelse(neo1$PostW == 1, neo1$logW, neo1$logD) ## resets after Wall
neo1$logD7 <- with(neo1, ifelse(PostP == 0 & PostW == 0, logD, ifelse(PostP == 1 & PostW ==0, logP, ifelse(PostP ==0 & PostW == 1, logW, ifelse(PreF > Wall, logP, logW))))) ## Resets after any break


##### Try fitting various functions to all data #####

## Start with very basic models
L0.0 <- lm(AbBiom ~ 1, data = neo1)
E0.0 <- lm(log(AbBiom) ~ 1, data = neo1)
	AIC0.1 <- c(AIC(L0.0), trueAIC(E0.0)[1])
	## Taking the log of the effect is unquestionably better (delta AIC ~ 7700)
	
## Add in a random intercept (SiteBank)
E0.1 <- lmer(log(AbBiom) ~ 1 + (1 | SiteBank), data = neo1, REML = FALSE)
	AIC0.2 <- AIC(E0.0, E0.1)
	## Adding the random effect unquestionably improves the fit (delta AIC ~ 100)
	
## Add distance
E0.2 <- lmer(log(AbBiom) ~ Dist + (1 | SiteBank), data = neo1, REML = FALSE)
P0.2 <- lmer(log(AbBiom) ~ logD + (1 | SiteBank), data = neo1, REML = FALSE)
	AIC0.3 <- AIC(E0.0, E0.1, E0.2, P0.2)
	## Adding distance may or may not help (delta AIC ~ 0-2)
	## Power better than exponential (delta AIC ~ 2)

## Add a random slope for distance
P0.3 <- lmer(log(AbBiom) ~ logD + (logD | SiteBank), data = neo1, REML = FALSE)
	AIC0.4 <- AIC(E0.1, P0.2, P0.3)
	## Adding a distance random slope unquestionably improves fit (delta AIC ~ 8)


##### Filter out bad sites #####

## Identify sites with positive initial slopes
SlopesP0.3 <- coef(P0.3)$SiteBank[order(coef(P0.3)$SiteBank$logD),]
MaybeBadP0.3 <- rownames(coef(P0.3)$SiteBank[coef(P0.3)$SiteBank$logD >= 0, ])
BadP0.3 <- MaybeBadP0.3[is.na(breaks[breaks$SiteBank %in% MaybeBadP0.3, 'predpreForest']) & is.na(breaks[breaks$SiteBank %in% MaybeBadP0.3, 'predWall'])]
	## There are 6 SiteBanks with positive slopes and no pre-Forest or Wall break. Stream signatures are going to be junk for these, so they must be removed.

## Make new dataset based only on good sites
neo2 <- neo1[!(neo1$SiteBank %in% BadP0.3), ]
neo2 <- droplevels(neo2)


##### Restart model fitting using new data #####

## Try linear, exponential, power models again
L0 <- lm(AbBiom ~ 1, data = neo2)
E0 <- lm(log(AbBiom) ~ 1, data = neo2)
E1 <- lmer(log(AbBiom) ~ 1 + (1 | SiteBank), data = neo2, REML = FALSE)
E2 <- lmer(log(AbBiom) ~ Dist + (1 | SiteBank), data = neo2, REML = FALSE)
P2 <- lmer(log(AbBiom) ~ logD + (1 | SiteBank), data = neo2, REML = FALSE)
	AIC1.0 <- c(AIC(L0), trueAIC(E0)[1])
	AIC1 <- AIC(E0, E1, E2, P2)
	## Same pattern. Distance helps a ton. Random intercept helps a ton. Power much better than exponential, which is much better than linear.

## Add a random slope for distance
P3 <- lmer(log(AbBiom) ~ logD + (logD | SiteBank), data = neo2, REML = FALSE)
	AIC2 <- AIC(E1, P2, P3)
	## Adding a distance random slope no longer improves model fit (now that the bad sites pulling the model are all gone!)

## Add intercepts for breaks
P4 <- lmer(log(AbBiom) ~ PostP + logD + (1 | SiteBank), data = neo2, REML = FALSE)
P5 <- lmer(log(AbBiom) ~ PostW + logD + (1 | SiteBank), data = neo2, REML = FALSE)
P6 <- lmer(log(AbBiom) ~ PostP + PostW + logD + (1 | SiteBank), data = neo2, REML = FALSE)
P7 <- lmer(log(AbBiom) ~ PostP * PostW + logD + (1 | SiteBank), data = neo2, REML = FALSE)
	AIC3 <- AIC(P2, P4, P5, P6, P7)
	## Adding a pre-Forest intercept helps, adding a Wall intercept does not.

## Try different slopes to account for breaks
P8 <- lmer(log(AbBiom) ~ PostP + logD2 + logP + (1 | SiteBank), data = neo2, REML = FALSE)
P9 <- lmer(log(AbBiom) ~ PostP + logD3 + logW + (1 | SiteBank), data = neo2, REML = FALSE)
P10 <- lmer(log(AbBiom) ~ PostP + logD4 + logP + (1 | SiteBank), data = neo2, REML = FALSE)
P11 <- lmer(log(AbBiom) ~ PostP + logD4 + logW + (1 | SiteBank), data = neo2, REML = FALSE)
P12 <- lmer(log(AbBiom) ~ PostP + logD4 + logP + logW + (1 | SiteBank), data = neo2, REML = FALSE)
P13 <- lmer(log(AbBiom) ~ PostP + logD4 + logA + (1 | SiteBank), data = neo2, REML = FALSE)
P14 <- lmer(log(AbBiom) ~ PostP + logD5 + (1 | SiteBank), data = neo2, REML = FALSE)
P15 <- lmer(log(AbBiom) ~ PostP + logD6 + (1 | SiteBank), data = neo2, REML = FALSE)
P16 <- lmer(log(AbBiom) ~ PostP + logD7 + (1 | SiteBank), data = neo2, REML = FALSE)
	AIC4 <- AIC(P4, P8, P9, P10, P11, P12, P13, P14, P15, P16)
	## Best option is to 0 out slope at breaks, with a new slope (logW) beyond Wall breaks only.

## Try random slopes one last time
P17 <- lmer(log(AbBiom) ~ PostP + logD4 + logW + (PostP | SiteBank), data = neo2, REML = FALSE)
P18 <- lmer(log(AbBiom) ~ PostP + logD4 + logW + (logD4 | SiteBank), data = neo2, REML = FALSE)
P19 <- lmer(log(AbBiom) ~ PostP + logD4 + logW + (logW | SiteBank), data = neo2, REML = FALSE)
	AIC5 <- AIC(P11, P17, P18, P19)
	## None of these help.
	
## Make a table leading logically from base to best fit
AICfit <- round(AIC(E2, P2, P3, P4, P5, P6, P10, P11, P12), 1)
	AICfit$Delta <- AICfit$AIC-min(AICfit$AIC)
	AICfitweight <- exp(-.5 * AICfit$Delta)
	AICfit$Weight <- round(AICfitweight / sum(AICfitweight), 4)
	AICfit$Model <- c('Base exponential', 'Base power', '+ Random slope', '+ Forest intercept', '+ Wall intercept', '+ Forest and Wall intercepts', '+ Forest intercept and slope', '+ Forest intercept and Wall slope', '+ Forest intercept and Forest and Wall slopes')
	## Best model includes distance random slope and intercept, a pre-Forest coefficient that comes into play at the forest break, and a slope that varies pre/post-Wall, with distance resetting at the Wall and zeroing out after forest breaks.


##### Get fitted data for each SiteBank #####	

## Get max and min AbBiom values for each SiteBank
maxs <- tapply(neo2$AbBiom, neo2$SiteBank, max)
mins <- tapply(neo2$AbBiom, neo2$SiteBank, min)

## Get fitted values for each distance, in 0.1 m increments, plus signatures
fits <- list()
sigs <- as.data.frame(names(maxs))
	colnames(sigs) <- 'SiteBank'
fitsbreaks <- breaks
fitsbreaks[is.na(fitsbreaks)] <- 1e100
disp <- (summary(P11)$sigma)^2
for(i in 1:length(maxs)){
	name <- names(maxs)[i]
	preF <- fitsbreaks[fitsbreaks$SiteBank == name, 'predpreForest']
	Wall <- fitsbreaks[fitsbreaks$SiteBank == name, 'predWall']
	LastBreak <- ifelse((preF != 1e100 & Wall != 1e100 & preF > Wall) | (preF != 1e100 & Wall == 1e100), preF, 0)
	coefs <- coef(P11)$SiteBank[rownames(coef(P11)$SiteBank) == name,]
	t1 <- as.data.frame(seq(0, 10000, .1))
		colnames(t1) <- 'Dist'
		t1$PostPreF <- ifelse(t1$Dist >= preF, 1, 0)
		t1$PostWall <- ifelse(t1$Dist >= Wall, 1, 0)
		t1$Wall <- rep(Wall, dim(t1)[1])
		t1$Last <- rep(LastBreak, dim(t1)[1])
		t1$Int <- as.numeric(rep(coefs[1], dim(t1)[1]))
		t1$FInt <- as.numeric(rep(coefs[2], dim(t1)[1]))
		t1$Slp <- as.numeric(rep(coefs[3], dim(t1)[1]))
		t1$WSlp <- as.numeric(rep(coefs[4], dim(t1)[1]))		
		t1$logD4 <- ifelse(t1$PostPreF == 1 | t1$PostW ==1, 0, ifelse(t1$Dist == 0, -3, log(t1$Dist)))
		t1$logW <- suppressWarnings(ifelse(t1$PostWall == 0, 0, ifelse(t1$Dist == t1$Wall, -3, log(t1$Dist - t1$Wall))))
		t1$logL <- suppressWarnings(ifelse((t1$Last == 0 | t1$PostPreF == 0), 0, ifelse(t1$Dist == t1$Last, -3, log(t1$Dist - t1$Last))))
		t1$logFitAB <- with(t1, ifelse(Last == 0, Int + FInt * PostPreF + Slp * logD4 + WSlp * logW, Int + FInt * PostPreF + Slp * logD4 + WSlp * logW + WSlp * logL))
		t1$FitAB <- round(exp(t1$logFitAB + 0.5 * disp), 2)
			## Note that need to take exp(log(y) + Dispersion/2) to get true estimates (not just exp(log(y))). Dispersion parameter = residual variance, or sigma^2. Can also get from VarCorr
		t1$FitEf <- round((t1$FitAB - min(t1$FitAB))/ (max(t1$FitAB) - min(t1$FitAB)), 4)		
	fits[[i]] <- t1[, c('Dist', 'FitAB', 'FitEf')]
	sigs$sig50[i] <- fits[[i]][min(which(fits[[i]]$FitEf <= 0.501)), 'Dist']
	sigs$sig10L[i] <- fits[[i]][min(which(fits[[i]]$FitEf <= 0.101)), 'Dist']
	sigs$sig10U[i] <- fits[[i]][max(which(fits[[i]]$FitEf >= 0.099)), 'Dist']
}
names(fits) <- names(maxs)

###STOPPED HERE.
	## Adjust individual site signatures to have Wall slopes post-Forest, assuming the Forest is the last break encountered. 
	## Something is wrong with how sig10L and sig10U are calculated (look at plot of ARZO1LB, for instance--it doesn't match with what sigs is saying). Suggests something is also wrong with how the fit is calculated.
	## Need to start adding factors.
	## Need to figure out what to do with getting signatures post-Forest break. Probably just fit a model that uses the Wall slope after those forest breaks. Or take numbers that would have been given for pre-break, and just add the break to them?
	## Need to incorporate code to get to mega8, above.
	## Need consolidate "From P drive" folder.
	

##### Plot curves for each SiteBank #####

## Try a test plot
bolzplot <- function(){
bolz <- neo2[neo2$SiteBank == 'BOLZ1RB',]
bolzcoefs <- coef(P11)$SiteBank[rownames(coef(P11)$SiteBank) == 'BOLZ1RB',]
bolzbreak <- breaks[breaks$SiteBank == 'BOLZ1RB', 'predpreForest']
bolzfit <- fits[['BOLZ1RB']]
	bolzfit$PostPreF <- ifelse(bolzfit$Dist >= 128, 1, 0)
	bolzfit$Dist5 <- ifelse(bolzfit$Dist < 128, bolzfit$Dist, bolzfit$Dist - 128)
	bolzfit$logD5 <- ifelse(bolzfit$Dist5 == 0, -3, log(bolzfit$Dist5))
	bolzfit$Int <- rep(as.numeric(bolzcoefs[1]), dim(bolzfit)[1])
	bolzfit$FInt <- rep(as.numeric(bolzcoefs[2]), dim(bolzfit)[1])
	bolzfit$Slp <- rep(as.numeric(bolzcoefs[3]), dim(bolzfit)[1])
	bolzfit$logFitAB2 <- with(bolzfit, bolzfit$Int + bolzfit$FInt * bolzfit$PostPreF + bolzfit$Slp * bolzfit$logD5)
	bolzfit$FitAB2 <- round(exp(bolzfit$logFitAB2 + 0.5 * disp), 2)
	bolzfit$FitEf2 <- round((bolzfit$FitAB2 - min(bolzfit$FitAB))/ (max(bolzfit$FitAB) - min(bolzfit$FitAB)), 4)		
plot(bolz$Effect ~ bolz$Dist, xlab = 'Distance', ylab = 'Effect', main = 'BOLZ1RB')
	with(bolzfit, points(Dist, FitEf2, type = 'l', col = 2))
	with(bolzfit, points(Dist, FitEf, type = 'l', col = 4))
	legend('top', legend = c('Real data', 'Fitted curve', 'Fitted curve with forced\n resetting after break'), lty = c(NA, 1, 1), pch = c(1, NA, NA), bty = 'n', col = c(1, 4, 2))
}
#bolzplot()

## Plot every SiteBank
plotting(Method = 'Predator', Data = neo2, DistVar = 'Dist', EffectVar = 'Effect', SiteVar = 'SiteBank', Ylab = '% Aquatic', NameAppend = 'wCurves', Type = 'pdf')

	
##### Figure out why the best models don't account for distance after the pre-Forest break #####

## Make unique dataframes for ante and post pre-Forest
neo2P <- neo2[neo2$PostP == 1, ]
	neo2P <- droplevels(neo2P)
neo2A <- neo2[neo2$AnteP == 1, ]
	neo2A <- droplevels(neo2A)

## Fit models only beyond pre-Forest break
L0P <- lm(AbBiom ~ 1, data = neo2P)
E0P <- lm(log(AbBiom) ~ 1, data = neo2P)
E1P <- lmer(log(AbBiom) ~ 1 + (1 | SiteBank), data = neo2P, REML = FALSE)
E2P <- lmer(log(AbBiom) ~ DistP + (1 | SiteBank), data = neo2P, REML = FALSE)
P2P <- lmer(log(AbBiom) ~ logP + (1 | SiteBank), data = neo2P, REML = FALSE)
P3P <- suppressWarnings(lmer(log(AbBiom) ~ logP + (logP | SiteBank), data = neo2P, REML = FALSE))
	## Does not converge.
P4P <- suppressWarnings(lmer(log(AbBiom) ~ PostW + logP + (logP | SiteBank), data = neo2P, REML = FALSE))
	## Does not converge.
P5P <- lmer(log(AbBiom) ~ logP + logW + (logP | SiteBank), data = neo2P, REML = FALSE)
P6P <- lmer(log(AbBiom) ~ PostW + logP + logW + (logP | SiteBank), data = neo2P, REML = FALSE)
	AICP <- AIC(L0P, E0P, E1P, E2P, P2P, P3P, P4P, P5P, P6P)
	## Confirmation that distance after pre-Forest break does not matter (best after pre-Forest model does not account for distance)

## Get typical distances sampled beyond pre-Forest break
PFmax <- tapply(neo2P$DistP, neo2P$SiteBank, function(x) max(x, na.rm = TRUE))
PFn <- tapply(neo2P$DistP, neo2P$SiteBank, length)
preVpost <- round(as.data.frame(c(median(PFmax), mean(PFmax), median(PFn), mean(PFn), length(PFn))), 1)
	colnames(preVpost) <- 'PostForest'
	rownames(preVpost) <- c('Median max distance', 'Mean max distance', 'Median samples per site', 'Mean samples per site', 'Number of sites')
AFmax <- tapply(neo2A$Dist, neo2A$SiteBank, function(x) max(x, na.rm = TRUE))
AFn <- tapply(neo2A$Dist, neo2A$SiteBank, length)
preVpost$AnteForest <- round(c(median(AFmax), mean(AFmax), median(AFn), mean(AFn), length(AFn)), 1)
	## Roughly 2-3x more samples per site and 2x more sites pre-Forest break than post-Forest break. Median max distance sampled post-Forest break also only 16 m. 
	## Probably a major reason is just not enough data post-break to justify fitting the distance decay, and many sites don't go far enough to show a decay.
	## Additionally, however, slope, when modeled above, is very low (near 0). Suggests perhaps decay not as rapid? Can't get further at this in a robust fashion though.
	

##### Compute overall signatures #####

## Create and fill dataframe of fitted values
Overall <- as.data.frame(seq(0, 10000, .1))
	colnames(Overall) <- 'Dist'
	Overall$logD4 <- ifelse(Overall$Dist == 0, -3, log(Overall$Dist))
	Overall$logFitAB <- fixef(P11)[1] + fixef(P11)[3] * Overall$logD4
	Overall$FitAB <- round(exp(Overall$logFitAB + 0.5 * disp), 2)
	Overall$FitEf <- round((Overall$FitAB - min(Overall$FitAB))/ (max(Overall$FitAB) - min(Overall$FitAB)), 4)		
	Overall$FitEf2 <- round((Overall$FitAB - min(Overall$FitAB))/ (max(Overall$FitAB) - min(Overall$FitAB)), 4)
	
## Get confidence intervals
OverallCIs <- suppressMessages(confint(P11, level = 0.95))
Overall$logFitABLower <- OverallCIs[3, 1] + OverallCIs[5, 1] * Overall$logD4
	Overall$logFitABUpper <- OverallCIs[3, 2] + OverallCIs[5, 2] * Overall$logD4
Overall$FitABLower <- round(exp(Overall$logFitABLower + 0.5 * disp), 2)
	Overall$FitABUpper <- round(exp(Overall$logFitABUpper + 0.5 * disp), 2)
Overall$FitEfLower <- round((Overall$FitABLower - min(Overall$FitABLower))/ (max(Overall$FitABLower) - min(Overall$FitABLower)), 4)	
	Overall$FitEfUpper <- round((Overall$FitABUpper - min(Overall$FitABUpper))/ (max(Overall$FitABUpper) - min(Overall$FitABUpper)), 4)	

## Get signatures
OverallSigs <- as.data.frame(round(c(Overall[min(which(Overall$FitEf <= 0.501)), 'Dist'], Overall[max(which(Overall$FitEf >= 0.099)), 'Dist']), 1))
	colnames(OverallSigs) <- 'Estimate'
	rownames(OverallSigs) <- c('Sig50' , 'Sig10')
	OverallSigs$LowerCI <- round(c(Overall[min(which(Overall$FitEfLower <= 0.501)), 'Dist'], Overall[max(which(Overall$FitEfLower >= 0.099)), 'Dist']), 1)
	OverallSigs$UpperCI <- round(c(Overall[min(which(Overall$FitEfUpper <= 0.501)), 'Dist'], Overall[max(which(Overall$FitEfUpper >= 0.099)), 'Dist']), 1)
	## Comparable (same order of magnitude) to estimates from meta-analysis.
OverallSigs
