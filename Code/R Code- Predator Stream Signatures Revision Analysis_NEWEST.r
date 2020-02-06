##### New analyses based on Ecological Monographs Reviews #####


##### Set working directory, read in data #####

## Working directory
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis')

## Create function to check if a package is installed. Install it if not and load it regardless.
pkgload <- function(pkg){
	newpkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
	if(length(newpkg) > 0){
		install.packages(newpkg, repos = 'http://cran.cnr.berkeley.edu/')
	}
	loaded <- lapply(pkg, require, character.only = TRUE, quietly = TRUE)
}

## Run the function for each package desired
pkgload(c('betareg', 'quantreg', 'lme4', 'nlme'))

## Load functions
source('R Code & Input Data/AICforLogResponse.r')

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


##### Try fitting various functions to all data #####

## Start with very basic models
L0 <- lm(AbBiom ~ 1, data = neo1)
E0 <- lm(log(AbBiom) ~ 1, data = neo1)
	AIC0 <- c(AIC(L0), trueAIC(E0)[1])
	## Taking the log of the effect is unquestionably better (delta AIC ~ 7700)
	
## Add in a random intercept (SiteBank)
E1 <- lmer(log(AbBiom) ~ 1 + (1 | SiteBank), data = neo1, REML = FALSE)
	AIC1 <- AIC(E0, E1)
	## Adding the random effect unquestionably improves the fit (delta AIC ~ 100)
	
## Add distance
E2 <- lmer(log(AbBiom) ~ Dist + (1 | SiteBank), data = neo1, REML = FALSE)
P2 <- lmer(log(AbBiom) ~ logD + (1 | SiteBank), data = neo1, REML = FALSE)
	AIC2 <- AIC(E0, E1, E2, P2)
	## Adding distance may or may not help (delta AIC ~ 0-2)
	## Power better than exponential (delta AIC ~ 2)

## Add a random slope for distance
P3 <- lmer(log(AbBiom) ~ logD + (logD | SiteBank), data = neo1, REML = FALSE)
	AIC3 <- AIC(E1, P2, P3)
	## Adding a distance random slope unquestionably improves fit (delta AIC ~ 8)
	
## Add a pre-Forest dummy variable, in isolation
P4 <- lmer(log(AbBiom) ~ PostP + (1 | SiteBank), data = neo1, REML = FALSE)
P5 <- lmer(log(AbBiom) ~ 1 + PostP + (PostP | SiteBank), data = neo1, REML = FALSE)
	AIC4 <- AIC(E1, P2, P3, P4, P5)
	## PostP alone does not help (delta AIC ~ 0 over base)

## Add pre-Forest dummy variable to distance model
P6 <- lmer(log(AbBiom) ~ logD + PostP + (logD | SiteBank), data = neo1, REML = FALSE)
P7 <- lmer(log(AbBiom) ~ PostP + logD:PostP + (logD | SiteBank), data = neo1, REML = FALSE)
P8 <- lmer(log(AbBiom) ~ logD + PostP + logD:PostP + (logD | SiteBank), data = neo1, REML = FALSE)
	AIC5 <- AIC(P3, P4, P5, P6, P7, P8)
	## Adding a pre-Forest intercept, and accounting for distance, unquestionably improves the model (delta AIC ~ 5)

## Add a random slope for pre-Forest
P9 <- lmer(log(AbBiom) ~ logD + PostP + (logD + PostP | SiteBank), data = neo1, REML = FALSE)
	AIC6 <- AIC(P6, P9)
	## Adding a random slope for pre-Forest does not improve model fit (delta AIC ~ -2)

## Add distance past the pre-Forest
P10 <- lmer(log(AbBiom) ~ logD + PostP + logP + (logD | SiteBank), data = neo1, REML = FALSE)
P11 <- lmer(log(AbBiom) ~ logD + PostP:logP + (logD | SiteBank), data = neo1, REML = FALSE)
P12 <- lmer(log(AbBiom) ~ logD + PostP + PostP:logP + (logD | SiteBank), data = neo1, REML = FALSE)
#P13 <- lmer(log(AbBiom) ~ logD + PostP + logP + PostP:logP + (logD | SiteBank), data = neo1, REML = FALSE)
	## Note: Would seem appealing to run this model, but that creates rank deficiencies due to perfect correlation, and model is identical to P10.
P14 <- lmer(log(AbBiom) ~ logD + logP + (logD | SiteBank), data = neo1, REML = FALSE)
	AIC7 <- AIC(P6, P10, P11, P12, P14)
	## No iteration of adding distance past the pre-Forest improves model fit (delta AIC ~ -2-10)

## Fit a curve before pre-Forest, and another after pre-Forest
P15 <- lmer(log(AbBiom) ~ logD + PostP + (logD | SiteBank), data = neo1, REML = FALSE)
	## Note: Would seem appealing to add a "+ AnteP" term to this model, but that creates rank deficiencies due to perfect correlation with PostP, and model is identical to this one.
P16 <- lmer(log(AbBiom) ~ logD + PostP + logD:AnteP + (logD | SiteBank), data = neo1, REML = FALSE)
#P17 <- lmer(log(AbBiom) ~ logD + PostP + AnteP + logD:AnteP + (logD | SiteBank), data = neo1, REML = FALSE)
	## Note: Would seem appealing to run this model, but that creates rank deficiencies due to perfect correlation with AnteP and PostP, and model is identical to P16.
P18 <- lmer(log(AbBiom) ~ logD:AnteP + logP:PostP + (logD | SiteBank), data = neo1, REML = FALSE)
P19 <- lmer(log(AbBiom) ~ logD + logD:AnteP + logP:PostP + (logD | SiteBank), data = neo1, REML = FALSE)
P20 <- lmer(log(AbBiom) ~ logD + AnteP + logD:AnteP + logP:PostP + (logD | SiteBank), data = neo1, REML = FALSE)
#P21 <- lmer(log(AbBiom) ~ logD + PostP + AnteP + logD:AnteP + logP:PostP + (logD | SiteBank), data = neo1, REML = FALSE)
	## Note: Would seem appealing to run this model, but that creates rank deficiencies due to perfect correlation with AnteP and PostP, and model is identical to P20.
	AIC8 <- AIC(P6, P15, P16, P18, P19, P20)
	## No iteration of fitting two curves improves model fit (delta AIC ~ -2-4)

## Add Wall dummy variable to the distance + pre-Forest model
P22 <- lmer(log(AbBiom) ~ logD + PostP + PostW + (logD | SiteBank), data = neo1, REML = FALSE)
P23 <- lmer(log(AbBiom) ~ logD + PostP + PostW + PostW:logW + (logD | SiteBank), data = neo1, REML = FALSE)
P24 <- lmer(log(AbBiom) ~ logD + PostP + PostW:logW + (logD | SiteBank), data = neo1, REML = FALSE)
P25 <- lmer(log(AbBiom) ~ AnteW:logD + PostP + PostW:logW + (logD | SiteBank), data = neo1, REML = FALSE)
	AIC9 <- AIC(P6, P22, P23, P24, P25)
	## Splitting the curve into separate power curves before and after the wall works best

## Add Wall dummy variable to distance model (no pre-Forest yet)
P26 <- lmer(log(AbBiom) ~ logD + PostW + (logD | SiteBank), data = neo1, REML = FALSE)
P27 <- lmer(log(AbBiom) ~ logD + PostW + PostW:logW + (logD | SiteBank), data = neo1, REML = FALSE)
P28 <- lmer(log(AbBiom) ~ logD + PostW:logW + (logD | SiteBank), data = neo1, REML = FALSE)
P29 <- lmer(log(AbBiom) ~ AnteW:logD + PostW:logW + (logD | SiteBank), data = neo1, REML = FALSE)
	AIC10 <- AIC(P3, P6, P25, P26, P27, P28, P29)
	## More important to add pre-Forest than to add Wall first (delta AIC ~ -4 to -10)
	
## Make a table leading logically from base to best fit
AICfit <- round(AIC(E0, E1, E2, P2, P3, P6, P23, P25), 1)
	AICfit$Delta <- AICfit$AIC-min(AICfit$AIC)
	AICfitweight <- exp(-.5 * AICfit$Delta)
	AICfit$Weight <- round(AICfitweight / sum(AICfitweight), 4)
	AICfit$Model <- c('Base', '+ Random intercept', '+ Distance', '+ log(Distance)', '+ Random slope', '+ Forest coefficient', '+ Wall coefficient & curve', '+ before-Wall curve')
	## Best model includes distance random slope and intercept, a pre-Forest coefficient that comes into play after the forest break, and a slope that varies pre/post-Wall.


##### Try fitting various functions using pre-pre-Forest only data #####

## Create pre-pre-Forest only dataframe
neo1P <- subset(neo1, Dist < PreF | is.na(PreF))

## Start with very basic models
L0P <- lm(AbBiom ~ 1, data = neo1P)
E0P <- lm(log(AbBiom) ~ 1, data = neo1P)
	AIC0P <- c(AIC(L0P), trueAIC(E0P)[1])
	## Taking the log of the effect is unquestionably better (delta AIC ~ 6600)
	
## Add in a random intercept (SiteBank)
E1P <- lmer(log(AbBiom) ~ 1 + (1 | SiteBank), data = neo1P, REML = FALSE)
	AIC1P <- AIC(E0P, E1P)
	## Adding the random effect unquestionably improves the fit (delta AIC ~ 90)
	
## Add distance
E2P <- lmer(log(AbBiom) ~ Dist + (1 | SiteBank), data = neo1P, REML = FALSE)
P2P <- lmer(log(AbBiom) ~ logD + (1 | SiteBank), data = neo1P, REML = FALSE)
	AIC2P <- AIC(E0P, E1P, E2P, P2P)
	## Adding distance helps (delta AIC ~ 1-4)
	## Power better than exponential (delta AIC ~ 3)

## Add a random slope for distance
P3P <- lmer(log(AbBiom) ~ logD + (logD | SiteBank), data = neo1P, REML = FALSE)
	AIC3P <- AIC(E1P, P2P, P3P)
	## Adding a distance random intercept unquestionably improves fit (delta AIC ~ 14)
	## Confirmation of the above that accounting for distance in the before the pre-forest break is a good idea.

	
##### Try fitting various functions using post-pre-Forest only data #####

## Create post-pre-Forest only dataframe
neo1P2 <- subset(neo1, Dist >= PreF)

## Start with very basic models
L0P2 <- lm(AbBiom ~ 1, data = neo1P2)
E0P2 <- lm(log(AbBiom) ~ 1, data = neo1P2)
	AIC0P2 <- c(AIC(L0P2), trueAIC(E0P2)[1])
	## Taking the log of the effect is unquestionably better (delta AIC ~ 752)
	
## Add in a random intercept (SiteBank)
E1P2 <- lmer(log(AbBiom) ~ 1 + (1 | SiteBank), data = neo1P2, REML = FALSE)
	AIC1P2 <- AIC(E0P2, E1P2)
	## Adding the random effect improves the fit (delta AIC ~ 4)
	
## Add distance
E2P2 <- lmer(log(AbBiom) ~ DistP + (1 | SiteBank), data = neo1P2, REML = FALSE)
P2P2 <- lmer(log(AbBiom) ~ logP + (1 | SiteBank), data = neo1P2, REML = FALSE)
	AIC2P2 <- AIC(E0P2, E1P2, E2P2, P2P2)
	## Adding distance does not improve fit (delta AIC ~ -1-2)
	## Power equivalent to exponential (delta AIC ~ 0.5)

## Add a random slope for distance
#P3P2 <- lmer(log(AbBiom) ~ logP + (logP | SiteBank), data = neo1P2, REML = FALSE)
	## Fails to converge
	## Confirmation that distance post-pre-Forest break does not matter.
	## Most likely just not enough data (not enough distances sampled) past those breaks to effectively model.
	## Model P6 (and subsequently P23, and finally P25) from above therefore truly are best.


##### Fit lines to each SiteBank #####	

## Get max and min AbBiom values for each SiteBank
maxs <- tapply(neo1$AbBiom, neo1$SiteBank, max)
mins <- tapply(neo1$AbBiom, neo1$SiteBank, min)

## Get fitted values for each distance, in 0.1 m increments
fits <- list()
sigs <- as.data.frame(names(maxs))
	colnames(sigs) <- 'SiteBank'
	sigs$sig10 <- sigs$sig50 <- NA
sigs2 <- sigs
fitsbreaks <- breaks
fitsbreaks[is.na(fitsbreaks)] <- 1e100
disp <- (summary(P25)$sigma)^2
for(i in 1:length(maxs)){
	name <- names(maxs)[i]
	preF <- fitsbreaks[fitsbreaks$SiteBank == name, 'predpreForest']
	Wall <- fitsbreaks[fitsbreaks$SiteBank == name, 'predWall']
	coefs <- coef(P25)$SiteBank[rownames(coef(P25)$SiteBank) == name,]
	t1 <- as.data.frame(seq(0, 10000, .1))
		colnames(t1) <- 'Dist'
		t1$PostpreF <- ifelse(t1$Dist >= preF, 1, 0)
		t1$AnteWall <- ifelse(t1$Dist < Wall, 1, 0)
		t1$PostWall <- ifelse(t1$Dist >= Wall, 1, 0)		
		t1$Slp <- as.numeric(rep(coefs[1], dim(t1)[1]))
		t1$Int <- as.numeric(rep(coefs[2], dim(t1)[1]))
		t1$FInt <- as.numeric(rep(coefs[3], dim(t1)[1]))
		t1$ASlp <- as.numeric(rep(coefs[4], dim(t1)[1]))
		t1$PSlp <- as.numeric(rep(coefs[5], dim(t1)[1]))
		t1$logD <- ifelse(t1$Dist == 0, -3.5, log(t1$Dist))
		t1$logW <- log(t1$Dist - Wall) * t1$PostWall
			t1$logW[is.na(t1$logW)] <- 0
		t1$logFitAB <- with(t1, Int + Slp * logD + PostpreF * FInt + AnteWall * logD * ASlp + PostWall * logW * PSlp)
		t1$FitAB <- round(exp(t1$logFitAB + 0.5 * disp), 2)
			## Note that need to take exp(log(y) + Dispersion/2) to get true estimates (not just exp(log(y))). Dispersion parameter = residual variance, or sigma^2. Can also get from VarCorr
		t1$FitEf <- round((t1$FitAB - mins[i])/ (maxs[i] - mins[i]), 4)
		t1$FitEf2 <- round((t1$FitAB - min(t1$FitAB))/ (max(t1$FitAB) - min(t1$FitAB)), 4)		
	fits[[i]] <- t1[, c('Dist', 'FitAB', 'FitEf', 'FitEf2')]
	sigs$sig50[i] <- fits[[i]][min(which(fits[[i]]$FitEf<=.501)), 'Dist']
	sigs$sig10[i] <- fits[[i]][max(which(fits[[i]]$FitEf>=.099)), 'Dist']
	sigs2$sig50[i] <- fits[[i]][min(which(fits[[i]]$FitEf2<=.501)), 'Dist']
	sigs2$sig10[i] <- fits[[i]][max(which(fits[[i]]$FitEf2>=.099)), 'Dist']
}
names(fits) <- names(maxs)

## Try a test plot
bolz <- neo1[neo1$SiteBank == 'BOLZ1RB',]
plot(bolz$Effect ~ bolz$Dist)
	with(fits[['BOLZ1RB']], points(Dist, FitEf, type = 'l', col = 2))
	with(fits[['BOLZ1RB']], points(Dist, FitEf2, type = 'l', col = 4))
	## Ef2 seems more reasonable (Effect size based on top and bottom of fitted, rather than real, values)

###STOPPED HERE.
	## Need to think critically about how to model slope past a break, and how the random slope factors into these models.
	## Need to incorporate code to get to mega8, above.
	## Need consolidate "From P drive" folder.
	
	
	
##### OLD BEYOND THIS POINT #####

## Basic Linear regressions
L0 <- lmer(AbBiom ~ Dist + (1|SiteBank), data = neo1, REML = FALSE)
E0 <- lmer(log(AbBiom) ~ Dist + (1|SiteBank), data = neo1, REML = FALSE)
P0 <- lmer(log(AbBiom) ~ logD + (1|SiteBank), data = neo1, REML = FALSE)

## Regressions with Charles' equation 1
L1 <- lmer(AbBiom ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP + (1|SiteBank), data = neo1, REML = FALSE)
E1 <- lmer(log(AbBiom) ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP + (1|SiteBank), data = neo1, REML = FALSE)
P1 <- lmer(log(AbBiom) ~ 0 + AnteP + PostP + AnteP:logD + PostP:logP + (1|SiteBank), data = neo1, REML = FALSE)

## Regressions with Charles' equation 2
L2 <- lmer(AbBiom ~ PostP + Dist + PostP:DistP + (1|SiteBank), data = neo1, REML = FALSE)
E2 <- lmer(log(AbBiom) ~ PostP + Dist + PostP:DistP + (1|SiteBank), data = neo1, REML = FALSE)
P2 <- lmer(log(AbBiom) ~ PostP + logD + PostP:logP + (1|SiteBank), data = neo1, REML = FALSE)

## Compare AICs
AICs <- data.frame(Model = c('Linear', 'Exponential', 'Power'))
	AICs$Eq2 <- AICs$Eq1 <- AICs$Basic <- NA
	AICs[ , 2:4] <- round(sapply(c(L0, E0, P0, L1, E1, P1, L2, E2, P2), AIC))
	
## Try power models with different forms
P4 <- lmer(log(AbBiom) ~ PostP + logD + (1|SiteBank), data = neo1, REML = FALSE)
P5 <- lmer(log(AbBiom) ~ logD + PostP:logP + (1|SiteBank), data = neo1, REML = FALSE)
P6 <- lmer(log(AbBiom) ~ 0 + AnteP + PostP + logD + logP + (1|SiteBank), data = neo1, REML = FALSE)
P7 <- lmer(log(AbBiom) ~ 0 + AnteP:logD + PostP:logP + (1|SiteBank), data = neo1, REML = FALSE)
P8 <- lmer(log(AbBiom) ~ AnteP:logD + PostP:logP + (1|SiteBank), data = neo1, REML = FALSE)
P9 <- lmer(log(AbBiom) ~ logD + (1 + logD|SiteBank), data = neo1, REML = FALSE)
P10 <- lmer(log(AbBiom) ~ logD + (1|SiteBank), data = neo1, REML = FALSE)
P11 <- lmer(log(AbBiom) ~ 1 + (1|SiteBank), data = neo1, REML = FALSE)
P12 <- lm(log(AbBiom) ~ logD, data = neo1)
P13 <- lmer(log(AbBiom) ~ 1 + (1 + logD|SiteBank), data = neo1, REML = FALSE)
P14 <- lmer(log(AbBiom) ~ PostP + (1 + logD|SiteBank), data = neo1, REML = FALSE)
P15 <- lmer(log(AbBiom) ~ PostP + (1 + logD + PostP|SiteBank), data = neo1, REML = FALSE)
P16 <- lmer(log(AbBiom) ~ PostP:logD + (1 + logD|SiteBank), data = neo1, REML = FALSE)
P17 <- lmer(log(AbBiom) ~ PostP*logD + (1 + logD|SiteBank), data = neo1, REML = FALSE)
P18 <- lmer(log(AbBiom) ~ PostP + PostP:logD + (1 + logD|SiteBank), data = neo1, REML = FALSE)
##Include logP somehow?
	AICs2 <- data.frame(sapply(list(P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19), AIC))
		colnames(AICs2) <- 'AIC'
		rownames(AICs2) <- 4:19
##### Subset data for a characteristic site for ease of viewing #####

## Subset only BOLZ1RB site
bolz<-neo1[neo1$SiteBank=='BOLZ1RB',]

## Subset data to before forest break
bolz2 <- bolz[bolz$Dist < bolz$PreF, ]


##### Try fitting various functions for BOLZ1RB pre-forest data #####

## Basic Linear regressions
L0.2 <- glm(AbBiom ~ Dist, data = bolz2)
E0.2 <- glm(log(AbBiom) ~ Dist, data = bolz2)
P0.2 <- glm(log(AbBiom) ~ logD, data = bolz2)

## Compare AICs
AICs <- data.frame(Model = c('Linear0', 'Exponential0', 'Power0', 'Linear1', 'Exponential1', 'Power1', 'Linear2', 'Exponential2', 'Power2'))
	AICs$AICpreForest <- round(c(AIC(L0.2), trueAIC(E0.2)[1], trueAIC(P0.2)[2]))
	## Power and exponential basically equivalent for BOLZ1RB pre-Forest data.
	## Both much better than linear though, assuming conversion is trusted.


##### Plot BOLZ1RB pre-Forest data	

## Set some plotting parameters
par(mar = c(4, 4, .5, .5))
cols <- c(2, 3, 4)

## Plot pre-Forest data, basic lognormal or linear models on AbBiom
plot(bolz$Dist, bolz$AbBiom, ylab = 'Abundance * biomass (#*mg/trap)', xlab = 'Distance from river (m)')
	curve(coef(L0.2)[1] + coef(L0.2)[2] * x, add = TRUE, col = cols[1], lty = 1)
	curve(exp(coef(E0.2)[1] + coef(E0.2)[2] * x), add = TRUE, col = cols[2], lty = 1)
	curve(exp(coef(E0.2)[1] + coef(E0.2)[2] * log(x+.05)), add = TRUE, col = cols[3], lty = 1)
	legend(30, 750, legend = paste(AICs$Model[1:3], ': AIC = ', AICs$AICpreForest[1:3], sep = ''), title = 'Model 0', col = cols, lty = 1, bty = 'n')
### THESE ARE STILL ATROCIOUS. ###
	
	
##### Try fitting various functions for BOLZ data #####

## Basic Linear regressions
L0 <- glm(AbBiom ~ Dist, data = bolz)
E0 <- glm(log(AbBiom) ~ Dist, data = bolz)
P0 <- glm(log(AbBiom) ~ logD, data = bolz)

## Regressions with Charles' equation 1
L1 <- glm(AbBiom ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = bolz)
E1 <- glm(log(AbBiom) ~ 0 + AnteP + PostP + AnteP:Dist + PostP:DistP, data = bolz)
P1 <- glm(log(AbBiom) ~ 0 + AnteP + PostP + AnteP:logD + PostP:logP, data = bolz)

## Regressions with Charles' equation 2
L2 <- glm(AbBiom ~ PostP + Dist + PostP:DistP, data = bolz)
E2 <- glm(log(AbBiom) ~ PostP + Dist + PostP:DistP, data = bolz)
P2 <- glm(log(AbBiom) ~ PostP + logD + PostP:logP, data = bolz)

## Compare AICs
AICs <- data.frame(Model = c('Linear0', 'Exponential0', 'Power0', 'Linear1', 'Exponential1', 'Power1', 'Linear2', 'Exponential2', 'Power2'))
	AICs$AICoverall <- round(c(AIC(L0), trueAIC(E0)[1], trueAIC(P0)[2], AIC(L1), trueAIC(E1)[1], trueAIC(P1)[2], AIC(L1), trueAIC(E1)[1], trueAIC(P1)[2]))
	## Using Charles' equations 1 and 2 do not improve AIC, at least for this 1 site.	
