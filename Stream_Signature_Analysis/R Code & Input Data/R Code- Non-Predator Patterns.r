##### Set working directory, read in raw data #####
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis')
raw1<-read.csv('R Code & Input Data/CSVBugs.csv',header=T,na.strings=c('n/a',' ','#N/A'))
env00<-read.csv('R Code & Input Data/All Sites Data.csv',header=T,na.strings=c('n/a','N/A','#N/A'))
breaks<-read.csv('R Code & Input Data/Breaks.csv')
	
##### Load some useful functions and libraries #####
library(MuMIn)
library(nlme)
source('C:/Users/jmuehlbauer/Documents/Misc/Trainings & Info/Stats & Data Analysis/R/AICforLogResponse.r')
source('R Code & Input Data/StreamSignatureFunction.r')
source('R Code & Input Data/PlottingFunction.r')

##### Create bugs dataframe #####

## Remove ground sampling sites without specific locations
bugs1<-raw1[raw1$Bank!='ALL'&raw1$Group!='--',]
	bugs1<-droplevels(bugs1)

## Convert Aq distances to NAs
bugs1$Dist<-suppressWarnings(as.numeric(as.character(bugs1$Dist)))

## Add trophic level
bugs2<-bugs1
	bugs2$Trophic<-rep(4,dim(bugs2)[1])
bugs2$Trophic<-ifelse(substr(bugs2$Group,1,2)=='Pr',2,ifelse(substr(bugs2$Group,1,2)=='He',1,ifelse(substr(bugs2$Group,3,6)=='Detr',0,ifelse(bugs2$Group=='AqAlga',0,1.5))))

			
##### Clean up env datasheet #####
env00$Order<-as.factor(env00$Order)
env.sites2<-env00[env00$Site=='BALL1'|env00$Site=='COWE1'|env00$Site=='LTEN1',]
	env.sites2$Site<-c('BALL2','BALL2','COWE2','COWE2','LTEN2')
env0<-rbind(env00,env.sites2)


##### Combine env data to abundance table #####
bugs3<-bugs2
m1<-match(substr(bugs3[,2],1,5),env0[,1])
	bugs3$Region<-env0$Region[m1]
	bugs3$Width<-env0$Channel.width[m1]
	bugs3$WidthClass<-env0$Width.class[m1]		
	bugs3$Order<-env0$Order[m1]
	bugs3$OrderClass<-env0$Order.class[m1]
	bugs3$Geomorph<-env0$Geomorphology[m1]
m2<-match(paste(substr(bugs3$Site,1,5),substr(bugs3$Bank,1,1)),paste(env0$Site,substr(env0$Bank,1,1)))	
	bugs3$Banks<-env0$Bank.type[m2]
	bugs3$VegFld<-env0$Floodplain.vegetation[m2]
	bugs3$VegUp<-env0$Upland.vegetation[m2]
	bugs3$VegShift<-env0$Vegetation.shift[m2]


##### Create datasheet for plotting and analysis #####
	
## Get rid of NAs
mega1<-droplevels(bugs3[!is.na(bugs3$Dist)&!is.na(bugs3$Bank),])

## Get rid of all sites with islands
mega2<-mega1[!substr(mega1$Site,1,5)=='CORN1'&!substr(mega1$Site,1,5)=='DANU1'&!substr(mega1$Site,1,5)=='ELBE1'&!substr(mega1$Site,1,5)=='FELL1'&!substr(mega1$Site,1,5)=='ISLA1'&!substr(mega1$Site,1,5)=='RESI1'&!substr(mega1$Site,1,5)=='RESI5',]

## Add true banks of island sites back in
mega2.1<-rbind(mega2,mega1[
	substr(mega1$Code,1,10)=='CORN1.1-RB'|
		substr(mega1$Code,1,10)=='CORN1.6-LB'|
	substr(mega1$Code,1,10)=='DANU1.2-RB'|
	substr(mega1$Code,1,10)=='ELBE1.2-LB'|
	substr(mega1$Code,1,10)=='FELL1.1-RB'|
	substr(mega1$Code,1,10)=='ISLA1.1-RB'|
		substr(mega1$Code,1,10)=='ISLA1.7-LB'|
	substr(mega1$Code,1,10)=='RESI1.1-RB'|
	substr(mega1$Code,1,10)=='RESI1.2-LB'|
	substr(mega1$Code,1,10)=='RESI5.1-RB'|
		substr(mega1$Code,1,10)=='RESI5.2-LB',])

## Rearrange alphabetically for ease
mega2.1<-droplevels(mega2.1[order(mega2.1$Code),])

## Get rid of predator and aquatic groups and night/malaise/ground sampling
mega3<-droplevels(mega2.1[mega2.1$Trophic!=2&substr(mega2.1$Group,1,2)!='Aq'&substr(mega2.1$Group,3,4)!='Aq'&mega2.1$Method!='Night Sampling'&mega2.1$Method!='Malaise'&mega2.1$Method!='Ground',])

## Add sweeps to each pitfall/bowl replicate
mega3Swp<-droplevels(mega3[mega3$Method=='Sweep'&mega3$Abund>0,])
	sumsAbund<-tapply(mega3Swp$Abund,mega3Swp$Code,sum)
	sumsBiomass<-tapply(mega3Swp$Biomass,mega3Swp$Code,sum)
mega4PitBowl<-mega3[mega3$Method=='Pitfall'|mega3$Method=='Bowl',]
mega4Swp<-mega3Swp[mega3Swp$Code%in%mega4PitBowl$Code==FALSE,]
mega4<-rbind(mega4PitBowl,mega4Swp)
	mega4$Abund<-ifelse(is.na(match(mega4$Code,names(sumsAbund))),mega4$Abund,mega4$Abund+sumsAbund[match(mega4$Code,names(sumsAbund))])
	mega4$Biomass<-ifelse(is.na(match(mega4$Code,names(sumsBiomass))),mega4$Biomass,mega4$Biomass+sumsBiomass[match(mega4$Code,names(sumsBiomass))])
	mega4<-droplevels(subset(mega4,select=-Method))
		rownames(mega4)<-c(1:dim(mega4)[1])
		
## Create some data columns for future use
mega4$shortSite<-as.factor(substr(mega4$Site,1,5))
mega4$SiteBank<-as.factor(paste(mega4$shortSite,mega4$Bank,sep=''))
mega4$AbBiom<-ifelse(mega4$Site=='BALL2'|mega4$Site=='COWE2'|mega4$Site=='LTEN2',mega4$Abund,mega4$Abund*mega4$Biomass)

## Remove NA rows in biomass and HawR, MudC, LTEN2, COWE2, BALL2
mega5<-droplevels(mega4[!is.na(mega4$Biomass)&!is.na(mega4$AbBiom)&mega4$Site!='HAWR1'&mega4$Site!='MUDC1'&mega4$Site!='BALL2'&mega4$Site!='COWE2'&mega4$Site!='LTEN2',])


##### Plot abundance, biomass, and abundance*biomass #####

## Plot various figures by SiteBank
#plotting('Non-Predator',mega5,'Dist','Abund','SiteBank','Abundance',Type='pdf')	
#plotting('Non-Predator',mega5,'Dist','Biomass','SiteBank','Biomass',Type='pdf')
#plotting('Non-Predator',mega5,'Dist','AbBiom','SiteBank','Abundance x Biomass',Type='pdf')


##### Remove outlier data #####
	
	## Based on Non-PredatorAbBiomSiteBank.pdf
	
## Remove extremely high values. TopCuts saved in Breaks.csv file
mega5$AbBiomTopCut<-mega5$AbBiom
for(i in 1:dim(mega5)[1]){
	bksub<-breaks[match(mega5$SiteBank[i],breaks$SiteBank),'nonpredTopCut']
	if(mega5$AbBiom[i]>bksub&is.na(bksub)==FALSE){mega5$AbBiomTopCut[i]<-NA}
}

## Plot AbBiom figures with these TopCut reductions
plotting('Non-Predator',mega5,'Dist','AbBiomTopCut','SiteBank','Abundance x Biomass',Type='pdf')


##### Inspect sites with forest boundaries and steep banks/levees #####

## Pinpoint locations of these breaks	
	## Breaks saved in Breaks.csv file.
mega5$preForest<-breaks[match(mega5$SiteBank,breaks$SiteBank),'predpreForest']
mega5$Forest<-breaks[match(mega5$SiteBank,breaks$SiteBank),'predForest']
mega5$Wall<-breaks[match(mega5$SiteBank,breaks$SiteBank),'predWall']

## Shift datapoints for sites where first few meters get flooded
mega5$DistShift<-mega5$Dist+breaks[match(mega5$SiteBank,breaks$SiteBank),'predDistShift']
mega6<-droplevels(mega5[!is.na(mega5$AbBiomTopCut)&mega5$DistShift>=0,])
#plotting('Non-Predator',mega6,'DistShift','AbBiomTopCut','SiteBank','Abundance x Biomass','DistShift',Type='pdf')


##### Convert Abundance x Biomass to effect size #####

## Percent of SiteBank max value
mega8<-mega6
maxs<-tapply(mega8$AbBiomTopCut,mega8$SiteBank,max)
mega8$Effect<-as.numeric(mega8$AbBiomTopCut/maxs[mega8$SiteBank])
plotting('Non-Predator',mega8,'DistShift','Effect','SiteBank','Relative Abundance x Biomass','DistShift')		

## Create some data columns for future use
mega8$logEffect<-log(mega8$Effect+.05)
mega8$logDist<-log(mega8$DistShift+.05)		

mega8preForest<-droplevels(rbind(subset(mega8,Dist<preForest),subset(mega8,is.na(preForest))))
	mega8preForest<-mega8preForest[order(mega8preForest$Code),]
	mega8preForest<-subset(mega8preForest,SiteBank!='DANU1RB')
mega8Forest<-droplevels(subset(mega8,Dist>=Forest))
mega8Wall<-droplevels(subset(mega8,Dist>=Wall))
mega8preWall<-droplevels(rbind(subset(mega8,Dist<Wall),subset(mega8,is.na(Wall))))
funkies<-c('DANU1RB','DFOR1RB','DRAV1LB','DSRB1LB','RESI5LB')	
mega8Reasonable<-rbind(
	mega8preForest[mega8preForest$SiteBank%in%funkies==FALSE,],
	mega8preWall[mega8preWall$SiteBank=='DFOR1RB',],
	mega8Wall[mega8Wall$SiteBank%in%'RESI5LB',],
	mega8Forest[mega8Forest$SiteBank%in%c('DANU1RB','DRAV1LB'),],
	mega8Forest[mega8Forest$SiteBank=='DSRB1LB'&mega8Forest$Dist<breaks[breaks$SiteBank=='DSRB1LB','predWall'],])	
	
	
##### Find SiteBank data subsets with positive (non-decaying) curves #####

## Make Regs.All function
Regs.All<-function(Data){
	Regs.All.data<-list()
	nas<-as.numeric(rep(NA,length(unique(Data$SiteBank))))
	for(i in 1:length(unique(Data$SiteBank))){
		subs<-subset(Data,SiteBank==unique(Data$SiteBank)[i])
		if(dim(subs)[1]>2){
			subs.Lin<-lm((Effect+.05)~DistShift,data=subs)
				results.Lin<-round(c(AIC(subs.Lin)[1],summary(subs.Lin)$adj.r.squared,summary(subs.Lin)$coefficients[2,3],summary(subs.Lin)$coefficients[2,4]),4)
			subs.Exp<-lm(logEffect~DistShift,data=subs)
				results.Exp<-round(c(trueAIC(subs.Exp)[1],summary(subs.Exp)$adj.r.squared,summary(subs.Exp)$coefficients[2,3],summary(subs.Exp)$coefficients[2,4]),4)
			subs.Pow<-lm(logEffect~logDist,data=subs)
				results.Pow<-round(c(trueAIC(subs.Pow)[1],summary(subs.Pow)$adj.r.squared,summary(subs.Pow)$coefficients[2,3],summary(subs.Pow)$coefficients[2,4]),4)
			Regs.All.data[[i]]<-as.data.frame(rbind(results.Lin,results.Exp,results.Pow))
				colnames(Regs.All.data[[i]])<-c('AIC','AdjR^2','Dist.t.val','Dist.p.val')
				rownames(Regs.All.data[[i]])<-c('Linear','Exp.', 'Power')
			}
		else{Regs.All.data[[i]]<-NA
			nas[i]<-i
			}
		}	
		names(Regs.All.data)<-unique(Data$SiteBank)
		Regs.All.data<-Regs.All.data[is.na(nas)]
		Regs.All.data
	}
	
## Using Reasonable data
	## Pos if curve of best fit (either linear, exponential, or power) show a positive curve based on AIC
Regs.All.Reasonable<-Regs.All(mega8Reasonable)
	Regs.Neg.Reasonable<-Regs.All.Reasonable[sapply(Regs.All.Reasonable, function(x) x[which(x$AIC==min(x$AIC,na.rm=T))[1],'Dist.t.val']<0)]
	Regs.Pos.Reasonable<-Regs.All.Reasonable[sapply(Regs.All.Reasonable, function(x) x[which(x$AIC==min(x$AIC,na.rm=T))[1],'Dist.t.val']>=0)]
	## Only 12/46 (26.09%) of SiteBanks have a decay with Distance; 34/46 SiteBanks (73.91%) are positively related to Distance.
	Regs.sigp.Reasonable<-Regs.All.Reasonable[sapply(Regs.All.Reasonable, function(x) x[which(x$AIC==min(x$AIC,na.rm=T))[1],'Dist.p.val']<.05)]
	## Only 12/46 (26.09%) of SiteBanks regression lines have significant slopes
	## Of these, only one is negative with Distance.
	
	
##### Try mixed effects models on data #####

## Use Site+Bank as the structural variable
lmeSB.lin<-lme((Effect+.05)~DistShift,random=~1|SiteBank,data=mega8,method='ML')
lmeSB.exp<-lme(logEffect~DistShift,random=~1|SiteBank,data=mega8,method='ML')
lmeSB.pow<-lme(logEffect~logDist,random=~1|SiteBank,data=mega8,method='ML')
rbind(AIC(lmeSB.lin),trueAIC(lmeSB.exp),trueAIC(lmeSB.pow));AIC(lmeSB.exp,lmeSB.pow)
	## All models are pretty poor, looking at values of slopes and intercepts (essentially flat lines). But all are definitely positive and non-decaying.
	## In light of this, not appropriate to look at observed vs. expected differences at breaks (don't trust the models).