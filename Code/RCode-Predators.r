##### Set working directory, read in raw data #####
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis')
raw1<-read.csv('R Code & Input Data/CSVBugs.csv',header=T,na.strings=c('n/a',' ','#N/A'))
env00<-read.csv('R Code & Input Data/All Sites Data.csv',header=T,na.strings=c('n/a','N/A','#N/A'))
islandsXS<-read.csv('R Code & Input Data/Islands Cross Section Data.csv',header=T,na.strings=c(' '))
isos.mega4<-read.csv('R Code & Input Data/IsotopeMega4forGoodBadFig.csv',row.names=1)
breaks<-read.csv('R Code & Input Data/Breaks.csv')
	
##### Load some useful functions and libraries #####
library(MuMIn)
library(nlme)
library(png)
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

## Get rid of non-predator groups and night/malaise/ground sampling
mega3<-droplevels(mega2.1[mega2.1$Trophic==2&mega2.1$Method!='Night Sampling'&mega2.1$Method!='Malaise'&mega2.1$Method!='Ground',])

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

## Plot various figures by shortSite and SiteBank
#plotting('Predator',mega5,'Dist','Abund','shortSite','Abundance',Type='pdf')	
#plotting('Predator',mega5,'Dist','Abund','SiteBank','Abundance',Type='pdf')	
#plotting('Predator',mega5,'Dist','Biomass','shortSite','Biomass',Type='pdf')	
#plotting('Predator',mega5,'Dist','Biomass','SiteBank','Biomass',Type='pdf')
#plotting('Predator',mega5,'Dist','AbBiom','shortSite','Abundance x Biomass',Type='pdf')	
#plotting('Predator',mega5,'Dist','AbBiom','SiteBank','Abundance x Biomass',Type='pdf')


##### Remove outlier data #####
	
	## Based on PredatorAbBiomSiteBank.pdf
	
## Remove extremely high values. TopCuts saved in Breaks.csv file
mega5$AbBiomTopCut<-mega5$AbBiom
for(i in 1:dim(mega5)[1]){
	bksub<-breaks[match(mega5$SiteBank[i],breaks$SiteBank),'predTopCut']
	if(mega5$AbBiom[i]>bksub&is.na(bksub)==FALSE){mega5$AbBiomTopCut[i]<-NA}
}

## Plot AbBiom figures with these TopCut reductions
#plotting('Predator',mega5,'Dist','AbBiomTopCut','shortSite','Abundance x Biomass',Type='pdf')
#plotting('Predator',mega5,'Dist','AbBiomTopCut','SiteBank','Abundance x Biomass',Type='pdf')
	
	
##### Inspect sites with forest boundaries and steep banks/levees #####

## Pinpoint locations of these breaks	
	## Breaks saved in Breaks.csv file.
mega5$preForest<-breaks[match(mega5$SiteBank,breaks$SiteBank),'predpreForest']
mega5$Forest<-breaks[match(mega5$SiteBank,breaks$SiteBank),'predForest']
mega5$Wall<-breaks[match(mega5$SiteBank,breaks$SiteBank),'predWall']

## Shift datapoints for sites where first few meters get flooded
mega5$DistShift<-mega5$Dist+breaks[match(mega5$SiteBank,breaks$SiteBank),'predDistShift']
mega6<-droplevels(mega5[!is.na(mega5$AbBiomTopCut)&mega5$DistShift>=0,])
#plotting('Predator',mega6,'DistShift','AbBiomTopCut','SiteBank','Abundance x Biomass','DistShift',Type='pdf')

## Remove other inexplicable outlier datapoints
mega6$AbBiomTopCutNoWeird<-
	ifelse(mega6$SiteBank=='ARZO1LB'&mega6$AbBiomTopCut>150&mega6$Dist==8,NA,
		ifelse(mega6$SiteBank=='ARZO1LB'&mega6$AbBiom<25&mega6$Dist<5,NA,
	ifelse(mega6$SiteBank=='ARZO1RB'&mega6$AbBiomTopCut>95&mega6$Dist<10,NA,
		ifelse(mega6$SiteBank=='ARZO1RB'&mega6$AbBiomTopCut>50&mega6$Dist==16,NA,
	ifelse(mega6$SiteBank=='BALL1LB'&mega6$AbBiomTopCut>25&mega6$Dist>38&mega6$Dist<60,NA,
	ifelse(mega6$SiteBank=='BECH1RB'&mega6$AbBiomTopCut>200&mega6$Dist>100,NA,
	ifelse(mega6$SiteBank=='BECH2RB'&mega6$AbBiomTopCut>50&mega6$Dist>100,NA,
	ifelse(mega6$SiteBank=='DRAV1LB'&mega6$AbBiomTopCut>50&mega6$Dist>40,NA,
	ifelse(mega6$SiteBank=='DSRB2LB'&mega6$Dist>250,NA,
		ifelse(mega6$SiteBank=='DSRB2LB'&mega6$AbBiomTopCut>50&mega6$Dist>50&mega6$Dist<150,NA,
	ifelse(mega6$SiteBank=='ELBE1LB'&mega6$AbBiomTopCut>30&mega6$Dist>50,NA,
		mega6$AbBiomTopCut)))))))))))
#plotting('Predator',mega6,'DistShift','AbBiomTopCutNoWeird','SiteBank','Abundance x Biomass','DistShift',Type='pdf')
		

##### Convert Abundance x Biomass to effect size #####

## Drop very high values at far distances to the max value near the bank
mega7<-droplevels(mega6[!is.na(mega6$AbBiomTopCutNoWeird),])
veryhighs<-c('ARZO1RB','BECH2RB','BOLZ1RB','DRAV1LB','FELL1RB')
maxnears<-vector()
for(i in 1:length(unique(mega7$SiteBank))){
	if(unique(mega7$SiteBank)[i]%in%veryhighs){
		maxnears[i]<-max(mega7[mega7$SiteBank==unique(mega7$SiteBank)[i]&mega7$Dist<=8,'AbBiomTopCutNoWeird'])}
	else{maxnears[i]<-max(mega7[mega7$SiteBank==unique(mega7$SiteBank)[i],'AbBiomTopCutNoWeird'])}
	}
	names(maxnears)<-unique(mega7$SiteBank)
mega7$AbBiomTopCutNoWeirdHighDrop<-ifelse(mega7$AbBiomTopCutNoWeird>maxnears[as.character(mega7$SiteBank)],maxnears[as.character(mega7$SiteBank)],mega7$AbBiomTopCutNoWeird)

## Percent of SiteBank max value
mega8<-mega7
maxs<-tapply(mega8$AbBiomTopCutNoWeirdHighDrop,mega8$SiteBank,max)
mega8$Effect<-as.numeric(mega8$AbBiomTopCutNoWeirdHighDrop/maxs[mega8$SiteBank])
#plotting('Predator',mega8,'DistShift','Effect','SiteBank','% Aquatic','DistShift')
		
## Create some data columns for future use
mega8$logEffect<-log(mega8$Effect)
mega8$logDist<-ifelse(mega8$DistShift == 0, -3, log(mega8$DistShift))
			## Note: -3 substitution for log(0) is ~ equal to log(0.05)
#write.csv(mega8,'R Code & Input Data/PredatorMega8.csv',row.names=FALSE)


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
	
	
##### Make figure to show effect of forest boundaries #####

## Use BOLZ1RB for break, BELG1RB for no-break
belg<-mega8[mega8$SiteBank=='BELG1LB',]
bolz<-mega8[mega8$SiteBank=='BOLZ1RB',]

## Read in Bolzano picture
bolzimg<-suppressWarnings(readPNG('Figures & Tables/Raw Figures/BOLZ1RB.png'))

## Plot figure panels
f1<-function(Type='pdf'){
if(Type=='png'){png('Figures & Tables/Predator Forest Break.png',width=6.8,height=9,res=300,bg=0,units='in')}
else{pdf('Figures & Tables/Predator Forest Break.pdf',width=6.8,height=9)}
oldpar<-par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(3.3,3.8,0.5,0.5),cex=1)
plot(belg$DistShift,belg$Effect,bty='l',xaxt='n',yaxt='n',ylab='',xlab='',type='n')
	axis(1,at=seq(0,80,20),lwd=0,lwd.tick=1)
	axis(2,las=2,lwd=0,lwd.tick=1)
	points(belg$DistShift,belg$Effect,cex=.75,pch=as.numeric(belg$Group))
	legend(70,1.15,legend='A',bty='n',cex=2)
	legend(25,1.1,legend='BELG1 left bank',bty='n')
	legend(50,1.1,legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3),bty='n',pt.cex=.75)
plot(bolz$DistShift,bolz$Effect,bty='l',xaxt='n',yaxt='n',ylab='',xlab='',type='n')
	axis(1,at=seq(0,160,40),lwd=0,lwd.tick=1)
	axis(2,las=2,lwd=0,lwd.tick=1)
		title(xlab="Distance from water's edge (m)",line=2.3)
	points(bolz$DistShift,bolz$Effect,cex=.75,pch=as.numeric(bolz$Group))
	legend(140,1.15,legend='B',bty='n',cex=2)
	legend(50,1.1,legend='BOLZ1 right bank',bty='n')	
	abline(v=130,lty=5)
mtext('Portion aquatic',side=2,line=2.8,at=1.21)
par(mar=c(.5,.5,.5,.5))
plot(bolz$DistShift,bolz$Effect, type='n',xlab='',ylab='',axes=F)
lim <- par()
	rasterImage(bolzimg, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
	legend(145,1.1,legend='C',bty='n',cex=2)
	legend(60,1,legend='BOLZ1 right bank',bty='n')		
par(xpd=NA)
	arrows(19,1.33,54,.1,col='grey85',lwd=4)
	arrows(136,1.33,160,.4,col='grey85',lwd=4)
dev.off()
}
#f1();f1(Type='png')


##### Make figure to show effect of wall boundaries #####

## Use RESI5LB for preds and DSRB2LB for isotopes 
res5<-mega8[mega8$SiteBank=='RESI5LB',]
res5XS<-islandsXS[substr(islandsXS$Site,1,5)=='RESI5'&!is.na(islandsXS$Height),]
dsr2<-isos.mega4[isos.mega4$SiteBank=='DSRB2LB',] 
prednames<-as.factor(c('PrBeet','PrSpHu','PrSpWe'))

## Plot figure panels of predator and isotope wall break data
f2<-function(Type='pdf'){
if(Type=='png'){png('Figures & Tables/Wall Break.png',width=6.5,height=8,res=300,bg=0,units='in')}
else{pdf('Figures & Tables/Wall Break.pdf',width=6.5,height=8)}
oldpar2<-par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(3.3,3.6,0.5,0.5),cex=1)
plot(res5$DistShift,res5$Effect,bty='l',xaxt='n',yaxt='n',ylab='',xlab='',type='n')
	axis(1,at=seq(0,15,3),lwd=0,lwd.tick=1)
	axis(2,las=2,lwd=0,lwd.tick=1)
	points(res5$DistShift,res5$Effect,pch=match(res5$Group,prednames))
	lines(c(7,7),c(-1,1),lty=3,lwd=2)
	legend(1,1.075,legend='RESI5 left bank',bty='n')
	legend(13,1.13,legend='A',bty='n',cex=2)
	legend(9,1.075,legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3),bty='n',pt.cex=.75)
plot(c(0,max(dsr2$Dist)),c(0,max(dsr2$PercAqRel)),bty='l',xaxt='n',yaxt='n',ylab='',xlab='',type='n')
	axis(1,lwd=0,lwd.tick=1)
	axis(2,las=2,lwd=0,lwd.tick=1)
	points(dsr2$Dist,dsr2$PercAqRel,pch=match(dsr2$Group,prednames))
	lines(c(128,128),c(-1,1),lty=3,lwd=2)
	lines(c(240,240),c(-1,1),lty=5)
	legend(20,1.07,legend='DSRB2 left bank',bty='n')
	legend(222,1.13,legend='B',bty='n',cex=2)
mtext("Distance from water's edge (m)",side=1,line=2.2)
mtext('Portion aquatic',side=2,line=-1,outer=T,adj=.53)
par(oldpar2)	
dev.off()
}
#f2();f2(Type='png')


##### Find SiteBank data subsets with positive (non-decaying) curves #####

## Make Regs.All function
Regs.All<-function(Data){
	Regs.All.data<-list()
	nas<-as.numeric(rep(NA,length(unique(Data$SiteBank))))
	for(i in 1:length(unique(Data$SiteBank))){
		subs<-subset(Data,SiteBank==unique(Data$SiteBank)[i])
		if(dim(subs)[1]>2){
			subs.Lin<-lm((Effect+.05)~DistShift,data=subs)
				results.Lin<-round(c(trueAIC(subs.Lin)[1],summary(subs.Lin)$adj.r.squared,summary(subs.Lin)$coefficients[2,3],summary(subs.Lin)$coefficients[2,4]),4)
			subs.Exp<-lm(logEffect~DistShift,data=subs)
				results.Exp<-round(c(trueAIC(subs.Exp)[1],summary(subs.Exp)$adj.r.squared,summary(subs.Exp)$coefficients[2,3],summary(subs.Exp)$coefficients[2,4]),4)
			subs.Pow<-lm(logEffect~logDist,data=subs)
				results.Pow<-round(c(trueAIC(subs.Pow)[1],summary(subs.Pow)$adj.r.squared,summary(subs.Pow)$coefficients[2,3],summary(subs.Pow)$coefficients[2,4]),4)
			Regs.All.data[[i]]<-as.data.frame(rbind(results.Lin,results.Exp,results.Pow))
				colnames(Regs.All.data[[i]])<-c('AIC','AdjR^2','Dist.t.val','Dist.p.val')
				rownames(Regs.All.data[[i]])<-c('Linear','Neg. Exp.', 'Inv. Power')
			}
		else{Regs.All.data[[i]]<-NA
			nas[i]<-i
			}
		}	
		names(Regs.All.data)<-unique(Data$SiteBank)
		Regs.All.data<-Regs.All.data[is.na(nas)]
		Regs.All.data
	}
	
## Using preForest and Reasonable data
	## Toss if any of lin, exp, or pow show a positive curve
	## Also toss COWE1RB, which is technically negative, but only because of a distinction between webbed spiders and other groups		
Regs.All.preForest<-Regs.All(mega8preForest)
	Regs.Keep.preForest<-Regs.All.preForest[sapply(Regs.All.preForest, function(x) max(x$Dist.t.val)<0)]
		Regs.Keep.preForest[['COWE1RB']]<-NULL
	Regs.Toss.preForest<-Regs.All.preForest[sapply(Regs.All.preForest, function(x) max(x$Dist.t.val)>0)]	
		Regs.Toss.preForest[['COWE1RB']]<-Regs.All.preForest[['COWE1RB']]
	mega9preForest<-mega8preForest[mega8preForest$SiteBank %in% names(Regs.Keep.preForest),]
Regs.All.Reasonable<-Regs.All(mega8Reasonable)
	Regs.Keep.Reasonable<-Regs.All.Reasonable[sapply(Regs.All.Reasonable, function(x) max(x$Dist.t.val)<0)]
		Regs.Keep.Reasonable[['COWE1RB']]<-NULL
	Regs.Toss.Reasonable<-Regs.All.Reasonable[sapply(Regs.All.Reasonable, function(x) max(x$Dist.t.val)>0)]	
		Regs.Toss.Reasonable[['COWE1RB']]<-Regs.All.Reasonable[['COWE1RB']]
	mega9Reasonable<-mega8Reasonable[mega8Reasonable$SiteBank %in% names(Regs.Keep.Reasonable),]

## Make table of tossed data
toss<-data.frame(levels(mega8Reasonable$OrderClass))
	colnames(toss)<-'OrderClass'
	toss$Percent.not.fit<-toss$Number.not.fit<-toss$Number.of.sites<-0
for(i in 1:dim(toss)[1]){
	toss$Number.of.sites[i]<-length(unique(mega8Reasonable[mega8Reasonable$OrderClass==toss$OrderClass[i],'SiteBank']))
	toss$Number.not.fit[i]<-length(unique(mega8Reasonable[mega8Reasonable$OrderClass==toss$OrderClass[i]&mega8Reasonable$SiteBank%in%names(Regs.Toss.Reasonable),'SiteBank']))
	toss$Percent.not.fit[i]<-round(length(unique(mega8Reasonable[mega8Reasonable$OrderClass==toss$OrderClass[i]&mega8Reasonable$SiteBank%in%names(Regs.Toss.Reasonable),'SiteBank']))/toss$Number.of.sites[i],4)
	}
toss.total<-c('Total',as.numeric(sum(toss$Number.of.sites)),as.numeric(sum(toss$Number.not.fit)),round(sum(toss$Number.not.fit)/sum(toss$Number.of.sites),4))
	toss$OrderClass<-as.character(toss$OrderClass)
toss<-rbind(toss,toss.total)
toss$Percent.of.all.non.fit.sites<-round(as.numeric(toss$Number.not.fit)/as.numeric(toss.total[3]),4)
#write.csv(toss,'Figures & Tables/R CSV Exports/PredatorNonFits.csv',row.names=F)

## Plot some good and bad examples
goodbadData<-isos.mega4[isos.mega4$SiteBank=='RESI4RB'|isos.mega4$SiteBank=='RESI2LB',c('Code','Site','Bank','Dist','Group','shortSite','SiteBank','PercAqRel')]
	colnames(goodbadData)<-c('Code','Site','Bank','DistShift','Group','shortSite','SiteBank','Effect')
	gb2<-rbind(mega8[mega8$SiteBank=='DSRB4LB'|mega8$SiteBank=='SAVA1LB',c('Code','Site','Bank','DistShift','Group','shortSite','SiteBank','Effect')],goodbadData)
goodbads<-c('RESI4RB','RESI2LB','DSRB4LB','SAVA1LB')
legnames<-c('RESI4 right bank','RESI 2 left bank','DSRB4 left bank','SAVA1 left bank')
f3<-function(Type='pdf'){
if(Type=='png'){png('Figures & Tables/GoodBad.png',width=6.5,height=6.5,units='in',bg=0,res=300)}
else{pdf('Figures & Tables/GoodBad.pdf',width=6.5,height=6.5)}
oldpar3<-par(mfrow=c(2,2),mar=c(3,3,.1,.1),oma=c(.8,.8,0,0),cex=1)
for(i in 1:4){
	Data<-gb2[gb2$SiteBank==goodbads[i],]
	plot(c(0,max(Data$DistShift)),c(0,max(Data$Effect)),xaxt='n',yaxt='n',axes=FALSE,type='n',xlab='',ylab='')
		points(Data$DistShift,Data$Effect,pch=as.numeric(Data$Group))
		box(bty='l')
		axis(1,lwd=0,lwd.tick=1)
		if(i%%2==0){axis(2,lwd=0,lwd.tick=1,las=2,labels=FALSE)}
		else{axis(2,lwd=0,lwd.tick=1,las=2)}
		if(i==4){legend(18,.95,legend=c('Beetles', 'Hunting spiders', 'Webbed spiders'),pch=c(1:3),bty='n')}
		legend('topright',legend=LETTERS[i],bty='n',cex=2,adj=c(0,-.6))
		legend('top',legend=legnames[i],bty='n',adj=c(-.05,-.6))
}
mtext("Distance from water's edge (m)",side=1,line=-.5,outer=TRUE,adj=.6)
mtext("Portion aquatic",side=2,line=-.2,outer=TRUE,adj=.55)
dev.off()
}
#f3();f3(Type='png')

##### Get stream signatures for each site subset #####

## Make sitesigs function
sitesigs<-function(Data){
	sitesigs.data<-data.frame(matrix(,nrow=1,ncol=2))
	for(i in 1:length(unique(Data$SiteBank))){
		sitesubs<-subset(Data,SiteBank==unique(Data$SiteBank)[i])
		site.Pow<-lm(logEffect~logDist,data=sitesubs)
			sig50<-round(exp((log(.5)-site.Pow$coefficients[1])/(site.Pow$coefficients[2])),2)
			sig90<-round(exp((log(.1)-site.Pow$coefficients[1])/(site.Pow$coefficients[2])),2)
		sitesigs.data[i,]<-cbind(sig50,sig90)
		}
		rownames(sitesigs.data)<-unique(Data$SiteBank)
		colnames(sitesigs.data)<-c('50%','90%')
		sitesigs.data<-sitesigs.data[order(rownames(sitesigs.data)),]	
	sitesigs.data
	}

## For preForest and Reasonable data
sitesigs.preForest<-sitesigs(mega9preForest)
sitesigs.Reasonable<-sitesigs(mega9Reasonable)
	## Reasonable results make more sense than preForest and are more readily justifiable as well
	## Includes all walls if going from 0 (i.e., walls not removed).
mega9<-mega9Reasonable<-droplevels(mega9Reasonable[order(mega9Reasonable$Code),])
mega9preForest<-droplevels(mega9preForest[order(mega9preForest$Code),])

## Subset only "normal" (non-braided, non-wall) sites for comparison with meta-analysis
normalsites<-c('ARZO1LB','BALL1LB','BELG1LB','BOLZ1RB','DANU1RB','DFOR1RB','DMDW1RB','DRAV1LB','DSRB2LB','ELBE1LB','ELBE2LB','ELBE3LB','LTEN1RB','RESI1RB','RESI2LB','RESI2RB','SAVA1LB')
mega9Normal<-mega9[mega9$SiteBank%in%normalsites,]
sitesigs.Normal<-sitesigs(mega9Normal)


##### Test model fits for data subsets #####

## Subset data by stream, bank
mega9subs<-mega9
	#mega9subs$subs<-mega9subs$SiteBank
	mega9subs$subs<-paste(mega9subs$SiteBank,mega9subs$Group)
subs.list<-list()
bads<-as.numeric(NA)
for(i in 1:length(unique(mega9subs$subs))){
	subs<-subset(mega9subs,subs==unique(mega9subs$subs)[i])
	if(dim(subs)[1]>2){
		subs.Lin<-lm((Effect+.05)~Dist,data=subs)
			results.Lin<-round(c(trueAIC(subs.Lin)[1],summary(subs.Lin)$adj.r.squared,summary(subs.Lin)$coefficients[2,3],summary(subs.Lin)$coefficients[2,4]),4)
		subs.Exp<-lm(logEffect~Dist,data=subs)
			results.Exp<-round(c(trueAIC(subs.Exp)[1],summary(subs.Exp)$adj.r.squared,summary(subs.Exp)$coefficients[2,3],summary(subs.Exp)$coefficients[2,4]),4)
		subs.Pow<-lm(logEffect~logDist,data=subs)
			results.Pow<-round(c(trueAIC(subs.Pow)[1],summary(subs.Pow)$adj.r.squared,summary(subs.Pow)$coefficients[2,3],summary(subs.Pow)$coefficients[2,4]),4)
		subs.list[[i]]<-as.data.frame(rbind(results.Lin,results.Exp,results.Pow))
			colnames(subs.list[[i]])<-c('AIC','AdjR^2','Dist.t.val','Dist.p.val')
			rownames(subs.list[[i]])<-c('Linear','Neg. Exp.', 'Inv. Power')
		}
	else{subs.list[[i]]<-NA
		bads[i]<-i
		}
	}	
	names(subs.list)<-unique(mega9subs$subs)
	subs.list<-subs.list[is.na(bads)]	
	
## See how often Inv power is the best, vs. neg. exp and linear
best1<-vector()
for(i in 1:length(subs.list)){best1[i]<-order(subs.list[[i]]$AIC)[1]}
best2<-ifelse(best1==1,'Linear',ifelse(best1==2,'Neg. Exp','Inv. Power'))
best<-cbind(table(best2),round(table(best2)/length(subs.list),4))
	colnames(best)<-c('Freq','Percent')
	## For SiteBank, 16/30(53.33%) pow, 14/30 (46.67%) exp.
	## For SiteBankGroup, 35/71(47.95%) pow, 37/73 (50.68%) exp, 1/73 (1.37%) lin.	

## Create a dataframe of the individual subset study characteristics
studies1<-unlist(strsplit(names(subs.list),' '))
studies<-as.data.frame(cbind(studies1[seq(1,2*length(subs.list),2)],studies1[seq(2,2*length(subs.list),2)]))
colnames(studies)<-c('Site','Group')

## See if there is anything special about the ones that were not best fit by inv.power	
non.Pows<-subs.list[which(best1!=3)]
non.Pows.list<-list()
non.Pows.list$Group<-round(table(studies$Group[which(best1!=3)])/table(studies$Group),4)
non.Pows.list$Site<-round(table(studies$Site[which(best1!=3)])/table(studies$Site),4)
	## Not much of interest here.

	
##### Do basic regressions on data #####
lm.lin<-lm((Effect+.05)~DistShift,data=mega9)
lm.exp<-lm(logEffect~DistShift,data=mega9)
lm.pow<-lm(logEffect~logDist,data=mega9)
#rbind(AIC(lm.lin),trueAIC(lm.exp),trueAIC(lm.pow));AIC(lm.exp,lm.pow)
	## Pow the winner here by ~40 AIC points (but both exp and pow much better than lin), but all R2s very low.

	
##### Try mixed effects models on data #####

## Use Site as the structural variable
lme.lin<-lme((Effect+.05)~DistShift,random=~1|shortSite,data=mega9,method='ML')
lme.exp<-lme(logEffect~DistShift,random=~1|shortSite,data=mega9,method='ML')
lme.pow<-lme(logEffect~logDist,random=~1|shortSite,data=mega9,method='ML')
#rbind(AIC(lme.lin),trueAIC(lme.exp),trueAIC(lme.pow));AIC(lme.exp,lme.pow)
	## Pow ~40 AIC point improvement over exp, and both a very big improvement over linear.	
	
## Use Site+Bank as the structural variable
lmeSB.lin<-lme((Effect+.05)~DistShift,random=~1|SiteBank,data=mega9,method='ML')
lmeSB.exp<-lme(logEffect~DistShift,random=~1|SiteBank,data=mega9,method='ML')
lmeSB.pow<-lme(logEffect~logDist,random=~1|SiteBank,data=mega9,method='ML')
#rbind(AIC(lmeSB.lin),trueAIC(lmeSB.exp),trueAIC(lmeSB.pow));AIC(lmeSB.exp,lmeSB.pow)
	## Probably the most sensible to use SiteBank because it is the truer experimental unit (and bank type is of interest), but Site technically better based on AIC.
aiccs<-as.data.frame(round(rbind(AICc(lmeSB.lin),trueAIC(lmeSB.exp)[2],trueAIC(lmeSB.pow)[2]),2))
	colnames(aiccs)<-'AICc'
	rownames(aiccs)<-c('Linear','Exponential','Power')
	aiccs$delta<-round(c(aiccs$AICc-min(aiccs$AICc)),2)
	aiccs$weight<-round(exp(-.5*aiccs$delta)/sum(exp(-.5*aiccs$delta)),4)
	aiccs<-aiccs[order(-aiccs$AICc),]
#write.csv(aiccs,'Figures & Tables/R CSV Exports/PredatorAICs.csv')

	
##### Compute single variable effects and stream signatures #####

## Make list of environmental variables of interest
vars<-c('Region','Width','Order','OrderClass','Geomorph','Banks','Group','VegFld','VegUp','VegShift')

## Run signature function
sigs.Reasonable<-signatures(mega9Reasonable,'logEffect','logDist','SiteBank',vars)	
sigs.preForest<-signatures(mega9preForest,'logEffect','logDist','SiteBank',vars)
sigs.Normal<-signatures(mega9Normal,'logEffect','logDist','SiteBank',vars)
	## Reasonable still makes the most sense, Normal and preForest aren't too different. 

## Write files for single model and signature results
w1<-function(){	
write.csv(sigs.Reasonable$fits,'Figures & Tables/R CSV Exports/PredatorEffectSinglesModels.csv')
out_file<-file('Figures & Tables/R CSV Exports/PredatorSignatures.csv', open='a')
	for(i in seq_along(sigs.Reasonable$sigs)){
		write.table(names(sigs.Reasonable$sigs)[i],file=out_file,sep=',',dec='.',quote=FALSE,col.names=FALSE,row.names=FALSE)  
		write.table(sigs.Reasonable$sigs[[i]],file=out_file,sep=',',dec='.',quote=FALSE,col.names=NA,row.names=TRUE)
		}
	close(out_file)
}
#w1()

##### Build overall multivariate model #####

## Decide on whether it is better to add stream order as Order, Width class, or Channel width (just using exp model for ease)
lmeWid.pow<-lme(logEffect~logDist+Width,random=~1|SiteBank,data=mega9,method='ML')
lmeOrd.pow<-lme(logEffect~logDist+Order,random=~1|SiteBank,data=mega9,method='ML')
lmeOC.pow<-lme(logEffect~logDist+OrderClass,random=~1|SiteBank,data=mega9,method='ML')
#AIC(lmeWid.pow,lmeOrd.pow,lmeOC.pow)
	## OC ~9 AIC points better than Order and 4 better than Width. Go with Order Class because interpretation is easier and results from singles.lme make more sense.

## Add variables en masse using dredge
	
## Start with additive terms only
#modtablep<-dredge(lme(logEffect~logDist+Group+OrderClass+Geomorph+Banks+Region+VegFld,random=~1|SiteBank,data=mega9,method='ML'),fixed='logDist',extra='adjR^2');subset(modtablep,delta<4)
	## Focusing on VegFld instead of VegShift or VegUp because it seems to improve these models the most and adding all of them is redundant, plus it's results from singles.lme are the most logical.
	## Also, VegShift can't be added with VegUp or VegFld due to singularities.
	
## Now try interaction terms
#modtablex<-dredge(lme(logEffect~logDist+logDist*Group+logDist*OrderClass+logDist*Geomorph+logDist*Banks+logDist*Region+logDist*VegFld,random=~1|SiteBank,data=mega9,method='ML'),fixed='logDist',extra='adjR^2');subset(modtablex,delta<4)
	## Improvement over additives alone.

## Export models
	## When subsetting best models, use delta<4.  This means the difference between AIC(i) and AIC(min) is <4.  "Substantial model support" according to Burnham & Anderson (2002) pp70.
#write.csv(modtablex,'Figures & Tables/R CSV Exports/PredatorEffectModels.csv')	
	

##### Compare predicted vs. observed values at breaks #####

## Add predicted (expected) and observed values to breaks dataframe
breaks$Wall.Exp<-breaks$Forest.Exp<-as.numeric(rep(NA,dim(breaks)[1]))
for(i in 1:dim(breaks)[1]){
	subsite.lm<-lm(logEffect~logDist,data=mega8[mega8$SiteBank==breaks$SiteBank[i],])
	breaks$Forest.Exp[i]<-exp(subsite.lm$coefficients[1]+subsite.lm$coefficients[2]*log(breaks$predpreForest[i]+.05))
	breaks$Wall.Exp[i]<-exp(subsite.lm$coefficients[1]+subsite.lm$coefficients[2]*log(breaks$predWall[i]+.05))
	}
forestbreaks<-mega8[!is.na(mega8$preForest)&mega8$Dist==mega8$preForest&mega8$SiteBank%in%mega9$SiteBank,c('SiteBank','Effect')]
	forestmeans<-data.frame(tapply(forestbreaks$Effect,forestbreaks$SiteBank,mean))	
	breaks$Forest.Obs<-forestmeans[match(breaks$SiteBank,rownames(forestmeans)),1]
wallbreaks<-mega8[!is.na(mega8$Wall)&mega8$Dist==mega8$Wall&mega8$SiteBank%in%mega9$SiteBank,c('SiteBank','Effect')]
	wallmeans<-data.frame(tapply(wallbreaks$Effect,wallbreaks$SiteBank,mean))	
	breaks$Wall.Obs<-wallmeans[match(breaks$SiteBank,rownames(wallmeans)),1]
## Run paired t-test to see if observed PercAqRel are higher at breaks than predicted
tFor<-t.test(breaks$Forest.Obs,breaks$Forest.Exp,alternative='greater',paired=T,var.equal=F,na.action='na.omit')
tWall<-t.test(breaks$Wall.Obs,breaks$Wall.Exp,alternative='greater',paired=T,var.equal=F,na.action='na.omit')
diffFor<-round(mean(breaks$Forest.Obs,na.rm=T)-mean(breaks$Forest.Exp,na.rm=T),4)
diffWall<-round(mean(breaks$Wall.Obs,na.rm=T)-mean(breaks$Wall.Exp,na.rm=T),4)
incFor<-round(100*diffFor/mean(breaks$Forest.Exp,na.rm=T),2)
incWall<-round(100*diffWall/mean(breaks$Wall.Exp,na.rm=T),2)
	tabFor<-cbind(diffFor,incFor,round(tFor$statistic,4),round(tFor$p.value,4))
	tabWall<-cbind(diffWall,incWall,round(tWall$statistic,4),round(tWall$p.value,4))
		tabForWall<-rbind(tabFor,tabWall)
			rownames(tabForWall)<-c('Forest','Wall')
			colnames(tabForWall)<-c('PercAqDiff','RelPercChange','t','p')
#write.csv(tabForWall,'Figures & Tables/R CSV Exports/PredatorForestWallObsVsExp.csv')
	## Both significant, with large increases in observed vs. model expected.


##### See if islands differ from true banks #####
	
## Make an island site dataframe
islands2<-droplevels(mega1[nchar(as.character(mega1$Site))>5,])

## Rearrange alphabetically for ease
islands2.1<-droplevels(islands2[order(as.character(islands2$Code)),])

## Get rid of non-predator groups and night/malaise/ground sampling
islands3<-droplevels(islands2.1[islands2.1$Trophic==2&islands2.1$Method!='Night Sampling'&islands2.1$Method!='Malaise'&islands2.1$Method!='Ground',])

## Add sweeps to each pitfall/bowl replicate
islands3Swp<-droplevels(islands3[islands3$Method=='Sweep'&islands3$Abund>0,])
	sumsAbund<-tapply(islands3Swp$Abund,islands3Swp$Code,sum)
	sumsBiomass<-tapply(islands3Swp$Biomass,islands3Swp$Code,sum)
islands4PitBowl<-islands3[islands3$Method=='Pitfall'|islands3$Method=='Bowl',]
islands4Swp<-islands3Swp[islands3Swp$Code%in%islands4PitBowl$Code==FALSE,]
islands4<-rbind(islands4PitBowl,islands4Swp)
	islands4$Abund<-ifelse(is.na(match(islands4$Code,names(sumsAbund))),islands4$Abund,islands4$Abund+sumsAbund[match(islands4$Code,names(sumsAbund))])
	islands4$Biomass<-ifelse(is.na(match(islands4$Code,names(sumsBiomass))),islands4$Biomass,islands4$Biomass+sumsBiomass[match(islands4$Code,names(sumsBiomass))])
	islands4<-droplevels(subset(islands4,select=-Method))
		rownames(islands4)<-c(1:dim(islands4)[1])
		
## Create some data columns for future use
islands4$shortSite<-as.factor(substr(islands4$Site,1,5))
islands4$SiteBank<-as.factor(paste(islands4$Site,islands4$Bank,sep=''))
islands4$AbBiom<-islands4$Abund*islands4$Biomass

## Percent of SiteBank max value
islands5<-islands4
maxs.islands<-tapply(islands5$AbBiom,islands5$shortSite,max)
islands5$Effect<-as.numeric(islands5$AbBiom/maxs.islands[islands5$shortSite])

## Add cross section distances to the dataframe
m3<-match(paste(islands5$SiteBank,islands5$Dist,sep=''),paste(islandsXS$Site,islandsXS$Bank,islandsXS$Dist,sep=''))
islands5$DistXS<-islandsXS$DistXS[m3]
islands5$BankIsland<-as.factor(ifelse(paste(islands5$Site,islands5$Bank)=='CORN1.1 RB'|paste(islands5$Site,islands5$Bank)=='CORN1.6 LB'|paste(islands5$Site,islands5$Bank)=='FELL1.1 RB'|paste(islands5$Site,islands5$Bank)=='ISLA1.1 RB'|paste(islands5$Site,islands5$Bank)=='ISLA1.7 LB'|paste(islands5$Site,islands5$Bank)=='RESI1.1 RB'|paste(islands5$Site,islands5$Bank)=='RESI1.2 LB'|paste(islands5$Site,islands5$Bank)=='RESI5.1 RB'|paste(islands5$Site,islands5$Bank)=='RESI5.2 LB'|paste(islands5$Site,islands5$Bank)=='ELBE1.2 LB','Bank','Island'))

## Run model to see if PercAq is higher on islands than banks
islands.mod<-lme(Effect~BankIsland,random=~1|shortSite,data=islands5)
	## No significant difference at this coarse level of analysis. 
	
## Plot
isla<-islands5[islands5$shortSite=='ISLA1'&islands5$Effect<1000,]
corn<-islands5[islands5$shortSite=='CORN1'&islands5$Effect<1000,]	
fell<-islands5[islands5$shortSite=='FELL1',]	
resi1<-islands5[islands5$shortSite=='RESI1'&islands5$Effect<1000,]	
resi5<-islands5[islands5$shortSite=='RESI5',]	
elbe<-islands5[islands5$shortSite=='ELBE1'&islands5$Effect<50,]
danu<-islands5[islands5$shortSite=='DANU1',]
#plot(isla$DistXS,isla$Effect);plot(corn$DistXS,corn$Effect);plot(fell$DistXS,fell$Effect);plot(resi1$DistXS,resi1$Effect);plot(resi5$DistXS,resi5$Effect);plot(elbe$DistXS,elbe$Effect);plot(danu$DistXS,danu$Effect)


##### Look at night samples #####

## Get night-only, aquatic adult-only data
night1<-mega1[mega1$Method=='Night Sampling'&substr(mega1$Group,1,4)=='AeAq',]

## Add a site-distance and Abundance*Biomass variables
night1$SiteDist<-as.factor(paste(substr(night1$Site,1,5),night1$Dist,sep='-'))
night1$AbBiom<-ifelse(is.na(night1$Biomass),night1$Abund,night1$Abund*night1$Biomass)

## Get site-dist sums
breaks<-c('ARZO2','BELG1','BOLZ1','FELL1','ISLA1','RESI1','RESI5')
night2<-data.frame(tapply(night1$AbBiom,night1$SiteDist,sum))
night3<-rbind(night2,night2[rownames(night2)=='ARZO1-1',1],0,0)
	rownames(night3)[(dim(night3)[1]-2):dim(night3)[1]]<-c('ARZO2-1','BATE1-32','COWE1-32')
	colnames(night3)<-'AbBiom'
	night3$Site<-substr(rownames(night3),1,5)
	night3$Dist<-as.numeric(substr(rownames(night3),7,9))
	night3$Loc<-as.factor(ifelse(night3$Dist<8,'Bank','Upland'))
	night3$Break<-as.factor(ifelse(night3$Site%in%breaks,'Break','None'))
night4<-night3[order(rownames(night3)),]
night5<-night4[rownames(night4)!='ARZO1-40',]

## Get difference between upland and bank values
night6<-data.frame(unique(night5$Site))
	colnames(night6)<-'Site'
	night6$Break<-night5[match(night6$Site,night5$Site),'Break']
	night6$AbBiomBank<-night5[night5$Loc=='Bank','AbBiom']
	night6$AbBiomUpland<-night5[night5$Loc=='Upland','AbBiom']
	night6$Difference<-night6$AbBiomUpland-night6$AbBiomBank
	night6$DiffRel<-night6$Difference/night6$AbBiomBank
	
## See if there is a significant decrease from bank to upland
night.all<-t.test(night6$DiffRel)
night.none<-t.test(night6[night6$Break=='None','DiffRel',],alternative='less')
night.break<-t.test(night6[night6$Break=='Break','DiffRel'],alternative='greater')
	## No overall significant difference (t = -.8931, p = .3869).
	## But for sites without breaks, there was a significant, 79.72% decrease (t = -9.6309, p = 1.37e-5).
	## 39.87% INCREASE from banks to upland for sites with breaks, but difference not significant (t = .8457, p = .2151).

## Plot night figures
f4<-function(Type='pdf'){
if(Type=='png'){png('Figures & Tables/Night.png',width=6.5,height=6.5,units='in',bg=0,res=300)}
else{pdf('Figures & Tables/Night.pdf',width=6.5,height=6.5)}
oldpar4<-par(mar=c(3.5,4,.4,.1))
boxplot(night6$DiffRel~as.factor(-as.numeric(night6$Break)),xaxt='n',yaxt='n',axes=FALSE,ylab='Difference from 0 m sample (%)')
axis(1,at=c(1,2),labels=c('Gradual/none','Sharp'),lwd=0,lwd.tick=1)
axis(2,at=seq(-1,2.5,.5),labels=seq(-1,2.5,.5)*100,las=2,lwd=0,lwd.tick=1)
mtext('Forest transition',side=1,line=2.5)
box(bty='l')
par(oldpar4)
dev.off()
}
#f4();f4(Type='png')


##### Compare isotope PercAq to AbBiom PercAq

## Create dataframe
isopred<-isos.mega4[,c('Code','Site','Bank','Dist','Group','shortSite','SiteBank','PercAq','PercAqAdj','PercAqRel')]
	isopred$PredRel<-isopred$AbBiomTopCut<-isopred$AbBiom<-NA

## Add predator abundance*biomass data to dataframe
for(i in 1:dim(isopred)[1]){
	t1<-mega8[which(mega8$Code%in%isopred[i,'Code']),c('AbBiom','AbBiomTopCut','Effect')]
	t2<-colMeans(t1)
	isopred$AbBiom[i]<-t2[1]
	isopred$AbBiomTopCut[i]<-t2[2]
	isopred$PredRel[i]<-t2[3]
}

## Test correlation with linear model
ip1<-lm(PercAqRel~PredRel,data=isopred)
plot(isopred$PercAqRel~isopred$PredRel)
abline(ip1)
summary(ip1)
	## Positive correlation (t = 2.975, df = 423, p = 0.0031, R^2 = 0.0181, slope = 0.2036)
	
## Try mixed-effects model with SiteBank as the random effect
ip2<-lme(PercAqRel~PredRel,random=~1|SiteBank,data=isopred,method='ML',na.action=na.omit)
	## Very similar result (t = 2.9829, df = 387, p = 0.003)