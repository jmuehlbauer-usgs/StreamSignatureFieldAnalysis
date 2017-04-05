##### Set working directory, read in raw data #####

## CSVs originating from larger Stream Signature Excel sheet, plus ISLA site cross section
setwd('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis')
raw0<-read.csv('R Code & Input Data/CSVExport.csv',header=T,na.strings=c('n/a',' ','#N/A'))
env00<-read.csv('R Code & Input Data/All Sites Data.csv',header=T,na.strings=c('n/a','N/A','#N/A'))
islandsXS<-read.csv('R Code & Input Data/Islands Cross Section Data.csv',header=T,na.strings=c(' '))
breaks<-read.csv('R Code & Input Data/Breaks.csv')


##### Load some useful functions and libraries #####

## Libraries
library(MuMIn)
library(nlme)
library(lme4)

## Functions
source('C:/Users/jmuehlbauer/Documents/Misc/Trainings & Info/Stats & Data Analysis/R/AICforLogResponse.r')
source('R Code & Input Data/StreamSignatureFunction.r')
source('R Code & Input Data/PercentAquaticFunction.r')
source('R Code & Input Data/PlottingFunction.r')


##### Clean up raw dataframe #####

## Delete any spurious columns, name columns
raw1<-raw0[,c(1:11)]
colnames(raw1)<-c('Code','Site','Bank','Dist','Group','Abund','Biomass','PercC','PercN','d13C','d15N')

## Add trophic level and site names
raw2<-raw1
	raw2$Trophic<-rep(4,dim(raw2)[1])
raw2$Trophic<-ifelse(substr(raw2$Group,1,2)=='Pr',2,ifelse(substr(raw2$Group,1,2)=='He',1,ifelse(substr(raw2$Group,3,6)=='Detr',0,ifelse(raw2$Group=='AqAlga',0,1.5))))
raw2$shortSite<-as.factor(substr(raw2$Site,1,5))
raw2$SiteBank<-as.factor(paste(raw2$shortSite,raw2$Bank,sep=''))

## Convert Aq distances to NAs, and distances to numeric
raw2$Dist<-suppressWarnings(as.numeric(as.character(raw2$Dist)))

## Remove ground sampling sites without specific locations
raw3<-droplevels(subset(raw2,Bank!='ALL'))
	

##### Clean up env datasheet #####

## Convert columns to factors as relevant
env00$Order<-as.factor(env00$Order)

## Add site rows for sites repeated sampled over 2 years (Coweeta Area "2" sites)
env.sites2<-env00[env00$Site=='BALL1'|env00$Site=='COWE1'|env00$Site=='LTEN1',]
	env.sites2$Site<-c('BALL2','BALL2','COWE2','COWE2','LTEN2')
env0<-rbind(env00,env.sites2)

	
##### Create isotopes dataframe #####

## Remove samples with no isotope results
iso<-subset(raw3,PercC!='NA' & PercN!='NA' & PercC!=0 & PercN!=0)


##### Compute percent aquatic values #####

isoSites<-unique(substr(iso$Site,1,5))
isoSites<-as.factor(isoSites[order(isoSites)])
isoRel<-list()
for(i in 1:length(isoSites)){
	isoRel[[i]]<-PercAq(iso,isoSites[i])
		names(isoRel)[i]<-isoSites[i]
	}
names(isoRel)<-isoSites
	
## Plot
prednames<-as.factor(c('PrBeet','PrSpHu','PrSpWe'))
f1<-function(){
pdf('Figures & Tables/Raw Figures/IsotopeRaw.pdf',width=6.5,height=6.5)
for(i in 1:length(isoSites)){
	preds<-isoRel[[i]][isoRel[[i]]$Trophic==2,]
	preds<-droplevels(preds)
	oldpar<-par(mfcol=c(2,2))
	if(length(na.omit(preds$PercCaqVegRel))>0){
	plot(preds$Dist,preds$PercCaqVegRel,col='white',xlab='Distance from bank (m)',ylab='% Aquatic C (plant baseline)')
		points(preds$Dist,preds$PercCaqVegRel,pch=match(preds$Group,prednames))
		legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
	title(names(isoRel[i]),outer=T,line=-2)
	plot(preds$Dist,preds$PercNaqVegRel,col='white',xlab='Distance from bank (m)',ylab='% Aquatic N (plant baseline)')
		points(preds$Dist,preds$PercNaqVegRel,pch=match(preds$Group,prednames))
		legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))}
	if(length(na.omit(preds$PercCaqAnRel))>0){
	plot(preds$Dist,preds$PercCaqAnRel,col='white',xlab='Distance from bank (m)',ylab='% Aquatic C (animal baseline)')
		points(preds$Dist,preds$PercCaqAnRel,pch=match(preds$Group,prednames))
		legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
	plot(preds$Dist,preds$PercNaqAnRel,col='white',xlab='Distance from bank (m)',ylab='% Aquatic N (animal baseline)')
		points(preds$Dist,preds$PercNaqAnRel,pch=match(preds$Group,prednames))
		legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
	title(names(isoRel[i]),outer=T,line=-2)}
	else{plot(0,0,col='white')
	title(names(isoRel[i]),outer=T,line=-2)
	text(0,0,'No Data')}
	par(oldpar)
	}
	dev.off()
	}
#f1()	
	## Note that BALL2, HAWR1, MUDC1, and LTEN2 don't have enough isotope data to plot anything. Consider removing from all analysis? Also, COWE2 doesn't exist as an actual site (from Envirothon sampling instead).

##### Combine isotope results with all data in big table #####

## Combine isotope list elements into 1 datasheet
comb<-do.call('rbind',isoRel)
comb$PercCaq<-ifelse(is.na(comb$PercCaqAnRel),comb$PercCaqVegRel,comb$PercCaqAnRel)
comb$PercNaq<-ifelse(is.na(comb$PercNaqAnRel),comb$PercNaqVegRel,comb$PercNaqAnRel)	

## Combine env and isotope data to abundance table	
mega<-raw3
mega$PercC<-round(mega$PercC/100,4)
mega$PercN<-round(mega$PercN/100,4)
m1<-match(mega[,1],comb[,1])
	mega$d13Cadj<-round(comb$d13Cadj[m1],2)
	mega$d15Nadj<-round(comb$d15Nadj[m1],2)
	mega$PercCaq<-round(comb$PercCaq[m1],4)
	mega$PercNaq<-round(comb$PercNaq[m1],4)
m2<-match(substr(mega[,2],1,5),env0[,1])
	mega$Region<-env0$Region[m2]
	mega$Width<-env0$Channel.width[m2]
	mega$Order<-env0$Order[m2]
	mega$OrderClass<-env0$Order.class[m2]
	mega$Geomorph<-env0$Geomorphology[m2]
	mega$Banks<-as.factor(ifelse(mega$Bank=='LB',as.character(env0$Left.bank[m2]),ifelse(mega$Bank=='RB',as.character(env0$Right.bank[m2]),NA)))
m3<-match(paste(substr(mega$Site,1,5),substr(mega$Bank,1,1)),paste(env0$Site,substr(env0$Bank,1,1)))	
	mega$Banks<-env0$Bank.type[m3]
	mega$VegFld<-env0$Floodplain.vegetation[m3]
	mega$VegUp<-env0$Upland.vegetation[m3]
	mega$VegShift<-env0$Vegetation.shift[m3]	

	
##### Compute percent aquatic #####

## Drop aquatic samples, non-predators, and extemporaneous sites
drops<-c('BALL2','HAWR1','LTEN2','MUDC1','COWE2')
mega1<-droplevels(mega[!is.na(mega$Dist)&mega$Trophic==2&mega$Site%in%drops==FALSE,])

## Drop middle island/braid locations, leaving only true banks
mega2.0<-droplevels(mega1[nchar(as.character(mega1$Site))==5,])
mega2<-rbind(mega2.0,mega1[
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
mega2<-droplevels(mega2[order(as.character(mega2$SiteBank)),])

## Drop SiteBanks with no predator isotope data
samplenums<-vector()
for(i in 1:length(unique(mega2$SiteBank))){
	sitebank<-mega2[mega2$SiteBank==unique(mega2$SiteBank)[i],]
	samplenums[i]<-dim(sitebank[!is.na(sitebank$PercCaq),])[1]
	}
	names(samplenums)<-unique(mega2$SiteBank)
mega3<-droplevels(mega2[mega2$SiteBank%in%names(samplenums[samplenums>0]),])

## Plot percent aquatic figures for C and N by SiteBank
f2<-function(){
pdf('Figures & Tables/Raw Figures/IsotopePercCN-Aq.pdf',width=6.5,height=6.5)
oldpar<-par(mfrow=c(2,2))
for(i in 1:length(unique(mega3$SiteBank))){
	preds<-droplevels(mega3[mega3$SiteBank==unique(mega3$SiteBank)[i],])
	plot(preds$Dist,preds$PercCaq,col='white',xlab='Distance from bank (m)',ylab='% Aquatic C')
		points(preds$Dist,preds$PercCaq,pch=match(preds$Group,prednames))
		legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
	title(unique(mega3$SiteBank)[i])
	plot(preds$Dist,preds$PercNaq,col='white',xlab='Distance from bank (m)',ylab='% Aquatic N')
		points(preds$Dist,preds$PercNaq,pch=match(preds$Group,prednames))
		legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
	title(unique(mega3$SiteBank)[i])	
	}
	par(oldpar)
	dev.off()
}
#f2()
	
##### Clean up percent aquatic data #####
	
## Based on IsotopePercCN-Aq.pdf outputted above, decide whether each site's percent aquatic values should be based on C vs. N.
	## Decisions saved in Breaks.csv file.
mega4.1<-mega3
mega4.1$PercAq<-ifelse(breaks[match(mega3$SiteBank,breaks$SiteBank),'CvsN']=='N',mega3$PercNaq,mega3$PercCaq)

## Toss NAs and SiteBanks with too few data to make inferences or where data are nonsense (e.g., due to too little fractionation)
mega4.2<-droplevels(mega4.1[!is.na(mega4.1$PercAq),])
mega4.3<-mega4.2[mega4.2$SiteBank %in% names(which(table(mega4.2$SiteBank)>3)),]
tsct<-tapply(mega4.3$Dist,mega4.3$SiteBank,max)
	mega4.4<-mega4.3[mega4.3$SiteBank %in% names(tsct[tsct>8]),]
frac<-tapply(mega4.4$PercAq,mega4.4$SiteBank,function(x) (max(x)-min(x))>.3)
	mega4.5<-mega4.4[mega4.4$SiteBank %in% names(frac[frac>.3]),]
mega4<-droplevels(mega4.5)
mega4<-droplevels(mega4.4)

## Plot percent aquatic figures for best of C or N by SiteBank
#plotting('Isotope',mega4,'Dist','PercAq','SiteBank','% Aquatic',Type='pdf')


##### Relativize percent aquatic data #####
	
## Adjust percent aquatic to clean up baselines, based on isotopes figures.pdf
	## Baselines and toplines saved in Breaks.csv file.
mods<-as.data.frame(breaks$SiteBank[breaks$SiteBank%in%mega4$SiteBank])
	colnames(mods)<-'SiteBank'
	mods$baselines<-breaks$isoBaselines[match(mods$SiteBank,breaks$SiteBank)]
	mods$toplines<-breaks$isoToplines[match(mods$SiteBank,breaks$SiteBank)]
	mods$maxs<-mods$mins<-NA
	mods<-droplevels(mods)
for(i in 1:dim(mods)[1]){
	mods$mins[i]<-min(mega4[mega4$SiteBank==mods$SiteBank[i]&mega4$PercAq>=mods$baselines[i],'PercAq'])
	mods$maxs[i]<-max(mega4[mega4$SiteBank==mods$SiteBank[i]&mega4$PercAq<=mods$toplines[i],'PercAq'])
	}
mega4$PercAqAdj<-ifelse(mods[match(mega4$SiteBank,mods$SiteBank),'mins']>mega4$PercAq,mods[match(mega4$SiteBank,mods$SiteBank),'mins'],
	ifelse(mods[match(mega4$SiteBank,mods$SiteBank),'maxs']<mega4$PercAq,mods[match(mega4$SiteBank,mods$SiteBank),'maxs'],mega4$PercAq))

## Relativize data based on these shifts
mega4$PercAqRel<-(mega4$PercAqAdj-mods[match(mega4$SiteBank,mods$SiteBank),'mins'])/(mods[match(mega4$SiteBank,mods$SiteBank),'maxs']-mods[match(mega4$SiteBank,mods$SiteBank),'mins'])

## Plot relativized percent aquatic figures by SiteBank
#plotting('Isotope',mega4,'Dist','PercAqRel','SiteBank','% Aquatic')

## Plot DSRB2LB
dsr2<-mega4[mega4$SiteBank=='DSRB2LB',] 
f3<-function(Type='pdf'){
	if(Type=='png'){png('Figures & Tables/DSRB2.png',width=6.5,height=4.5,units='in',bg=0,res=300)}
	else{pdf('Figures & Tables/DSRB2.pdf',width=6.5,height=4.5)}
oldpar3<-par(oma=c(0,0,0,0),mar=c(3.3,3.6,0.5,0.5),cex=1)
plot(c(0,max(dsr2$Dist)),c(0,max(dsr2$PercAqAdj)),bty='l',xaxt='n',yaxt='n',ylab='',xlab='',type='n')
	axis(1,lwd=0,lwd.tick=1)
	axis(2,las=2,lwd=0,lwd.tick=1)
		title(xlab="Distance from water's edge (m)",line=2.2)
		title(ylab='% Aquatic',line=2.6)
	points(dsr2$Dist,dsr2$PercAqAdj,pch=match(dsr2$Group,prednames))
	legend('topright',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3),bty='n',pt.cex=.75)
	lines(c(128,128),c(-1,.8),lty=3,lwd=2)
	lines(c(240,240),c(-1,.8),lty=5)
par(oldpar3)	
dev.off()
}
#f3(Type='pdf');f3(Type='png')	

## Write data to csv for use in making GoodBad figure (see Predator Stream Signature code)
#write.csv(mega4,'R Code & Input Data/IsotopeMega4forGoodBadFig.csv') 	

## Account for presence of vegetation shifts and walls
	## Decisions saved in Breaks.csv file.
	## Note that unlike Predator code, preForest here is preForest and Reasonable breaks.
mega4$preForest<-breaks[match(mega4$SiteBank,breaks$SiteBank),'isopreForest']
mega4$Forest<-breaks[match(mega4$SiteBank,breaks$SiteBank),'isoForest']
mega4$Wall<-breaks[match(mega4$SiteBank,breaks$SiteBank),'isoWall']

## Create some data columns for future use
mega4$Effect<-mega4$PercAqRel
mega4$logEffect<-log(mega4$Effect+.05)
mega4$logDist<-log(mega4$Dist+.05)	
mega4$AntepreForest<-as.factor(ifelse(mega4$Dist<mega4$preForest|mega4$preForest %in% NA,1,0))
mega4$PostpreForest<-as.factor(ifelse(mega4$AntepreForest==0,1,0))
mega4$AnteForest<-as.factor(ifelse(mega4$Dist<mega4$Forest|mega4$Forest %in% NA,1,0))
mega4$PostForest<-as.factor(ifelse(mega4$AnteForest==0,1,0))
mega4$AnteWall<-as.factor(ifelse(mega4$Dist<mega4$Wall|mega4$Wall %in% NA,1,0))
mega4$PostWall<-as.factor(ifelse(mega4$AnteWall==0,1,0))
mega4<-transform(mega4,closestBreak=pmin(preForest,Forest,Wall,na.rm=TRUE))
mega4$AnteBreak<-as.factor(ifelse(mega4$Dist<mega4$closestBreak|mega4$closestBreak %in% NA,1,0))
mega4$PostBreak<-as.factor(ifelse(mega4$AnteBreak==0,1,0))
mega4$EffectNo0<-mega4$Effect+.05
mega4$logDistPpF<-log(.05+ifelse(mega4$Dist<mega4$preForest|mega4$preForest %in% NA,0,mega4$Dist-mega4$preForest))
mega4$logDistPF<-log(.05+ifelse(mega4$Dist>=mega4$Forest,mega4$Dist-mega4$Forest,NA))
mega4$logDistPW<-log(.05+ifelse(mega4$Dist>=mega4$Wall,mega4$Dist-mega4$Wall,NA))
mega4$logDist2<-log(.05+ifelse(mega4$Dist<mega4$preForest|mega4$preForest %in% NA,mega4$Dist,mega4$Dist-mega4$preForest))

## Create subset dataframes
mega4preForest<-droplevels(rbind(subset(mega4,Dist<preForest),subset(mega4,is.na(preForest))))
mega4Forest<-droplevels(subset(mega4,Dist>=Forest))
mega4Wall<-droplevels(subset(mega4,Dist>=Wall))
mega4preWall<-droplevels(rbind(subset(mega4,Dist<Wall),subset(mega4,is.na(Wall))))


##### Find SiteBank data subsets with positive (non-decaying) curves #####

## Make Regs.All function
Regs.All<-function(Data){
	Regs.All.data<-list()
	nas<-as.numeric(rep(NA,length(unique(Data$SiteBank))))
	for(i in 1:length(unique(Data$SiteBank))){
		subs<-subset(Data,SiteBank==unique(Data$SiteBank)[i])
		if(dim(subs)[1]>2){
			subs.Lin<-lm((EffectNo0)~Dist,data=subs)
				results.Lin<-round(c(trueAIC(subs.Lin)[1],summary(subs.Lin)$adj.r.squared,summary(subs.Lin)$coefficients[2,3],summary(subs.Lin)$coefficients[2,4]),4)
			subs.Exp<-lm(logEffect~Dist,data=subs)
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
	
## Using preForest data
	## Toss if any of lin, exp, or pow show a positive curve
Regs.All.preForest<-Regs.All(mega4preForest)
	Regs.Keep.preForest<-Regs.All.preForest[sapply(Regs.All.preForest, function(x) max(x$Dist.t.val)<0)]
	Regs.Toss.preForest<-Regs.All.preForest[sapply(Regs.All.preForest, function(x) max(x$Dist.t.val)>0)]	
	mega5<-mega5preForest<-mega4preForest[mega4preForest$SiteBank %in% names(Regs.Keep.preForest),]
	
## Make table of tossed data
toss<-data.frame(levels(mega4preForest$OrderClass))
	colnames(toss)<-'OrderClass'
	toss$Percent.not.fit<-toss$Number.not.fit<-toss$Number.of.sites<-0
for(i in 1:dim(toss)[1]){
	toss$Number.of.sites[i]<-length(unique(mega4preForest[mega4preForest$OrderClass==toss$OrderClass[i],'SiteBank']))
	toss$Number.not.fit[i]<-length(unique(mega4preForest[mega4preForest$OrderClass==toss$OrderClass[i]&mega4preForest$SiteBank%in%names(Regs.Toss.preForest),'SiteBank']))
	toss$Percent.not.fit[i]<-round(length(unique(mega4preForest[mega4preForest$OrderClass==toss$OrderClass[i]&mega4preForest$SiteBank%in%names(Regs.Toss.preForest),'SiteBank']))/toss$Number.of.sites[i],4)
	}
toss.total<-c('Total',as.numeric(sum(toss$Number.of.sites)),as.numeric(sum(toss$Number.not.fit)),round(sum(toss$Number.not.fit)/sum(toss$Number.of.sites),4))
	toss$OrderClass<-as.character(toss$OrderClass)
toss<-rbind(toss,toss.total)
toss$Percent.of.all.non.fit.sites<-round(as.numeric(toss$Number.not.fit)/as.numeric(toss.total[3]),4)
#write.csv(toss,'Figures & Tables/R CSV Exports/IsotopeNonFits.csv',row.names=F)


##### Get stream signatures for each site subset #####

## Make sitesigs function
sitesigs<-function(Data){
	sitesigs.data<-data.frame(matrix(,nrow=1,ncol=2))
	for(i in 1:length(unique(Data$SiteBank))){
		sitesubs<-subset(Data,SiteBank==unique(Data$SiteBank)[i])
		site.Pow<-lm(logEffect~logDist,data=sitesubs)
			sig50<-round(exp((log(.5)-site.Pow$coefficients[1])/(site.Pow$coefficients[2])),2)
			sig10<-round(exp((log(.1)-site.Pow$coefficients[1])/(site.Pow$coefficients[2])),2)
		sitesigs.data[i,]<-cbind(sig50,sig10)
		}
		rownames(sitesigs.data)<-unique(Data$SiteBank)
		colnames(sitesigs.data)<-c('50%','10%')
		sitesigs.data<-sitesigs.data[order(rownames(sitesigs.data)),]	
	sitesigs.data
	}

## For preForest data
sitesigs.preForest<-sitesigs(mega5preForest)


##### Test model fits for data subsets #####

## Subset data by stream, bank
mega5subs<-mega5
	#mega5subs$subs<-mega5subs$SiteBank
	mega5subs$subs<-paste(mega5subs$SiteBank,mega5subs$Group)
subs.list<-list()
bads<-as.numeric(NA)
for(i in 1:length(unique(mega5subs$subs))){
	subs<-subset(mega5subs,subs==unique(mega5subs$subs)[i])
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
	names(subs.list)<-unique(mega5subs$subs)
	subs.list<-subs.list[is.na(bads)]	
	
## See how often Inv power is the best, vs. neg. exp and linear
best1<-vector()
for(i in 1:length(subs.list)){best1[i]<-order(subs.list[[i]]$AIC)[1]}
best2<-ifelse(best1==1,'Linear',ifelse(best1==2,'Neg. Exp','Inv. Power'))
best<-cbind(table(best2),round(table(best2)/length(subs.list),4))
	colnames(best)<-c('Freq','Percent')
	## For SiteBank, 6/23 (26.09%) pow, 16/23 (69.57%) exp, 1/23 (4.35%) lin.
	## For SiteBankGroup, 15/41 (36.59%) pow, 26/41 (63.41%) exp.	

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
lm.lin<-lm((Effect+.05)~Dist,data=mega5)
lm.exp<-lm(logEffect~Dist,data=mega5)
lm.pow<-lm(logEffect~logDist,data=mega5)
#rbind(AIC(lm.lin),trueAIC(lm.exp)[1],trueAIC(lm.pow)[1]);AIC(lm.exp,lm.pow)
	## Pow ~15 AIC point improvement over exp, and ~55 improvement over linear. R2s ~ 25%.


##### Try mixed effects models on data #####

## Use Site as the structural variable
lme.lin<-lme((Effect+.05)~Dist,random=~1|shortSite,data=mega5,method='ML')
lme.exp<-lme(logEffect~Dist,random=~1|shortSite,data=mega5,method='ML')
lme.pow<-lme(logEffect~logDist,random=~1|shortSite,data=mega5,method='ML')
#rbind(AIC(lme.lin),trueAIC(lme.exp)[1],trueAIC(lme.pow)[1]);AIC(lme.exp,lme.pow)
	## Pow still clearly the best.
	
## Use Site+Bank as the structural variable
lmeSB.lin<-lme((Effect+.05)~Dist,random=~1|SiteBank,data=mega5,method='ML')
lmeSB.exp<-lme(logEffect~Dist,random=~1|SiteBank,data=mega5,method='ML')
lmeSB.pow<-lme(logEffect~logDist,random=~1|SiteBank,data=mega5,method='ML')
#rbind(AIC(lmeSB.lin),trueAIC(lmeSB.exp)[1],trueAIC(lmeSB.pow)[1]);AIC(lmeSB.exp,lmeSB.pow)
	## Probably the most sensible to use SiteBank because it is the truer experimental unit (and bank type is of interest). Also improves AIC of pow model a bit. Power substantially better than exp.	
aiccs<-as.data.frame(round(rbind(AICc(lmeSB.lin),trueAIC(lmeSB.exp)[2],trueAIC(lmeSB.pow)[2]),2))
	colnames(aiccs)<-'AICc'
	rownames(aiccs)<-c('Linear','Exponential','Power')
	aiccs$delta<-round(c(aiccs$AICc-min(aiccs$AICc)),2)
	aiccs$weight<-round(exp(-.5*aiccs$delta)/sum(exp(-.5*aiccs$delta)),4)
	aiccs<-aiccs[order(-aiccs$AICc),]
#write.csv(aiccs,'Figures & Tables/R CSV Exports/IsotopeAICs.csv')


##### Compare more sophisticated models using all data (not just preForest) #####

## Fixed effects only
lm.lin0<-lm((EffectNo0)~Dist,data=mega4)
lm.exp0<-lm(logEffect~Dist,data=mega4)
lm.pow0<-lm(logEffect~logDist,data=mega4)
lm.pow1<-lm(logEffect~0+AntepreForest+PostpreForest+AntepreForest:logDist+PostpreForest:logDistPpF,data=mega4);summary(lm.pow1);AIC(lm.pow1)
lm.pow2<-lm(logEffect~PostpreForest+logDist+PostpreForest:logDistPpF,data=mega4);summary(lm.pow2);AIC(lm.pow1)

## Random intercepts
lmeSB.lin0<-lmer(EffectNo0~Dist + (1|SiteBank),data=mega4,REML=FALSE)
lmeSB.exp0<-lmer(logEffect~Dist + (1|SiteBank),data=mega4,REML=FALSE)
lmeSB.pow0<-lmer(logEffect~logDist,random=~1|SiteBank,data=mega4,method='ML');AIC(lmeSB.pow0)
lmeSB.pow1<-lmer(logEffect~0+AntepreForest+PostpreForest+AntepreForest:logDist+PostpreForest:logDistPpF +(1|SiteBank),data=mega4,REML=FALSE);AIC(lmeSB.pow1)
lmeSB.pow2<-lmer(logEffect~PostpreForest+logDist+PostpreForest:logDistPpF+(1|SiteBank),data=mega4,REML=FALSE);AIC(lmeSB.pow2)
mega4$PredLin0<-coef(lmeSB.lin0)[[1]][1,2]*mega4$Dist+coef(lmeSB.lin0)[[1]][mega4$SiteBank,1]

##### Compute single variable effects and stream signatures #####

## Make list of environmental variables of interest
vars<-c('Region','Width','Order','OrderClass','Geomorph','Banks','Group','VegFld','VegUp','VegShift')

## Run signature function
sigs.preForest<-signatures(mega5preForest,'logEffect','logDist','SiteBank',vars)

## Write files for single model and signature results	
w1<-function(){
write.csv(sigs.preForest$fits,'Figures & Tables/R CSV Exports/IsotopeSinglesModels.csv')
out_file<-file('Figures & Tables/R CSV Exports/IsotopeSignatures.csv', open='a')
	for(i in seq_along(sigs.preForest$sigs)){
		write.table(names(sigs.preForest$sigs)[i],file=out_file,sep=',',dec='.',quote=FALSE,col.names=FALSE,row.names=FALSE)  
		write.table(sigs.preForest$sigs[[i]],file=out_file,sep=',',dec='.',quote=FALSE,col.names=NA,row.names=TRUE)
		}
	close(out_file)
}
#w1()


##### Build overall multivariate model #####

## Decide on whether it is better to add stream order as Order, Order class, or Channel width (just using power model for ease)
lmeWid.pow<-lme(logEffect~logDist+Width,random=~1|SiteBank,data=mega5,method='ML')
lmeOrd.pow<-lme(logEffect~logDist+Order,random=~1|SiteBank,data=mega5,method='ML')
lmeOC.pow<-lme(logEffect~logDist+OrderClass,random=~1|SiteBank,data=mega5,method='ML')
#AIC(lmeWid.pow,lmeOrd.pow,lmeOC.pow)
	## OC beats Order by ~7, and ~1 worse than Width. Go with Order Class also because interpretation is easier and results from singles.lme make more sense.

## Add variables en masse using dredge

## Start with additive terms only
	## Can't do all at once due to singularities. So have to do a couple iterations
runmodsp<-function(){
modtablep1<-dredge(lme(logEffect~logDist+Group+Geomorph+Banks+VegFld+OrderClass,random=~1|SiteBank,data=mega5,method='ML'),fixed='logDist',extra='adjR^2');subset(modtablep1,delta<4)
modtablep2<-dredge(lme(logEffect~logDist+Group+Geomorph+Banks+VegFld+Region,random=~1|SiteBank,data=mega5,method='ML'),fixed='logDist',extra='adjR^2');subset(modtablep2,delta<4)
modtablep1$Region<-modtablep2$OrderClass<-NA
modtablep<-rbind(modtablep1,modtablep2[!is.na(modtablep2$Region),])
	modtablep$Region<-as.factor(modtablep$Region)
	modtablep<-modtablep[order(modtablep$delta),c(1:6,14,7:13)]
	m1<-exp(-.5*modtablep$delta)
	modtablep$weight<-m1/sum(m1)
print(subset(modtablep,delta<4))
assign('modtablep',modtablep,envir=.GlobalEnv)
}
#runmodsp()
	## Focusing on VegFld instead of VegShift or VegUp because it seems to improve these models the most and adding all of them is redundant, plus its results from singles.lme are the most logical.
	## Also, VegShift can't be added with VegUp or VegFld due to singularities.

## Now try interaction terms
runmodsx<-function(){
modtablex1<-dredge(lme(logEffect~logDist+logDist*Group+logDist*Geomorph+logDist*Banks+logDist*VegFld+logDist*OrderClass,random=~1|SiteBank,data=mega5,method='ML'),fixed='logDist',extra='adjR^2');subset(modtablex1,delta<4)
modtablex2<-dredge(lme(logEffect~logDist+logDist*Group+logDist*Geomorph+logDist*Banks+logDist*VegFld+logDist*Region,random=~1|SiteBank,data=mega5,method='ML'),fixed='logDist',extra='adjR^2');subset(modtablex2,delta<4)
modtablex1$logDistxRegion<-modtablex1$Region<-modtablex2$logDistxOrderClass<-modtablex2$OrderClass<-NA
	colnames(modtablex1)[11]<-'logDistxOrderClass'
	colnames(modtablex2)[11]<-'logDistxRegion'
modtablex<-rbind(modtablex1,modtablex2[!is.na(modtablex2$Region),])
	modtablex$Region<-as.factor(modtablex$Region)
	modtablex$logDistxRegion<-as.factor(modtablex$logDistxRegion)
	modtablex<-modtablex[order(modtablex$delta),c(1:6,19,7:11,20,12:18)]
	m1<-exp(-.5*modtablex$delta)
	modtablex$weight<-m1/sum(m1)
print(subset(modtablex,delta<4))
assign('modtablex',modtablex,envir=.GlobalEnv)
}
#runmodsx()
	## Best models are still only additive, but this does add a few new ones to the fold.

## Export models
	## When subsetting best models, use delta<4.  This means the difference between AIC(i) and AIC(min) is <4.  "Substantial model support" according to Burnham & Anderson (2002) pp70.
#write.csv(modtablex,'Figures & Tables/R CSV Exports/IsotopeModels.csv')	
	

##### Compare predicted vs. observed values at breaks #####

## Add predicted values to mods dataframe
mods$preForest<-breaks[match(mods$SiteBank,breaks$SiteBank),'isopreForest']
mods$Forest<-breaks[match(mods$SiteBank,breaks$SiteBank),'isoForest']
mods$Wall<-breaks[match(mods$SiteBank,breaks$SiteBank),'isoWall']
mods$Wall.Exp<-mods$Forest.Exp<-as.numeric(rep(NA,dim(mods)[1]))
for(i in 1:dim(mods)[1]){
	subsite.lm<-lm(logEffect~logDist,data=mega4[mega4$SiteBank==mods$SiteBank[i],])
	mods$Forest.Exp[i]<-exp(subsite.lm$coefficients[1]+subsite.lm$coefficients[2]*log(mods$preForest[i]+.05))
	mods$Wall.Exp[i]<-exp(subsite.lm$coefficients[1]+subsite.lm$coefficients[2]*log(mods$Wall[i]+.05))	
	}
forestbreaks<-mega4[!is.na(mega4$preForest)&mega4$Dist==mega4$preForest&mega4$SiteBank%in%mega5$SiteBank,c('SiteBank','Effect')]
	forestmeans<-data.frame(tapply(forestbreaks$Effect,forestbreaks$SiteBank,mean))	
	mods$Forest.Obs<-forestmeans[match(mods$SiteBank,rownames(forestmeans)),1]
wallbreaks<-mega4[!is.na(mega4$Wall)&mega4$Dist==mega4$Wall&mega4$SiteBank%in%mega5$SiteBank,c('SiteBank','Effect')]
	wallmeans<-data.frame(tapply(wallbreaks$Effect,wallbreaks$SiteBank,mean))	
	mods$Wall.Obs<-wallmeans[match(mods$SiteBank,rownames(wallmeans)),1]	

## Run model to see if observed PercAqRel are higher at breaks than predicted
tFor<-t.test(mods$Forest.Obs,mods$Forest.Exp,alternative='greater',paired=T,var.equal=F,na.action='na.omit')
tWall<-t.test(mods$Wall.Obs,mods$Wall.Exp,alternative='greater',paired=T,var.equal=F,na.action='na.omit')
	## Both significant (t=2.4028,2.762, p=.01857,.02537 for forest and wall, respectively).
	## For forest, a 16.26 PercAq increase from expected to observed. 20.00 for wall.
tFor<-t.test(mods$Forest.Obs,mods$Forest.Exp,alternative='greater',paired=T,var.equal=F,na.action='na.omit')
tWall<-t.test(mods$Wall.Obs,mods$Wall.Exp,alternative='greater',paired=T,var.equal=F,na.action='na.omit')
diffFor<-round(mean(mods$Forest.Obs,na.rm=T)-mean(mods$Forest.Exp,na.rm=T),4)
diffWall<-round(mean(mods$Wall.Obs,na.rm=T)-mean(mods$Wall.Exp,na.rm=T),4)
incFor<-round(100*diffFor/mean(mods$Forest.Exp,na.rm=T),2)
incWall<-round(100*diffWall/mean(mods$Wall.Exp,na.rm=T),2)
	tabFor<-cbind(diffFor,incFor,round(tFor$statistic,4),round(tFor$p.value,4))
	tabWall<-cbind(diffWall,incWall,round(tWall$statistic,4),round(tWall$p.value,4))
		tabForWall<-rbind(tabFor,tabWall)
			rownames(tabForWall)<-c('Forest','Wall')
			colnames(tabForWall)<-c('PercAqDiff','RelPercChange','t','p')
#write.csv(tabForWall,'Figures & Tables/R CSV Exports/IsotopeForestWallObsVsExp.csv')
	## Both significant, with large increases in observed vs. model expected.	
##### Make islands-only dataframe #####

## Get islands only
islands1<-droplevels(mega1[nchar(as.character(mega1$Site))>5,])

## Rearrange alphabetically for ease
islands2<-droplevels(islands1[order(as.character(islands1$Code)),])

## Make Percent Aquatic effect
islands2$PercAq<-ifelse(breaks[match(islands2$SiteBank,breaks$SiteBank),'CvsN']=='N',islands2$PercNaq,islands2$PercCaq)
islands3<-droplevels(islands2[!is.na(islands2$PercAq),])

## Add cross section distances to the dataframe
m3<-match(paste(islands3$Site,islands3$Bank,islands3$Dist,sep=''),paste(islandsXS$Site,islandsXS$Bank,islandsXS$Dist,sep=''))
islands3$DistXS<-islandsXS$DistXS[m3]
islands3$BankIsland<-as.factor(ifelse(paste(islands3$Site,islands3$Bank)=='CORN1.1 RB'|paste(islands3$Site,islands3$Bank)=='CORN1.6 LB'|paste(islands3$Site,islands3$Bank)=='FELL1.1 RB'|paste(islands3$Site,islands3$Bank)=='ISLA1.1 RB'|paste(islands3$Site,islands3$Bank)=='ISLA1.7 LB'|paste(islands3$Site,islands3$Bank)=='RESI1.1 RB'|paste(islands3$Site,islands3$Bank)=='RESI1.2 LB'|paste(islands3$Site,islands3$Bank)=='RESI5.1 RB'|paste(islands3$Site,islands3$Bank)=='RESI5.2 LB'|paste(islands3$Site,islands3$Bank)=='ELBE1.2 LB','Bank','Island'))

## Relativize PercAq
mods.islands<-data.frame(unique(islands3$shortSite))
	colnames(mods.islands)<-'shortSite'
	mods.islands$maxs<-mods.islands$mins<-c(rep(0,length(unique(islands3$shortSite))))
for(i in 1:dim(mods.islands)[1]){
	mods.islands$mins[i]<-min(islands3[islands3$shortSite==mods.islands$shortSite[i],'PercAq'])
	mods.islands$maxs[i]<-max(islands3[islands3$shortSite==mods.islands$shortSite[i],'PercAq'])
	}
islands3$PercAqRel<-(islands3$PercAq-mods.islands[match(islands3$shortSite,mods.islands$shortSite),'mins'])/(mods.islands[match(islands3$shortSite,mods.islands$shortSite),'maxs']-mods.islands[match(islands3$shortSite,mods.islands$shortSite),'mins'])
	
## Run model to see if PercAq is higher on islands than banks
islands.mod<-lme(PercAqRel~BankIsland,random=~1|shortSite,data=islands3)
	## Yes. Looking at intercept, Islands 16.08 AqPerc higher than banks (which are 49.07%). t = 4.192481, p = 0.00002 (p-value is half the reported value, as what I really want is a 1-sided test). 

## Plot
isla<-islands3[islands3$shortSite=='ISLA1'&islands3$PercAq,]
corn<-islands3[islands3$shortSite=='CORN1'&islands3$PercAq,]
fell<-islands3[islands3$shortSite=='FELL1',]
resi1<-islands3[islands3$shortSite=='RESI1'&islands3$PercAq,]
resi5<-islands3[islands3$shortSite=='RESI5',]
elbe<-islands3[islands3$shortSite=='ELBE1'&islands3$PercAq,]
danu<-islands3[islands3$shortSite=='DANU1',]
#plot(isla$DistXS,isla$PercAq);plot(corn$DistXS,corn$PercAq);plot(fell$DistXS,fell$PercAq);plot(resi1$DistXS,resi1$PercAq);plot(resi5$DistXS,resi5$PercAq);plot(elbe$DistXS,elbe$PercAq);plot(danu$DistXS,danu$PercAq)

## Make pdf of CORN1 as an example
cornXS<-islandsXS[substr(islandsXS$Site,1,5)=='CORN1'&!is.na(islandsXS$Height),]
f5<-function(Type='pdf'){
if(Type=='png'){png('Figures & Tables/Isotope Braids.png',width=6.5,height=8,res=300,bg=0,units='in')}
else{pdf('Figures & Tables/Isotope Braids.pdf',width=6.5,height=8)}
oldpar2<-par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(1,3.6,0.5,0.5))
plot(c(-5,365),c(0,1.03),type='n',xaxs='i',yaxs='i',xlab='',ylab='',xaxt='n',yaxt='n',bty='l')
	axis(1,at=seq(0,360,60),labels=F,lwd=0,lwd.ticks=1)
	axis(2,las=2,lwd=0,lwd.ticks=1)
		title(ylab='Portion aquatic',line=2.5)
	points(corn$DistXS,corn$PercAq,pch=as.numeric(corn$Group))
	lines(lowess(corn$DistXS,corn$PercAq,f=.5))
	abline(v=c(15,335),lty=2,xpd=TRUE)
	legend(220,.215,legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3),bty='n')
	legend(142,.1,legend='CORN1',bty='n')
	legend(-20,1.1,legend='A',cex=2,bty='n')
par(mar=c(3.8,3.6,.1,.5))
plot(c(-5,365),c(min(cornXS$Height),max(cornXS$Height)),type='n',xaxs='i',yaxs='i',xlab='',ylab='',xaxt='n',yaxt='n',bty='l')
	axis(1,at=seq(0,360,60),lwd=0,lwd.ticks=1)
	axis(2,las=2,lwd=0,lwd.ticks=1)
		title(xlab='Distance from left bank-most point (m)',line=2.5)
		title(ylab='Height above water level (m)',line=2.5)
	#polygon(c(3600,cornXS$DistXS,0),c(-1,cornXS$Height,-1),col='gray80')
	abline(v=c(15,335),lty=2,xpd=FALSE)
	lines(cornXS$DistXS,cornXS$Height)
	abline(h=0,lty=3,lwd=2)
	legend(-20,3.8,legend='B',cex=2,bty='n')
par(oldpar2)
dev.off()
}
#f5();f5(Type='png')	


##### Get percent of watershed area encompassed by these stream signatures #####

## Load libraries
library(sp)
library(rgdal)

## Read in shapefiles
basin<-readOGR('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis/Site Data, Physical Data & Maps/GIS Data & Maps/Coweeta Map Files/StreamSigs','coweeta_basin',verbose=FALSE)
pred10<-readOGR('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis/Site Data, Physical Data & Maps/GIS Data & Maps/Coweeta Map Files/StreamSigs/Buffers','BuffPred10',verbose=FALSE)
#pred50<-readOGR('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis/Site Data, Physical Data & Maps/GIS Data & Maps/Coweeta Map Files/StreamSigs/Buffers','BuffPred50',verbose=FALSE)
	## Note: This one doesn't work because all values (stream signatures) are 0.
iso10<-readOGR('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis/Site Data, Physical Data & Maps/GIS Data & Maps/Coweeta Map Files/StreamSigs/Buffers','BuffIso10',verbose=FALSE)
iso50<-readOGR('C:/Users/jmuehlbauer/Documents/Projects/Stream Signature/Stream Signature Analysis/Site Data, Physical Data & Maps/GIS Data & Maps/Coweeta Map Files/StreamSigs/Buffers','BuffIso50',verbose=FALSE)

## Create data frame
GIS<-as.data.frame(rep(c('Predator','Isotope'),c(2,2)))
	colnames(GIS)<-'Metric'
	GIS$Signature<-rep(c('50%','10%'),2)
## Get areas
GIS$Area<-c(0,pred10$Shape_Area,iso50$Shape_Area,iso10$Shape_Area)
GIS$PercentBasinArea<-round(GIS$Area/basin[[1]],4)

## Export
#write.csv(GIS,'Figures & Tables/R CSV Exports/GIS Percent of Basin Area.csv',row.names=F)