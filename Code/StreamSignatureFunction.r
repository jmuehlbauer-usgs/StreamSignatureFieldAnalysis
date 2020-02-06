signatures<-function(Data,EffectVar,DistanceVar,SiteVar,varnames){

## Load required libraries
library(MuMIn)
library(nlme)

## Make variables
perc<-Data[,EffectVar]
dstn<-Data[,DistanceVar]
strm<-Data[,SiteVar]

## Add predictors to intercept and slope & intercept
lme.Stream<-lme(perc~dstn,random=~1|strm,method='ML')
singles.lme.list<-list()
for(i in 1:length(varnames)){
	variable<-Data[,varnames[i]]
	singles.lme.list[[i]]<-lme(perc~dstn+variable,random=~1|strm,method='ML')
	singles.lme.list[[i+length(varnames)]]<-lme(perc~dstn*variable,random=~1|strm,method='ML')
	}
	names(singles.lme.list)<-c(paste('p',varnames,sep=''),paste('x',varnames,sep=''))
singles.lme<-round(as.data.frame(sapply(c(list(lme.Stream),singles.lme.list),function(x) ((as.numeric(VarCorr(lme.Stream)[1,1]) - as.numeric(VarCorr(x)[1,1])) / as.numeric(VarCorr(lme.Stream)[1,1])))),4)
	rownames(singles.lme)<-c('Distance',names(singles.lme.list))
	colnames(singles.lme)<-'Pseudo R^2'
	singles.lme$LikRatio<-c('NA',round(sapply(singles.lme.list,function(x) anova(lme.Stream,x)[2,8]),4))
	singles.lme$p<-c('NA',round(sapply(singles.lme.list,function(x) anova(lme.Stream,x)[2,9]),4))
		## Anova p is from likelihood ratio test vs. the distance only model.
	singles.lme$DF<-c(4,round(sapply(singles.lme.list,function(x) anova(lme.Stream,x)[2,3]),4))
	singles.lme$AICc<-round(sapply(c(list(lme.Stream),singles.lme.list),function(x) AICc(x)),2)
	singles.lme$delta<-c(singles.lme$AICc-min(singles.lme$AICc))
singles.lme$weight<-round(exp(-.5*singles.lme$delta)/sum(exp(-.5*singles.lme$delta)),4)
	singles.lme<-singles.lme[order(-singles.lme$AICc),]


##### Get 50% and 10% distances and confidence intervals from singles.lme results #####

## Subset only the variables of interest
lme.AICcs<-as.data.frame(1,nrow=1,ncol=1)
	for(i in 1:length(singles.lme.list)){
		lme.AICcs[i,1]<-AICc(singles.lme.list[[i]])}
	colnames(lme.AICcs)<-c('AICc')
	rownames(lme.AICcs)<-names(singles.lme.list)
	lme.AICcs$Keepers<-ifelse(lme.AICcs$AIC<AICc(lme.Stream),'Keep','Toss')
sub.vars<-rownames(subset(lme.AICcs,Keepers=='Keep'))
sub.singles.list1<-list()
	sub.singles.list1$Overall<-lme.Stream
sub.singles.list<-c(sub.singles.list1,singles.lme.list[sub.vars])

## Get confidence intervals
sub.singles.CI<-list()
for(i in 1:length(sub.singles.list)){
sub.singles.CI[[i]]<-intervals(sub.singles.list[[i]],level=.5)
	}
names(sub.singles.CI)<-names(sub.singles.list)

## Compute signatures for these variable conditions
percents<-c(.5,.1)
addon<-c(0,3)
sigs<-list()
for(i in 1:length(sub.singles.CI)){
	sigs[[i]]<-as.data.frame(list())
		for(k in 1:2){for (l in 1:3){
			sigs[[i]][1,addon[k]+l]<-exp((log(percents[k])-sub.singles.CI[[i]]$fixed[1,'est.'])/(sub.singles.CI[[i]]$fixed[2,l]))
			}}
		if(substr(names(sub.singles.CI)[i],1,1)=='x'){
			if(is.numeric(Data[,substr(names(sub.singles.CI)[[i]],2,100)])){
				for(k in 1:2){for (l in 1:3){
					sigs[[i]][1,addon[k]+l]<-exp((log(percents[k])-sub.singles.CI[[i]]$fixed[1,'est.']-sub.singles.CI[[i]]$fixed[3,'est.'])/(sub.singles.CI[[i]]$fixed[2,'est.']+sub.singles.CI[[i]]$fixed[4,l]))
					}}}			
			else{
				for(j in 1:(dim(sub.singles.CI[[i]]$fixed)[1]/2-1)){for(k in 1:2){for (l in 1:3){
					sigs[[i]][j+1,addon[k]+l]<-exp((log(percents[k])-sub.singles.CI[[i]]$fixed[1,'est.']-sub.singles.CI[[i]]$fixed[j+2,'est.'])/(sub.singles.CI[[i]]$fixed[2,'est.']+sub.singles.CI[[i]]$fixed[j+2+(length(sub.singles.CI[[i]]$fixed[,l])-2)/2,l]))
					}}}}}
		if(substr(names(sub.singles.CI)[i],1,1)=='p'){
			if(is.numeric(Data[,substr(names(sub.singles.CI)[[i]],2,100)])){
				for(k in 1:2){for (l in 1:3){
					sigs[[i]][1,addon[k]+l]<-exp((log(percents[k])-sub.singles.CI[[i]]$fixed[1,'est.']-sub.singles.CI[[i]]$fixed[3,'est.'])/sub.singles.CI[[i]]$fixed[2,l])
					}}}
			else{
				for(j in 1:(dim(sub.singles.CI[[i]]$fixed)[1]-2)){for(k in 1:2){for (l in 1:3){			
					sigs[[i]][j+1,addon[k]+l]<-exp((log(percents[k])-sub.singles.CI[[i]]$fixed[1,'est.']-sub.singles.CI[[i]]$fixed[j+2,'est.'])/sub.singles.CI[[i]]$fixed[2,l])
					}}}}}
		colnames(sigs[[i]])<-c('Lower CI 50%','50%','Upper CI 50%','Lower CI 10%','10%','Upper CI 10%')
		if(names(sub.singles.CI)[i]=='Overall'||names(sub.singles.CI)[i]=='pWidth'||names(sub.singles.CI)[i]=='xWidth'){rownames(sigs[[i]])<-''}
		else{rownames(sigs[[i]])<-levels(droplevels(Data[,substr(names(sub.singles.CI),2,100)[i]]))}
	sigs[[i]]<-signif(round(sigs[[i]],2))
		}
	names(sigs)<-c('Overall',names(sub.singles.CI[-1]))
final<-list()
	final$fits<-singles.lme
	final$sigs<-sigs
final
}