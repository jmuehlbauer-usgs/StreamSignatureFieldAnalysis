PercAq<-function(Data,SiteName){
##### Subset data by sites for analysis #####
site<-Data[substr(Data$Site,1,5)==SiteName,]


##### Adjust to trophic baseline #####
## Get average d13C and d15N values for all groups
avgsN<-tapply(site$d15N,site$Group,mean)
avgsC<-tapply(site$d13C,site$Group,mean)

## Adjust all data based on d13C and d15N fractionation between hemiptera herbivores and terrestrial detritus
site$d15Nadj<-site$d15N-site$Trophic*(avgsN['HeHemi']-avgsN['TeDetr'])
site$d13Cadj<-site$d13C-site$Trophic*(avgsC['HeHemi']-avgsC['TeDetr'])
## If no hemipterans are available for comparison, use Post (2002) values of 3.4 and .4 per mil.
if(is.na(avgsC['HeHemi'])){
	site$d15Nadj<-site$d15N-site$Trophic*(3.4)
	site$d13Cadj<-site$d13C-site$Trophic*(.4)
	}
avgsNnew<-tapply(site$d15Nadj,site$Group,mean)
avgsCnew<-tapply(site$d13Cadj,site$Group,mean)

	
##### Calculate Percent Aquatic #####
## Calculate using terrestrial detritus and aquatic algae as baselines
site$PercNaqVeg<-round((site$d15Nadj-avgsNnew['TeDetr'])/(avgsNnew['AqAlga']-avgsNnew['TeDetr']),4)
site$PercCaqVeg<-round((site$d13Cadj-avgsCnew['TeDetr'])/(avgsCnew['AqAlga']-avgsCnew['TeDetr']),4)

## If no aquatic algae is available, use aquatic detritus instead
if(is.na(avgsNnew['AqAlga'])){
	site$PercNaqVeg<-round((site$d15Nadj-avgsNnew['TeDetr'])/(avgsNnew['AqDetr']-avgsNnew['TeDetr']),4)
	site$PercCaqVeg<-round((site$d13Cadj-avgsCnew['TeDetr'])/(avgsCnew['AqDetr']-avgsCnew['TeDetr']),4)
	}

## Calculate using hemiptera herbivores and aquatic insects as baselines
site$PercNaqAn<-round((site$d15Nadj-avgsNnew['HeHemi'])/((sum(site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',]$d15Nadj*site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',]$Biomass)/sum(site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',]$Biomass))-avgsNnew['HeHemi']),4)
site$PercCaqAn<-round((site$d13Cadj-avgsCnew['HeHemi'])/((sum(site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',]$d13Cadj*site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',]$Biomass)/sum(site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',]$Biomass))-avgsCnew['HeHemi']),4)

## If no aquatic insects are available, use aquatic misc
if(dim(site[is.na(site$Dist)&site$Trophic>0&site$Group!='AqMisc',])[1]==0){
	site$PercNaqAn<-round((site$d15Nadj-avgsNnew['HeHemi'])/((sum(site[is.na(site$Dist)&site$Trophic>0,]$d15Nadj*site[is.na(site$Dist)&site$Trophic>0,]$Biomass)/sum(site[is.na(site$Dist)&site$Trophic>0,]$Biomass))-avgsNnew['HeHemi']),4)
	site$PercCaqAn<-round((site$d13Cadj-avgsCnew['HeHemi'])/((sum(site[is.na(site$Dist)&site$Trophic>0,]$d13Cadj*site[is.na(site$Dist)&site$Trophic>0,]$Biomass)/sum(site[is.na(site$Dist)&site$Trophic>0,]$Biomass))-avgsCnew['HeHemi']),4)
	}
	
##### Relativize data to between 0% and 100% #####
## Translate down or up, respectively, the top 5% and bottom 5% of data values to reduce outlier effect
## Shift values so smallest value is 0%
## Relativize so highest point is 100%
PercRed<-list()
Percgrps<-c('PercNaqVeg','PercCaqVeg','PercNaqAn','PercCaqAn')
PercgrpsRel<-paste(Percgrps,'Rel',sep='')
for(i in 1:4){
	PercRed[[i]]<-site[,Percgrps[i]]
		PercRed[[i]][c(order(site[,Percgrps[i]])[1:ceiling(dim(site)[1]*.05)])]<-PercRed[[i]][order(site[,Percgrps[i]])[ceiling(dim(site)[1]*.05)+1]]
		PercRed[[i]][c(order(-site[,Percgrps[i]])[1:ceiling(dim(site)[1]*.05)])]<-PercRed[[i]][order(-site[,Percgrps[i]])[ceiling(dim(site)[1]*.05)+1]]
	PercRed[[i+4]]<-PercRed[[i]]-min(PercRed[[i]])
	site[,PercgrpsRel[i]]<-round(PercRed[[i+4]]/max(PercRed[[i+4]]),4)
	}

## Flip data about y-axis (% aquatic) if necessary)
trop<-site[site$Trophic==2,]
if(length(na.omit(trop$PercCaqVegRel))>0){
	if(suppressWarnings(lm(log(PercCaqVegRel+.05)~log(Dist+.05),data=trop)$coefficients[2]>0)){
		site$PercCaqVegRel<-1-site$PercCaqVegRel}
	if(suppressWarnings(lm(log(PercNaqVegRel+.05)~log(Dist+.05),data=trop)$coefficients[2]>0)){
		site$PercNaqVegRel<-1-site$PercNaqVegRel}
	}
if(length(na.omit(trop$PercCaqAnRel))>0){
	if(suppressWarnings(lm(log(PercCaqAnRel+.05)~log(Dist+.05),data=trop)$coefficients[2]>0)){
		site$PercCaqAnRel<-1-site$PercCaqAnRel}
	if(suppressWarnings(lm(log(PercNaqAnRel+.05)~log(Dist+.05),data=trop)$coefficients[2]>0)){
		site$PercNaqAnRel<-1-site$PercNaqAnRel}
	}

site
}