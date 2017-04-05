library(lme4)
mega4<-read.csv('P:/BIOLOGICAL/FLyco/Jeff/mega4.csv')
dat<-as.data.frame(mega4[,c('SiteBank')])
	colnames(dat)<-'SiteBank'
dat$Break<-breaks[match(dat$SiteBank,breaks$SiteBank),'isopreForest']
dat$Pre<-ifelse(mega4$Dist<dat$Break|dat$Break %in% NA,1,0)
dat$Post<-ifelse(dat$Pre==0,1,0)
dat$factPre<-ifelse(mega4$Dist<dat$Break|dat$Break %in% NA,'Yes','No')
dat$factPost<-ifelse(dat$Pre=='No','Yes','No')
dat$Effect<-round(mega4$PercAqRel+.05,4)
dat$Dist<-mega4$Dist+.05
dat$DistB<-ifelse(dat$Dist<dat$Break|dat$Break %in% NA,NA,dat$Dist-dat$Break)
dat$DistBr<-ifelse(dat$Dist<dat$Break|dat$Break %in% NA,0,dat$Dist-dat$Break)
dat$DistC<-ifelse(dat$Dist<dat$Break|dat$Break %in% NA,dat$Dist,dat$Dist-dat$Break)
dat$logEffect<-round(log(dat$Effect),4)
dat$logDist<-round(log(dat$Dist),4)
dat$logDistB<-round(log(dat$DistB),4)
dat$logDistBr<-ifelse(dat$DistBr==0,0,round(log(dat$DistBr),4))
dat$logDistC<-round(log(dat$DistC),4)
write.csv(dat,'P:/Biological/Flyco/Jeff/dat.csv',row.names=FALSE)
dat<-read.csv('P:/BIOLOGICAL/FLyco/Jeff/dat.csv')


