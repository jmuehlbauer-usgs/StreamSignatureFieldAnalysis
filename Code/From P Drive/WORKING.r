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
mega4$Dist2<-ifelse(mega4$Dist<mega4$preForest|mega4$preForest %in% NA,mega4$Dist,mega4$Dist-mega4$preForest)
mega4$logDist2<-log(.05+mega4$Dist2)

## Create subset dataframes
mega4preForest<-droplevels(rbind(subset(mega4,Dist<preForest),subset(mega4,is.na(preForest))))
mega4Forest<-droplevels(subset(mega4,Dist>=Forest))
mega4Wall<-droplevels(subset(mega4,Dist>=Wall))
mega4preWall<-droplevels(rbind(subset(mega4,Dist<Wall),subset(mega4,is.na(Wall))))


## Fixed effects only
lm.lin0<-lm((EffectNo0)~Dist,data=mega4);AIC(lm.lin0)
lm.exp0<-lm(logEffect~Dist,data=mega4);AIC(lm.exp0)
lm.pow0<-lm(logEffect~logDist,data=mega4);AIC(lm.pow0)
lm.pow1<-lm(logEffect~0+AntepreForest+PostpreForest+AntepreForest:logDist+PostpreForest:logDistPpF,data=mega4);AIC(lm.pow1)
lm.pow2<-lm(logEffect~PostpreForest+logDist+PostpreForest:logDistPpF,data=mega4);AIC(lm.pow2)
round(sapply(list(lm.lin0,lm.exp0,lm.pow0,lm.pow1,lm.pow2),AIC),1)

lm.d2lin0<-lm((EffectNo0)~Dist2,data=mega4);AIC(lm.d2lin0)
lm.d2exp0<-lm(logEffect~Dist2,data=mega4);AIC(lm.d2exp0)
lm.d2pow0<-lm(logEffect~logDist2,data=mega4);AIC(lm.d2pow0)
lm.d2pow1<-lm(logEffect~0+AntepreForest+PostpreForest+AntepreForest:logDist2+PostpreForest:logDistPpF,data=mega4);AIC(lm.d2pow1)
lm.d2pow2<-lm(logEffect~PostpreForest+logDist2+PostpreForest:logDistPpF,data=mega4);AIC(lm.d2pow2)
round(sapply(list(lm.lin0,lm.exp0,lm.pow0,lm.pow1,lm.pow2),AIC),1)

## Random intercepts
lmeSB.lin0<-lmer(EffectNo0~Dist + (1|SiteBank),data=mega4,REML=FALSE)
lmeSB.exp0<-lmer(logEffect~Dist + (1|SiteBank),data=mega4,REML=FALSE)
lmeSB.pow0<-lmer(logEffect~logDist+(1|SiteBank),data=mega4,REML=FALSE)
lmeSB.pow1<-lmer(logEffect~0+AntepreForest+PostpreForest+AntepreForest:logDist+PostpreForest:logDistPpF +(1|SiteBank),data=mega4,REML=FALSE)
lmeSB.pow2<-lmer(logEffect~PostpreForest+logDist+PostpreForest:logDistPpF+(1|SiteBank),data=mega4,REML=FALSE)
round(sapply(list(lmeSB.lin0,lmeSB.exp0,lmeSB.pow0,lmeSB.pow1,lmeSB.pow2),AIC),1)



mega4$PredLin0<-coef(lmeSB.lin0)[[1]][1,2]*mega4$Dist+coef(lmeSB.lin0)[[1]][mega4$SiteBank,1]
