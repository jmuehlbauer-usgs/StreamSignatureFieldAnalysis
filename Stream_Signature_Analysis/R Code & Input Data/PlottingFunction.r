plotting<-function(Method,Data,DistVar,EffectVar,SiteVar,Ylab,NameAppend='',Type='pdf&png'){
if(Method!='Non-Predator'){
	pages<-seq(1,length(unique(Data[,SiteVar])),6)
	pages.end<-c(pages[-1]-1,length(unique(Data[,SiteVar])))
		names(pages)<-LETTERS[1:length(pages)]
	prednames<-as.factor(c('PrBeet','PrSpHu','PrSpWe'))
	if(Type=='png'||Type=='pdf&png'){
		for(j in 1:length(pages)){
			png(paste('Figures & Tables/Raw Figures/',Method,EffectVar,SiteVar,NameAppend,'_',names(pages)[j],'.png',sep=''),width=6.5,height=9.75,res=300,bg=0,units='in')
			oldpar<-par(mfrow=c(3,2),mar=c(2.5,1.5,0,0),oma=c(2,2.5,.1,.1),cex=1)
			for(i in pages[j]:pages.end[j]){
				preds<-droplevels(Data[Data[,SiteVar]==unique(Data[,SiteVar])[i],])
				plot(c(0,max(preds[,DistVar])),c(0,1),col='white',xlab='',ylab='',yaxt='n')
					if(i%%2==1){axis(2,las=2,lwd=0,lwd.ticks=1)}
					else{axis(2,labels=FALSE,lwd=0,lwd.ticks=1)}
					points(preds[,DistVar],preds[,EffectVar],pch=match(preds$Group,prednames))
					legend('topright',legend=paste(substr(unique(Data[,SiteVar])[i],1,5),'-',substr(unique(Data[,SiteVar])[i],6,7),sep=''),bty='n')
					if(i%in%seq(2,200,6)){					
						legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
					}
					if(pages[j]==max(pages)&&max(pages)==max(pages.end)){
						legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
					}					
				if(pages.end[j]%%6==0){
					mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=.5)
					mtext(Ylab,side=2,outer=TRUE,line =1.5)	
				}
				else{
					if(pages[j]==max(pages)&&pages.end[j]-pages[j]==4){
						mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-.5)
						mtext(Ylab,side=2,outer=TRUE,line =1.5)	
					}
					else{
						if(pages[j]==max(pages)&&pages.end[j]-pages[j]==3){
							mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-15)
							mtext(Ylab,side=2,outer=TRUE,at=.69,line =1.5)	
						}
						else{
							if(pages[j]==max(pages)&&pages.end[j]-pages[j]==2){
								mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-15)
								mtext(Ylab,side=2,outer=TRUE,at=.69,line =1.5)	
							}
							else{
								if(pages[j]==max(pages)&&pages.end[j]-pages[j]==1){
									mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-31)
									mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)	
								}	
								else{
									mtext("Distance from water's edge (m)",side=1,outer=TRUE,at=.28,line=-31)
									mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)								
				}}}}}				
			}
		par(oldpar)
		dev.off()
		}
	}
	if(Type=='pdf'||Type=='pdf&png'){
		pdf(paste('Figures & Tables/Raw Figures/',Method,EffectVar,SiteVar,NameAppend,'.pdf',sep=''),width=6.5,height=9.75)
		oldpar<-par(mfrow=c(3,2),mar=c(2.5,1.5,0,0),oma=c(2,2.5,.1,.1),cex=1)
		for(i in 1:length(unique(Data[,SiteVar]))){
			preds<-droplevels(Data[Data[,SiteVar]==unique(Data[,SiteVar])[i],])
			plot(c(0,max(preds[,DistVar])),c(0,1),col='white',xlab='',ylab='',yaxt='n')
				if(i%%2==1){axis(2,las=2,lwd=0,lwd.ticks=1)}
				else{axis(2,labels=FALSE,lwd=0,lwd.ticks=1)}
				points(preds[,DistVar],preds[,EffectVar],pch=match(preds$Group,prednames))
				legend('topright',legend=paste(substr(unique(Data[,SiteVar])[i],1,5),'-',substr(unique(Data[,SiteVar])[i],6,7),sep=''),bty='n')
				if(i%in%seq(2,200,6)){					
					legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
				}
				if(i==max(pages)&&max(pages)==max(pages.end)){
					legend('top',bty='n',legend=c('Beetles','Hunting spiders','Webbed spiders'),pch=c(1:3))
				}	
				if(i%in%pages&&i!=max(pages)){
					mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=.5)
					mtext(Ylab,side=2,outer=TRUE,line =1.5)	
				}
				if(i%in%pages&&i==max(pages)){
					if(max(pages.end)-max(pages)==5||max(pages.end)-max(pages)==4){
						mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-.5)
						mtext(Ylab,side=2,outer=TRUE,line =1.5)						
					}
					else{
						if(max(pages.end)-max(pages)==3||max(pages.end)-max(pages)==2){
							mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-15)
							mtext(Ylab,side=2,outer=TRUE,at=.69,line =1.5)						
						}
						else{
							if(max(pages.end)-max(pages)==1){
								mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-31)
								mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)								
							}
							else{
								mtext("Distance from water's edge (m)",side=1,outer=TRUE,at=.28,line=-31)
								mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)								
					}}}
				}				
		}
		par(oldpar)
		dev.off()
	}
}
if(Method=='Non-Predator'){
	pages<-seq(1,length(unique(Data[,SiteVar])),6)
	pages.end<-c(pages[-1]-1,length(unique(Data[,SiteVar])))
		names(pages)<-LETTERS[1:length(pages)]
	prednames<-as.factor(c('HeHemi','HeLepi','AeLepi','TeMisc'))
	if(Type=='png'||Type=='pdf&png'){
		for(j in 1:length(pages)){
			png(paste('Figures & Tables/Raw Figures/',Method,EffectVar,SiteVar,NameAppend,'_',names(pages)[j],'.png',sep=''),width=6.5,height=9.75,res=300,bg=0,units='in')
			oldpar<-par(mfrow=c(3,2),mar=c(2.5,1.5,0,0),oma=c(2,2.5,.1,.1),cex=1)
			for(i in pages[j]:pages.end[j]){
				preds<-droplevels(Data[Data[,SiteVar]==unique(Data[,SiteVar])[i],])
				plot(c(0,max(preds[,DistVar])),c(0,max(preds[,EffectVar],na.rm=T)),col='white',xlab='',ylab='',yaxt='n')
					if(i%%2==1){axis(2,las=2,lwd=0,lwd.ticks=1)}
					else{axis(2,labels=FALSE,lwd=0,lwd.ticks=1)}
					points(preds[,DistVar],preds[,EffectVar],pch=match(preds$Group,prednames))
					legend('topright',legend=paste(substr(unique(Data[,SiteVar])[i],1,5),'-',substr(unique(Data[,SiteVar])[i],6,7),sep=''),bty='n')
					if(i%in%seq(2,200,6)){					
						legend('top',bty='n',legend=c('Hemipteran herbivores','Lepidoptera herbivores','Other lepidoptera','Miscellaneous other insects'),pch=c(1:4))
					}
					if(pages[j]==max(pages)&&max(pages)==max(pages.end)){
						legend('top',bty='n',legend=c('Hemipteran herbivores','Lepidoptera herbivores','Other lepidoptera','Miscellaneous other insects'),pch=c(1:4))
					}					
				if(pages.end[j]%%6==0){
					mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=.5)
					mtext(Ylab,side=2,outer=TRUE,line =1.5)	
				}
				else{
					if(pages[j]==max(pages)&&pages.end[j]-pages[j]==4){
						mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-.5)
						mtext(Ylab,side=2,outer=TRUE,line =1.5)	
					}
					else{
						if(pages[j]==max(pages)&&pages.end[j]-pages[j]==3){
							mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-15)
							mtext(Ylab,side=2,outer=TRUE,at=.69,line =1.5)	
						}
						else{
							if(pages[j]==max(pages)&&pages.end[j]-pages[j]==2){
								mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-15)
								mtext(Ylab,side=2,outer=TRUE,at=.69,line =1.5)	
							}
							else{
								if(pages[j]==max(pages)&&pages.end[j]-pages[j]==1){
									mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-31)
									mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)	
								}	
								else{
									mtext("Distance from water's edge (m)",side=1,outer=TRUE,at=.28,line=-31)
									mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)								
				}}}}}				
			}
		par(oldpar)
		dev.off()
		}
	}
	if(Type=='pdf'||Type=='pdf&png'){
		pdf(paste('Figures & Tables/Raw Figures/',Method,EffectVar,SiteVar,NameAppend,'.pdf',sep=''),width=6.5,height=9.75)
		oldpar<-par(mfrow=c(3,2),mar=c(2.5,1.5,0,0),oma=c(2,2.5,.1,.1),cex=1)
		for(i in 1:length(unique(Data[,SiteVar]))){
			preds<-droplevels(Data[Data[,SiteVar]==unique(Data[,SiteVar])[i],])
			plot(c(0,max(preds[,DistVar])),c(0,max(preds[,EffectVar],na.rm=T)),col='white',xlab='',ylab='',yaxt='n')
				if(i%%2==1){axis(2,las=2,lwd=0,lwd.ticks=1)}
				else{axis(2,labels=FALSE,lwd=0,lwd.ticks=1)}
				points(preds[,DistVar],preds[,EffectVar],pch=match(preds$Group,prednames))
				legend('topright',legend=paste(substr(unique(Data[,SiteVar])[i],1,5),'-',substr(unique(Data[,SiteVar])[i],6,7),sep=''),bty='n')
				if(i%in%seq(2,200,6)){					
					legend('top',bty='n',legend=c('Hemipteran herbivores','Lepidoptera herbivores','Other lepidoptera','Miscellaneous other insects'),pch=c(1:4))
				}
				if(i==max(pages)&&max(pages)==max(pages.end)){
					legend('top',bty='n',legend=c('Hemipteran herbivores','Lepidoptera herbivores','Other lepidoptera','Miscellaneous other insects'),pch=c(1:4))
				}	
				if(i%in%pages&&i!=max(pages)){
					mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=.5)
					mtext(Ylab,side=2,outer=TRUE,line =1.5)	
				}
				if(i%in%pages&&i==max(pages)){
					if(max(pages.end)-max(pages)==5||max(pages.end)-max(pages)==4){
						mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-.5)
						mtext(Ylab,side=2,outer=TRUE,line =1.5)						
					}
					else{
						if(max(pages.end)-max(pages)==3||max(pages.end)-max(pages)==2){
							mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-15)
							mtext(Ylab,side=2,outer=TRUE,at=.69,line =1.5)						
						}
						else{
							if(max(pages.end)-max(pages)==1){
								mtext("Distance from water's edge (m)",side=1,outer=TRUE,line=-31)
								mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)								
							}
							else{
								mtext("Distance from water's edge (m)",side=1,outer=TRUE,at=.28,line=-31)
								mtext(Ylab,side=2,outer=TRUE,at=.86,line =1.5)								
					}}}
				}				
		}
		par(oldpar)
		dev.off()
	}
}	
}
