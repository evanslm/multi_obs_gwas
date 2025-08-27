### Plot patterns as functions of parameter values

d<-NULL
pwr<-NULL
param.table<-read.table("ParamCombs.txt",header=T)
files<-list.files(path="results/",pattern=paste(paste0("^sims",1:24,".RDat"),collapse="|")) ## the "general" simulations. Simulation parameter combinations 25-32 vary N or missingness, and are plotted separately.
for(rr in 1:length(files)){
	file<-paste0("results/",files[rr])
	load(file,verbose=T)
	pcombos<-param.table[which(param.table$combination==simresults$params['combination']),]
	d<-rbind(d,cbind(simresults$U,pcombos,nrow(simresults$U)))
	# pwr<-rbind(pwr,cbind(rr,"single",as.vector(simresults$U[,1]),pcombos))
	# pwr<-rbind(pwr,cbind(rr,"mean",as.vector(simresults$U[,2]),pcombos))
}

apply(d[,grep("null",colnames(d))],2,table) ### Almost no false positives, even for the single observation.

pwr.single<-d$nsig_single/d$M
pwr.mean<-d$nsig_mean/d$M
pwr.med<-d$nsig_median/d$M
pwr.lowess<-d$nsig_lowess/d$M

pwr<-c(pwr.single,pwr.mean,pwr.med,pwr.lowess)
corbeta<-c(d$single_corb,d$mean_corb,d$median_corb,d$lowess_corb)
prsr2<-c(d$single_prsr2,d$mean_prsr2,d$mean_prsr2,d$lowess_prsr2)
trait<-c(rep("1single",nrow(d)),rep("2mean",nrow(d)),rep("3median",nrow(d)),rep("4lowess",nrow(d)))
n<-as.factor(c(d$n,d$n,d$n,d$n))
h2<-as.factor(c(d$h2,d$h2,d$h2,d$h2))
h2ge<-as.factor(c(d$h2ge,d$h2ge,d$h2ge,d$h2ge))
mean.slope<-as.factor(c(d$mean.slope,d$mean.slope,d$mean.slope,d$mean.slope))
het.slope<-as.factor(c(d$het.slope,d$het.slope,d$het.slope,d$het.slope))
d.model<-data.frame(pwr,corbeta,trait,n,h2,h2ge,mean.slope,het.slope)

touse<-1:nrow(d.model) #which(d.model$mean.slope!=.1)
cols<-rep(c(gray(0.1,0.2),rgb(1,0,0,.2),rgb(0,1,0,0.2),rgb(0,0,1,0.2)),3)
box.rel.loc<-(1:4)*0.5/3-median((1:4)*0.5/3)
box.loc<-c(1+box.rel.loc,2+box.rel.loc,3+box.rel.loc,4+box.rel.loc)

### models to test parameter effects
fit <- lm( pwr ~ trait + h2 + h2ge + trait:h2 + trait:h2ge + h2:h2ge + het.slope + mean.slope, data=d.model[touse,] )
summary(fit)
anova(fit)

fit <- lm( corbeta ~ trait + h2 + h2ge + trait:h2 + trait:h2ge + h2:h2ge + het.slope + mean.slope, data=d.model[touse,] )
summary(fit)
anova(fit)

fit <- lm( prsr2 ~ trait + h2 + h2ge + trait:h2 + trait:h2ge + h2:h2ge + het.slope + mean.slope, data=d.model[touse,] )
summary(fit)
anova(fit)

## Plot the results:
jpeg("Fig.1.jpg",height=12,width=6,units='in',res=300)
par(mfcol=c(3,1),mar=c(3,3,3,.2),mgp=c(2,.75,0),oma=c(1,1,1,1))
bp<-boxplot(pwr~h2+h2ge+trait,cex.axis=0.5,data=d.model[touse,],boxwex=0.15,at=box.loc,cex=0.2,col=cols,border=1:4,xaxt='n',ylab="Power to Detect Associated Loci",xlab='Trait Definition')
axis(1,at=1:4,labels=c("single","mean","median","lowess"))
legend("topleft",c("0.1,0","0.5,0","0.1,0.05","0.5,0.05"),col=1:4,pch=15,title="h2, h2gxe",bty='n')

bp<-boxplot(corbeta~h2+h2ge+trait,cex.axis=0.5,data=d.model[touse,],boxwex=0.15,at=box.loc,cex=0.2,col=cols,border=1:4,xaxt='n',ylab="Correlation of True and Estiamted Betas",xlab='Trait Definition')
axis(1,at=1:4,labels=c("single","mean","median","lowess"))
legend("topleft",c("0.1, 0","0.5, 0","0.1, 0.05","0.5, 0.05"),col=1:4,pch=15,title="h2, h2gxe",bty='n')

bp<-boxplot(prsr2 ~h2+h2ge+trait,cex.axis=0.5,data=d.model[touse,],boxwex=0.15,at=box.loc,cex=0.2,col=cols,border=1:4,xaxt='n',ylab="Out-Of-Sample PRS r^2",xlab='Trait Definition')
axis(1,at=1:4,labels=c("single","mean","median","lowess"))
legend("topleft",c("0.1, 0","0.5, 0","0.1, 0.05","0.5, 0.05"),col=1:4,pch=15,title="h2, h2gxe",bty='n')
dev.off()

rbc<-rainbow(24,start=0,end=1)

jpeg("Fig.S1.jpg",height=8,width=12,units='in',res=300)
par(mar=c(3,3,3,.2),mgp=c(2,.75,0),oma=c(1,1,1,1))
bp<-boxplot(pwr~mean.slope + het.slope + h2 + h2ge+ trait,cex.axis=0.35,data=d.model[touse,],boxwex=0.4,ylab="Power to Detect Associated Loci",cex=0.25,las=2,pch=20,col=rbc)
abline(v=c(24.5,48.5,72.5),lty=3,lwd=0.5,col=gray(0.2,0.7))
legend("bottomright",gsub('[.]1single','',bp$names[1:24]),col=rbc,pch=15,bty='n',cex=0.75)
dev.off()

jpeg("Fig.S2.jpg",height=8,width=12,units='in',res=300)
par(mar=c(3,3,3,.2),mgp=c(2,.75,0),oma=c(1,1,1,1))
bp<-boxplot(corbeta~mean.slope + het.slope + h2 + h2ge +  trait,cex.axis=0.35,data=d.model[touse,],boxwex=0.4,ylab="Correlation of True and Estiamted Betas",cex=0.25,las=2,pch=20,col=rbc)
abline(v=c(24.5,48.5,72.5),lty=3,lwd=0.5,col=gray(0.2,0.7))
legend("bottomright",gsub('[.]1single','',bp$names[1:24]),col=rbc,pch=15,bty='n',cex=0.75)
dev.off()

jpeg("Fig.S3.jpg",height=8,width=12,units='in',res=300)
par(mar=c(3,3,3,.2),mgp=c(2,.75,0),oma=c(1,1,1,1))
bp<-boxplot(prsr2~mean.slope + het.slope + h2 + h2ge +  trait,cex.axis=0.35,data=d.model[touse,],boxwex=0.4,ylab="Out-Of-Sample PRS r^2",cex=0.25,las=2,pch=20,col=rbc)
abline(v=c(24.5,48.5,72.5),lty=3,lwd=0.5,col=gray(0.2,0.7))
# legend("bottomright",gsub('[.]1single','',bp$names[1:24]),col=rbc,pch=15,bty='n',cex=0.75)
dev.off()
