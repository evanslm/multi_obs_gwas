### Plot patterns as functions of parameter values

d<-NULL
pwr<-NULL
param.table<-read.table("ParamCombs.txt",header=T)
files<-list.files(path="results/",pattern=paste(paste0("^sims",c(1,7,13,19,25:28),".RDat"),collapse="|"))
for(rr in 1:length(files)){
	file<-paste0("results/",files[rr])
	load(file,verbose=T)
	pcombos<-param.table[which(param.table$combination==simresults$params['combination']),]
	d<-rbind(d,cbind(simresults$U,pcombos))
	# pwr<-rbind(pwr,cbind(rr,"single",as.vector(simresults$U[,1]),pcombos))
	# pwr<-rbind(pwr,cbind(rr,"mean",as.vector(simresults$U[,2]),pcombos))
}
unique(d[,23:33]) ## Parameter combinations included here
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
cols<-rep(c(gray(0.1,0.2),gray(0.1,0.8),rgb(1,0,0,.2),rgb(1,0,0,.8),rgb(0,1,0,0.2),rgb(0,1,0,0.8),rgb(0,0,1,0.2),rgb(0,0,1,0.8)),4)
bor.cols<-rep(c(1,1,2,2,3,3,4,4),4)

box.rel.loc<-(1:8)*0.3/7-median((1:8)*0.3/7)
box.loc<-c(1+box.rel.loc,2+box.rel.loc,3+box.rel.loc,4+box.rel.loc)

### models to test parameter effects
fit <- lm( pwr ~ trait + h2 + h2ge + n + trait:h2 + trait:h2ge + trait:n + h2:n + h2:h2ge + h2ge:n, data=d.model[touse,] )
summary(fit)
anova(fit)

fit <- lm( corbeta ~ trait + h2 + h2ge + n + trait:h2 + trait:h2ge + trait:n + h2:n + h2:h2ge + h2ge:n, data=d.model[touse,] )
summary(fit)
anova(fit)

fit <- lm( prsr2 ~ trait + h2 + h2ge + n + trait:h2 + trait:h2ge + trait:n + h2:n + h2:h2ge + h2ge:n, data=d.model[touse,] )
summary(fit)
anova(fit)

## Plot the results:
jpeg("Fig.S4.jpg",height=12,width=8,units='in',res=300)
par(mfcol=c(3,1),mar=c(3,3,3,.2),mgp=c(2,.75,0),oma=c(1,1,1,1))
bp<-boxplot(pwr~n+h2+h2ge+trait,cex.axis=0.5,data=d.model[touse,],boxwex=0.05,at=box.loc,cex=0.2,col=cols,border=bor.cols,xaxt='n',ylab="Power to Detect Associated Loci",xlab='Trait Definition')
axis(1,at=1:4,labels=c("single","mean","median","lowess"))
legend("topleft",c("0.1, 0","0.5, 0","0.1, 0.05","0.5, 0.05","10K","20K"),col=c(1:4,gray(0.1,0.2),gray(0.1,0.8)),pch=15,title="h2, h2gxe and N",bty='n')

bp<-boxplot(corbeta~n+h2+h2ge+trait,cex.axis=0.5,data=d.model[touse,],boxwex= 0.05,at=box.loc,cex=0.2,col=cols,border=bor.cols,xaxt='n',ylab="Correlation of True and Estiamted Betas",xlab='Trait Definition')
axis(1,at=1:4,labels=c("single","mean","median","lowess"))
legend("bottomright",c("0.1, 0","0.5, 0","0.1, 0.05","0.5, 0.05","10K","20K"),col=c(1:4,gray(0.1,0.2),gray(0.1,0.8)),pch=15,title="h2, h2gxe and N",bty='n')

bp<-boxplot(prsr2 ~ n+h2+h2ge+trait,cex.axis=0.5,data=d.model[touse,],boxwex= 0.05,at=box.loc,cex=0.2,col=cols,border=bor.cols,xaxt='n',ylab="Out-Of-Sample PRS r^2",xlab='Trait Definition',ylim=c(0,0.505))
axis(1,at=1:4,labels=c("single","mean","median","lowess"))
legend("left",c("0.1, 0","0.5, 0","0.1, 0.05","0.5, 0.05","10K","20K"),col=c(1:4,gray(0.1,0.2),gray(0.1,0.8)),pch=15,title="h2, h2gxe and N",bty='n')
abline(h=c(0.1,0.1,0.5,0.5),col=c(1,3,2,4),lty=c(1,2,1,2),lwd=0.6)
dev.off()

