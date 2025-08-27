#!/curc/tools/x_86_64/rh6/R/3.1.3/bin/Rscript
library(foreach) # setting up environment to use multiple cores
library(iterators)
library(doMC)
registerDoMC(cores=6)
source("functions.R")
library(MASS)

args = commandArgs(trailingOnly=TRUE)

params<-as.numeric(args[1]) ### Parameter combination number. Input from command line (e.g., in slurm bash job)
load("ParamCombs.RDat",verbose=T)

simresults<-list() ### Object to store everything

param.values<-sapply(param.combins,"[[",params)
simresults$params<-param.values

cat("parameter values: ",param.values,"\n")
simresults$params<-param.values
nrep<-param.combins[['nrep']][params]
M<-param.combins[['M']][params]
n<-param.combins[['n']][params]
h2<-param.combins[['h2']][params]
nwaves<-param.combins[['nwaves']][params]
het.slope<-param.combins[['het.slope']][params]
mean.slope<-param.combins[['mean.slope']][params]
h2ge<-param.combins[['h2ge']][params]
gxe_maf<-TRUE ### This forces a MAF-BetaGxE relationship. If false, no MAF relationship with effect sizes

waves<-1:nwaves
waves.s<-(waves-mean(waves))/sd(waves)
midwave<-1-waves.s[floor(nwaves/2)]*het.slope
wave.var.fact<-het.slope*waves.s+midwave ### This makes the middle wave be the set h2 below. Also allows for heteroskedasticity. If het.slope=0, no heteroskedasticity
wave.mean<-mean.slope*waves.s ### making there be a trend over time. Could make it ^2 or something
meanvar<-cbind(wave.var.fact,wave.mean)
alpha<-0.05/M
traits<-c('single','mean','median','lowess')
ntraits<-length(traits)
trait.pairs<-combn(1:ntraits,2)
ntpairs<-ncol(trait.pairs)


U<-foreach(nn=1:1000,.combine=rbind,.errorhandling="remove")%dopar%{
	p<-runif(M,0.01,0.5)
	g1gen<-function(rep){rbinom(M,2,p)}
	g<-t(cbind(sapply(rep(0,n),g1gen)))
	g2<-t(cbind(sapply(rep(0,n),g1gen)))
	g.null<-t(cbind(sapply(rep(0,n),g1gen)))

	if(h2ge>0) {
		if(gxe_maf){
			bgvar<-1/M*h2/(2*p*1-p)
			bivar<-1/M*h2ge/(2*p*1-p)
			rho_gi<-0.2*sqrt(bgvar*bivar)
			b<-t(sapply(1:M,bgxe.samp,bgvar=bgvar,bivar=bivar,rho_gi=rho_gi))
		} else {
			rho_gi<-0.5*sqrt(h2*h2ge)
			sigma<-1/M*matrix(c(h2,rho_gi,rho_gi,h2ge),2,2)
			b<-mvrnorm(M,mu=rep(0,2),Sigma=sigma) ## Forcing it to have correlated main and GxE effects
		}
		bi<-b[,2]
		b<-b[,1]
		gp<-as.vector(g%*%b)
		gp2<-as.vector(g2%*%b)
		Vg<-var(gp)
		gxep<-g%*%bi%*%matrix(waves.s,nrow=1) ##Giving the Gx(scaled)Wave effect here
		gxep2<-g2%*%bi%*%matrix(waves.s,nrow=1) ##Giving the Gx(scaled)Wave effect here
		Vge<-var(as.vector(gxep)) #apply(gxep,2,var) ## Doing total gxep variance vs. per-wave. Probably former?
		VG<-Vg+Vge
		VE<-VG/(h2+h2ge)-VG
		ep<-matrix(rnorm(n*nwaves,0,sqrt(VE)),n,nwaves) ## Just doing total variance here
		ep2<-matrix(rnorm(n*nwaves,0,sqrt(VE)),n,nwaves) ## Just doing total variance here
		if(mean.slope>0) for(i in waves) {
			ep[,i]<-ep[,i]+wave.mean[i]
			ep2[,i]<-ep2[,i]+wave.mean[i]
		}
		Ve<-apply(ep,2,var)
		y<-gp+gxep+ep
		y2<-gp2+gxep2+ep2
		Vt<-apply(y,2,var)
	}else { 
		# b<-rnorm(M,0,sqrt(h2/M)) ## If wanting no MAF-beta relationship
		b<-rnorm(M,0,sqrt(h2/(M*2*p*(1-p))))
		gp<-as.vector(g%*%b)
		gp2<-as.vector(g2%*%b)
		Vg<-var(gp)
		VE<-Vg/h2-Vg
		ep<-matrix(rnorm(n*nwaves,0,sqrt(VE)),n,nwaves) ## Just doing total variance here
		ep2<-matrix(rnorm(n*nwaves,0,sqrt(VE)),n,nwaves)
		# ep<-sapply(wave.var.fact*sqrt(VE),rnorm,n=n,mean=0) ### This makes the middle wave be the set h2.
		if(mean.slope>0) for(i in waves) {
			ep[,i]<-ep[,i]+wave.mean[i]
			ep2[,i]<-ep2[,i]+wave.mean[i]
		}
		Ve<-apply(ep,2,var)
		y<-gp+ep
		y2<-gp2+ep2
		Vt<-apply(y,2,var)
	}
	# plot(waves,y[1,],ylim=range(y),type='l')
	# for(i in 1:10) points(waves,y[i,],type='l',col=i)
	# points(waves,apply(y,2,mean),type='l',col=1,lwd=3)
	# Vg/Vt		
	# Ve/Vt
	low.y<-sapply(1:n,loess.integ,y=y,x=waves)
	med.y<-apply(y,1,median,na.rm=T)
	mean.y<-apply(y,1,mean,na.rm=T)	
	gwas<-list()
	gwas[["single"]]<-t(sapply(1:M,gwas.funct,y=y[,nwaves],g=g))
	gwas[["mean"]]<-t(sapply(1:M,gwas.funct,y=mean.y,g=g))
	gwas[["median"]]<-t(sapply(1:M,gwas.funct,y=med.y,g=g))
	gwas[["low"]]<-t(sapply(1:M,gwas.funct,y=low.y,g=g))
	gwas[["single.null"]]<-t(sapply(1:M,gwas.funct,y=y[,nwaves],g=g.null))
	gwas[["mean.null"]]<-t(sapply(1:M,gwas.funct,y=mean.y,g=g.null))
	gwas[["median.null"]]<-t(sapply(1:M,gwas.funct,y=med.y,g=g.null))
	gwas[["low.null"]]<-t(sapply(1:M,gwas.funct,y=low.y,g=g.null))
	nsig<-unlist(lapply(c("single","mean","median","low","single.null","mean.null","median.null","low.null"),count.sig,gwas=gwas,alpha=alpha))
	bcor<-unlist(lapply(1:ntraits,beta.corr,gwas=gwas,b=b))
	bcor.pairs<-sapply(1:ntpairs,beta.corr.pairs,gwas=gwas,trait.pairs=trait.pairs)
	
	### Test PRS in independent sample (y2) with a single observation
	prs2.single<-g2%*%gwas[['single']][,1]
	prs2.mean<-g2%*%gwas[['mean']][,1]
	prs2.med<-g2%*%gwas[['median']][,1]	
	prs2.loess<-g2%*%gwas[['low']][,1]	
	y2.single.obs<-t(sapply(1:n,samp.ind.obs,y=y2,nwaves=nwaves))
	r2single<-summary(lm(y2.single.obs[,1]~prs2.single+y2.single.obs[,2]))$adj.r.squared-	summary(lm(y2.single.obs[,1]~y2.single.obs[,2]))$adj.r.squared
	r2mean<-summary(lm(y2.single.obs[,1]~prs2.mean+y2.single.obs[,2]))$adj.r.squared-	summary(lm(y2.single.obs[,1]~y2.single.obs[,2]))$adj.r.squared
	r2med<-summary(lm(y2.single.obs[,1]~prs2.med+y2.single.obs[,2]))$adj.r.squared-	summary(lm(y2.single.obs[,1]~y2.single.obs[,2]))$adj.r.squared
	r2lowess<-summary(lm(y2.single.obs[,1]~prs2.loess+y2.single.obs[,2]))$adj.r.squared-	summary(lm(y2.single.obs[,1]~y2.single.obs[,2]))$adj.r.squared
		
	c(nsig,bcor,bcor.pairs,r2single,r2mean,r2med,r2lowess)
}
colnames(U)<-c(paste0("nsig_",traits),paste0("nsig_",traits,".null"),paste0(traits[1],"_corb"),paste0(traits[2],"_corb"),paste0(traits[3],"_corb"),paste0(traits[4],"_corb"),paste0(traits[trait.pairs[,1]],collapse="_"),paste0(traits[trait.pairs[,2]],collapse="_"),paste0(traits[trait.pairs[,3]],collapse="_"),paste0(traits[trait.pairs[,4]],collapse="_"),paste0(traits[trait.pairs[,5]],collapse="_"),paste0(traits[trait.pairs[,6]],collapse="_"),paste0(traits,"_prsr2"))

summary(U)
simresults$U<-U
save(simresults,file=paste0("sims",params,".RDat"))


