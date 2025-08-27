### Functions for defining longitudinal phenotypes and also to run GWAS sims

make.missing<-function(ind,y,nwaves,nmiss){
	make.NA<-sample(nwaves,nmiss)
	y[ind,make.NA]<-NA
	return(y[ind,])
}

samp.ind.obs<-function(ind,y,nwaves){
	sw<-sample(nwaves,1,replace=F)
	c(y[ind,sw],sw)
}

samp.last.obs<-function(ind,y){
	y[ind,max(which(!is.na(y[ind,])))]
}

bgxe.samp<-function(l,bgvar,bivar,rho_gi,mu=c(0,0)){
	sigma<-matrix(c(bgvar[l],rho_gi[l],rho_gi[l],bivar[l]),2,2)
	mvrnorm(1,mu=mu,Sigma=sigma)
}

gwas.funct<-function(y=y,g=g,loc=loc){
	fit<-summary(lm(y~g[,loc]))
	(coef(fit)[2,])
}

gwas.funct.cov<-function(y=y,g=g,loc=loc,cov=cov){
	fit<-summary(lm(y~g[,loc]+cov))
	(coef(fit)[2,])
}

count.sig<-function(l=l,gwas=gwas,alpha=0.05,pcol=4){
	length(which(gwas[[l]][,pcol]<=alpha))
}

beta.corr<-function(l=l,gwas=gwas,b=b,bcol=1,pval=FALSE){
	be.corr<-cor.test(b,gwas[[l]][,bcol])
	if(pval) {
		return(c(be.corr$estimate,be.corr$p.value))
	}else{
		return(be.corr$estimate)
	}
}

beta.corr.pairs<-function(gwas=gwas,trait.pairs=trait.pairs,tp=tp){
	cor(gwas[[trait.pairs[1,tp]]][,1],gwas[[trait.pairs[2,tp]]][,1],use='pairwise')
}

harm.mean<-function(ind,y){
	harm.n<-length(na.omit(y[ind,]))
	harm.n/sum(1/y[ind,],na.rm=T)
}

geom.mean<-function(ind,y){
	geom.n<-length(na.omit(y[ind,]))
	exp(1/geom.n * sum(log(y[ind,]),na.rm=T))
}

loess.integ<-function(y=y,x=waves,ind=ind){
	l <- loess(y[ind,] ~ x, span=0.5, degree=2, family="symmetric", iterations=10)
	f <- function(x) predict(l,newdata=x)
	x.range<-diff(range(na.omit(x)))
	# perform integration
	(integrate(f, lower = min(x), upper = max(x)))$value/x.range
}
# plx<-predict(l, se=T)
# plot(x,y[ind,])
# points(x,plx$fit,type='l',col=2)
# abline(h=mean(y[ind,]))
# abline(h=(integrate(f, lower = min(x), upper = max(x)))$value/x.range,col=3)