#########################################################################################
## LOAD PACKAGES AND USEFUL FUNCTIONS
#########################################################################################

library(huge)
library(qgraph)
library(Matrix)
library(mvtnorm)
#use setwd to the directory where sp500.RData is stored in your computer
load("sp500.RData")



#Function to extract EBIC from 'huge' object
getEBIC= function(fit, n, d, ebic.gamma=1) {
  if (class(fit) != 'huge') stop("fit must be of class 'huge'")
  ebic= -n * fit$loglik + fit$df * (log(n) + 4 * ebic.gamma * log(d)) #formula used in huge package
  sel.model= which(ebic==min(ebic))
  K= forceSymmetric(fit$icov[[sel.model]])
  ans= list(ebic=ebic, opt.icov= K)
  return(ans)
}


#Obtain partial correlations given inverse covariance matrix
partialcor= function(K) {
    D= diag(1/sqrt(diag(K)))
    ans= -D %*% K %*% D
    diag(ans)=1
    return(ans)
}


#Non-parametric estimate of partial correlation between functions (from partial residuals)
partialcorf= function(x, f1=function(z) z, f2= f1) {
    idx= expand.grid(1:ncol(x),1:ncol(x))
    idx= idx[idx[,1]<idx[,2],]
    foo= function(i) {
      idx1= idx[i,1]; idx2= idx[i,2]
      r1= residuals( lm(x[,idx1] ~ x[,-c(idx1,idx2)]) ) 
      r2= residuals( lm(x[,idx2] ~ x[,-c(idx1,idx2)]) ) 
      if ((i %% (nrow(idx)/10)) == 0) cat(".")
      return( cor(f1(r1),f2(r2)) )
    }
    ans= unlist(mclapply(1:nrow(idx), foo, mc.cores=4))
    ans= data.frame(idx, partialcor=ans)
    return(ans)
}



#########################################################################################
## SP500 ANALYSIS
#########################################################################################


#Format data
y= log(sp500[-1,-1] / sp500[-nrow(sp500),-1])  #daily log-returns
z= scale(y)                                    #center to 0 mean, scale to variance 1
z.npn= huge.npn(z[,1:100])                     # marginal transformations to improve normality
S= cov(z.npn)                                  #Pearson covariance on transformed data (elliptical model)
S.skep= huge.npn(z[,1:100],npn.func="skeptic") # Covariance estimate from Kendall's tau (transelliptical model)
znorm= scale(rmvnorm(n=nrow(z), sigma=S))
znorm.skep= scale(rmvnorm(n=nrow(z), sigma=S.skep))


# Figure: Check Normality via Mahalanobis distances
d.npn= mahalanobis(z.npn, center=rep(0,ncol(z.npn)), cov=S)
xlim= c(0, max(d.npn))
dseq= seq(0,1.5*max(d.npn),length=1000)
plot(dseq, dchisq(dseq, df=ncol(S)), xlim=xlim, col='gray', lwd=2, lty=2, type='l', ylab='Density', xlab='Mahalanobis distance', cex.lab=1.2, cex.axis=1.3)
lines(density(d.npn))
legend('topright',c('Observed data','Chi-square'), lty=c(1,2), lwd=2, col=c(1,'gray'), cex=1.4)


# Figure: Marginal tail dependence
r2= cor(z.npn^2); r2norm= cor(znorm^2); r2expected= S^2
r2dif= r2[upper.tri(r2)] - r2expected[upper.tri(r2expected)]
r2difnorm= r2norm[upper.tri(r2norm)] - r2expected[upper.tri(r2expected)]
xlim= 1.2*c(min(min(r2dif),min(r2difnorm)), max(max(r2dif),max(r2difnorm)))

plot(density(r2difnorm), col='gray', lty=2, xlim=xlim, main='', xlab='Marginal tail dependence', cex.lab=1.3, cex.axis=1.3)
lines(density(r2dif))
legend('topright', c('Observed data','Simulated normal data'), lty=1:2, col=c(1,'gray'), cex=1.3)


#Fit model (elliptical on transformed data)
zfit= huge(z.npn, method="glasso", nlambda=200, lambda.min.ratio=0.01)
ebic= getEBIC(zfit, n=nrow(z.npn), d= ncol(z.npn), ebic.gamma=1)
icov.sel= ebic$opt.icov[upper.tri(ebic$opt.icov)]
sign.sel= sign(icov.sel)

table(sign.sel)


#Fit model (transelliptical on original data)
zfit.skep = huge(S.skep, method="glasso", nlambda=200, lambda.min.ratio=0.01)
ebicskep= getEBIC(zfit.skep, n=nrow(z.npn), d= ncol(z.npn), ebic.gamma=1)
icovskep.sel= ebicskep$opt.icov[upper.tri(ebicskep$opt.icov)]
signskep.sel= sign(icovskep.sel)

table(signskep.sel)


#Compare elliptical and trans-elliptical models: edges, Spearman rho between estimates
table(sign.sel,signskep.sel)

cor(icov.sel, icovskep.sel, method='spearman')  


#Fit model to simulated normal data
zfitnorm= huge(znorm, method="glasso", nlambda=200, lambda.min.ratio=0.01)
ebicnorm= getEBIC(zfitnorm, n=nrow(znorm), d= ncol(znorm), ebic.gamma=1)
signnorm.sel= sign(ebicnorm$opt.icov[upper.tri(ebicnorm$opt.icov)])


#Figure: Conditional tail dependence
icovsel= ebic$opt.icov
icovselskep= ebicskep$opt.icov
icovselnorm= ebicnorm$opt.icov

rpartial= partialcor(K=icovsel)                    #partial correlations (model-based)
r2partial= partialcorf(z.npn, f1=function(z) z^2)  #partial correlations between squares
r2partialexpected= rpartial[upper.tri(rpartial)]^2 #expected r2partial under Normality

#tail dependence for simulated normal data
rpartialnorm= partialcor(K=icovselnorm)  #partial correlations (model-based)
r2partialnorm= partialcorf(znorm, f1=function(z) z^2) #partial correlations between squares
r2partialnormexpected= rpartialnorm[upper.tri(rpartialnorm)]^2 #expected r2partial under Normality

r2dif= r2partial$partialcor - r2partialexpected
r2difnorm= r2partialnorm$partialcor - r2partialnormexpected


xlim= 1.2*c(min(min(r2dif),min(r2difnorm)), max(max(r2dif),max(r2difnorm)))
plot(density(r2difnorm), col='gray', lty=2, xlim=xlim, main='', xlab='Conditional tail dependence', cex.lab=1.3, cex.axis=1.3)
lines(density(r2dif))
legend('topright', c('Observed data','Simulated normal data'), lty=1:2, col=c(1,'gray'), cex=1.3)


#Figure: Conditional tail dependence vs squared partial correlations
xlim= c(0,0.4); ylim= c(-0.1,0.4)
plot(r2partial$partialcor ~ r2partialexpected, xlim=xlim, ylim=ylim, xlab='Squared partial correlation', ylab='Correlation between squares', cex.lab=1.3, cex.axis=1.3)
lines(lowess(r2partial$partialcor ~ r2partialexpected))
lines(lowess(r2partialnorm$partialcor ~ r2partialnormexpected), col='gray', lty=2, lwd=2)
legend('bottomright', c('Observed data','Simulated normal data'), lty=1:2, col=c(1,'gray'), cex=1.3)




#########################################################################################
## ANALYSIS UNDER POSITIVE PARTIAL CORRELATIONS RESTRICTION
#########################################################################################

library(devtools)
install_github("pzwiernik/golazo")
library(golazo)
S= huge.npn(z[,1:100],npn.func="skeptic")
# use the GOLAZO package
resM <- positive.golazo(S,rho=Inf)
# this is the estimated inverse covariance matrix
KM <- resM$K
# this is the number of estimated edges
thres <- 1e-6
(sum(abs(cov2cor(KM))>thres)-100)/2
# cmopute the log-likelihood
n=nrow(z)

#compare log-likelihood to GLASSO solution
Sglasso= as.matrix(solve(ebicskep$opt.icov))
loglglasso= sum(dmvnorm(z.npn, sigma= Sglasso, log=TRUE))
# log-likelihood of the MTP2 estimate
Smtp2= resM$Sig
loglmtp2= sum(dmvnorm(z.npn, sigma= Smtp2, log=TRUE))
# refitted MLE for the GLASSO MODEL (we can also do it using the GOLAZO package)
U <- Inf*(abs(cov2cor(solve(Sglasso)))<thres)
U[is.na(U)] <- 0
diag(U) <- 0
resREFIT <- golazo(S,L=-U,U=U)
loglrefit= sum(dmvnorm(z.npn, sigma= resREFIT$Sig, log=TRUE))

# compare the EBIC criteria
gamma <- 0
dglasso <- (sum(abs(cov2cor(solve(Sglasso)))>thres)-100)/2 # 1600 edges
EBICglasso <- -2*loglglasso+dglasso*(log(n) +4*gamma*log(100))
dMTP2 <- (sum(abs(cov2cor(KM))>thres)-100)/2
EBICMTP2 <- -2*loglmtp2+dMTP2*(log(n) +4*gamma*log(100))
drefit <- (sum(abs(cov2cor(resREFIT$K))>thres)-100)/2 # 1600 edges
EBICrefit <- -2*loglrefit+drefit*(log(n) +4*gamma*log(100))

c(EBICglasso,EBICrefit,EBICMTP2)

