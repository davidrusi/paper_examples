###############################################################################
##
## THIS SCRIPT REPRODUCES FIG 1 AND RELATED SUPPLEMENTARY FIGURES
## SEE simulation_iid.R FOR THE FULLER SIMULATION STUDY
##
###############################################################################

library(mombf)
library(tidyverse)


truemean= function(x,z) {
    ans= double(nrow(x))
    group1= (x[,1]==1)
    ans[group1]= ifelse(z[group1] <=0, cos(z[group1]), 1)
    ans[!group1]= ifelse(z[!group1]<=0, cos(z[!group1]), 1/(z[!group1]+1)^2)
    return(ans)
}


###############################################################################################
## 1. SIMULATED DATA WITH n=1000
###############################################################################################

set.seed(1)
n= 1000
x1= rep(0:1,c(n/2,n/2))
x2= x1 + rnorm(n)
x= cbind(x1,x2)
z= rep(seq(-3,3,length=n/2), 2)
m= truemean(x,z)
y= truemean(x,z) + rnorm(n, 0, .25)
znew= matrix(rep(seq(-2.9,2.9,length=50),2),ncol=1)
xnew= cbind(rep(0:1,each=nrow(znew)/2), mean(x2))
newdata= list(x=xnew, z=znew)

# Plot simulated data

par(mar=c(4.2,4.2,.1,.1))
sel= x[,1]==1
plot(z[sel], m[sel], ylim=range(y), xlab='z', ylab='y', type='l', col='gray', lwd=5, cex.axis=1.5, cex.lab=1.5)
lines(z[!sel], m[!sel], col='gray', lwd=5)
sel= x[,1]==1; points(z[sel], y[sel]); points(z[!sel], y[!sel], pch=2, col='gray')
text(x=3, y=c(-0.5,1.5), c('Group 0','Group 1'), cex=1.5, pos=2)
legend('topleft', c('Group 0','Group 1'), pch=2:1, col=c('gray',1), cex=1.5)


## 1.1 ANALYSIS BASED ON 20 KNOTS, 9 LOCAL KNOTS ##
###################################################

#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=20, nlocalknots=9, localgridsize=500)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=20, nlocalknots=9, localgridsize=500)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis


#Projection of the true mean on the spline basis
w0= mombf:::bspline(z, degree=3, knots=seq(-4,4,length=20))
w1= x[,1] * w0
w= cbind(w0, w1)
mproj.uncut= predict(lm(m ~ w)) #uncut basis


#Plot true means & uncut cubic spline projection
par(mar=c(4.2,4.2,.1,.1))
ylim= range(m)
sel= x[,1]==1
plot(z[sel], m[sel], ylim=ylim, xlab='z', ylab='y', type='l', col='gray', lwd=5, cex.axis=1.5, cex.lab=1.5)
lines(z[!sel], m[!sel], col='gray', lwd=5)
points(k, rep(min(m), length(k)), pch=16, cex=2)
sel1= ((x[,1]==1) & (abs(z) < 3))
sel0= ((x[,1]==0) & (abs(z) < 3))
lines(z[sel1], mproj.uncut[sel1], type='l')
lines(z[sel0], mproj.uncut[sel0])
text(x=3, y=c(0,1.02), c('Group 0','Group 1'), cex=1.5, pos=2)
segments(x0=k[5], x1=k[6], y0=ylim[1]-0.02, lwd=2)
segments(x0=k[5], x1=k[6], y0= ylim[1]+0.02, lwd=2)
segments(x0=k[5], y0= ylim[1]-0.02, y1= ylim[1]+0.02, lwd=2)
segments(x0=k[6], y0= ylim[1]-0.02, y1= ylim[1]+0.02, lwd=2)
abline(v=k[5:6], lty=2)


#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend('topleft', c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))




#Plot cubic spline basis
par(mar=c(4.2,4.2,.1,.1))
k= fit0$regionbounds[[1]][[1]]
zseq= seq(-3,3,length=200)
w= mombf:::bspline(zseq, degree=3, knots=k)
col= c('grey','black','black','black','black','grey','grey')
lty= ifelse(col=='grey', 2, 1)
plot(zseq, w[,1], type='l', xlab='z', ylab='Cubic B-spline', col=col[1], lty=lty[1], cex.lab=1.5, cex.axis=1.5)
for (i in 2:ncol(w)) lines(zseq, w[,i], col=col[i], lty=lty[i])
points(k, rep(0, length(k)), pch=16, cex=2)
segments(x0=k[5], x1=k[6], y0=-0.01, lwd=2)
segments(x0=k[5], x1=k[6], y0= 0.01, lwd=2)
segments(x0=k[5], y0= -0.01, y1= 0.01, lwd=2)
segments(x0=k[6], y0= -0.01, y1= 0.01, lwd=2)


#Plot cut cubic spline basis
par(mar=c(4.2,4.2,.1,.1))
wcut= cutbasisuniv(zseq, degree=3, regionbounds=k)$basis
col= c(rep('grey',9),rep('black',4),rep('grey',18))
lty= ifelse(col=='grey', 2, 1)
plot(zseq, wcut[,1], type='l', xlab='z', ylab='Cut cubic B-spline', col=col[1], lty=lty[1], cex.lab=1.5, cex.axis=1.5)
for (i in 2:ncol(wcut)) lines(zseq, wcut[,i], col=col[i], lty=lty[i])
points(k, rep(0, length(k)), pch=16, cex=2)
segments(x0=k[5], x1=k[6], y0=-0.01, lwd=2)
segments(x0=k[5], x1=k[6], y0= 0.01, lwd=2)
segments(x0=k[5], y0= -0.01, y1= 0.01, lwd=2)
segments(x0=k[6], y0= -0.01, y1= 0.01, lwd=2)



## 1.2 ANALYSIS BASED ON 20 KNOTS, 15 LOCAL KNOTS ##
####################################################

nbaseknots= 20; nlocalknots= 15

#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis

#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend('topleft', c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))




## 1.3 ANALYSIS BASED ON 30 KNOTS, 30 LOCAL KNOTS ##
####################################################

nbaseknots= 30; nlocalknots= 30

#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis


#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend('topleft', c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))



## 1.4 ANALYSIS BASED ON 12 KNOTS, 8 LOCAL KNOTS ##
###################################################

nbaseknots= 12; nlocalknots= 8
#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=nbaseknots, nlocalknots=nlocalknots)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=nbaseknots, nlocalknots=nlocalknots)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis


#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', ylim=0:1, cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend(0.8,.2, c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))
rect(xleft= k[4], xright= k[5], ybottom=-0.5, ytop=1.5, density=2, border=0, col='gray')





###############################################################################################
## 2. SIMULATED DATA WITH n=2000
###############################################################################################


## 2.1 ANALYSIS BASED ON 20 KNOTS, 15 LOCAL KNOTS. LARGER n ##
##############################################################

nbaseknots= 20; nlocalknots= 15

#Simulate data with n=2000
set.seed(1)
n= 2000
x1= rep(0:1,c(n/2,n/2))
x2= x1 + rnorm(n)
x= cbind(x1,x2)
z= rep(seq(-3,3,length=n/2), 2)
m= truemean(x,z)
y= truemean(x,z) + rnorm(n, 0, .25)
znew= matrix(rep(seq(-2.9,2.9,length=50),2),ncol=1)
xnew= cbind(rep(0:1,each=nrow(znew)/2), mean(x2))
newdata= list(x=xnew, z=znew)

#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis



#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend('topleft', c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))
rect(xleft= k[8], xright= k[9], ybottom=-0.5, ytop=1.5, density=3, border=0, col='gray')




## 2.2 ANALYSIS BASED ON 30 KNOTS, 30 LOCAL KNOTS ##
####################################################

nbaseknots= 30; nlocalknots= 30

#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=nbaseknots, nlocalknots=nlocalknots, localgridsize=500)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis


#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend('topleft', c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))
rect(xleft= k[15:16], xright= k[16:17], ybottom=-0.5, ytop=1.5, density=3, border=0, col='gray')



## 2.3 ANALYSIS BASED ON 12 KNOTS, 8 LOCAL KNOTS, LARGER n=2000 ##
##################################################################

nbaseknots= 12; nlocalknots= 8
#Perform local test with cut splines and standard splines
fit0= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, nbaseknots=nbaseknots, nlocalknots=nlocalknots)
fit0.uncut= localnulltest(y, x=x[,1,drop=FALSE], z=z, cutdegree=3, usecutbasis=FALSE, nbaseknots=nbaseknots, nlocalknots=nlocalknots)
b0= coef(fit0)
b0.uncut= coef(fit0.uncut)
k= fit0$regionbounds[[1]][[1]] #knots of the local test basis


#Plot local test posterior probabilities
par(mar=c(4.2,4.2,.1,.1))
plot(b0.uncut[,'z1'], b0.uncut[,'margpp'], type='l', lty=1, xlab='z', ylab='Posterior probability of a local effect', ylim=0:1, cex.axis=1.5, cex.lab=1.5, lwd=3, col='gray')
lines(b0[,'z1'], b0[,'margpp'], lty=2)
points(k, rep(0, length(k)), pch=16, cex=2)
legend(0.5, 0.2, c('Cubic spline','Cut spline'), lty=1:2, cex=1.5, col=c('grey','black'), lwd=c(3,1))
rect(xleft= k[3:4], xright= k[4:5], ybottom=-0.5, ytop=1.5, density=3, border=0, col='gray')




