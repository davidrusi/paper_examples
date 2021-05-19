#######################################################################
##
## APPROXIMATE LAPLACE APPROXIMATION. R CODE TO REPRODUCE EXAMPLES
##
#######################################################################


## TABLE OF CONTENTS
##
## 1. LOGISTIC REGRESSION SIMULATIONS
##   1.1 FIXED p, INCREASING n
##   1.2 FIXED n, GROWING p
## 2. POISSON REGRESSION SIMULATIONS (IN SUPPLEMENTARY MATERIAL)
##   2.1 FIXED p, INCREASING n
##   2.2 FIXED n, GROWING p
## 3. ALA ACCURACY FOR NON-LOCAL PRIORS IN GAUSSIAN REGRESSION
##   3.1 n=50, correlated covariates
##   3.2 n=200, correlated covariates
##   3.3 n=500, correlated covariates
##   3.4 n=50, uncorrelated covariates
##   3.5 n=200, uncorrelated covariates
##   3.6 n=500, uncorrelated covariates
##   3.7 FIGURE IN MAIN PAPER SUMMARIZING SIMULATION RESULTS FROM SECTIONS 3.1-3.6
## 4. POVERTY LINE EXAMPLE
##   4.1 MAIN EFFECTS ONLY
##   4.2 MAIN EFFECTS + INTERACTIONS
## 5. SURVIVAL DATA EXAMPLES
##   5.1 SIMULATION SCENARIO 1 (TRULY AFT MODEL)
##   5.2 SIMULATION SCENARIO 5 (TRULY PH MODEL)
##   5.3 COLON CANCER DATA (IN SUPPLEMENTARY MATERIAL). ORIGINAL VARIABLES
##   5.4 COLON CANCER DATA (IN SUPPLEMENTARY MATERIAL). ORIGINAL VARIABLES + FAKE VARIABLES
##   5.5 COLON CANCER DATA (IN SUPPLEMENTARY MATERIAL). SUMMARIZE RESULTS
##   5.6 COLON CANCER DATA (IN SUPPLEMENTARY MATERIAL). CONCORDANCE INDEX
## 6. IMPORTANCE SAMPLING AND VARIABLE SCREENING
##   6.1 LOGISTIC REGRESSION
##   6.2 POISSON REGRESSION
## 7. GROUP CONSTRAINTS IN NON-LINEAR GAUSSIAN REGRESSION
## 8. GROUP CONSTRAINTS IN LINEAR GAUSSIAN REGRESSION WITH CATEGORICAL PREDICTORS (IN SUPPLEMENTARY MATERIAL)


# Load packages
library(mombf)
library(mvtnorm)


#######################################################################
## 1. LOGISTIC REGRESSION SIMULATIONS
#######################################################################

#Set priors
prCoef= groupzellnerprior(taustd=1)
prDelta= modelbbprior()

#Function to simulate data in logistic and poisson models
simdata= function(n,beta,rho,model) {
  Sigma= diag(length(beta)-1)
  Sigma[upper.tri(Sigma)]= Sigma[lower.tri(Sigma)]= rho
  X= cbind(rep(1,n), rmvnorm(n, sigma=Sigma))
  linpred= X %*% beta
  if (model== 'logistic') {
    y= rbinom(n, size=1, prob=1/(1+exp(-linpred)))
  } else if (model== 'poisson') {
    y= rpois(n, exp(linpred))
  } else { stop("Wrong model") }
  return(list(X=X,y=y))
}

## 1.1 FIXED p, INCREASING n ##
###############################

#Run simulations
rho=0.5; p=10; nseq= c(100,500,1000, 5000)
beta= c(2,rep(0,p-3),.5,1)
model='logistic'; nsims= 5
t.la.newton= t.ala= t.ala1= t.ala2= matrix(NA,nrow=nsims,ncol=length(nseq))
margpp.la= margpp.ala= margpp.ala1= margpp.ala2= matrix(NA,nrow=length(nseq),ncol=2)
colnames(margpp.la)= colnames(margpp.ala)= colnames(margpp.ala1)= colnames(margpp.ala2)= c('Truly zero','Truly non-zero')
for (i in 1:length(nseq)) {
    tmp.la= tmp.ala= tmp.ala1= tmp.ala2= c(0,0)
    for (j in 1:nsims) {
        data= simdata(n=nseq[i],beta=beta,rho=rho,model=model)
        #
        t.la.newton[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optimMethod='Newton', enumerate=TRUE, verbose=FALSE))['elapsed']
        tmp.la= tmp.la + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        #
        t.ala[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='ALA', enumerate=TRUE, verbose=FALSE))['elapsed']
        tmp.ala= tmp.ala + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        t.ala1[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optimMethod='Newton', optim_maxit=1, enumerate=TRUE, verbose=FALSE))['elapsed']
        tmp.ala1= tmp.ala1 + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        t.ala2[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optimMethod='Newton', optim_maxit=2, enumerate=TRUE, verbose=FALSE))['elapsed']
        tmp.ala2= tmp.ala2 + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        cat(".")
    }
    margpp.la[i,]= tmp.la/nsims
    margpp.ala[i,]= tmp.ala/nsims
    margpp.ala1[i,]= tmp.ala1/nsims
    margpp.ala2[i,]= tmp.ala2/nsims
    cat("\n")
}

meantime= cbind(LA=colMeans(t.la.newton), ALA=colMeans(t.ala), ALA1=colMeans(t.ala1), ALA2=colMeans(t.ala2))


#CPU time figure
plot(nseq, meantime[,'LA'], type='l', lty=2, ylab='Mean CPU time (seconds)', ylim= range(meantime), xlab='n', cex.lab=1.5, cex.axis=1.4)
points(nseq, meantime[,'LA'])
lines(nseq, meantime[,'ALA']); points(nseq, meantime[,'ALA'], pch=16)
lines(nseq, meantime[,'ALA1'], col='gray'); points(nseq, meantime[,'ALA1'], pch=16, col='gray')
lines(nseq, meantime[,'ALA2'], col='gray', lty=2); points(nseq, meantime[,'ALA2'], col='gray')
legend('topleft', c('Laplace','ALA - 2 iter','ALA - 1 iter','ALA'), lty=c(2,2,1,1), pch=c(1,1,16,16), col=c('black','grey','grey','black'), cex=1.5)


#Posterior inclusion probabilities figure
plot(nseq, margpp.la[,'Truly zero'], type='l', lty=2, ylim=c(.03,1), ylab='Marginal posterior inclusion probability', xlab='n', cex.lab=1.5, cex.axis=1.4, log='y')
points(nseq, margpp.la[,'Truly zero'])
lines(nseq, margpp.la[,'Truly non-zero'], lty=2)
points(nseq, margpp.la[,'Truly non-zero'])
text(nseq[length(nseq)], margpp.la[length(nseq),1]+.01, 'LA, truly zero', pos=2, cex=1.3)
text(nseq[1], margpp.la[1,2]+.07, 'LA, truly non-zero', pos=4, cex=1.3)
#
lines(nseq, margpp.ala2[,'Truly zero'], col='darkgray',lty=2)
points(nseq, margpp.ala2[,'Truly zero'], col='darkgray')
lines(nseq, margpp.ala2[,'Truly non-zero'], col='darkgray',lty=2)
points(nseq, margpp.ala2[,'Truly non-zero'], col='darkgray')
text(nseq[length(nseq)], margpp.ala2[length(nseq),1], 'ALA 2 iter, truly zero', pos=2, cex=1.3)
text(nseq[1], margpp.ala2[1,2], 'ALA 2 iter, truly non-zero', pos=4, cex=1.3)
#
lines(nseq, margpp.ala1[,'Truly zero'], col='darkgray')
points(nseq, margpp.ala1[,'Truly zero'], col='darkgray', pch=16)
lines(nseq, margpp.ala1[,'Truly non-zero'], col='darkgray')
points(nseq, margpp.ala1[,'Truly non-zero'], col='darkgray', pch=16)
text(nseq[length(nseq)], margpp.ala1[length(nseq),1], 'ALA 1 iter, truly zero', pos=2, cex=1.3)
text(nseq[1], margpp.ala1[1,2]-.03, 'ALA 1 iter, truly non-zero', pos=4, cex=1.3)
#
lines(nseq, margpp.ala[,'Truly zero'], col='black', )
points(nseq, margpp.ala[,'Truly zero'], col='black', pch=16)
lines(nseq, margpp.ala[,'Truly non-zero'], col='black')
points(nseq, margpp.ala[,'Truly non-zero'], col='black', pch=16)
text(nseq[length(nseq)], margpp.ala[length(nseq),1], 'ALA, truly zero', pos=2, cex=1.3)
text(nseq[1], margpp.ala[1,2]-.03, 'ALA, truly non-zero', pos=4, cex=1.3)



## 1.2 FIXED n, GROWING p ##
############################


#Run simulations
rho=0.5; n=500; pseq= c(5, 10, 25, 50)
model= 'logistic'; nsims= 10
t.la.newton= t.ala= t.ala1= t.ala2= matrix(NA,nrow=nsims,ncol=length(pseq))
margpp.la= margpp.ala= margpp.ala1= margpp.ala2= matrix(NA,nrow=length(pseq),ncol=2)
colnames(margpp.la)= colnames(margpp.ala)= colnames(margpp.ala1)= colnames(margpp.ala2)= c('Truly zero','Truly non-zero')
for (i in 1:length(pseq)) {
    tmp.la= tmp.ala= tmp.ala1= tmp.ala2= c(0,0)
    for (j in 1:nsims) {
        beta= c(rep(0,pseq[i]-2),.5,1)
        data= simdata(n=n,beta=beta,rho=rho,model=model)
        #
        t.la.newton[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta= prDelta, method='Laplace', optimMethod='Newton', verbose=FALSE))['elapsed']
        tmp.la= tmp.la + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        #
        t.ala[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta= prDelta, method='ALA', verbose=FALSE))['elapsed']
        tmp.ala= tmp.ala + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        t.ala1[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta= prDelta, method='Laplace', optimMethod='Newton', optim_maxit=1, verbose=FALSE))['elapsed']
        tmp.ala1= tmp.ala1 + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        t.ala2[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta= prDelta, method='Laplace', optimMethod='Newton', optim_maxit=2, verbose=FALSE))['elapsed']
        tmp.ala2= tmp.ala2 + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        cat(".")
    }
    margpp.la[i,]= tmp.la/nsims
    margpp.ala[i,]= tmp.ala/nsims
    margpp.ala1[i,]= tmp.ala1/nsims
    margpp.ala2[i,]= tmp.ala2/nsims
    cat("\n")
}


meantime= cbind(LA=colMeans(t.la.newton), ALA=colMeans(t.ala), ALA1= colMeans(t.ala1), ALA2=colMeans(t.ala2))
mintime= apply(meantime[,1:2], 1, 'min')

#CPU time figure
plot(pseq, meantime[,'LA'], type='l', lty=2, ylab='Mean CPU time (seconds)', ylim= range(meantime), xlab='p', cex.lab=1.5, cex.axis=1.4)
points(pseq, meantime[,'LA'])
lines(pseq, meantime[,'ALA']); points(pseq, meantime[,'ALA'], pch=16)
lines(pseq, meantime[,'ALA1'], col='gray'); points(pseq, meantime[,'ALA1'], pch=16, col='gray')
lines(pseq, meantime[,'ALA2'], col='gray', lty=2); points(pseq, meantime[,'ALA2'], col='gray')
legend('topleft', c('Laplace','ALA - 2 iter','ALA - 1 iter','ALA'), lty=c(2,2,1,1), pch=c(1,1,16,16), col=c('black','grey','grey','black'), cex=1.5)


#Posterior inclusion probabilities figure
plot(pseq, margpp.la[,'Truly zero'], type='l', lty=2, ylim=c(.01,1), ylab='Marginal posterior inclusion probability', xlab='p', cex.lab=1.5, cex.axis=1.4, log='y')
points(pseq, margpp.la[,'Truly zero'])
lines(pseq, margpp.la[,'Truly non-zero'], lty=2)
points(pseq, margpp.la[,'Truly non-zero'])
text(pseq[length(pseq)], margpp.la[length(pseq),1]+.01, 'LA, truly zero', pos=2, cex=1.3)
text(pseq[1], margpp.la[1,2]+.1, 'LA, truly non-zero', pos=4, cex=1.3)
#
lines(pseq, margpp.ala2[,'Truly zero'], col='darkgray',lty=2)
points(pseq, margpp.ala2[,'Truly zero'], col='darkgray')
lines(pseq, margpp.ala2[,'Truly non-zero'], col='darkgray',lty=2)
points(pseq, margpp.ala2[,'Truly non-zero'], col='darkgray')
text(pseq[length(pseq)]+.01, margpp.ala2[length(pseq),1]+.005, 'ALA 2 iter, truly zero', pos=2, cex=1.3)
text(pseq[1], margpp.ala2[1,2]-.05, 'ALA 2 iter, truly non-zero', pos=4, cex=1.3)
#
lines(pseq, margpp.ala1[,'Truly zero'], col='darkgray')
points(pseq, margpp.ala1[,'Truly zero'], col='darkgray', pch=16)
lines(pseq, margpp.ala1[,'Truly non-zero'], col='darkgray')
points(pseq, margpp.ala1[,'Truly non-zero'], col='darkgray', pch=16)
text(pseq[length(pseq)], margpp.ala1[length(pseq),1]-0.001, 'ALA 1 iter, truly zero', pos=2, cex=1.3)
text(pseq[1], margpp.ala1[1,2]-.2, 'ALA 1 iter, truly non-zero', pos=4, cex=1.3)
#
lines(pseq, margpp.ala[,'Truly zero'], col='black', )
points(pseq, margpp.ala[,'Truly zero'], col='black', pch=16)
lines(pseq, margpp.ala[,'Truly non-zero'], col='black')
points(pseq, margpp.ala[,'Truly non-zero'], col='black', pch=16)
text(pseq[length(pseq)], margpp.ala[length(pseq),1], 'ALA, truly zero', pos=2, cex=1.3)
text(pseq[1], margpp.ala[1,2]-.3, 'ALA, truly non-zero', pos=4, cex=1.3)





#######################################################################
## 2. POISSON REGRESSION
#######################################################################

## 2.1 FIXED p, INCREASING n ##
###############################

#Run simulations
set.seed(1)
rho=0.5; p=10; nseq= c(100,500,1000,5000)
beta= c(rep(0,p-2),.5,1)
model='poisson'; nsims= 50
t.la.cda= t.la.newton= t.ala= matrix(NA,nrow=nsims,ncol=length(nseq))
margpp.la= margpp.ala= margpp.ala.unadj= matrix(NA,nrow=length(nseq),ncol=2)
colnames(margpp.la)= colnames(margpp.ala)= colnames(margpp.ala.unadj)= c('Truly zero','Truly non-zero')
for (i in 1:length(nseq)) {
    tmp.la= tmp.ala= tmp.ala.unadj= c(0,0)
    for (j in 1:nsims) {
        data= simdata(n=nseq[i],beta=beta,rho=rho,model=model)
        t.la.cda[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optimMethod='CDA', enumerate=TRUE, verbose=FALSE))['elapsed']
        #
        t.la.newton[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optimMethod='Newton', enumerate=TRUE, verbose=FALSE))['elapsed']
        tmp.la= tmp.la + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        #
        t.ala[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta=prDelta, method='ALA', enumerate=TRUE, verbose=FALSE))['elapsed']
        tmp.ala= tmp.ala + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        #
        ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta=prDelta, method='ALA', adj.overdisp='none', enumerate=TRUE, verbose=FALSE)
        tmp.ala.unadj= tmp.ala.unadj + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        cat(".")
    }
    margpp.la[i,]= tmp.la/nsims
    margpp.ala[i,]= tmp.ala/nsims
    margpp.ala.unadj[i,]= tmp.ala.unadj/nsims
    cat("\n")
}

meantime= cbind(LA.CDA=colMeans(t.la.cda), LA.Newton=colMeans(t.la.newton), ALA=colMeans(t.ala))
mintime= apply(meantime[,1:2], 1, 'min')

#Run times figure
plot(nseq, mintime, type='l', lty=2, ylab='Mean CPU time (seconds)', ylim= range(c(mintime,meantime[,'ALA'])), xlab='n', cex.lab=1.5, cex.axis=1.4)
points(nseq, mintime)
lines(nseq, meantime[,'ALA'])
points(nseq, meantime[,'ALA'], pch=16)
legend('topleft', c('Laplace','ALA'), lty=c(2,1), pch=c(1,16), cex=1.5)

#Posterior inclusion probabilities figure
plot(nseq, margpp.la[,'Truly zero'], type='l', lty=2, ylim=0:1, ylab='Marginal posterior inclusion probability', xlab='n', cex.lab=1.5, cex.axis=1.4)
points(nseq, margpp.la[,'Truly zero'])
lines(nseq, margpp.la[,'Truly non-zero'], lty=2)
points(nseq, margpp.la[,'Truly non-zero'])
text(nseq[1], margpp.la[1,], c('LA, truly zero','LA, truly non-zero'), pos=4, cex=1.3)
#
lines(nseq, margpp.ala[,'Truly zero'], col='darkgray')
points(nseq, margpp.ala[,'Truly zero'], col='darkgray')
lines(nseq, margpp.ala[,'Truly non-zero'], col='darkgray')
points(nseq, margpp.ala[,'Truly non-zero'], col='darkgray')
text(nseq[1], margpp.ala[1,], c('ALA adjusted, truly zero','ALA adjusted, truly non-zero'), pos=4, cex=1.3)
#
lines(nseq, margpp.ala.unadj[,'Truly zero'], lty=3, lwd=2, col='darkgray')
points(nseq, margpp.ala.unadj[,'Truly zero'], col='darkgray')
lines(nseq, margpp.ala.unadj[,'Truly non-zero'], lty=3, lwd=2, col='darkgray')
points(nseq, margpp.ala.unadj[,'Truly non-zero'], col='darkgray')
text(nseq[1], margpp.ala.unadj[1,]-.02, c('ALA unadjusted, truly zero','ALA unadjusted, truly non-zero'), pos=4, cex=1.3)



## 2.2 FIXED n, GROWING p ##
############################

#Run simulations
set.seed(1)
rho=0.5; n=500; pseq= c(5, 10, 25, 50)
model= 'poisson'; nsims= 50
t.la.cda = t.la.newton= t.ala= matrix(NA,nrow=nsims,ncol=length(pseq))
margpp.la= margpp.ala= margpp.ala.unadj= matrix(NA,nrow=length(pseq),ncol=2)
colnames(margpp.la)= colnames(margpp.ala)= colnames(margpp.ala.unadj)= c('Truly zero','Truly non-zero')
for (i in 1:length(pseq)) {
    tmp.la= tmp.ala= tmp.ala.unadj= c(0,0)
    for (j in 1:nsims) {
        beta= c(rep(0,pseq[i]-2),.5,1)
        data= simdata(n=n,beta=beta,rho=rho,model=model)
        t.la.cda[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta= prDelta, method='Laplace', optimMethod='CDA', verbose=FALSE))['elapsed']
        #
        t.la.newton[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta= prDelta, method='Laplace', optimMethod='Newton', verbose=FALSE))['elapsed']
        tmp.la= tmp.la + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        #
        t.ala[j,i]= system.time(ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta= prDelta, method='ALA', verbose=FALSE))['elapsed']
        tmp.ala= tmp.ala + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        #
        ms <- modelSelection(y=data$y, x=data$X, family='poisson', priorCoef=prCoef, priorDelta=prDelta, method='ALA', adj.overdisp='none', verbose=FALSE)
        tmp.ala.unadj= tmp.ala.unadj + c(mean(ms$margpp[-1][ beta[-1] == 0]), mean(ms$margpp[-1][ beta[-1] != 0]))
        cat(".")
    }
    margpp.la[i,]= tmp.la/nsims
    margpp.ala[i,]= tmp.ala/nsims
    margpp.ala.unadj[i,]= tmp.ala.unadj/nsims
    cat("\n")
}


meantime= cbind(LA.CDA=colMeans(t.la.cda), LA.Newton=colMeans(t.la.newton), ALA=colMeans(t.ala))
mintime= apply(meantime[,1:2], 1, 'min')

#Run times figure
plot(pseq, mintime, type='l', lty=2, ylab='Mean CPU time (seconds)', xlab='p', cex.lab=1.5, cex.axis=1.4)
points(pseq, mintime)
lines(pseq, meantime[,'ALA'])
points(pseq, meantime[,'ALA'], pch=16)
legend('topleft', c('Laplace','ALA'), lty=c(2,1), pch=c(1,16), cex=1.5)


#Posterior inclusion probabilities figure
sel= rep(TRUE,length(pseq))
plot(pseq[sel], margpp.la[sel,'Truly zero'], type='l', lty=2, ylim=0:1, ylab='Marginal posterior inclusion probability', xlab='p', cex.lab=1.5, cex.axis=1.4)
points(pseq[sel], margpp.la[sel,'Truly zero'])
lines(pseq[sel], margpp.la[sel,'Truly non-zero'], lty=2)
points(pseq[sel], margpp.la[sel,'Truly non-zero'])
text(pseq[1], margpp.la[1,]+.02, c('LA, truly zero','LA, truly non-zero'), pos=4, cex=1.3)
#
lines(pseq[sel], margpp.ala[sel,'Truly zero'], col='darkgray')
points(pseq[sel], margpp.ala[sel,'Truly zero'], col='darkgray')
lines(pseq[sel], margpp.ala[sel,'Truly non-zero'], col='darkgray')
points(pseq[sel], margpp.ala[sel,'Truly non-zero'], col='darkgray')
text(pseq[1], margpp.ala[1,]-.02, c('ALA adjusted, truly zero','ALA adjusted, truly non-zero'), pos=4, cex=1.3)
#
lines(pseq[sel], margpp.ala.unadj[sel,'Truly zero'], lty=3, lwd=2, col='darkgray')
points(pseq[sel], margpp.ala.unadj[sel,'Truly zero'], col='darkgray')
lines(pseq[sel], margpp.ala.unadj[sel,'Truly non-zero'], lty=3, lwd=2, col='darkgray')
points(pseq[sel], margpp.ala.unadj[sel,'Truly non-zero'], col='darkgray')
text(pseq[1], margpp.ala.unadj[1,]-.01, c('ALA unadjusted, truly zero','ALA unadjusted, truly non-zero'), pos=4, cex=1.3)





###############################################################################################################
# 3. ALA ACCURACY FOR NON-LOCAL PRIORS IN GAUSSIAN REGRESSION
###############################################################################################################


simlm <- function(seed,n,p,theta,sigma,correlated) {
  # Simulate data from linear model
  set.seed(seed)
  z     <- rnorm(n,0,sigma)
  if (correlated) {
      b = matrix(rnorm(p*p,0,1),p,p)
      v = cov2cor(t(b)%*%b)
      x = mvrnorm(n,rep(0,p),v)
  } else {
      x = matrix(rnorm(n*p),nrow=n,ncol=p)
  }
  p1  <- length(theta)
  th1 <- c(theta,rep(0,(p - p1)))
  y   <- x%*%th1 + z
  return(list(y=y,x=x,th=th1))
}


# --------------------------------------------------
# 3.1 n=50, correlated covariates
# --------------------------------------------------

n= 50; p= 10
theta= c(0.4,0.6,1.2,0.8) #non-zero parameter values
priorCoef= momprior(tau=1/3); priorVar= igprior(alpha=.01, lambda=.01)

nsims= 100 #number of simulations
nvars= rep(1:p, nsims)
ml = matrix(NA, nrow=length(nvars), ncol=7)
colnames(ml)= c('nvars','Exact','LA','ALA','time.Exact','time.LA','time.ALA')
for(i in 1:length(nvars)){
    d  <- simlm(seed=i, n=n, p=p, theta=theta, sigma=1.0, correlated=TRUE)
    y= d$y; x= d$x
    ml[i,1]= nvars[i]
    ml[i,5]= system.time(ml[i,2] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='auto'))['elapsed'] # exact
    ml[i,6]= system.time(ml[i,3] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='Laplace'))['elapsed'] # LA
    ml[i,7]= system.time(ml[i,4] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='plugin'))['elapsed'] #ALA
    cat('.')
}


e= rbind(data.frame(method='LA',  nvars=ml[,'nvars'], error= ml[,'LA'] - ml[,'Exact']), data.frame(method='ALA', nvars=ml[,'nvars'], error= ml[,'ALA'] - ml[,'Exact']))
e$method= factor(e$method); e$nvars= factor(e$nvars)
save(ml, e, file="la_ala_error_n50_corr.RData")



# --------------------------------------------------
# 3.2 n=200, correlated covariates
# --------------------------------------------------

n= 200; p= 10
theta= c(0.4,0.6,1.2,0.8) #non-zero parameter values
priorCoef= momprior(tau=1/3); priorVar= igprior(alpha=.01, lambda=.01)

nsims= 100 #number of simulations
nvars= rep(1:p, nsims)
ml = matrix(NA, nrow=length(nvars), ncol=7)
colnames(ml)= c('nvars','Exact','LA','ALA','time.Exact','time.LA','time.ALA')
for(i in 1:length(nvars)){
    d  <- simlm(seed=i, n=n, p=p, theta=theta, sigma=1.0, correlated=TRUE)
    y= d$y; x= d$x
    ml[i,1]= nvars[i]
    ml[i,5]= system.time(ml[i,2] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='auto'))['elapsed'] # exact
    ml[i,6]= system.time(ml[i,3] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='Laplace'))['elapsed'] # LA
    ml[i,7]= system.time(ml[i,4] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='plugin'))['elapsed'] #ALA
    cat('.')
}

e= rbind(data.frame(method='LA',  nvars=ml[,'nvars'], error= ml[,'LA'] - ml[,'Exact']), data.frame(method='ALA', nvars=ml[,'nvars'], error= ml[,'ALA'] - ml[,'Exact']))
e$method= factor(e$method); e$nvars= factor(e$nvars)
save(ml, e, file="la_ala_error_n200_corr.RData")


#Supplementary figure. Error for n=200
ggplot(aes(y = error, x = nvars, colour = method, dodge=method), data = e) + theme_bw() + geom_boxplot() + coord_cartesian(ylim= c(-15,1)) +  scale_x_discrete(name="Number of covariates") + scale_y_continuous("Estimated - exact log marginal likelihood") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.position=c(0.1,0.1), legend.title=element_text(size=20), legend.text=element_text(size=20)) + scale_colour_manual(values=c("black", "grey"))



# --------------------------------------------------
# 3.3 n=500, correlated covariates
# --------------------------------------------------

n= 500; p= 10
theta= c(0.4,0.6,1.2,0.8) #non-zero parameter values
priorCoef= momprior(tau=1/3); priorVar= igprior(alpha=.01, lambda=.01)

nsims= 100 #number of simulations
nvars= rep(1:p, nsims)
ml = matrix(NA, nrow=length(nvars), ncol=7)
colnames(ml)= c('nvars','Exact','LA','ALA','time.Exact','time.LA','time.ALA')
for(i in 1:length(nvars)){
    d  <- simlm(seed=i, n=n, p=p, theta=theta, sigma=1.0, correlated=TRUE)
    y= d$y; x= d$x
    ml[i,1]= nvars[i]
    ml[i,5]= system.time(ml[i,2] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='auto'))['elapsed'] # exact
    ml[i,6]= system.time(ml[i,3] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='Laplace'))['elapsed'] # LA
    ml[i,7]= system.time(ml[i,4] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='plugin'))['elapsed'] #ALA
    cat('.')
}

e= rbind(data.frame(method='LA',  nvars=ml[,'nvars'], error= ml[,'LA'] - ml[,'Exact']), data.frame(method='ALA', nvars=ml[,'nvars'], error= ml[,'ALA'] - ml[,'Exact']))
e$method= factor(e$method); e$nvars= factor(e$nvars)
save(ml, e, file="la_ala_error_n500_corr.RData")

#Supplementary figure. Error for n=500
ggplot(aes(y = error, x = nvars, colour = method, dodge=method), data = e) + theme_bw() + geom_boxplot() + coord_cartesian(ylim= c(-15,1)) +  scale_x_discrete(name="Number of covariates") + scale_y_continuous("Estimated - exact log marginal likelihood") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.position=c(0.1,0.1), legend.title=element_text(size=20), legend.text=element_text(size=20)) + scale_colour_manual(values=c("black", "grey"))



# --------------------------------------------------
# 3.4 n=50, uncorrelated covariates
# --------------------------------------------------

n= 50; p= 10
theta= c(0.4,0.6,1.2,0.8) #non-zero parameter values
priorCoef= momprior(tau=1/3); priorVar= igprior(alpha=.01, lambda=.01)

nsims= 100 #number of simulations
nvars= rep(1:p, nsims)
ml = matrix(NA, nrow=length(nvars), ncol=7)
colnames(ml)= c('nvars','Exact','LA','ALA','time.Exact','time.LA','time.ALA')
for(i in 1:length(nvars)){
    d  <- simlm(seed=i, n=n, p=p, theta=theta, sigma=1.0, correlated=FALSE)
    y= d$y; x= d$x
    ml[i,1]= nvars[i]
    ml[i,5]= system.time(ml[i,2] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='auto'))['elapsed'] # exact
    ml[i,6]= system.time(ml[i,3] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='Laplace'))['elapsed'] # LA
    ml[i,7]= system.time(ml[i,4] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='plugin'))['elapsed'] #ALA
    cat('.')
}

e= rbind(data.frame(method='LA',  nvars=ml[,'nvars'], error= ml[,'LA'] - ml[,'Exact']), data.frame(method='ALA', nvars=ml[,'nvars'], error= ml[,'ALA'] - ml[,'Exact']))
e$method= factor(e$method); e$nvars= factor(e$nvars)
save(ml, e, file="la_ala_error_n50_uncorr.RData")


#Supplementary figure. Error for n=50
ggplot(aes(y = error, x = nvars, colour = method, dodge=method), data=e) + theme_bw() + geom_boxplot() + coord_cartesian(ylim= c(-15,1)) +  scale_x_discrete(name="Number of covariates") + scale_y_continuous("Estimated - exact log marginal likelihood") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.position=c(0.1,0.1), legend.title=element_text(size=20), legend.text=element_text(size=20)) + scale_colour_manual(values=c("black", "grey"))



# --------------------------------------------------
# 3.5 n=200, uncorrelated covariates
# --------------------------------------------------

n= 200; p= 10
theta= c(0.4,0.6,1.2,0.8) #non-zero parameter values
priorCoef= momprior(tau=1/3); priorVar= igprior(alpha=.01, lambda=.01)

nsims= 100 #number of simulations
nvars= rep(1:p, nsims)
ml = matrix(NA, nrow=length(nvars), ncol=7)
colnames(ml)= c('nvars','Exact','LA','ALA','time.Exact','time.LA','time.ALA')
for(i in 1:length(nvars)){
    d  <- simlm(seed=i, n=n, p=p, theta=theta, sigma=1.0, correlated=FALSE)
    y= d$y; x= d$x
    ml[i,1]= nvars[i]
    ml[i,5]= system.time(ml[i,2] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='auto'))['elapsed'] # exact
    ml[i,6]= system.time(ml[i,3] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='Laplace'))['elapsed'] # LA
    ml[i,7]= system.time(ml[i,4] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='plugin'))['elapsed'] #ALA
    cat('.')
}

e= rbind(data.frame(method='LA',  nvars=ml[,'nvars'], error= ml[,'LA'] - ml[,'Exact']), data.frame(method='ALA', nvars=ml[,'nvars'], error= ml[,'ALA'] - ml[,'Exact']))
e$method= factor(e$method); e$nvars= factor(e$nvars)
save(e, file="la_ala_error_n200_uncorr.RData")


#Supplementary figure. Error for n=200
ggplot(aes(y = error, x = nvars, colour = method, dodge=method), data = e) + theme_bw() + geom_boxplot() + coord_cartesian(ylim= c(-15,1)) +  scale_x_discrete(name="Number of covariates") + scale_y_continuous("Estimated - exact log marginal likelihood") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.position=c(0.1,0.1), legend.title=element_text(size=20), legend.text=element_text(size=20)) + scale_colour_manual(values=c("black", "grey"))



# --------------------------------------------------
# 3.6 n=500, uncorrelated covariates
# --------------------------------------------------

n= 500; p= 10
theta= c(0.4,0.6,1.2,0.8) #non-zero parameter values
priorCoef= momprior(tau=1/3); priorVar= igprior(alpha=.01, lambda=.01)

nsims= 100 #number of simulations
nvars= rep(1:p, nsims)
ml = matrix(NA, nrow=length(nvars), ncol=7)
colnames(ml)= c('nvars','Exact','LA','ALA','time.Exact','time.LA','time.ALA')
for(i in 1:length(nvars)){
    d  <- simlm(seed=i, n=n, p=p, theta=theta, sigma=1.0, correlated=FALSE)
    y= d$y; x= d$x
    ml[i,1]= nvars[i]
    ml[i,5]= system.time(ml[i,2] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='auto'))['elapsed'] # exact
    ml[i,6]= system.time(ml[i,3] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='Laplace'))['elapsed'] # LA
    ml[i,7]= system.time(ml[i,4] <- nlpMarginal(sel=1:nvars[i], y=y, x=x, priorCoef=priorCoef, method='plugin'))['elapsed'] #ALA
    cat('.')
}

e= rbind(data.frame(method='LA',  nvars=ml[,'nvars'], error= ml[,'LA'] - ml[,'Exact']), data.frame(method='ALA', nvars=ml[,'nvars'], error= ml[,'ALA'] - ml[,'Exact']))
e$method= factor(e$method); e$nvars= factor(e$nvars)
save(e, file="la_ala_error_n500_uncorr.RData")

#Supplementary figure. Error for n=500
ggplot(aes(y = error, x = nvars, colour = method, dodge=method), data = e) + theme_bw() + geom_boxplot() + coord_cartesian(ylim= c(-15,1)) +  scale_x_discrete(name="Number of covariates") + scale_y_continuous("Estimated - exact log marginal likelihood") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.position=c(0.1,0.1), legend.title=element_text(size=20), legend.text=element_text(size=20)) + scale_colour_manual(values=c("black", "grey"))



# --------------------------------------------------
# 3.7 FIGURE IN MAIN PAPER SUMMARIZING SIMULATION RESULTS FROM SECTIONS 3.1-3.6
# --------------------------------------------------

#Compute mean error
library(plyr)
load("la_ala_error_n50_uncorr.RData")
e50u= ddply(.data=e, .variables=c('method','nvars'), .fun=summarize, error=mean(error))
e50u$nvars= as.numeric(as.character(e50u$nvars))
load("la_ala_error_n50_corr.RData")
e50c= ddply(.data=e, .variables=c('method','nvars'), .fun=summarize, error=mean(error))
e50c$nvars= as.numeric(as.character(e50c$nvars))



par(mar=c(5,5,.1,.1))
sel= e50u$method=='LA'
plot(e50u$nvars[sel], e50u$error[sel], type='l', xlab="Number of covariates", ylab="Estimated - exact log marginal likelihood", cex.lab=1.4, cex.axis=1.4, xlim=c(0,11.5), ylim=c(-5.5,0), xaxt='n')
axis(1,at=1:10, cex.axis=1.4)
text(10, tail(e50u$error[sel],1), "LA, I", cex=1.25, pos=4)
lines(e50u$nvars[!sel], e50u$error[!sel], col='gray', lwd=2)
text(10, tail(e50u$error[!sel],1), "ALA, I", cex=1.25, pos=4)
#
lines(e50c$nvars[sel], e50c$error[sel], lty=2, lwd=2)
text(10, tail(e50c$error[sel],1), "LA, V", cex=1.25, pos=4)
lines(e50c$nvars[!sel], e50c$error[!sel], col='gray', lty=2, lwd=2)
text(10, tail(e50c$error[!sel],1), "ALA, V", cex=1.25, pos=4)
#
text(4.5, -5.5, "No spurious covariates", cex=1.25, pos=2)
text(4.5, -5.5, "Spurious covariates", cex=1.25, pos=4)
segments(x0=4.5, y0=-7, y1=-5.25); arrows(x0=4.5, x1=2.5, y0= -5.25); arrows(x0=4.5, x1=6.5, y0= -5.25)




#############################################################################################
# 4. POVERTY LINE EXAMPLE
#############################################################################################

library(mombf)
library(grplasso)
library(xtable)
load("salary_2010_2019.RData")

#SELECT CASES. AGED 18-65, SINGLE RACE, NON-MILITARY EMPLOYED 35-40H/WEEK
salary= salary[(salary$age>=18) & (salary$age<=65),]
salary= salary[salary$employment== "At work",]
salary= salary[!(salary$occ %in% c("unemployed/neverworked","military")), ]
salary= salary[(salary$hoursworked >=35) & (salary$hoursworked <=40), ]
salary= salary[salary$wkstat == "Full-time", ]
salary$occ= factor(salary$occ)

#SELECT VARIABLES NEEDED FOR THE ANALYSIS
xnames= c("female","hispanic","marital","edu","age","race","citizen","nativity","occ","firmsize","difficulty","classworker","movedstate","hoursworked")
data= salary[,c("poverty",xnames)]


#BASIC DESCRIPTIVE STATISTICS
apply(data[,xnames], 2, table)



# --------------------------------------------------
# 4.1. MAIN EFFECTS ONLY
# --------------------------------------------------

prCoef= groupzellnerprior(taustd=1); prDelta= modelbbprior()
ncol(model.matrix( ~ ., data=data[,-1])) #p=60

#Model selection via ALA and LA
system.time(ms.ala <- modelSelection(poverty ~ ., data=data, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='ALA'))
system.time(ms.la <- modelSelection(poverty ~ ., data=data, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='Laplace')) #takes long!

#Marginal posterior inclusion probabilities
cbind(ALA=ms.ala$margpp, LA=ms.la$margpp)
cor(ms.ala$margpp, ms.la$margpp) #>0.999


#Model selection via group LASSO
system.time(lambda <- lambdamax(poverty ~ ., data=data, model=LogReg())) #1 sec
lambda <- exp(seq(log(lambdamax), log(lambdamax/2^7), length=50))
system.time(glasso <- grplasso(poverty ~ ., data=data, model=LogReg(), lambda=lambda)) #41.7 sec

#Compute BIC for models in group LASSO path
logl= apply(glasso$fitted, 2, function(z) sum(dbinom(data$poverty, size=1, prob=z, log=TRUE)))
npars= colSums(glasso$coef != 0)
bic= -2*logl + npars * log(nrow(data))
bic= cbind(lambda,bic,logl,npars)
rownames(bic)= NULL

b.glasso= glasso$coef[,which.min(bic[,'bic'])] #coefficient estimates
table(b.glasso != 0) #number of selected variables


#Get P-values for each term based on the MLE fit under the full model
getPvalues= function(data, xnames, family='binomial') {
    pvalue= double(length(xnames))
    glmfit= glm(poverty ~ ., data=data, family=family)
    for (i in 1:length(xnames)) {
        datasel= data[,c('poverty',xnames[-i])]
        glmfit2= glm(poverty ~ ., data=datasel, family=family)
        pvalue[i]= anova(glmfit2,glmfit,test='Chisq')[2,'Pr(>Chi)']
        cat('.')
    }
    names(pvalue)= xnames
    return(pvalue)
}

pvals= getPvalues(data, xnames, family='binomial')
pvals



#get point estimates and 95% CI
sel.ala= ms.ala$margpp > 0.5
sel.ala[is.na(sel.ala)]= FALSE
xsel= ms.ala$xstd[,sel.ala]
glm1 <- glm(ms.ala$ystd ~ -1 + xsel, family='binomial')

tt= summary(glm1)$coef
ci= cbind(tt[,'Estimate'], tt[,'Estimate'] - 1.96 * tt[,'Std. Error'], tt[,'Estimate'] + 1.96 * tt[,'Std. Error'])
ci= data.frame(MLE=round(ci[,1],3),CI=paste("(",round(ci[,2],3),",",round(ci[,3],3),")",sep=""))
rownames(ci)= sub("xsel","",rownames(ci))
rownames(ci)= sub("classworker","class",rownames(ci))
rownames(ci)= sub("Government","Gov",rownames(ci))
round(ci, 3)



# --------------------------------------------------
# 4.2 MAIN EFFECTS + INTERACTIONS
# --------------------------------------------------

prCoef= groupzellnerprior(taustd=1); prDelta= modelbbprior()
data= salary[,c("poverty",xnames)]

#remove columns with zero sd (no observations in the interaction category). 
colsd= apply(z,2,'sd')
which(colsd==0) #edu:nativity; edu:occ; race:occ; citizen:nativity; nativity:occ; firmsize:classworker

#Create formula object with all main effects and pairwise interactions
varpairs= expand.grid(1:length(xnames),1:length(xnames))
varpairs= varpairs[varpairs[,1] < varpairs[,2],]
inter= paste(xnames[varpairs[,1]], xnames[varpairs[,2]], sep=":")
inter= inter[!(inter %in% c("edu:nativity","edu:occ","race:occ","citizen:nativity","nativity:occ","firmsize:classworker"))]
f= paste(c(xnames,inter),collapse="+")
f= as.formula(paste("poverty ~",f))


#Model selection via ALA and LA
system.time(ms.ala <- modelSelection(f, data=data, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='ALA', initSearch='none'))
system.time(ms.la <- modelSelection(f, data=data, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='Laplace', initSearch='none')) #takes long!

#Posterior marginal inclusion probabilities
margpp= cbind(group=ms.la$groups, LA=ms.la$margpp, ALA=ms.ala$margpp)
margpp[is.na(margpp[,'ALA']), 'ALA']= margpp[is.na(margpp[,'ALA']), 'ALA2']
margpp= unique(margpp)
margpp

colSums(margpp[,c('LA','ALA')]>0.5) #number of selected terms

cor(margpp[,'LA'],margpp[,'ALA']) #correlation between marginal posterior inclusion probabilities



#Model selection via group LASSO
library(grplasso)
system.time(lambdamax <- lambdamax(f, data=data, model=LogReg()))
lambda <- exp(seq(log(lambdamax), log(lambdamax/2^7), length=50))
system.time(glasso <- grplasso(f, data=data, model=LogReg(), lambda=lambda))

#Obtain BIC for all models in group LASSO path
logl= apply(glasso$fitted, 2, function(z) sum(dbinom(data$poverty, size=1, prob=z, log=TRUE)))
npars= colSums(glasso$coef != 0)
bic= -2*logl + npars * log(nrow(data))
bic= cbind(lambda,bic,logl,npars)
rownames(bic)= NULL

b.glasso= glasso$coef[,which.min(bic[,'bic'])] #coefficient estimates
table(b.glasso != 0) #number of selected variables



#Get point estimates and 95% CI from MLE fit for model selected by ALA
sel.ala= ms.ala$margpp > 0.5
sel.ala[is.na(sel.ala)]= FALSE
xsel= ms.ala$xstd[,sel.ala]
glm1 <- glm(ms.ala$ystd ~ -1 + xsel, family='binomial')

tt= summary(glm1)$coef
ci= cbind(tt[,'Estimate'], tt[,'Estimate'] - 1.96 * tt[,'Std. Error'], tt[,'Estimate'] + 1.96 * tt[,'Std. Error'])
ci= data.frame(MLE=round(ci[,1],3),CI=paste("(",round(ci[,2],3),",",round(ci[,3],3),")",sep=""))
rownames(ci)= sub("xsel","",rownames(ci))
rownames(ci)= sub("classworker","class",rownames(ci))
rownames(ci)= sub("Government","Gov",rownames(ci))
round(ci, 3)


#Hierarchical LASSO
library(hierNet)
fhier= formula(paste(c('~ -1',xnames),collapse="+"))

des= model.matrix(fhier, data=data)
system.time(hierpath <- hierNet.logistic.path(des,data$poverty)) #XX sec

print(hierpath)

system.time(hnetfit <- hierNet.logistic(x=des, y=data$poverty, lam=hierpath$lamhat)) #fails to run, unsufficient memory with 16Gb RAM



#Compare ALA & LA posterior probabilities for selected model subset
models= unique(rbind(unique(ms.ala$postSample), unique(ms.la$postSample)))==TRUE  #select models saved in ALA or LA MCMC
mssel.ala <- modelSelection(ms.ala$y, x=ms.ala$xstd, groups=ms.ala$groups, models=models, data=data, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='ALA')
mssel.la <- modelSelection(ms.ala$y, x=ms.ala$xstd, groups=ms.ala$groups, models=models, data=data, family='binomial', priorCoef=prCoef, priorDelta=prDelta, method='Laplace')

margppsel= data.frame(group=mssel.ala$group, ALA= mssel.ala$margpp, LA= mssel.la$margpp)
margppsel= unique(margppsel)
cor(margppsel[,2:3])  #0.972




#############################################################################################
# 5. SURVIVAL DATA EXAMPLES
#############################################################################################

library(mombf)
library(mvtnorm)
library(survival)
library(survcomp)
library(hgu133plus2.db)
library(parallel)
source('routines_AFT.R')


# --------------------------------------------------
# 5.1 SIMULATION SCENARIO 1 (TRULY AFT MODEL)
# --------------------------------------------------

beta= c(1/3,1/2,rep(0,48)); sigmae = 0.5; corr=0.5; censtimes1=.5; censtimes2= 1


## group MOM + group MOM prior ##
priorCoef= groupmomprior(taustd=0.192 * 3); priorGroup=groupmomprior(taustd= 0.192 * 3)

#n=100
sim1= vector("list",250); t1= double(length(sim1))
for (i in 1:length(sim1)) {
  t1[i]= system.time(sim1[[i]] <- simScenarioWithSpur(i,beta=beta,n=100,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA'))['elapsed']
}
save(sim1, t1, file='data/simAFT/sim1_gmomgmom_n100_p50.RData')

#n=200
f= function(i) simScenarioWithSpur(i,beta=beta,n=200,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim1= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim1, file='data/simAFT/sim1_gmomgmom_n200_p50.RData')

#n=500
sim1= vector("list",250); t1= double(length(sim1))
for (i in 1:length(sim1)) {
  t1[i]= system.time(sim1[[i]] <- simScenarioWithSpur(i,beta=beta,n=500,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA'))['elapsed']
}
save(sim1, t1, file='data/simAFT/sim1_gmomgmom_n500_p50.RData')


## group Zellner + group Zellner ##
priorCoef= groupzellnerprior(taustd=1); priorGroup=groupzellnerprior(taustd=1)

#n=100, ALA
sim1= vector("list",250); t1= double(length(sim1))
for (i in 1:length(sim1)) {
  t1[i]= system.time(sim1[[i]] <- simScenarioWithSpur(i,beta=beta,n=100,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA'))['elapsed']
}
save(sim1, t1, file='data/simAFT/sim1_gzellgzell_ala_n100_p50.RData')

#n=200, ALA
f= function(i) simScenarioWithSpur(i,beta=beta,n=200,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim1= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim1, file='data/simAFT/sim1_gzellgzell_ala_n200_p50.RData')

#n=500, ALA
sim1= vector("list",250); t1= double(length(sim1))
for (i in 1:length(sim1)) {
  t1[i]= system.time(sim1[[i]] <- simScenarioWithSpur(i,beta=beta,n=500,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA'))['elapsed']
}
save(sim1, t1, file='data/simAFT/sim1_gzellgzell_ala_n500_p50.RData')


#n=200, LA
f= function(i) simScenarioWithSpur(i,beta=beta,n=200,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='Laplace')
sim1= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim1, file='data/simAFT/sim1_gzellgzell_la_n200_p50.RData')

#n=500, LA
sim1= vector("list",250); t1= double(length(sim1))
for (i in 1:length(sim1)) {
  t1[i]= system.time(sim1[[i]] <- simScenarioWithSpur(i,beta=beta,n=500,sigmae=sigmae,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='Laplace'))['elapsed']
}
save(sim1, t1, file='data/simAFT/sim1_gzellgzell_la_n500_p50.RData')



## SUMMARIZE RESULTS. PROPORTION OF CORRECT MODEL SELECTIONS ##
load('data/simAFT/sim1_gmomgmom_n100_p50.RData')
tabs1.gmom.n100= data.frame(method='gMOM ALA',n=100,correctselmombf(sim1,p=50,active=1:2))
load('data/simAFT/sim1_gmomgmom_n500_p50.RData')
tabs1.gmom.n500= data.frame(method='gMOM ALA',n=500,correctselmombf(sim1,p=50,active=1:2))

load('data/simAFT/sim1_gzellgzell_ala_n100_p50.RData')
tabs1.gzell.ala.n100= data.frame(method='gZellner ALA',n=100,correctselmombf(sim1,p=50,active=1:2))
load('data/simAFT/sim1_gzellgzell_ala_n500_p50.RData')
tabs1.gzell.ala.n500= data.frame(method='gZellner ALA',n=500,correctselmombf(sim1,p=50,active=1:2))

load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/data/sim_scen12_gzell_n100_p50.RData')
tabs1.gzell.la.n100= data.frame(method='gZellner LA',n=100,correctselmombf(simout1,p=50,active=1:2,selcens=TRUE))
load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/data/sim_scen12_gzell_n500_p50.RData')
tabs1.gzell.la.n500= data.frame(method='gZellner LA',n=500,correctselmombf(simout1,p=50,active=1:2,selcens=TRUE))

load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/data/sim_aftlasso_scen12_n100_p50.RData')
tabs1.aftlasso.n100= data.frame(method='AFT LASSO',n=100,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))
load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/data/sim_aftlasso_scen12_n500_p50.RData')
tabs1.aftlasso.n500= data.frame(method='AFT LASSO',n=500,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))

load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/data/sim_coxlasso_scen12_n100_p50.RData')
tabs1.coxlasso.n100= data.frame(method='Cox LASSO',n=100,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))
load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/data/sim_coxlasso_scen12_n500_p50.RData')
tabs1.coxlasso.n500= data.frame(method='Cox LASSO',n=500,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))


tabs1= rbind(tabs1.gmom.n100,tabs1.gmom.n500,tabs1.gzell.ala.n100,tabs1.gzell.ala.n500,tabs1.gzell.la.n100,tabs1.gzell.la.n500,tabs1.aftlasso.n100,tabs1.aftlasso.n500,tabs1.coxlasso.n100,tabs1.coxlasso.n500)
tabs1


#Figure with proportion of correct model selections
par(mar=c(4,4.5,.3,.3))
tab= tabs1
col= c('black','black','black','gray','darkgray'); lty= c(1,3,2,1,3); pch= c(1,2,3,2,1)
sel= (tab$method=='gMOM ALA')
plot(tab$n[sel],tab[sel,'correct.model'],type='l',xlab='n',ylab='Proportion of correct model selections',cex.lab=1.3,cex.axis=1.3,ylim=c(0,1),col=col[1],lty=lty[1])
points(tab$n[sel],tab[sel,'correct.model'],col=col[1],pch=pch[1])
sel= (tab$method=='gZellner ALA'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[2],lty=lty[2]); points(tab$n[sel],tab[sel,'correct.model'],col=col[2],pch=pch[2])
sel= (tab$method=='gZellner LA'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[3],lty=lty[3]); points(tab$n[sel],tab[sel,'correct.model'],col=col[3],pch=pch[3])
sel= (tab$method=='AFT LASSO'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[4],lty=lty[4]); points(tab$n[sel],tab[sel,'correct.model'],col=col[4],pch=pch[4])
sel= (tab$method=='Cox LASSO'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[5],lty=lty[5]); points(tab$n[sel],tab[sel,'correct.model'],col=col[5],pch=pch[5])
legend('topleft',c('gMOM ALA','gZellner ALA','gZellner LA','AFT LASSO','Cox LASSO'),col=col,lty=lty,pch=pch,cex=1.3)



# --------------------------------------------------
# 5.2 SIMULATION SCENARIO 5 (TRULY PH MODEL)
# --------------------------------------------------

mu= 0; sigma= 0.5; parH= c(mu,sigma); beta= c(0.75,-1.25,rep(0,48)); betaH = c(0,0,rep(0,48)); corr=0.5; censtimes1=.55; censtimes2= 0.95


## group MOM prior ##
priorCoef= groupmomprior(taustd=0.192 * 3); priorGroup=groupmomprior(taustd= 0.192 * 3)

#n=100
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=100,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gmomgmom_n100_p50.RData')

#n=200
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=200,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gmomgmom_n200_p50.RData')

#n=500
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=500,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gmomgmom_n500_p50.RData')


## group Zellner prior ##
priorCoef= groupzellnerprior(taustd=1); priorGroup=groupzellnerprior(taustd=1)

#n=100
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=100,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gzellgzell_ala_n100_p50.RData')

#n=200
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=200,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gzellgzell_ala_n200_p50.RData')

#n=500
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=500,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='ALA')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gzellgzell_ala_n500_p50.RData')

#n=200, LA
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=200,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='Laplace')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gzellgzell_la_n200_p50.RData')

#n=500, LA
f= function(i) simScenarioWithSpur35(i,parH=c(mu,sigma),betaH=betaH,beta=beta,n=500,corr=corr,censtimes=censtimes1,priorCoef=priorCoef,priorGroup=priorGroup,method='Laplace')
sim5= mclapply(1:250, f, mc.cores=4, mc.preschedule= FALSE)
save(sim5, file='data/simAFT/sim5_gzellgzell_la_n500_p50.RData')


## SUMMARIZE RESULTS. PROPORTION OF CORRECT MODEL SELECTIONS ##
load('data/simAFT/sim5_gmomgmom_n100_p50.RData')
tab.gmom.n100= data.frame(method='gMOM ALA',n=100,correctselmombf(sim5,p=50,active=1:2))
load('data/simAFT/sim5_gmomgmom_n500_p50.RData')
tab.gmom.n500= data.frame(method='gMOM ALA',n=500,correctselmombf(sim5,p=50,active=1:2))

load('data/simAFT/sim5_gzellgzell_ala_n100_p50.RData')
tab.gzell.ala.n100= data.frame(method='gZellner ALA',n=100,correctselmombf(sim5,p=50,active=1:2))
load('data/simAFT/sim5_gzellgzell_ala_n500_p50.RData')
tab.gzell.ala.n500= data.frame(method='gZellner ALA',n=500,correctselmombf(sim5,p=50,active=1:2))

load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/Sim_Javier/PH7/sim_scen12_gzell_n100_p50.RData')
tab.gzell.la.n100= data.frame(method='gZellner LA',n=100,correctselmombf(simout1,p=50,active=1:2,selcens=TRUE))
load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/Sim_Javier/PH7/sim_scen12_gzell_n500_p50.RData')
tab.gzell.la.n500= data.frame(method='gZellner LA',n=500,correctselmombf(simout1,p=50,active=1:2,selcens=TRUE))

load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/Sim_Javier/PH7/sim_aftlasso_scen12_n100_p50.RData')
tab.aftlasso.n100= data.frame(method='AFT LASSO',n=100,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))
load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/Sim_Javier/PH7/sim_aftlasso_scen12_n500_p50.RData')
tab.aftlasso.n500= data.frame(method='AFT LASSO',n=500,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))

load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/Sim_Javier/PH7/sim_coxlasso_scen12_n100_p50.RData')
tab.coxlasso.n100= data.frame(method='Cox LASSO',n=100,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))
load('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/Sim_Javier/PH7/sim_coxlasso_scen12_n500_p50.RData')
tab.coxlasso.n500= data.frame(method='Cox LASSO',n=500,correctselLASSO(simout1,p=50,active=1:2,uncens=FALSE))


tab= rbind(tab.gmom.n100,tab.gmom.n500,tab.gzell.ala.n100,tab.gzell.ala.n500,tab.gzell.la.n100,tab.gzell.la.n500,tab.aftlasso.n100,tab.aftlasso.n500,tab.coxlasso.n100,tab.coxlasso.n500)
tab


#Figure with proportion of correct model selections
par(mar=c(4,4.5,.3,.3))
col= c('black','black','black','gray','darkgray'); lty= c(1,3,2,1,3); pch= c(1,2,3,2,1)
sel= (tab$method=='gMOM ALA')
plot(tab$n[sel],tab[sel,'correct.model'],type='l',xlab='n',ylab='Proportion of correct model selections',cex.lab=1.3,cex.axis=1.3,ylim=c(0,1),col=col[1],lty=lty[1])
points(tab$n[sel],tab[sel,'correct.model'],col=col[1],pch=pch[1])
sel= (tab$method=='gZellner ALA'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[2],lty=lty[2]); points(tab$n[sel],tab[sel,'correct.model'],col=col[2],pch=pch[2])
sel= (tab$method=='gZellner LA'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[3],lty=lty[3]); points(tab$n[sel],tab[sel,'correct.model'],col=col[3],pch=pch[3])
sel= (tab$method=='AFT LASSO'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[4],lty=lty[4]); points(tab$n[sel],tab[sel,'correct.model'],col=col[4],pch=pch[4])
sel= (tab$method=='Cox LASSO'); lines(tab$n[sel],tab[sel,'correct.model'],col=col[5],lty=lty[5]); points(tab$n[sel],tab[sel,'correct.model'],col=col[5],pch=pch[5])
legend('topleft',c('gMOM ALA','gZellner ALA','gZellner LA','AFT LASSO','Cox LASSO'),col=col,lty=lty,pch=pch,cex=1.3)



# --------------------------------------------------
# 5.3 COLON CANCER DATA. ORIGINAL VARIABLES
# --------------------------------------------------

load("colon_p173.RData") #load supplementary file with data
y= Surv(log(survdata$Timerecurrence), event=survdata$recurrence)
X= survdata[,-1:-3]; X$stage= factor(X$stage)
Xfake= X[,2:(1+pfake)] + rnorm(nrow(X)*pfake,0,sd=1)
Xdes= cbind(intercept=1, stage2= ifelse(X$stage==2,1,0), stage3= ifelse(X$stage==3,1,0), X[,-1])

f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))

priorCoef= groupmomprior(taustd=0.192 * 3); priorGroup=groupmomprior(taustd= 0.192 * 3); priorDelta= modelunifprior()
system.time(ms.gmom1 <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^3, method='ALA'))
system.time(ms.gmom2 <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^4, method='ALA'))
system.time(ms.gmom3 <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=2 * 10^4, method='ALA'))

priorCoef= groupzellnerprior(taustd=1); priorGroup=groupzellnerprior(taustd=1); priorDelta= modelunifprior()
system.time(ms.gzell1 <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^3, method='ALA'))
system.time(ms.gzell2 <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^4, method='ALA'))
system.time(ms.gzell3 <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=2 * 10^4, method='ALA'))

system.time(ms.gzell1.la <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^3, method='Laplace'))
system.time(ms.gzell2.la <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^4, method='Laplace'))
system.time(ms.gzell3.la <- modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=2 * 10^4, method='Laplace'))

ms.gzell1.la= ms.gzell2.la= ms.gzell3.la= NA

aa <- modelSelection(f, priorCoef=priorCoef, priorDelta= priorDelta, priorGroup=priorGroup, niter=0, method='ALA')
system.time(b.coxlasso <- fitCoxLASSOCV(y=aa$ystd, x=aa$xstd, maxit=20000)) #Cox LASSO

system.time(b.aftlasso <-  fitAFTLASSOCV(y=aa$ystd, aa$xstd[,-1])[1:ncol(aa$xstd)]) #AFT LASSO

save(ms.gmom1, ms.gmom2, ms.gmom3, ms.gzell1, ms.gzell2, ms.gzell3, ms.gzell1.la, ms.gzell2.la, ms.gzell3.la, b.coxlasso, b.aftlasso, file='ms_p173.RData')


# --------------------------------------------------
# 5.4 COLON CANCER DATA. ORIGINAL VARIABLES + FAKE VARIABLES
# --------------------------------------------------

load("data/colon/colon_p173.RData")
y= Surv(log(survdata$Timerecurrence), event=survdata$recurrence)
X= survdata[,-1:-3]; X$stage= factor(X$stage)

set.seed(1); pfake=50
Xfake= X[,2:(1+pfake)] + rnorm(nrow(X)*pfake,0,sd=1)

r= sapply(1:ncol(Xfake), function(i) cor(X[,i+1],Xfake[,i]))
mean(r) #average correlation between X and fake X's is 0.70

f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
Z= cbind(X, Xfake)
ffake= formula(paste('y ~ ',paste('Z[,',1:ncol(Z),']',sep='',collapse="+"),sep=''))


#Run Bayesian model selection under MOM ALA, Zellner ALA, Zellner LA
priorCoef= groupmomprior(taustd=0.192 * 3); priorGroup=groupmomprior(taustd= 0.192 * 3); priorDelta= modelunifprior()
system.time(ms.gmom1 <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^3, method='ALA'))
system.time(ms.gmom2 <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^4, method='ALA'))
system.time(ms.gmom3 <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=2 * 10^4, method='ALA'))

priorCoef= groupzellnerprior(taustd=1); priorGroup=groupzellnerprior(taustd=1); priorDelta= modelunifprior()
system.time(ms.gzell1 <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^3, method='ALA'))
system.time(ms.gzell2 <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^4, method='ALA'))
system.time(ms.gzell3 <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=2 * 10^4, method='ALA'))

system.time(ms.gzell1.la <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^3, method='Laplace'))
system.time(ms.gzell2.la <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=10^4, method='Laplace'))
system.time(ms.gzell3.la <- modelSelection(ffake, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=2 * 10^4, method='Laplace'))

ms.gzell1.la= ms.gzell2.la= ms.gzell3.la= NA


#Run Cox LASSO and AFT LASSO
aa <- modelSelection(ffake, priorCoef=priorCoef, priorDelta= priorDelta, priorGroup=priorGroup, niter=0, method='ALA')
system.time(b.coxlasso <- fitCoxLASSOCV(y=aa$ystd, x=aa$xstd, maxit=20000)) #Cox LASSO

system.time(b.aftlasso <-  fitAFTLASSOCV(y=aa$ystd, aa$xstd[,-1])[1:ncol(aa$xstd)]) #AFT LASSO

save(ms.gmom1, ms.gmom2, ms.gmom3, ms.gzell1, ms.gzell2, ms.gzell3, ms.gzell1.la, ms.gzell2.la, ms.gzell3.la, b.coxlasso, b.aftlasso, file='ms_fake_p173.RData')



# --------------------------------------------------
# 5.5 COLON CANCER DATA. SUMMARIZE RESULTS
# --------------------------------------------------

#Compare top model found with increasing MCMC iterations
head(postProb(ms.gmom1))
head(postProb(ms.gmom2))
head(postProb(ms.gmom3))

head(postProb(ms.gzell1))
head(postProb(ms.gzell2))
head(postProb(ms.gzell3))


#Compare marginal inclusion probabilities estimated by Rao-Blackwellization vs. re-normalizing posterior model probabilities
ms= ms.gmom2
plot(pp2margpp(ms)$margpp[-1], ms$margpp[-1])
cor(pp2margpp(ms)$margpp[-1], ms$margpp[-1])


#Selected variables
ms= modelSelection(y=y, x=Xdes, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, niter=0, method='ALA')
nn= colnames(ms$xstd)

pp= postProb(ms.gmom2)
sel.gmom= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))
nn.gmom= sub("X","",nn[sel.gmom]); nn.gmom= nn.gmom[nn.gmom != 'intercept']
unlist(mget(nn.gmom, hgu133plus2SYMBOL))     

pp= postProb(ms.gzell2)
sel.gzell= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))
nn.gzell= sub("X","",nn[sel.gzell]); nn.gzell= nn.gzell[nn.gzell != 'intercept']
unlist(mget(nn.gzell[-1:-2], hgu133plus2SYMBOL)) 

sel.coxlasso= which(b.coxlasso != 0)
nn.coxlasso= sub("X","",nn[sel.coxlasso]); nn.coxlasso= nn.coxlasso[nn.coxlasso != 'intercept']
unlist(mget(nn.coxlasso[-1], hgu133plus2SYMBOL)) 

sel.aftlasso= which(b.aftlasso != 0)
nn.aftlasso= sub("X","",nn[sel.aftlasso]); nn.aftlasso= nn.aftlasso[nn.aftlasso != 'intercept']
unlist(mget(nn.aftlasso, hgu133plus2SYMBOL)) 


sel.gmom
sel.gzell
sel.coxlasso
sel.aftlasso

#gMOM selected ESM1, GAS1, PDPN
#gZellner selected stage2, stage3, ESM1    FBXO32     KANK4    TAGLN3    IGFBP3     DACT1  CALB2     PGBD5   ANKRD44
#Cox LASSO selected stage3,  TRIB1 NET1      ESM1    FBXO32     KANK4     EFNB2      ETV6     FLT1    IGFBP3     CALB2 
#AFT LASSO selected  ESM1     ODAPH      RYBP    DNAJB5       LIF    RASL12 


#Number of selected variables in original data, original + fake data, and overlap
load('ms_p173.RData')
pp= postProb(ms.gmom2)
sel.gmom= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))

pp= postProb(ms.gzell2)
sel.gzell= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))

sel.coxlasso= which(b.coxlasso != 0)
sel.aftlasso= which(b.aftlasso != 0)

sel.gmom= sel.gmom[sel.gmom != 1]
sel.gzell= sel.gzell[sel.gzell != 1]
sel.coxlasso= sel.coxlasso[sel.coxlasso != 1]
sel.aftlasso= sel.aftlasso[sel.aftlasso != 1]


load('ms_fake_p173.RData')
pp= postProb(ms.gmom2)
sel.gmom.fake= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))

pp= postProb(ms.gzell2)
sel.gzell.fake= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))

sel.coxlasso.fake= which(b.coxlasso != 0)
sel.aftlasso.fake= which(b.aftlasso != 0)

sel.gmom.fake= sel.gmom.fake[sel.gmom.fake != 1]
sel.gzell.fake= sel.gzell.fake[sel.gzell.fake != 1]
sel.coxlasso.fake= sel.coxlasso.fake[sel.coxlasso.fake != 1]
sel.aftlasso.fake= sel.aftlasso.fake[sel.aftlasso.fake != 1]


tab= matrix(NA,nrow=4, ncol=3)
rownames(tab)= c('gMOM ALA', 'gZellner ALA', 'Cox LASSO', 'AFT LASSO')
colnames(tab)= c('Original','Original + fake','Overlap')
tab[1,]= c(length(sel.gmom), length(sel.gmom.fake), sum(sel.gmom %in% sel.gmom.fake))
tab[2,]= c(length(sel.gzell), length(sel.gzell.fake), sum(sel.gzell %in% sel.gzell.fake))
tab[3,]= c(length(sel.coxlasso), length(sel.coxlasso.fake), sum(sel.coxlasso %in% sel.coxlasso.fake))
tab[4,]= c(length(sel.aftlasso), length(sel.aftlasso.fake), sum(sel.aftlasso %in% sel.aftlasso.fake))
tab




# --------------------------------------------------
# 5.6 COLON CANCER DATA. CONCORDANCE INDEX
# --------------------------------------------------

load("colon_p173.RData")
y= Surv(log(survdata$Timerecurrence), event=survdata$recurrence)
X= survdata[,-1:-3]; X$stage= factor(X$stage)
f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
Xd= model.matrix(f)
nfolds= nrow(X)

#Within-sample concordance index
load('ms_p173.RData')
yexp= y; yexp[,1]= exp(yexp[,1])

pp= postProb(ms.gmom2)
sel.gmom= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))
idx=sel.gmom; tmp= survreg(yexp ~ -1 + Xd[,idx], dist="lognormal")
ypred.gmom= Xd[,idx,drop=FALSE] %*% matrix(coef(tmp),ncol=1)
ci.gmom= concordance.index(-ypred.gmom, surv.time=y[,1], surv.event=y[,2])$c.index

pp= postProb(ms.gzell2)
sel.gzell= as.numeric(unlist(strsplit(as.character(pp[1,'modelid']),split=',')))
idx=sel.gzell; tmp= survreg(yexp ~ -1 + Xd[,idx], dist="lognormal")
ypred.gzell= Xd[,idx,drop=FALSE] %*% matrix(coef(tmp),ncol=1)
ci.gzell= concordance.index(-ypred.gzell, surv.time=y[,1], surv.event=y[,2])$c.index

ypred.coxlasso= Xd %*% matrix(b.coxlasso,ncol=1)
ci.coxlasso= concordance.index(ypred.coxlasso, surv.time=y[,1], surv.event=y[,2])$c.index

ypred.aftlasso= Xd %*% matrix(b.aftlasso,ncol=1)
ci.aftlasso= concordance.index(-ypred.aftlasso, surv.time=y[,1], surv.event=y[,2])$c.index

round(c(ci.gmom, ci.gzell, ci.coxlasso, ci.aftlasso),2)


#Cross-validation
nfolds= 10

priorCoef= groupmomprior(taustd=0.192 * 3); priorGroup=groupmomprior(taustd= 0.192 * 3); priorDelta= modelunifprior()
pred.mom= cv.aftmom(y=y, X=Xd, nfolds=nfolds, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, method='ALA', niter=1000)

priorCoef= groupzellnerprior(taustd=1); priorGroup=groupzellnerprior(taustd=1); priorDelta= modelunifprior()
pred.gzell= cv.aftmom(y=y, X=Xd, nfolds=nfolds, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, method='ALA', niter=1000)

pred.coxlasso= cv.coxlasso(y=y, X=Xd[,-1], nfolds=nfolds)
pred.aftlasso= cv.aftlasso(y=y, X=Xd, nfolds=nfolds))

ci.mom= concordance.index(-pred.mom$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.gzell= concordance.index(-pred.gzell$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.coxlasso= concordance.index(pred.coxlasso$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.aftlasso= concordance.index(-pred.aftlasso$pred, surv.time=y[,1], surv.event=y[,2])$c.index

tab= rbind(c(ci.mom, mean(pred.mom$nsel)), c(ci.gzell, mean(pred.gzell$nsel)), c(ci.coxlasso, mean(pred.coxlasso$nsel)), c(ci.aftlasso, mean(pred.aftlasso$nsel)))
rownames(tab)= c('gMOM ALA','gZellner ALA','gZellner LA','Cox LASSO','AFT LASSO')
colnames(tab)= c('CI','nvars')
tab







#######################################################################
## 6. IMPORTANCE SAMPLING AND VARIABLE SCREENING
#######################################################################

library(mombf)
library(mvtnorm)
prCoef= groupzellnerprior(taustd=1)
prDelta= modelbbprior()

#Simulate data
simdata= function(n,beta,rho,model) {
  Sigma= diag(length(beta)-1)
  Sigma[upper.tri(Sigma)]= Sigma[lower.tri(Sigma)]= rho
  X= cbind(rep(1,n), rmvnorm(n, sigma=Sigma))
  linpred= X %*% beta
  if (model== 'logistic') {
    y= rbinom(n, size=1, prob=1/(1+exp(-linpred)))
  } else if (model== 'poisson') {
    y= rpois(n, exp(linpred))
  } else { stop("Wrong model") }
  return(list(X=X,y=y))
}

## 6.1 LOGISTIC REGRESSION ##
#############################

set.seed(1)
rho=0.5; p=10; n= 1000
beta= c(2,rep(0,p-3),.5,1)
model='logistic'; nsims= 5
data= simdata(n=n,beta=beta,rho=rho,model=model)

#Run model selection with ALA, ALA 1 Newton-Raphson iter, ALA 2 Newton-Raphson iter, LA
ms0 <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='ALA')
ms1 <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optim_maxit=1)
ms2 <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optim_maxit=2)
ms <- modelSelection(y=data$y, x=data$X, family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace')

#Supplementary table: marginal posterior probabilities
margpp= cbind(ALA= ms0$margpp, ALA1= ms1$margpp, ALA2= ms2$margpp, LA= ms$margpp)
margpp

#Marginal posterior probabilities when using ALA for screening
sel= c(TRUE, ms0$margpp[-1] > .5) #keep intercept
ms.screen <- modelSelection(y=data$y, x=data$X[,sel], family='binomial', priorCoef=prCoef, priorDelta=prDelta , method='Laplace')
pp.screen= postProb(ms.screen)
round(cbind(exact= ms$margpp[sel], exact.ALAscreened= ms.screen$margpp), 4)


#Supplementary figure: distribution of importance sampling weights
pp0= postProb(ms0)
pp1= postProb(ms1)
pp2= postProb(ms2)
pp= postProb(ms)
nn= rownames(pp)
pp= data.frame(modelid=pp[nn,'modelid'], ALA= pp0[nn,'pp'], ALA1= pp1[nn,'pp'], ALA2= pp2[nn,'pp'], LA= pp[nn,'pp'])
w= pp[,'LA']/pp[,2:4] #importance sampling weights
w[w >10^4] = 10^4
for (i in 1:ncol(w)) { w[is.nan(w[,i]),i]= 1 } #0 divided by 0


wc= apply(w, 2, function(z) table(cut(z, breaks=c(-.01,0,.1,.25,.5,1,2,4,10,10^4 - .001, 10^4))))
rownames(wc)[1]= '0'
rownames(wc)[nrow(wc)-1]= '>10'
rownames(wc)[nrow(wc)]= 'infinite'
barplot(t(wc), xlab='', ylab='Number of models', cex.lab=1.3, beside=TRUE, las=2, col=c('white','gray','black'), legend=TRUE, args.legend=list(x='topleft',cex=1.3,legend=c('ALA','ALA - 1 ter','ALA - 2 iter')))
title(xlab='Importance sampling weights', line=5, cex.lab=1.3)




## 6.2 POISSON REGRESSION ##
############################

set.seed(1)
rho=0.5; p=5
beta= c(0,rep(0,p-2),.5,1)
model='poisson'


#Run model selection with ALA, ALA 1 Newton-Raphson iter, ALA 2 Newton-Raphson iter, LA
data= simdata(n=1000,beta=beta,rho=rho,model=model)
x= data$X
x= cbind(x, x[,-1]^2)
includevars= c(TRUE,rep(FALSE,ncol(x)-1)) #include intercept
ms0 <- modelSelection(y=data$y, x=x, family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='ALA', includevars=includevars)
ms1 <- modelSelection(y=data$y, x=x, family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optim_maxit=1, includevars=includevars)
ms2 <- modelSelection(y=data$y, x=x, family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', optim_maxit=2, includevars=includevars)
ms <- modelSelection(y=data$y, x=x, family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', includevars=includevars)


#Supplementary table: marginal posterior probabilities
margpp= cbind(ALA= ms0$margpp, ALA1= ms1$margpp, ALA2= ms2$margpp, LA= ms$margpp)[-1,]
margpp

#Marginal posterior probabilities when using ALA for screening
sel= c(TRUE, ms0$margpp[-1] > .5) #keep intercept
includevars= c(TRUE, rep(FALSE, sum(sel)-1)) #include intercept
ms.screen= modelSelection(y=data$y, x=x[,sel], family='poisson', priorCoef=prCoef, priorDelta=prDelta , method='Laplace', includevars=includevars)
pp.screen= postProb(ms.screen)

round(cbind(exact= ms$margpp[sel], exact.ALAscreened= ms.screen$margpp), 4)


#Supplementary figure: distribution of importance sampling weights
pp0= postProb(ms0)
pp1= postProb(ms1)
pp2= postProb(ms2)
pp= postProb(ms)
nn= rownames(pp)
pp= data.frame(modelid=pp[nn,'modelid'], ALA= pp0[nn,'pp'], ALA1= pp1[nn,'pp'], ALA2= pp2[nn,'pp'], LA= pp[nn,'pp'])
w= pp[,'LA']/pp[,2:4] #importance sampling weights
w[w >10^4] = 10^4
for (i in 1:ncol(w)) { w[is.nan(w[,i]),i]= 1 } #0 divided by 0


wc= apply(w, 2, function(z) table(cut(z, breaks=c(-.01,0,.1,.25,.5,1,2,4,10,10^4 - .001, 10^4))))
rownames(wc)[1]= '0'
rownames(wc)[nrow(wc)-1]= '>10'
rownames(wc)[nrow(wc)]= 'infinite'
barplot(t(wc), xlab='', ylab='Number of models', cex.lab=1.3, beside=TRUE, las=2, col=c('white','gray','black'), legend=TRUE, args.legend=list(x='topleft',cex=1.3,legend=c('ALA','ALA - 1 ter','ALA - 2 iter')))
title(xlab='Importance sampling weights', line=5, cex.lab=1.3)


#######################################################################
## 7. GROUP CONSTRAINTS IN NON-LINEAR GAUSSIAN REGRESSION
#######################################################################

library("mombf")
library("grpreg")
library("dplyr")
library("ggplot2")
library("ggforce")

# --------------------------------------------------
# 7.1 RUN SIMULATIONS
# --------------------------------------------------

get_data <- function(n, p, seed, freq1, freq2) {
  set.seed(seed)
  sigma <- matrix(.5, nrow=p, ncol=p)
  diag(sigma) <- 1
  z <- rnorm(n)
  X <- get_corr_data(p, z, corr_mat=sigma)
  colnames(X) <- paste0("x", seq(ncol(X)))
  y <- cos(-freq1 * X[, 1] + .1) + cos(freq2*X[, 2]) + rnorm(n)
  return(data.frame(y, X))
}

get_first_spline_name <- function(xstd_names, x_var_names) {
    spline_names <- c()
    for (x_name in x_var_names) {
        spline_names <- c(spline_names, xstd_names[grepl(paste0(x_name, ".s"), xstd_names)][1])
    }
    return(spline_names)
}

n_vec <- c(50, 100, 200, 350, 500, 750, 1000)
n_rep <- 150
prior_names <- c(
  "mom_gmom",
  "normid_gzell"
)

p <- 25
n_knots <- 9
freq1 <- 1
freq2 <- 2.5
data0 <- get_data(500, p, 3, freq1, freq2)
x_names <- colnames(data0)[2:ncol(data0)]

form <- as.formula(paste("~", paste(x_names, collapse="+")))
form_y <- as.formula(paste("y~", paste(x_names, collapse="+")))

covariates_to_plot <- c(x_names, paste0(x_names, ".s"))

aux <- create_data.table_n_rep_prior(
  n_vec, n_rep, prior_names, covariates = covariates_to_plot
)
model_dt <- aux$model_dt
marg_dt <- aux$marg_dt

penalties <- c("grLasso", "grMCP", "grSCAD")

compare_fit_to_model <- function(fit, model_idxs) {
  p <- length(model_idxs)
  coef_bools <- !near(coef(fit)[-1], 0, tol=1e-3)
  if (p != length(coef_bools)) { print(fit);  stop("error in lasso real model check")}
  with_linear_terms <- all(coef_bools == model_idxs)
  no_linear_terms <- all(coef_bools[3:p] == model_idxs[3:p])
    return(as.integer(with_linear_terms | no_linear_terms))
}

n_col <- rep(n_vec, each=length(penalties)*n_rep)
rep_col <- rep(seq(n_rep), times=length(n_vec), each=length(penalties))
pen_col <- rep(penalties, times=length(n_vec)*n_rep)
lasso_dt <- data.table::data.table(
  n=n_col,
  penalty=pen_col,
  repetition=rep_col,
  model_top=0,
  elapsed_time=0
)

row_count <- 1
lasso_row <- 1
for (n in n_vec) {
    message(paste("n =", n))
    for (idx in 1:n_rep) {
        data <- get_data(n, p, idx, freq1, freq2)
        for (name in prior_names) {
            priors <- get_prior_pair(name, n)
            pCoef <- priors$pCoef; pGroup <- priors$pGroup
            log <- capture.output(
              time <- system.time(
                fit <- modelSelection(
                    form_y, smoothterms=form, data=data,
                    priorCoef=pCoef, priorGroup=pGroup,
                    priorDelta=modelbbprior(1,1),
                    includevars=c(1), nknots=n_knots
                )
              )["elapsed"]
            )
            covariate_names <- colnames(fit$xstd)
            bool_idxs <- grepl("x1.s", covariate_names) | grepl("x2.s", covariate_names)
            bool_idxs[c(1,2,3)] <- TRUE
            real_model_idxs <- paste0(which(bool_idxs), collapse=",")
            update_data.table(
                model_dt = model_dt,
                marg_dt = marg_dt,
                row = row_count,
                fit = fit,
                model_idx = real_model_idxs,
                covariates_out = covariates_to_plot,
                covariates_in = c(x_names, get_first_spline_name(covariate_names, x_names)),
                time=time
            )
            row_count <- row_count + 1
        }
        for (penalty in penalties) {
          time <- system.time(
            grpreg_fit <- cv.grpreg(fit$xstd[, -1], fit$ystd, fit$groups[-1], penalty=penalty)
          )["elapsed"]
          lasso_dt[lasso_row, "model_top"] <- compare_fit_to_model(grpreg_fit, bool_idxs[-1])
          lasso_dt[lasso_row, "elapsed_time"] <- time
          lasso_row <- lasso_row + 1
        }
    }
}

save_df <- cbind(marg_dt, subset(model_dt, select=c(model_top, model_prob, elapsed_time)))
readr::write_csv(save_df, paste0("results_p_",p,".csv"))

readr::write_csv(lasso_dt, paste0("lasso_p_",p,".csv"))


# --------------------------------------------------
# 7.2 PRODUCE FIGURE
# --------------------------------------------------

ps <- c(5, 25)
first <- TRUE
for (p in ps) {
  aux_i <- readr::read_csv(paste0("results_p_",p,".csv"))
  aux_i <- tibble::add_column(
      aux_i,
      p=rep(p, nrow(aux_i)),
      package=rep("mombf", nrow(aux_i)),
      .before=TRUE
  ) %>% rename(method=prior)
  lasso_i <- readr::read_csv(paste0("lasso_p_",p,".csv"))
  lasso_i <- tibble::add_column(
      lasso_i,
      p=rep(p, nrow(lasso_i)),
      package=rep("grpreg", nrow(lasso_i)),
      .before=TRUE
  ) %>% rename(method=penalty)
  if (first) {
      first <- FALSE
      aux <- bind_rows(aux_i, lasso_i)
  } else {
      aux <- bind_rows(aux, aux_i, lasso_i)
  }
}
aux$p = factor(aux$p)
aux$package = factor(aux$package, levels=c("mombf", "grpreg"))

marg_df = subset(aux, select=-c(model_top, model_prob, elapsed_time))
covariates_to_plot <- colnames(marg_df)
covariates_to_plot <- covariates_to_plot[6:length(covariates_to_plot)]
covariates_to_plot
marg_df <- tidyr::drop_na(tidyr::gather(marg_df, covariates_to_plot, key="covariate", value="margpp"), margpp)

model_df = subset(aux, select=setdiff(names(aux), covariates_to_plot))
model_df <- tidyr::drop_na(
    tidyr::gather(model_df, "model_top", "model_prob", "elapsed_time", key="quantity", value="value"),
    value
)

img_dir <- "/home/oriol/Public/bgam-ala/draft/images/"

par(mar=c(4,4.5,.3,.3))
col= c('black','black','gray','gray','gray'); lty= c(1,3,2,1,3); pch= c(1,2,3,2,1)
legend_labels <- c('gMOM ALA','gZellner ALA','grLASSO','grMCP','grSCAD')
sel_labels <- c("mom_gmom", "normid_gzell", "grLasso", "grMCP", "grSCAD")

dat <- model_df %>%
    filter(method %in% sel_labels) %>%
    group_by(quantity, n, method, package, p) %>%
    summarise(mean_summary = mean(value))

quantity_to_plot <- "elapsed_time"
filename_base <- "time"
ylim <- NULL
ylabel <- "Computation time (s)"
for (p in ps) {
    tab_aux <- dat %>% filter(quantity == !!quantity_to_plot & p == !!p)
    if (is.null(ylim)) {ylim_ <- c(0, max(tab_aux$mean_summary))} else {ylim_ <- ylim}
    pdf(paste0(img_dir, filename_base, "_nonlinear_p",p,".pdf"))
    # first method
    tab <- tab_aux %>% filter(method==sel_labels[1])
    plot(tab$n,tab$mean_summary,type='l',xlab='n',ylab=ylabel,cex.lab=1.3,cex.axis=1.3,ylim=ylim_,col=col[1],lty=lty[1])
    points(tab$n,tab$mean_summary,col=col[1],pch=pch[1])
    # loop over rest of methods
    for (i in 2:length(legend_labels)) {
      tab <- tab_aux %>% filter(method==sel_labels[i])
      lines(tab$n,tab$mean_summary,col=col[i],lty=lty[i])
      points(tab$n,tab$mean_summary,col=col[i],pch=pch[i])
    }
    legend('bottomright',legend_labels,col=col,lty=lty,pch=pch,cex=1.3)
    dev.off()
}


#######################################################################
## 8. GROUP CONSTRAINTS IN LINEAR GAUSSIAN REGRESSION WITH CATEGORICAL PREDICTORS (IN SUPPLEMENTARY MATERIAL)
#######################################################################

library("mombf")
library("grpreg")
library("ggplot2")
library("dplyr")
library("ggforce")


# --------------------------------------------------
# 8.1 RUN SIMULATIONS
# --------------------------------------------------

theme_set(
    theme_grey(base_size = 20)
)

colors <- c("#E69F00", "#56B4E9", "#D55E00", "#0072B2", "#000000", "#CC79A7", "#009E73", "#F0E442")

# run parameters
p <- 50
beta <- "small"
n_rep <- 150
n_vec <- c(50, 100, 200, 350, 500, 750, 1000)
penalties <- c("grLasso", "grMCP", "grSCAD")

get_data <- function(n, p, seed, beta) {
  set.seed(seed)
  corr_mat <- matrix(.5, nrow=p, ncol=p)
  diag(corr_mat) <- 1
  z <- rnorm(n)
  x11 <- as.integer(z < -1)
  x12 <- as.integer(z > 1)
  x1 <- rep(2, n)
  x1[which(x11 == 1)] <- 1
  x1[which(x12 == 1)] <- 3
  Xaux <- get_corr_data(p, z, corr_mat=corr_mat)
  X <- data.frame(x11, x12, Xaux[,2:p])
  var_names <- c("x1_2", "x1_3", paste0("x", seq(2, p)))
  colnames(X) <- var_names
  if (beta == "big") {
    y <- -.6 * x11 + .6 * x12 + .5 * Xaux[, 2] + rnorm(n)
  } else {
    y <- -.3 * x11 + .3 * x12 + .5 * Xaux[, 2] + rnorm(n)
  }

  return(list(y=y, X=X))
}

compare_fit_to_model <- function(fit, model_idxs) {
    return(as.integer(all((!near(coef(fit)[-1], 0, tol=1e-3)) == real_model_idxs)))
}

real_model_idxs <- c(rep(TRUE, 3), rep(FALSE, p-2))
data0 <- get_data(200, p, 1, beta)
covariate_names <- colnames(data0$X)
groups <- c("x1", "x1", covariate_names[3:length(covariate_names)])

fit <- cv.grpreg(data0$X, data0$y, groups, penalty="grLasso")

n_col <- rep(n_vec, each=length(penalties)*n_rep)
rep_col <- rep(seq(n_rep), times=length(n_vec), each=length(penalties))
penalty_col <- rep(penalties, times=length(n_vec)*n_rep)
model_dt <- data.table::data.table(
    n=n_col,
    penalty=penalty_col,
    repetition=rep_col,
    model_top=0,
    elapsed_time=0
)

row_count <- 1
for (n in n_vec) {
    message(paste("n =", n))
    for (idx in 1:n_rep) {
        data <- get_data(n, p, idx, beta)
        for (penalty in penalties) {
            time <- system.time(
                fit <- cv.grpreg(data$X, data$y, groups, penalty=penalty)
            )["elapsed"]
            model_dt[row_count, "model_top"] <- compare_fit_to_model(fit, real_model_idxs)
            model_dt[row_count, "elapsed_time"] <- time
            row_count <- row_count + 1
        }
    }
}

readr::write_csv(model_dt, paste0("lasso_p_",p,"_beta_",beta,".csv"))


# --------------------------------------------------
# 8.2 PRODUCE FIGURE
# --------------------------------------------------

ps <- c(5, 50)
betas <- c("small", "big")
beta_val <- list(small=0.3, big=0.6)
first <- TRUE
for (p in ps) {
    for (beta in betas) {
            aux_i <- readr::read_csv(paste0("results_p_",p,"_beta_", beta,".csv"))
            aux_i <- tibble::add_column(
                aux_i,
                p=rep(p, nrow(aux_i)),
                beta=rep(beta_val[[beta]], nrow(aux_i)),
                package=rep("mombf", nrow(aux_i)),
                .before=TRUE
            ) %>% rename(method=prior)
            lasso_i <- readr::read_csv(paste0("lasso_p_",p,"_beta_", beta,".csv"))
            lasso_i <- tibble::add_column(
                lasso_i,
                p=rep(p, nrow(lasso_i)),
                beta=rep(beta_val[[beta]], nrow(lasso_i)),
                package=rep("grpreg", nrow(lasso_i)),
                .before=TRUE
            ) %>% rename(method=penalty)
            if (first) {
                first <- FALSE
                aux <- bind_rows(aux_i, lasso_i)
            } else {
                aux <- bind_rows(aux, aux_i, lasso_i)
            }
    }
}
aux$p = factor(aux$p)
aux$beta = factor(aux$beta)
aux$package = factor(aux$package, levels=c("mombf", "grpreg"))

marg_df = subset(aux, select=-c(model_top, model_prob, elapsed_time))
covariates_to_plot <- colnames(marg_df)
covariates_to_plot <- covariates_to_plot[7:length(covariates_to_plot)]
marg_df <- tidyr::drop_na(tidyr::gather(marg_df, covariates_to_plot, key="covariate", value="margpp"), margpp)

model_df = subset(aux, select=setdiff(names(aux), covariates_to_plot))
model_df <- tidyr::drop_na(
    tidyr::gather(model_df, "model_top", "model_prob", "elapsed_time", key="quantity", value="value"),
    value
)

img_dir <- "/home/oriol/Public/bgam-ala/draft/images/"

par(mar=c(4,4.5,.3,.3))
col= c('black','black','gray','gray','gray'); lty= c(1,3,2,1,3); pch= c(1,2,3,2,1)
legend_labels <- c('gMOM ALA','gZellner ALA','grLASSO','grMCP','grSCAD')
sel_labels <- c("mom_gmom", "normid_gzell", "grLasso", "grMCP", "grSCAD")

dat <- model_df %>%
    filter(method %in% sel_labels) %>%
    group_by(quantity, n, method, package, p, beta) %>%
    summarise(mean_summary = mean(value))


quantity_to_plot <- "elapsed_time"
filename_base <- "time"
ylim <- NULL
ylabel <- "Computation time (s)"
for (p in ps) {
  for (beta_name in betas) {
    beta <- beta_val[[beta_name]]
    tab_aux <- dat %>% filter(quantity == !!quantity_to_plot & beta == !!beta & p == !!p)
    if (is.null(ylim)) {ylim_ <- c(0, max(tab_aux$mean_summary))} else {ylim_ <- ylim}
    pdf(paste0(img_dir, filename_base, "_discrete_p",p,"_beta_",beta,".pdf"))
    # first plot
    tab <- tab_aux %>% filter(method==sel_labels[1])
    plot(tab$n,tab$mean_summary,type='l',xlab='n',ylab=ylabel,cex.lab=1.3,cex.axis=1.3,ylim=ylim_,col=col[1],lty=lty[1])
    points(tab$n,tab$mean_summary,col=col[1],pch=pch[1])
    for (i in 2:length(legend_labels)) {
      tab <- tab_aux %>% filter(method==sel_labels[i])
      lines(tab$n,tab$mean_summary,col=col[i],lty=lty[i])
      points(tab$n,tab$mean_summary,col=col[i],pch=pch[i])
    }
    legend('bottomright',legend_labels,col=col,lty=lty,pch=pch,cex=1.3)
    dev.off()
  }
}
