####################################################################################################
##
## ANALYSIS OF COLON CANCER DATA
##
####################################################################################################

#Load data and R packages
library(mombf)
library(AdapEnetClass)
library(BVSNLP)
library(hgu133plus2.db)
library(mvtnorm)
library(survival)
library(survcomp)
library(xtable)
source('routines.R')
survdata= read.table("ftbrs_survdata.txt",header=TRUE) #data available at GitHub
y= Surv(log(survdata$Timerecurrence), event=survdata$recurrence)
X= survdata[,-1:-3]; X$stage= factor(X$stage)
priorGroup= groupzellnerprior(tau=nrow(X))


####################################################################################################
## 1. PERMUTATION EXERCISE
####################################################################################################

## PERMUTATION. ONLY STAGING AND TGFB AS COVARIATES ##
######################################################

B= 100
ppnull.aftmom= double(B); margpp.aftmom= matrix(NA,nrow=B,ncol=9)
ppnull.coxbvs= double(B); margpp.coxbvs= matrix(NA,nrow=B,ncol=9)
sel.aftmom= sel.coxbvs= sel.aftlasso= sel.coxlasso= matrix(0,nrow=B,ncol=9)
for (i in 1:B) {
    set.seed(i)
    yperm= y[sample(1:nrow(y),replace=TRUE),]
    #AFT pMOM
    msperm= modelSelection(yperm ~ X$stage + X$tgfb, smooth= ~ X$tgfb, priorCoef= momprior(0.192), priorDelta= modelbbprior(1,1), priorGroup=priorGroup, verbose=FALSE)
    pp.aftmom= postProb(msperm)
    ppnull.aftmom[i]= sum(pp.aftmom$pp[pp.aftmom$modelid=='' | pp.aftmom$modelid=='1'])
    margpp.aftmom[i,]= msperm$margpp
    sel= as.numeric(as.character(pp.aftmom$model[1]))
    if (length(sel)>0) sel.aftmom[i,sel]= 1
    #Cox iMOM
    fit= runCoxBVS(y=msperm$ystd, x=msperm$xstd[,-1], groups=msperm$groups[-1])
    ppnull.coxbvs[i]= sum(fit$ppmodels[fit$ppmodels$model == '',2])
    margpp.coxbvs[i,-1]= fit$margpp
    sel= as.numeric(strsplit(as.character(fit$topmodel$model[1]), split=',')[[1]])
    if (length(sel)>0) sel.coxbvs[i,sel+1]= 1
    #AFT LASSO
    fit= try(fitAFTLASSOCV(msperm$ystd,msperm$xstd[,-1]), silent=TRUE)
    if (class(fit) != "try-error") sel.aftlasso[i,fit[1:ncol(msperm$xstd)] !=0 ]= 1
    #Cox LASSO
    resp= cbind(time=exp(msperm$ystd[,1]),status=msperm$ystd[,2])
    cv.fit= try(cv.glmnet(msperm$xstd, y=resp, family="cox", maxit=10000, nfolds=10, alpha=1), silent=TRUE)
    fit= try(glmnet(msperm$xstd, y=resp, family = "cox", maxit=10000, alpha=1), silent=TRUE)
    bhat= as.double(coef(fit, s=cv.fit$lambda.min))
    sel.coxlasso[i, bhat!=0]= 1
    cat('.')
}
colnames(margpp.aftmom)= colnames(margpp.coxbvs)= colnames(msperm$xstd)
save(ppnull.aftmom, margpp.aftmom, ppnull.coxbvs, margpp.coxbvs, sel.aftmom, sel.coxbvs, sel.aftlasso, sel.coxlasso, file='colon_stagetgfb_perm.RData')


## PERMUTATION. ALL COVARIATES ##
#################################

B= 100
ppnull.aftmom= double(B); margpp.aftmom= matrix(NA,nrow=B,ncol=176)
ppnull.coxbvs= double(B); margpp.coxbvs= matrix(NA,nrow=B,ncol=176)
sel.aftmom= sel.coxbvs= sel.aftlasso= sel.coxlasso= matrix(0,nrow=B,ncol=176)
for (i in 1:B) {
    set.seed(i)
    yperm= y[sample(1:nrow(y),replace=TRUE),]
    #AFT pMOM
    f= formula(paste('yperm ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
    system.time(msperm <- modelSelection(f, priorCoef= momprior(0.192), priorDelta= modelbbprior(1,1), priorGroup=priorGroup, niter=5000, verbose=FALSE))
    pp.aftmom= postProb(msperm)
    ppnull.aftmom[i]= sum(pp.aftmom$pp[pp.aftmom$modelid=='' | pp.aftmom$modelid=='1'])
    margpp.aftmom[i,]= msperm$margpp
    sel= as.numeric(as.character(pp.aftmom$model[1]))
    if (length(sel)>0) sel.aftmom[i,sel]= 1
    #Cox iMOM
    fit= runCoxBVS(y=msperm$ystd, x=msperm$xstd[,-1], groups=msperm$groups[-1])
    ppnull.coxbvs[i]= sum(fit$ppmodels[fit$ppmodels$model == '',2])
    margpp.coxbvs[i,-1]= fit$margpp
    sel= as.numeric(strsplit(as.character(fit$topmodel$model[1]), split=',')[[1]])
    if (length(sel)>0) sel.coxbvs[i,sel+1]= 1
    #AFT LASSO
    fit= try(fitAFTLASSOCV(msperm$ystd,msperm$xstd[,-1]), silent=TRUE)
    if (class(fit) != "try-error") sel.aftlasso[i,fit[1:ncol(msperm$xstd)] !=0 ]= 1
    #Cox LASSO
    resp= cbind(time=exp(msperm$ystd[,1]),status=msperm$ystd[,2])
    cv.fit= try(cv.glmnet(msperm$xstd, y=resp, family="cox", maxit=10000, nfolds=10, alpha=1), silent=TRUE)
    fit= try(glmnet(msperm$xstd, y=resp, family = "cox", maxit=10000, alpha=1), silent=TRUE)
    bhat= as.double(coef(fit, s=cv.fit$lambda.min))
    sel.coxlasso[i, bhat!=0]= 1
    cat('.')
}
colnames(margpp.aftmom)= colnames(margpp.coxbvs)= colnames(msperm$xstd)
save(ppnull.aftmom, margpp.aftmom, ppnull.coxbvs, margpp.coxbvs, sel.aftmom, sel.coxbvs, sel.aftlasso, sel.coxlasso, file='colon_perm.RData')


## REPORT PERMUTATION RESULTS ##
################################

load('colon_stagetgfb_perm.RData')
mean(ppnull.aftmom)
mean(ppnull.coxbvs)

tab= matrix(0,nrow=4,ncol=11); colnames(tab)= c(0:9,'>=10')
z= table(rowSums(sel.aftmom[,-1])); tab[1,names(z)]= z
z= table(rowSums(sel.coxbvs[,-1])); tab[2,names(z)]= z
z= table(rowSums(sel.aftlasso[,-1])); tab[3,names(z)]= z
z= table(rowSums(sel.coxlasso[,-1])); tab[4,names(z)]= z
rownames(tab)= c('AFT-pMOMZ','Cox-piMOM','AFT-LASSO','Cox-LASSO')

par(mar=c(4.5,4.5,1,.1))
barplot(tab, beside=TRUE, xlab='Number of selected variables', ylab='Percentage of simulations', cex.lab=1.3, cex.names=1.3, cex.axis=1.3, legend=TRUE, args.legend= list(cex=1.25), col=gray(c(0,.3,.6,1)))


tab1= rbind(mean(sel.aftmom[,-1]), mean(sel.coxbvs[,-1]), mean(sel.aftlasso[,-1]>0,na.rm=TRUE), mean(sel.coxlasso[,-1]>0))
correctsel1= rbind(mean(sapply(1:nrow(sel.aftmom), function(i) all(sel.aftmom[i,-1]==0))), mean(sapply(1:nrow(sel.coxbvs), function(i) all(sel.coxbvs[i,-1]==0))), mean(sapply(1:nrow(sel.aftlasso), function(i) all(sel.aftlasso[i,-1]==0))), mean(sapply(1:nrow(sel.coxlasso), function(i) all(sel.coxlasso[i,-1]==0))))
tab1= cbind(tab1, correctsel1)
rownames(tab1)= c('AFT-pMOMZ','Cox-piMOM','AFT-LASSO','Cox-LASSO')
colnames(tab1)= c('False positives','% correct selections')


load('colon_perm.RData')
mean(ppnull.aftmom)
mean(ppnull.coxbvs)

tab= matrix(0,nrow=4,ncol=11); colnames(tab)= c(0:9,'>=10')
z= table(rowSums(sel.aftmom[,-1])); tab[1,names(z)]= z
z= table(rowSums(sel.coxbvs[,-1])); tab[2,names(z)]= z
z= table(rowSums(sel.aftlasso[,-1])); tab[3,]= c(z[1:10],sum(z[-1:-10]))
z= table(rowSums(sel.coxlasso[,-1])); tab[4,]= c(z[1:10],sum(z[-1:-10]))
rownames(tab)= c('AFT-pMOMZ','Cox-piMOM','AFT-LASSO','Cox-LASSO')

par(mar=c(4.5,4.5,1,.1))
barplot(tab, beside=TRUE, xlab='Number of selected variables', ylab='Percentage of simulations', cex.lab=1.3, cex.names=1.3, cex.axis=1.3, legend=TRUE, args.legend= list(cex=1.25), col=gray(c(0,.3,.6,1)))


tab2= rbind(mean(sel.aftmom[,-1]), mean(sel.coxbvs[,-1]), mean(sel.aftlasso[,-1]>0,na.rm=TRUE), mean(sel.coxlasso[,-1]>0))
correctsel2= rbind(mean(sapply(1:nrow(sel.aftmom), function(i) all(sel.aftmom[i,-1]==0))), mean(sapply(1:nrow(sel.coxbvs), function(i) all(sel.coxbvs[i,-1]==0))), mean(sapply(1:nrow(sel.aftlasso), function(i) all(sel.aftlasso[i,-1]==0))), mean(sapply(1:nrow(sel.coxlasso), function(i) all(sel.coxlasso[i,-1]==0))))
tab2= cbind(tab2, correctsel2)
rownames(tab2)= c('MOM','Cox iMOM','AFT LASSO','Cox LASSO')
colnames(tab2)= c('False positives','% correct selections')

tab= cbind(tab1,tab2)
xtable(100*tab, digits=c(0,1,1,1,1))




####################################################################################################
## 2. DATA ANALYSIS FOR MODEL ONLY WITH TGFB AND STAGE
####################################################################################################

#Run Bayesian variable selection
ms0= modelSelection(y ~ X$stage + X$tgfb, smooth= ~ X$tgfb, priorCoef= momprior(0.192), priorDelta= modelbbprior(1,1), priorGroup=priorGroup)
pp0= postProb(ms0)
pp0
ms0$margpp #TGFB gets 0.9977 inclusion prob for linear effect, 0.009036 for non-linear effect


#MLE for top model
yexp= y; yexp[,1]= exp(yexp[,1])
fittop= survreg(yexp ~ X[,1] + X[,2], dist = "lognormal")
summary(fittop)
ypred= predict(fittop,type='linear')
b= coef(fittop)

exp(-b[-1])


#Plot predicted survival
par(mar=c(4.5,4.5,.1,.1))
ylim= c(1/20,20); xlim= c(-2.5,2.5); x2plot= X$tgfb; o= order(x2plot)
des= cbind(1,X$stage==2,X$stage==3,X[,2])
sel= c(1,2,3)
y2plot= ypred - des[,sel] %*% matrix(b[sel],ncol=1)
plot(x2plot[o], exp(-y2plot[o]), xlab='TGFB expression (standardized)', ylab='Accelaration in time to recurrence', xlim=xlim, ylim=ylim, type='l',log='y', cex.axis=1.15, cex.lab=1.4)
y2plot= ypred - des[,sel] %*% matrix(b[sel],ncol=1) + b[2]
lines(x2plot[o], exp(-y2plot[o]), col=1, lty=2)
y2plot= ypred - des[,sel] %*% matrix(b[sel],ncol=1) + b[3]
lines(x2plot[o], exp(-y2plot[o]), col='gray', lwd=2)
legend('bottomright',c('Stage 1','Stage 2','Stage 3'),col=c(1,1,'gray'),lty=c(1,2,1),lwd=c(1,1,2),cex=1.4)





####################################################################################################
## 3. DATA ANALYSIS FOR MODEL WITH TGFB AND OTHER COVARIATES
####################################################################################################

## LINEAR EFFECTS ##
####################

f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
mslin= modelSelection(f, priorCoef= momprior(0.192), priorDelta= modelbbprior(1,1), priorGroup=priorGroup, enumerate=FALSE, niter=10000)

fit.coxbvs= runCoxBVS(y=mslin$ystd, x=mslin$xstd[,-1], groups=mslin$groups[-1]) #Cox iMOM
b.aftlasso= fitAFTLASSOCV(mslin$ystd,mslin$xstd[,-1])[1:ncol(mslin$xstd)] #AFT LASSO
resp= cbind(time=exp(mslin$ystd[,1]),status=mslin$ystd[,2]) #Cox LASSO
cv.fit= try(cv.glmnet(mslin$xstd, y=resp, family="cox", maxit=10000, nfolds=10, alpha=1), silent=TRUE)
fit= try(glmnet(mslin$xstd, y=resp, family = "cox", maxit=10000, alpha=1), silent=TRUE)
b.coxlasso= as.double(coef(fit, s=cv.fit$lambda.min))
save(mslin, fit.coxbvs, b.aftlasso, b.coxlasso, file='ftbrs_mslin.RData')

#Marginal posterior inclusion probabilities for Bayesian methods
margpp.aftmom= mslin$margpp
margpp.coxbvs= fit.coxbvs$margpp

#Selected variables by each method
pp.aftmom= postProb(mslin)
pp.coxbvs= fit.coxbvs$ppmodels; pp.coxbvs= pp.coxbvs[order(pp.coxbvs[,'pp'], decreasing=TRUE),]
sel.coxbvs= as.numeric(strsplit(as.character(fit.coxbvs$topmodel$model[1]), split=',')[[1]])
sel.aftlasso= which(b.aftlasso !=0)
sel.coxlasso= which(b.coxlasso !=0)


#Get gene symbols for genes selected by MOM-AFT
#Top model had post prob 0.110 and selected FLT1.
#Three top models (total 0.25 post prob) included stage, ESM1, FLT1 and GAS1 with marginal inclusion prob 0.329, 0.604, 0.276, 0.471.
#Five top models (total 0.297 post prob) included stage, ESM1, FBXO32, FLT1, EFNB2 and GAS1 with marginal inclusion prob 0.329, 0.604, 0.068, 0.276, 0.141, 0.471.
idx= which(margpp.aftmom>.1); data.frame(colnames(mslin$xstd)[idx], margpp= margpp.aftmom[idx])
nn= sub("X","",colnames(X)[c(13,14,43,29,57)])
unlist(mget(nn, hgu133plus2SYMBOL))

#The top model for Cox BVS had post prob 0.142 and selected stage and FLT1.
#Three top models (total 0.334 post prob) included stage, NET1 and FLT1 with marginal inclusion prob 0.573, 0.272, 0.884.
#Five top models (total 0.462 post prob) included stage, NET1, EFNB2 and FLT1 with marginal inclusion prob 0.573, 0.272, 0.192, 0.884.
idx= which(margpp.coxbvs>.15); data.frame(colnames(mslin$xstd[,-1])[idx], margpp= margpp.coxbvs[idx])
nn= sub("X","",colnames(X)[c(9,29,43)])
unlist(mget(nn, hgu133plus2SYMBOL))

#AFT LASSO selected stage and 6 genes. FBXO32, KANK4, GPR161, SETBP1, SNX30, DACT1. Fitting an AFT model via MLE on selected covars gave P-value<0.05 for FXO32, KANK4, SETBP1.
colnames(mslin$xstd[,-1])[which(b.aftlasso != 0)]
nn= sub("X","",colnames(X)[c(14,23,47,54,85,101)])
unlist(mget(nn, hgu133plus2SYMBOL))

yexp= mslin$ystd; yexp[,1]= exp(yexp[,1])
tmp= survreg(yexp ~ mslin$xstd[,1+which(b.aftlasso != 0)], dist = "lognormal")
summary(tmp)

#Cox LASSO selected stage and 10 genes. Columns 8, 9, 13, 14, 23, 29, 35, 43, 76, 142 in X
#Gene symbols: TRIB1 NET1 ESM1 FBXO32 KANK4 EFNB2 ETV6 FLT1 IGFBP3 CALB2
#Fitting a Cox model on the selected covariates gave P-value<0.05 for stage, NET1, ESM1, EFNB2, ETV6, FLT1. For P-value < 0.01 only gene ETV6 is selected.
colnames(mslin$xstd)[sel.coxlasso]
nn= sub("X","",colnames(X)[c(8,9,13,14,23,29,35,43,76,142)])
unlist(mget(nn, hgu133plus2SYMBOL))

resp= cbind(time=exp(mslin$ystd[,1]),status=mslin$ystd[,2])
tmp= coxph(Surv(resp) ~ mslin$xstd[,sel.coxlasso])





## NON-LINEAR EFFECTS ##
########################

f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
s= formula(paste('~ ',paste('X[,',2:ncol(X),']',sep='',collapse="+"),sep=''))
system.time(ms.smooth <- modelSelection(f, smooth=s, priorCoef= momprior(0.192), priorDelta= modelbbprior(1,1), priorGroup=priorGroup, enumerate=FALSE, niter=10000, nknots=8))
save(ms.smooth, file='ftbrs_mssmooth.RData')


#MLE with top variables
load('ftbrs_mslin.RData')
mslin= mslin; pplin= postProb(mslin)
head(pplin)
o= order(mslin$margpp,decreasing=TRUE)
head(mslin$margpp[o],n=8)

yexp= y; yexp[,1]= exp(yexp[,1])
fittop= survreg(yexp ~ X[,1] + X[,2] + X[,13] + X[,43] + X[,57], dist = "lognormal") #top vars + TGFB
summary(fittop)
ypred= predict(fittop,type='linear')
b= coef(fittop)


#Plot predicted survival
par(mar=c(4.5,4.5,.1,.1))
ylim= c(1/20,20); xlim= c(-2.5,2.5); x2plot= X$tgfb; o= order(x2plot)
des= cbind(1,X$stage==2,X$stage==3,X[,2],X[,13],X[,43],X[,57])
sel= c(1,2,3,5,6,7)
y2plot= ypred - des[,sel] %*% matrix(b[sel],ncol=1)
plot(x2plot[o], exp(-y2plot[o]), xlab='TGFB expression (standardized)', ylab='Accelaration in time to recurrence', xlim=xlim, ylim=ylim, type='l',log='y', cex.axis=1.15, cex.lab=1.4)
y2plot= ypred - des[,sel] %*% matrix(b[sel],ncol=1) + b[2]
lines(x2plot[o], exp(-y2plot[o]), col=1, lty=2)
y2plot= ypred - des[,sel] %*% matrix(b[sel],ncol=1) + b[3]
lines(x2plot[o], exp(-y2plot[o]), col='gray', lwd=2)
legend('bottomright',c('Stage 1','Stage 2','Stage 3'),col=c(1,1,'gray'),lty=c(1,2,1),lwd=c(1,1,2),cex=1.4)




####################################################################################################
## 4. CROSS-VALIDATION ANALYSIS OF CONCORDANCE INDEX
####################################################################################################

## 4.1 LINEAR TERMS ##
######################

f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
Xd= model.matrix(f)

nfolds= nrow(Xd) #LOO-CV
pred.aftmom= cv.aftmom(y=y, X=Xd, nfolds=nfolds, priorCoef=momprior(0.192), priorDelta=modelbbprior(1,1), priorGroup=groupzellnerprior(tau=nrow(Xd)), niter=1000) #identical results for niter=1000
pred.coxbvs= cv.coxbvs(y=y, X=Xd[,-1], nfolds=nfolds)
pred.coxlasso= cv.coxlasso(y=y, X=Xd[,-1], nfolds=nfolds)
pred.aftlasso= cv.aftlasso(y=y, X=Xd, nfolds=nfolds)

ci.aftmom= concordance.index(-pred.aftmom$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.coxbvs= concordance.index(pred.coxbvs$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.coxlasso= concordance.index(pred.coxlasso$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.aftlasso= concordance.index(-pred.aftlasso$pred, surv.time=y[,1], surv.event=y[,2])$c.index

tab= rbind(c(ci.aftmom, mean(pred.aftmom$nsel)), c(ci.coxbvs, mean(pred.coxbvs$nsel)), c(ci.coxlasso, mean(pred.coxlasso$nsel)), c(ci.aftlasso, mean(pred.aftlasso$nsel)))
rownames(tab)= c('AFT pMOMZ','Cox piMOM','Cox LASSO','AFT LASSO')
colnames(tab)= c('CI','nvars')
round(tab,2)



## 4.2 LINEAR AND NON-LINEAR TERMS ##
#####################################

#Create design matrix with spline terms
f= formula(paste('y ~ ',paste('X[,',1:ncol(X),']',sep='',collapse="+"),sep=''))
s= formula(paste('~ ',paste('X[,',2:ncol(X),']',sep='',collapse="+"),sep=''))
ms0= modelSelection(f, smooth=s, priorCoef= momprior(0.192), priorDelta= modelbbprior(1,1), priorGroup=priorGroup, enumerate=FALSE, niter=0, initSearch='none', verbose=FALSE)
Xd= ms0$xstd; groups= ms0$groups

#Define hierarchical constraints for modelSelection (inclusion of non-linear terms requires the presence of linear terms)
constraints= vector("list",length(unique(groups)))
names(constraints)= unique(groups)
for (k in 1:175) constraints[[k]]= integer(0)
for (k in 176:348) constraints[[k]]= as.integer(k - 174)

nfolds= nrow(X) #LOO-CV
pred.aftmom= cv.aftmom(y=y, X=Xd, nfolds=nfolds, priorCoef=momprior(0.192), priorDelta=modelbbprior(1,1), priorGroup=groupzellnerprior(nrow(X)), seed=1234, niter=1000, groups=groups, constraints=constraints)
pred.coxbvs= cv.coxbvs(y=y, X=Xd[,-1], nfolds=nfolds, seed=1234)
pred.coxlasso= cv.coxlasso(y=y, X=Xd[,-1], nfolds=nfolds, seed=1234)
pred.aftlasso= cv.aftlasso(y=y, X=Xd, nfolds=nfolds, seed=1234)
save(pred.aftmom, pred.coxbvs, pred.coxlasso, pred.aftlasso, file='colon_loocv_splines.RData')


load('colon_loocv_splines.RData')
ci.aftmom= concordance.index(-pred.aftmom$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.coxbvs= concordance.index(pred.coxbvs$pred, surv.time=y[,1], surv.event=y[,2])$c.index
ci.coxlasso= concordance.index(pred.coxlasso$pred, surv.time=y[,1], surv.event=y[,2])$c.index
sel= (pred.aftlasso$pred != 0)  #remove cases where method failed to run
ci.aftlasso= concordance.index(-pred.aftlasso$pred[sel], surv.time=y[sel,1], surv.event=y[sel,2])$c.index

tab= round(rbind(c(ci.aftmom, mean(pred.aftmom$nsel)),c(ci.coxbvs, mean(pred.coxbvs$nsel)),c(ci.coxlasso, mean(pred.coxlasso$nsel)),c(ci.aftlasso, mean(pred.aftlasso$nsel))),2)
rownames(tab)= c('AFT pMOMZ','Cox piMOM','Cox LASSO','AFT LASSO')
colnames(tab)= c('CI','nvars')
tab
