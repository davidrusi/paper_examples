#####################################################################################################
##
## AUXILIARY ROUTINES
##
## - ROUTINES TO PRODUCE CROSS-VALIDATED PREDICTIONS OF SURVIVAL TIME / HAZARD
##
## - ROUTINES TO SIMULATE DATA IN SCENARIOS 1-2. SCENARIO 1: TRUTH IS X1 beta1 + log(abs(X2)) beta2. SCENARIO 2: TRUTH IS X1 beta1 + log(1 + X2) beta2
##
## - ROUTINES TO RUN VARIABLE SELECTION METHODS ON SIMULATED DATA
##
## - ROUTINES TO POST-PROCESS SIMULATION RESULTS
##
#####################################################################################################



#####################################################################################################
## ROUTINES TO PRODUCE CROSS-VALIDATED PREDICTIONS OF SURVIVAL TIME / HAZARD
#####################################################################################################


cv.aftmom= function(y, X, nfolds=10, priorCoef=momprior(0.192), priorDelta=modelbbprior(1,1), priorGroup=groupzellnerprior(tau=nrow(X)), seed, ...) {
  #K-fold cross-validated predictions for variables selected by MOM-AFT (union of posterior mode and median prob model). Predictions are based on the MLE under a Normal AFT model.
  if (!missing(seed)) set.seed(seed)
  if (nfolds < nrow(X)) {
      subset= rep(1:nfolds,ceiling(nrow(X)/nfolds))[1:nrow(X)]
      subset= sample(subset,size=nrow(X),replace=FALSE)
  } else {
      subset= 1:nrow(X)
  }
  pred= double(nrow(X)); nsel= integer(nfolds)
  for (k in 1:nfolds) {
      sel= (subset==k)
      ytrain= y[!sel,]; Xtrain= X[!sel,,drop=FALSE]
      ms= modelSelection(ytrain, Xtrain, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, enumerate=FALSE, verbose=FALSE, ...)
      pp= postProb(ms)
      idx= which((ms$postMode==1) | (ms$margpp>0.5))
      nsel[k]= length(idx)
      yexp= ytrain; yexp[,1]= exp(yexp[,1])
      tmp= survreg(yexp ~ -1 + ms$xstd[,idx], dist="lognormal")
      pred[sel]= X[sel,idx,drop=FALSE] %*% matrix(coef(tmp),ncol=1)
      cat(".")
  }
  cat("\n")
  return(list(pred=pred,nsel=nsel))
}


cv.coxbvs= function(y, X, nfolds=10, seed) {
  #K-fold cross-validated predictions for variables selected by Cox iMOM (union of posterior mode and median prob model). Predictions are based on the MPLE under a Cox model.
  if (!missing(seed)) set.seed(seed)
  if (nfolds < nrow(X)) {
      subset= rep(1:nfolds,ceiling(nrow(X)/nfolds))[1:nrow(X)]
      subset= sample(subset,size=nrow(X),replace=FALSE)
  } else {
      subset= 1:nrow(X)
  }
  pred= double(nrow(X)); nsel= integer(nfolds)
  for (k in 1:nfolds) {
      sel= (subset==k)
      ytrain= y[!sel,]; Xtrain= X[!sel,,drop=FALSE]
      fit.coxbvs= runCoxBVS(y=ytrain, x=Xtrain)
      idx= as.numeric(strsplit(as.character(fit.coxbvs$topmodel$model), split=",")[[1]])
      idx= unique(c(which(fit.coxbvs$margpp>.5), idx))
      nsel[k]= length(idx)
      yexp= ytrain; yexp[,1]= exp(yexp[,1])
      tmp= coxph(yexp ~ Xtrain[,idx,drop=FALSE])
      pred[sel]= X[sel,idx,drop=FALSE] %*% matrix(coef(tmp),ncol=1)
      cat(".")
  }
  cat("\n")
  return(list(pred=pred,nsel=nsel))
}


cv.coxlasso= function(y, X, nfolds=10, seed) {
  #K-fold cross-validated predictions for variables selected by Cox LASSO
  if (!missing(seed)) set.seed(seed)
  if (nfolds < nrow(X)) {
      subset= rep(1:nfolds,ceiling(nrow(X)/nfolds))[1:nrow(X)]
      subset= sample(subset,size=nrow(X),replace=FALSE)
  } else {
      subset= 1:nrow(X)
  }
  pred= double(nrow(X)); nsel= integer(nfolds)
  for (k in 1:nfolds) {
      sel= (subset==k)
      ytrain= y[!sel,]; Xtrain= X[!sel,,drop=FALSE]
      resp= cbind(time=exp(ytrain[,1]),status=ytrain[,2])
      cv.fit= try(cv.glmnet(Xtrain, y=resp, family="cox", maxit=10000, nfolds=10, alpha=1), silent=TRUE)
      fit= try(glmnet(Xtrain, y=resp, family="cox", maxit=10000, alpha=1), silent=TRUE)
      b.coxlasso= matrix(as.double(coef(fit, s=cv.fit$lambda.min)), ncol=1)
      nsel[k]= sum(b.coxlasso!=0)
      pred[sel]= X[sel,] %*% b.coxlasso
      cat(".")
  }
  cat("\n")
  return(list(pred=pred,nsel=nsel))
}


cv.aftlasso= function(y, X, nfolds=10, seed) {
  #K-fold cross-validated predictions for variables selected by AFT LASSO
  if (!missing(seed)) set.seed(seed)
  if (nfolds < nrow(X)) {
      subset= rep(1:nfolds,ceiling(nrow(X)/nfolds))[1:nrow(X)]
      subset= sample(subset,size=nrow(X),replace=FALSE)
  } else {
      subset= 1:nrow(X)
  }
  pred= double(nrow(X)); nsel= integer(nfolds)
  for (k in 1:nfolds) {
      sel= (subset==k)
      ytrain= y[!sel,]; Xtrain= X[!sel,,drop=FALSE]
      fit= try(fitAFTLASSOCV(ytrain,Xtrain[,-1]))
      if (class(fit) != "try-error") {
          b.aftlasso= fit[1:ncol(Xtrain)]
          pred[sel]= b.aftlasso[1] + X[sel,-1] %*% matrix(b.aftlasso[-1], ncol=1)
          nsel[k]= sum(b.aftlasso != 0)
      }
      cat(".")
  }
  cat("\n")
  return(list(pred=pred,nsel=nsel))
}






#####################################################################################################
# ROUTINES TO SIMULATE DATA IN SCENARIOS 1-2
#####################################################################################################


simDataScen1= function(seed,beta,n,corr,sigmae,censtimes,errors='normal',asym=0) {
  #Simulate data for Scenario 1
  set.seed(seed)
  p= length(beta)
  X.Sigma <- matrix(corr, nrow=(p), ncol=(p)) + diag(p)*(1 - corr) # Correlation matrix
  X= rmvnorm(n, sigma=X.Sigma) # Design matrix
  if (errors=='normal') {
      e= rnorm(n,0,sigmae)
  } else if (errors=='laplace') {
      #Set scale such that var(e)= sigmae^2. Recall: var of alaplace = 2 * scale * (1+asym^2)
      scale= 0.5 * sigmae^2 / (1+asym^2)
      e= ralapl(n,0,scale=scale,alpha=asym)
      #e= ralapl(n,0,scale=0.5*sigmae^2,alpha=asym)
  } else { stop("Invalid error distribution") }
  lifestimes= exp(X[,1]*beta[1] + log(abs(X[,2]))*beta[2] + e) # times to event
  y= Surv(time=log(pmin(lifestimes, censtimes)), event= censtimes>lifestimes)
  yuncens= y; yuncens[,1]= log(lifestimes); yuncens[,2]=1
  return(list(y=y,yuncens=yuncens,X=X))
}


simDataScen2= function(seed,beta,n,corr,sigmae,censtimes,errors='normal',asym=0) {
  #Simulate data for Scenario 2
  set.seed(seed)
  p= length(beta)
  X.Sigma <- matrix(corr, nrow=(p), ncol=(p)) + diag(p)*(1 - corr) # Correlation matrix
  X= rmvnorm(n, sigma=X.Sigma) # Design matrix
  X[,2]= abs(X[,2])
  if (errors=='normal') {
      e= rnorm(n,0,sigmae)
  } else if (errors=='laplace') {
      e= ralapl(n,0,scale=0.5*sigmae^2,alpha=asym)
  } else { stop("Invalid error distribution") }
  lifestimes= exp(X[,1]*beta[1] + log(1+X[,2])*beta[2] + e) # times to event
  y= Surv(time=log(pmin(lifestimes, censtimes)), event= censtimes>lifestimes)
  yuncens= y; yuncens[,1]= log(lifestimes); yuncens[,2]=1
  return(list(y=y,yuncens=yuncens,X=X))
}



#####################################################################################################
# ROUTINES TO RUN VARIABLE SELECTION METHODS ON SIMULATED DATA
#####################################################################################################


simScenario= function(seed,beta,corr,sigmae,n=100,censtimes=1,scenario,errors='normal',asym=0,priorCoef= momprior(tau=0.192), priorDelta= modelbbprior(1,1), priorGroup=groupzellnerprior(tau=n), omitX2=FALSE, nknots=9) {
  #Simulate data (potentially with many covariates, but typically only use 2) and run modelSelection on the first 2 covariates
  # - seed: seed for random number generation
  # - beta: true AFT model parameters
  # - corr: each row in the covariate matrix arises from a Normal with mean 0, variance 1 and pairwise correlations = corr
  # - sigmae: residual standard deviation
  # - n: sample size
  # - censtimes: administrative censoring time (i.e. any observation > censtimes is censored)
  # - scenario: scenario==1 simulates from simDataScen1, scenario==2 from simDataScen2
  # - error: distribution of the errors, must be 'normal' or 'laplace'
  # - asym: if error=='laplace', asym is the asymmetry parameter of the asymmetric Laplace (asym=0 gives the standard Laplace)
  # - priorCoef, priorDelta, priorGroup: priors on individual coefficients, models, and on coefficient groups, passed to modelSelection
  # - omitX2: set to TRUE to only use X2 to simulate data, but omit it from the actual analysis
  #
  # OUTPUT
  # - pplin: posterior model probabilities when only considering linear effects
  # - ppnolin: posterior model probabilities when only considering non-linear effects (splines)
  # - ppsplines: posterior model probabilities when considering both linear and non-linear effects
  p= length(beta)
  if (scenario==1) {
      sim= simDataScen1(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else if (scenario==2) {
      sim= simDataScen2(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else { stop("Wrong scenario") }
  y= sim$y; yuncens= sim$yuncens; X= sim$X
  if (!omitX2) {
      #Fit model to uncensored data
      msuncenslin= modelSelection(yuncens ~ X[,1] + X[,2], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=FALSE)
      msuncens= modelSelection(yuncens ~ X[,1] + X[,2], smooth=~ X[,1] + X[,2], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, nknots=nknots, verbose=FALSE)
      ppuncenslin= postProb(msuncenslin)
      ppuncens= postProb(msuncens)
      #Fit model to censored data
      mslin= modelSelection(y ~ X[,1] + X[,2], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=FALSE)
      ms= modelSelection(y ~ X[,1] + X[,2], smooth=~ X[,1] + X[,2], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, nknots=nknots, verbose=FALSE)
      pplin= postProb(mslin)
      pp= postProb(ms)
      #Post-process output
      idx1= which(msuncens$groups==1 | msuncens$groups==3)
      idx2= which(msuncens$groups==2 | msuncens$groups==4)
      pplin= cbind(pplinWithoutIntercept(ppuncenslin),pplinWithoutIntercept(pplin))[,-3]
      ppnolin= cbind(ppWithoutLinear(ppuncens, idx1=idx1, idx2=idx2), ppWithoutLinear(pp, idx1=idx1, idx2=idx2))[,-3]
      S1= which(msuncens$groups==3)
      S2= which(msuncens$groups==4)
      ppsplines= data.frame(ppWithoutIntercept(ppuncens, L1=2, L2=3, S1=S1, S2=S2), ppWithoutIntercept(pp, L1=2, L2=3, S1=S1, S2=S2))[,-3]
      #S1= which(msuncens$groups==3); S2= which(msuncens$groups==4)
      #ppsplines= data.frame(ppWithoutIntercept(ppuncens,S1=S1,S2=S2), ppWithoutIntercept(pp,S1=S1,S2=S2))[,-3]
      names(pplin)= names(ppnolin)= names(ppsplines)= c('model','uncensored','censored')
  } else {
      #Fit model to uncensored data
      msuncenslin= modelSelection(yuncens ~ X[,1], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=FALSE)
      msuncens= modelSelection(yuncens ~ X[,1], smooth=~ X[,1], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, nknots=nknots, verbose=FALSE)
      ppuncens= postProb(msuncens)
      #Fit model to censored data
      mslin= modelSelection(y ~ X[,1], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=FALSE)
      ms= modelSelection(y ~ X[,1], smooth=~ X[,1], priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, nknots=nknots, verbose=FALSE)
      pp= postProb(ms)
      #Post-process output
      pplin= c(uncensored=msuncenslin$margpp[2], censored=mslin$margpp[2])
      tmp= ppuncens[!(as.character(ppuncens$modelid) %in% c('2','1,2')),]
      ppnolin1= sum(ppuncens$pp[grep('2,',as.character(tmp$modelid))])
      tmp= pp[!(as.character(pp$modelid) %in% c('2','1,2')),]
      ppnolin2= sum(pp$pp[grep('2,',as.character(tmp$modelid))])
      ppnolin= c(uncensored= ppnolin1, censored= ppnolin2)
      ppsplines= c(uncensored= msuncens$margpp[2], censored= ms$margpp[2])
  }
  ans= list(pplin=pplin,ppnolin=ppnolin,ppsplines=ppsplines)
  cat(".")
  return(ans)
}



simScenarioCoxiMOM= function(seed,beta,corr,sigmae,n=100,censtimes=1,scenario,errors='normal',asym=0,tau=0.25,a=1,b=1,linearfit=FALSE,nknots=9) {
  #Simulate data (potentially with many covariates, but typically only use 2) and run Nikooienejad, Wang & Johnson's iMOM-prior on Cox models on the first 2 covariates
  # - seed: seed for random number generation
  # - beta: true AFT model parameters
  # - corr: each row in the covariate matrix arises from a Normal with mean 0, variance 1 and pairwise correlations = corr
  # - sigmae: residual standard deviation
  # - n: sample size
  # - censtimes: administrative censoring time (i.e. any observation > censtimes is censored)
  # - scenario: scenario==1 simulates from simDataScen1, scenario==2 from simDataScen2
  # - error: distribution of the errors, must be 'normal' or 'laplace'
  # - asym: if error=='laplace', asym is the asymmetry parameter of the asymmetric Laplace (asym=0 gives the standard Laplace)
  # - tau: iMOM prior dispersion parameter (a default tau=0.25 is recommended by the authors, Section 2.3)
  # - a,b: prior on models is Beta-Binomial(a,b)
  # - linearfit: if linearfit==TRUE, perform an additional analysis where only linear covariate effects are considered
  #
  # OUTPUT
  # - margpp: marginal posterior inclusion probabilities for uncensored and censored data
  # - topmodeluncens: highest posterior probability model (uncensored data)
  # - topmodelcens: highest posterior probability model (censored data)
  # - margpplin: if linearfit=TRUE, marginal posterior inclusion probabilities for uncensored and censored data when attention is restricted to linear terms
  # - topmodeluncenslin: same as topmodeluncens when attention is restricted to linear terms
  # - topmodelcenslin: same as topmodelcens when attention is restricted to linear terms
  require(BVSNLP)
  p= length(beta)
  if (scenario==1) {
      sim= simDataScen1(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else if (scenario==2) {
      sim= simDataScen2(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else { stop("Wrong scenario") }
  y= sim$y; yuncens= sim$yuncens; X= sim$X
  des= designmatrixWithSpur(y=y,X=X,smoothterms=TRUE,nknots=nknots)
  desuncens= designmatrixWithSpur(y=yuncens,X=X,smoothterms=TRUE,nknots=nknots)
  L= 10; J= 10; d= 2 * ceiling(log(p)); temps= seq(3, 1, length.out = L); r= 1
  fituncens= cox_bvs(cbind(exp(desuncens$ystd[,1]),desuncens$ystd[,2],desuncens$xstd[,-1]),cur_cols=c(1,2),nf=0,tau=tau,r=1,nlptype=0,a=a,b=b,d=d,L=L,J=J,temps=temps)
  ppuncens= ppbvs(fituncens, groups=desuncens$groups[-1])
  fit= cox_bvs(cbind(exp(des$ystd[,1]),des$ystd[,2],des$xstd[,-1]),cur_cols=c(1,2),nf=0,tau=tau,r=1,nlptype=0,a=1,b=1,d=d,L=L,J=J,temps=temps)
  pp= ppbvs(fit, groups=des$groups[-1])
  margpp= cbind(ppuncens$margpp, pp$margpp)
  margppgroups= cbind(ppuncens$margppgroups, pp$margppgroups)
  groupnames=  c(paste('X',1:p,sep=''), paste('S',1:p,sep=''))
  rownames(margpp)= colnames(des$xstd)[-1]
  rownames(margppgroups)= groupnames; colnames(margpp)= colnames(margppgroups)= c('uncensored','censored')
  ans= list(margpp=margpp, margppgroups=margppgroups, topmodeluncens= ppuncens$topmodel, topmodel=pp$topmodel)
  if (linearfit) {  #Perform analysis that only considers linear terms
      deslin= designmatrixWithSpur(y=y,X=X,smoothterms=FALSE)
      desuncenslin= designmatrixWithSpur(y=yuncens,X=X,smoothterms=FALSE)
      fituncenslin= cox_bvs(cbind(exp(desuncenslin$ystd[,1]),desuncenslin$ystd[,2],desuncenslin$xstd[,-1]),cur_cols=c(1,2),nf=0,tau=tau,r=1,nlptype=0,a=1,b=1,d=d,L=L,J=J,temps=temps)
      ppuncenslin= ppbvs(fituncenslin, groups=desuncenslin$groups[-1])
      fitlin= cox_bvs(cbind(exp(deslin$ystd[,1]),deslin$ystd[,2],deslin$xstd[,-1]),cur_cols=c(1,2),nf=0,tau=tau,r=1,nlptype=0,a=1,b=1,d=d,L=L,J=J,temps=temps)
      pplin= ppbvs(fitlin, groups=deslin$groups[-1])
      margpplin= cbind(ppuncenslin$margpp, pplin$margpp)
      margpplingroups= cbind(ppuncenslin$margppgroups, pplin$margppgroups)
      rownames(margpplin)= colnames(deslin$xstd)[-1]
      rownames(margpplingroups)= paste('X',1:p,sep=''); colnames(margpplin)= colnames(margpplingroups)= c('uncensored','censored')
      ans$margpplin= margpplin; ans$margpplingroups= margpplingroups; ans$topmodeluncenslin= ppuncenslin$topmodel; ans$topmodellin= pplin$topmodel
  }
  cat(".")
  return(ans)
}


ppbvs= function(fit,groups,allmodels=FALSE) {
    #Extract top model, its post model prob and marginal inclusion post prob from an object returned by BVSNLP. If allmodels==TRUE, return posterior prob of all visited models
    models= sapply(fit$vis_covs_list, function(z) paste(z[,1]+1,collapse=','))
    pp= fit$all_probs - fit$max_prob
    pp= exp(pp) / sum(exp(pp))
    sel= which.max(pp)
    topmodel= data.frame(model=models[sel],pp=pp[sel])
    nvar= ncol(fit$max_model)
    modelsbin= t(sapply(fit$vis_covs_list, function(z) { ans= rep(0,nvar); ans[z[,1]+1]=1; return(ans) }))
    margpp= matrix(pp,nrow=1) %*% modelsbin
    groupsbin= sapply(unique(groups), function(g) rowSums(modelsbin[,groups==g,drop=FALSE])>0)
    margppgroups= matrix(pp,nrow=1) %*% groupsbin
    if (allmodels) { ppmodels= data.frame(model=models, pp=pp) } else { ppmodels= NULL }
    ans= list(margpp=as.vector(margpp), margppgroups=as.vector(margppgroups), topmodel=topmodel, ppmodels=ppmodels)
    return(ans)
}

runCoxBVS= function(y, x, groups= 1:ncol(x), tau=0.25) {
  #Run cox_bvs from package NLPBVS, return posterior marginal inclusion probabilities and posterior model probabilities
  L= 10; J= 10; d= 2 * ceiling(log(ncol(x))); temps= seq(3, 1, length.out = L); r= 1
  fit= cox_bvs(cbind(exp(y[,1]),y[,2],x),cur_cols=1,nf=0,tau=tau,r=1,nlptype=0,a=1,b=1,d=d,L=L,J=J,temps=temps)
  pp= ppbvs(fit, groups=groups, allmodels=TRUE)
  return(pp)
}



simScenarioCoxLASSO= function(seed,beta,corr,sigmae,n=100,censtimes=1,scenario,errors='normal',asym=0,smoothterms=TRUE,nknots=9,verbose=TRUE) {
  #Simulate data with many covariates and run Cox LASSO from glmnet package
  # - seed: seed for random number generation
  # - beta: true AFT model parameters
  # - corr: each row in the covariate matrix arises from a Normal with mean 0, variance 1 and pairwise correlations = corr
  # - sigmae: residual standard deviation
  # - n: sample size
  # - censtimes: administrative censoring time (i.e. any observation > censtimes is censored)
  # - scenario: scenario==1 simulates from simDataScen1, scenario==2 from simDataScen2
  # - error: distribution of the errors, must be 'normal' or 'laplace'
  # - asym: if error=='laplace', asym is the asymmetry parameter of the asymmetric Laplace (asym=0 gives the standard Laplace)
  # - smoothterms: if TRUE cubic spline terms are added to the design matrix
  #
  # OUTPUT: list with two elements
  # - beta: estimated parameters at LASSO regularization parameter given by 10-fold cross-validation, and group that each parameter belongs to (linear / spline term for each of the p variables)
  # - nselgroup: number of non-zero beta's for each group
  # - betauncens: same as beta, for uncensored data
  # - nselgroupuncens: same as nselgroup, for uncensored data
  p= length(beta)
  #Simulate data
  if (scenario==1) {
      sim= simDataScen1(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else if (scenario==2) {
      sim= simDataScen2(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else { stop("Wrong scenario") }
  y= sim$y; yuncens= sim$yuncens; X= sim$X
  des= designmatrixWithSpur(y=y,X=X,smoothterms=smoothterms,nknots=nknots)
  desuncens= designmatrixWithSpur(y=yuncens,X=X,smoothterms=smoothterms,nknots=nknots)
  if (smoothterms) {
      groupnames=  c('Intercept',paste('X',1:p,sep=''), paste('S',1:p,sep=''))
  } else {
      groupnames=  c('Intercept',paste('X',1:p,sep=''))
  }
  #Run for censored data
  resp= cbind(time=exp(des$ystd[,1]),status=des$ystd[,2])
  cv.fit= try(cv.glmnet(des$xstd, y=resp, family="cox", maxit=10000, nfolds=10, alpha=1), silent=TRUE)   # Cox LASSO (alpha = 1 in the elastic net)
  fit= try(glmnet(des$xstd, y=resp, family = "cox", maxit=10000, alpha=1), silent=TRUE)
  bhat= as.double(coef(fit, s=cv.fit$lambda.min))
  bhat= data.frame(estimate= bhat, group= groupnames[des$groups+1])
  rownames(bhat)= NULL
  nselgroup= tapply(bhat[,1],INDEX=bhat$group,FUN= function(z) sum(z!=0))
  nselgroup= nselgroup[groupnames]
  #Run for uncensored data
  resp= cbind(time=exp(desuncens$ystd[,1]),status=desuncens$ystd[,2])
  cv.fit= try(cv.glmnet(desuncens$xstd, y=resp, family="cox", maxit=10000, nfolds=10, alpha=1), silent=TRUE)   # Cox LASSO (alpha = 1 in the elastic net)
  fit= try(glmnet(desuncens$xstd, y=resp, family = "cox", maxit=10000, alpha=1), silent=TRUE)
  bhatuncens= as.double(coef(fit, s=cv.fit$lambda.min))
  bhatuncens= data.frame(estimate= bhatuncens, group= groupnames[desuncens$groups+1])
  rownames(bhatuncens)= NULL
  nselgroupuncens= tapply(bhatuncens[,1],INDEX=bhatuncens$group,FUN= function(z) sum(z!=0))
  nselgroupuncens= nselgroupuncens[groupnames]
  #Return output
  ans= list(beta=bhat, nselgroup=nselgroup, betauncens=bhatuncens, nselgroupuncens=nselgroupuncens)
  if (verbose) cat('.')
  return(ans)
}


simScenarioAFTLASSO= function(seed,beta,corr,sigmae,n=100,censtimes=1,scenario,errors='normal',asym=0,smoothterms=TRUE,nknots=9,verbose=TRUE) {
  #Simulate data with many covariates and run AFT LASSO from AdapEnetClass package
  # - seed: seed for random number generation
  # - beta: true AFT model parameters
  # - corr: each row in the covariate matrix arises from a Normal with mean 0, variance 1 and pairwise correlations = corr
  # - sigmae: residual standard deviation
  # - n: sample size
  # - censtimes: administrative censoring time (i.e. any observation > censtimes is censored)
  # - scenario: scenario==1 simulates from simDataScen1, scenario==2 from simDataScen2
  # - error: distribution of the errors, must be 'normal' or 'laplace'
  # - asym: if error=='laplace', asym is the asymmetry parameter of the asymmetric Laplace (asym=0 gives the standard Laplace)
  # - smoothterms: if TRUE cubic spline terms are added to the design matrix
  #
  # OUTPUT: list with elements
  # - beta: estimated parameters at LASSO regularization parameter given by 10-fold cross-validation, and group that each parameter belongs to (linear / spline term for each of the p variables)
  # - nselgroup: number of non-zero beta's for each group
  # - betauncens: same as beta, for uncensored data
  # - nselgroupuncens: same as nselgroup, for uncensored data
  p= length(beta)
  if (scenario==1) {
      sim= simDataScen1(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else if (scenario==2) {
      sim= simDataScen2(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else { stop("Wrong scenario") }
  y= sim$y; yuncens= sim$yuncens; X= sim$X
  des= designmatrixWithSpur(y=y,X=X,smoothterms=smoothterms,nknots=nknots)
  desuncens= designmatrixWithSpur(y=yuncens,X=X,smoothterms=smoothterms,nknots=nknots)
  if (smoothterms) {
      groupnames=  c('Intercept',paste('X',1:p,sep=''), paste('S',1:p,sep=''))
  } else {
      groupnames=  c('Intercept',paste('X',1:p,sep=''))
  }
  #Run for censored data
  fit= try(fitAFTLASSOCV(des$ystd,des$xstd[,-1]), silent=TRUE)
  #fitk= fitAFTLASSO(y=des$ystd,X=des$xstd[,-1])
  if (class(fit) != "try-error") {
      beta= data.frame(estimate= fit[1:ncol(des$xstd)], group= groupnames[des$groups+1])
      rownames(beta)= NULL
      nselgroup= tapply(beta[,1],INDEX=beta$group,FUN= function(z) sum(z!=0))
      nselgroup= nselgroup[groupnames]
  } else {
      beta= data.frame(estimate= rep(NA,ncol(des$xstd)), group= groupnames[des$groups+1])
      rownames(beta)= NULL
      nselgroup= rep(NA, length(groupnames)); names(nselgroup)= groupnames
  }
  #Run for uncensored data
  fituncens= try(fitAFTLASSOCV(desuncens$ystd,desuncens$xstd[,-1]), silent=TRUE)
  if (class(fituncens) != "try-error") {
      betauncens= data.frame(estimate= fituncens[1:ncol(desuncens$xstd)], group= groupnames[desuncens$groups+1])
      rownames(betauncens)= NULL
      nselgroupuncens= tapply(betauncens[,1],INDEX=betauncens$group,FUN= function(z) sum(z!=0))
      nselgroupuncens= nselgroupuncens[groupnames]
  } else {
      betauncens= data.frame(estimate= rep(NA,ncol(desuncens$xstd)), group= groupnames[desuncens$groups+1])
      rownames(betauncens)= NULL
      nselgroupuncens= rep(NA, length(groupnames)); names(nselgroupuncens)= groupnames
  }
  ans= list(beta=beta, nselgroup=nselgroup, betauncens=betauncens, nselgroupuncens=nselgroupuncens)
  if (verbose) cat('.')
  return(ans)
}


fitAFTLASSOCV= function(y,X,K=10) {
#Fit AFT LASSO via AEnet.aft, setting the penalization parameter via K-fold cross-validation
  subset= rep(1:K,ceiling(nrow(X)/K))[1:nrow(X)]
  subset= sample(subset,size=nrow(X),replace=FALSE)
  neglogl= vector("list",length=K)
  for (k in 1:K) {
      sel <- subset==k
      fitk= fitAFTLASSO(y=y[!sel],X=X[!sel,,drop=FALSE],ynew=y[sel],Xnew=X[sel,,drop=FALSE])
      neglogl[[k]]= fitk[,'neglogl']
  }
  nlambda= min(sapply(neglogl,length))
  neglogl= sapply(neglogl,function(z) z[1:nlambda])
  neglogl.cv= apply(neglogl,1,function(z) sum(z,na.rm=TRUE))
  fit= fitAFTLASSO(y=y,X=X)
  ans= fit[min(which.min(neglogl.cv),nrow(fit)),]
  return(ans)
}


fitAFTLASSO= function(y,X,ynew,Xnew) {
  #Fit AFT LASSO model via AEnet.aft. For all LASSO grid return estimated regression parameters, rho=log(residual precision) based on X, and cross-validated log-likelihood evaluated at (ynew,Xnew)
  require(AdapEnetClass)
  if (!missing(ynew) && !missing(Xnew)) { crossval= TRUE } else { crossval= FALSE }
  l= try(mrbj(y ~ X, mcsize=100, trace=FALSE, gehanonly=FALSE), silent=TRUE)
  if (class(l) != "try-error") {
    wt= l$enet
    ft= AEnet.aft(X=X,Y=y[,1], delta=y[,2], weight=wt, lambda2=0, maxit=1000)
    neglogl= rhoopt= rep(NA,nrow(ft$beta))
    for (i in 1:nrow(ft$beta)) {
      bhat= c(ft$mu,ft$beta[i,])
      opt= nlminb(start=0,negloglikAFT_rho,y=y,X=cbind(1,X),beta=bhat)
      rhoopt[i]= opt$par
      if (crossval) neglogl[i]= negloglikAFT_rho(rhoopt[i],beta=bhat,y=ynew,X=cbind(1,Xnew))
    }
    if (crossval) {
      ans= cbind(mu=ft$mu,ft$beta,rho=rhoopt,neglogl)
    } else {
      ans= cbind(mu=ft$mu,ft$beta,rho=rhoopt)
    }
    colnames(ans)[2:(ncol(X)+1)]= paste('beta',1:ncol(X),sep='')
  } else {
      ans= matrix(NA,nrow=ncol(X)+1,ncol=ncol(X)+2+crossval)
      nn= c('mu',paste('beta',1:ncol(X),sep=''),'rho')
      if (crossval) nn= c(nn,'neglogl')
      colnames(ans)= nn
  }
  return(ans)
}


negloglikAFT_rho <- function(rho,beta,y,X) {
    tau= exp(rho)
    alpha= beta * tau
    -loglikAFT(c(alpha,rho),y=y,X=X)
}


designmatrixWithSpur= function(y,X,smoothterms=TRUE,nknots=9,priorCoef= momprior(tau=0.192), priorDelta= modelbbprior(1,1), priorGroup=groupzellnerprior(tau=n)) {
  #Return design matrix with standardized outcome y, linear terms X, and if smoothterms==TRUE also spline basis with nknots
  # - y: outcome
  # - X: matrix with linear covariates
  # - smoothterms: if TRUE cubic spline terms are added to the basis. These terms are made orthogonal to the linear terms in X
  # - nknots: number of knots (for each covariate)
  # Output
  # - ystd: y standardized to mean 0, variance 1
  # - xstd: x standardized to mean 0, variance 1 (only for numeric variables, factors are not transformed)
  # - groups: groups for columns in xstd. Each linear term in X is in its own separate group. Spline terms corresponding to a single column in x are grouped
  p= ncol(X); n= nrow(X)
  f= formula(paste("y ~",paste(paste('X[,',1:p,']',sep=''),collapse='+')))
  s= formula(paste("~",paste(paste('X[,',1:p,']',sep=''),collapse='+')))
  if (smoothterms) {
      ms= modelSelection(f, smooth=s, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=FALSE, nknots=nknots, niter=0, initSearch='none')
  } else {
      ms= modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=FALSE, nknots=nknots, niter=0, initSearch='none')
  }
  return(list(ystd=ms$ystd, xstd=ms$xstd, groups=ms$groups))
}




simScenarioWithSpur= function(seed,beta,corr,sigmae,n=100,censtimes=1,scenario,errors='normal',asym=0,smoothterms=TRUE,nknots=9,verbose=FALSE,priorCoef= momprior(tau=0.192), priorDelta= modelbbprior(1,1), priorGroup=groupzellnerprior(tau=n),usevars=1:length(beta)) {
  #Simulate data with many covariates and run modelSelection on all of them
  # - seed: seed for random number generation
  # - beta: true AFT model parameters
  # - corr: each row in the covariate matrix arises from a Normal with mean 0, variance 1 and pairwise correlations = corr
  # - sigmae: residual standard deviation
  # - n: sample size
  # - censtimes: administrative censoring time (i.e. any observation > censtimes is censored)
  # - scenario: scenario==1 simulates from simDataScen1, scenario==2 from simDataScen2
  # - error: distribution of the errors, must be 'normal' or 'laplace'
  # - asym: if error=='laplace', asym is the asymmetry parameter of the asymmetric Laplace (asym=0 gives the standard Laplace)
  # - smoothterms: if FALSE consider only linear effects, if TRUE consider both linear and non-linear effects
  # - verbose: set to TRUE to print progress informations
  # - priorCoef, priorDelta, priorGroup: priors on individual coefficients, models, and on coefficient groups, passed to modelSelection
  # - usevars: covariates to be used in the model fit. Covariates not in usevars are used to generate the data, but are omitted in the actual analysis
  #
  # OUTPUT
  # - margpp: marginal posterior inclusion probabilities without and with censoring (1st and 2nd columns respectively)
  # - topmodel.uncens: highest posterior probability model (uncensored data)
  # - topmodel.cens: highest posterior probability model (censored data)
  p= length(beta)
  if (scenario==1) {
      sim= simDataScen1(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else if (scenario==2) {
      sim= simDataScen2(seed=seed,beta=beta,n=n,corr=corr,sigmae=sigmae,censtimes=censtimes,errors=errors,asym=asym)
  } else { stop("Wrong scenario") }
  y= sim$y; yuncens= sim$yuncens; X= sim$X
  #Fit model to uncensored data
  fu= formula(paste("yuncens ~",paste(paste('X[,',usevars,']',sep=''),collapse='+')))
  f= formula(paste("y ~",paste(paste('X[,',usevars,']',sep=''),collapse='+')))
  s= formula(paste("~",paste(paste('X[,',usevars,']',sep=''),collapse='+')))
  if (smoothterms) {
      msuncens= modelSelection(fu, smooth=s, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=verbose, nknots=nknots)
      ppuncens= pp2margpp(msuncens)
      ms= modelSelection(f, smooth=s, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=verbose, nknots=nknots)
      pp= pp2margpp(ms)
  } else {
      msuncens= modelSelection(fu, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=verbose, nknots=nknots)
      ppuncens= pp2margpp(msuncens)
      ms= modelSelection(f, priorCoef=priorCoef, priorDelta=priorDelta, priorGroup=priorGroup, verbose=verbose, nknots=nknots)
      pp= pp2margpp(ms)
  }
  #Marginal posterior inclusion probabilities
  margpp= cbind(ppuncens$margpp,pp$margpp)
  colnames(margpp)= c('No censoring','Censoring')
  smoothterms= sapply(1:length(beta),function(i) grep(paste("X\\[, ",i,"\\].s",sep=''),rownames(margpp))[1]) #select one margpp per smooth term
  margpp= margpp[c(1:(length(beta)+1),smoothterms),]
  rownames(margpp)= sub("\\.s.",".spline",rownames(margpp))
  #Highest posterior probability model (uncensored data)
  topmodel.uncens= ppuncens$pp[1,]
  varnames.uncens= rownames(ppuncens$margpp)[as.numeric(strsplit(as.character(topmodel.uncens$modelid),split=',')[[1]])]
  sel.smooth= sapply(1:length(beta),function(i) grep(paste("X\\[, ",i,"\\].s",sep=''),varnames.uncens)[1])
  sel.smooth= sel.smooth[!is.na(sel.smooth)]
  varnames.uncens= varnames.uncens[c(grep("\\]$",varnames.uncens), sel.smooth)]
  varnames.uncens= sub("\\.s.",".spline",varnames.uncens)
  topmodel.uncens= data.frame(varnames.uncens=paste(varnames.uncens, collapse=', '), pp=topmodel.uncens$pp)
  #Highest posterior probability model (censored data)
  topmodel.cens= pp$pp[1,]
  varnames.cens= rownames(pp$margpp)[as.numeric(strsplit(as.character(topmodel.cens$modelid),split=',')[[1]])]
  sel.smooth= sapply(1:length(beta),function(i) grep(paste("X\\[, ",i,"\\].s",sep=''),varnames.cens)[1])
  sel.smooth= sel.smooth[!is.na(sel.smooth)]
  varnames.cens= varnames.cens[c(grep("\\]$",varnames.cens), sel.smooth)]
  varnames.cens= sub("\\.s.",".spline",varnames.cens)
  topmodel.cens= data.frame(varnames.cens=paste(varnames.cens, collapse=', '), pp=topmodel.cens$pp)
  ans= list(margpp=margpp, topmodel.uncens=topmodel.uncens, topmodel.cens=topmodel.cens)
  cat(".")
  return(ans)
}


pp2margpp= function(msfit) {
#Obtain posterior model probabilities (renormalized) and corresponding marginal posterior probabilities
    pp= postProb(msfit)
    xnames= colnames(msfit$xstd)
    p= length(xnames)
    modelid= strsplit(as.character(pp[,1]), split=',')
    model= do.call(cbind,lapply(modelid, function(z) { ans= rep(0,p); ans[as.numeric(z)]= 1; return(ans) }))
    margpp= model %*% pp[,3]
    rownames(margpp)= xnames
    return(list(pp=pp,margpp=margpp))
}



pplinWithoutIntercept= function(pp) {
#Return post prob marginalizing over the intercept term, assuming the model contained only linear effects for X1,X2
    ans= pp
    ans$modelid=sub('1,','',ans$modelid)
    ans$modelid[ans$modelid==1]= ''
    ans= tapply(ans$pp,ans$modelid,sum)
    ans= ans[c(which(names(ans)==""),which(names(ans)=='2'),which(names(ans)=='3'),which(names(ans)=='2,3'))]
    ans= data.frame(model=c("","X1lin","X2lin","X1lin,X2lin"),pp=ans)
    return(ans)
}


ppWithoutIntercept= function(pp, L1=2, L2=3, S1= 4:13, S2= 14:23) {
#Return post prob marginalizing over the intercept term, assuming the model contained both linear effects and splines for X1,X2
    modelid= lapply(strsplit(as.character(pp$modelid),split=','),as.numeric)
    hasL1= sapply(modelid, function(z) L1 %in% z)
    hasL2= sapply(modelid, function(z) L2 %in% z)
    hasS1= sapply(modelid, function(z) any(S1 %in% z))
    hasS2= sapply(modelid, function(z) any(S2 %in% z))
    hasmat= cbind(hasL1,hasL2,hasS1,hasS2)
    ans= pp
    ans$modelid= apply(hasmat, 1, function(z) paste(c("1","2","S1","S2")[z], collapse=","))
    ans= tapply(ans$pp,INDEX=ans$modelid,FUN=sum)
    ans= data.frame(variables=names(ans),pp=ans)
    rownames(ans)[rownames(ans)==""]= "0"
    ans= ans[c("0","1","1,S1","2","2,S2","1,2","1,2,S1","1,2,S2","1,2,S1,S2"),]
    #ans= pp
    #ans$modelid= sub(paste(as.character(S1),collapse=','),'S1',ans$modelid)
    #ans$modelid= sub(paste(as.character(S2),collapse=','),'S2',ans$modelid)
    #ans$modelid= sub('1,2','2',ans$modelid)
    #ans$modelid= sub('1,3','3',ans$modelid)
    #ans$modelid[ans$modelid=='1']= ''
    #ans= tapply(ans$pp,INDEX=ans$modelid,FUN=sum)
    #ans= data.frame(variables=names(ans),pp=ans)
    #ans$pp= ans$pp/sum(ans$pp)
    #ans= ans[order(ans$pp,decreasing=TRUE),]
    #ans= ans[sapply(c("","2","2,S1","3","3,S2","2,3","2,3,S1","2,3,S2","2,3,S1,S2"), function(z) which(ans[,1]==z)),]
    #ans[,1]= c("","1","1,S1","2","2,S2","1,2","1,2,S1","1,2,S2","1,2,S1,S2")
    rownames(ans)= NULL
    return(ans)
}

ppWithoutLinear= function(pp, idx1= c(2,4:13), idx2= c(3,14:23)) {
    #Return posterior probabilities removing the models with only linear terms
    modelid= lapply(strsplit(as.character(pp$modelid),split=','),as.numeric)
    allidx1= sapply(modelid, function(z) all(idx1 %in% z))
    allidx2= sapply(modelid, function(z) all(idx2 %in% z))
    noidx1= sapply(modelid, function(z) all(!(idx1 %in% z)))
    noidx2= sapply(modelid, function(z) all(!(idx2 %in% z)))
    hasidx= cbind(allidx1,allidx2,noidx1,noidx2)
    sel= (noidx1 & noidx2) | (allidx1 & noidx2) | (noidx1 & allidx2) | (allidx1 & allidx2)
    ppsel= data.frame(hasidx[sel,],pp[sel,-2])
    nn= paste(ifelse(ppsel[,1],'X1',''),ifelse(ppsel[,2],'X2',''),sep='')
    ans= tapply(ppsel$pp,INDEX=nn,FUN=sum)
    ans= data.frame(variables=names(ans),pp=ans)
    ans$pp= ans$pp/sum(ans$pp)
    ans= ans[c(which(ans$variables==''),which(ans$variables=='X1'),which(ans$variables=='X2'),which(ans$variables=='X1X2')),]
    rownames(ans)= NULL
    return(ans)
}



#####################################################################################################
# ROUTINES TO POST-PROCESS SIMULATION RESULTS
#####################################################################################################

margpp= function(simout) {
    #Average marginal posterior inclusion probability for (X1,X2) under linear effects, non-linear effects, linear + non-linear effects
    ans= lapply(1:length(simout), function(z) matrix(NA,nrow=3,ncol=4))
    for (i in 1:length(simout)) {
        z= simout[[i]]
        ans[[i]][1,]= c(colSums(z$pplin[grep('X1',z$pplin$model),-1]), colSums(z$pplin[grep('X2',z$pplin$model),-1]))
        ans[[i]][2,]= c(colSums(z$ppnolin[grep('X1',z$ppnolin$model),-1]), colSums(z$ppnolin[grep('X2',z$ppnolin$model),-1]))
        ans[[i]][3,]= c(colSums(z$ppsplines[grep('1',z$ppsplines$model),-1]), colSums(z$ppsplines[grep('2',z$ppsplines$model),-1]))
        rownames(ans[[i]])= c('Linear','Non-linear','Linear and non-linear')
        colnames(ans[[i]])= c('X1, uncensored','X1, censored','X2, uncensored','X2, censored')
    }
    ans= ans[sapply(ans, function(z) any(is.na(z))) == FALSE]
    ans= Reduce("+",ans)/length(ans)
    return(ans)
}


averageSims= function(simout) {
    #Report average results over simulations
    pplin= simout[[1]]$pplin[,-1]; nsim=1
    for (i in 2:length(simout)) {
        notna= !any(is.na(simout[[i]]$pplin[,2])) & !any(is.na(simout[[i]]$pplin[,3]))
        if (notna) { pplin= pplin + simout[[i]]$pplin[,-1]; nsim= nsim+1 }
    }
    pplin= data.frame(model=simout[[1]]$pplin[,1],pplin/nsim)
    #
    ppnolin= simout[[1]]$ppnolin[,-1]; nsim=1
    for (i in 2:length(simout)) {
        notna= !any(is.na(simout[[i]]$ppnolin[,2])) & !any(is.na(simout[[i]]$ppnolin[,3]))
        if (notna) { ppnolin= ppnolin + simout[[i]]$ppnolin[,-1]; nsim= nsim+1 }
    }
    ppnolin= data.frame(model=simout[[1]]$ppnolin[,1],ppnolin/nsim)
    #
    ppsplines= simout[[1]]$ppsplines[,-1]; nsim=1
    for (i in 2:length(simout)) {
        notna= !any(is.na(simout[[i]]$ppsplines[,2])) & !any(is.na(simout[[i]]$ppsplines[,3]))
        if (notna) { ppsplines= ppsplines + simout[[i]]$ppsplines[,-1]; nsim= nsim+1 }
    }
    ppsplines= data.frame(model=simout[[1]]$ppsplines[,1],ppsplines/nsim)
    ans= list(pplin=pplin,ppnolin=ppnolin,ppsplines=ppsplines)
    return(ans)
}



correctselmombf= function(simout,p,activevars,uncens) {
    #Return proportion of correct model selections
    # - simout: list where each element is simulation output returned by simScenarioWithSpur
    # - p: number of variables
    # - activevars: indexes truly active variables, e.g. activevars=1:2
    # - uncens: if TRUE results are returned for uncensored data, if FALSE results are for censored data
    # Output: a vector with the following elements
    # - Proportion of correct model selections
    # - Average posterior probability of the correct model, conditional on it being selected
    # - Average number of truly active variables that were selected
    # - Average number of truly inactive variables that were selected
    inactivevars= setdiff(1:p, activevars)
    if (uncens) {
        topmodel= sapply(simout, function(z) as.character(z$topmodel.uncens[1,1]))
        pp= sapply(simout, function(z) z$topmodel.uncens[1,'pp'])
    } else {
        topmodel= sapply(simout, function(z) as.character(z$topmodel.cens[1,1]))
        pp= sapply(simout, function(z) z$topmodel.cens[1,'pp'])
    }
    #Number of truly active variables that were selected
    hasactive= matrix(0,nrow=length(topmodel),ncol=length(activevars))
    for (i in 1:length(activevars)) {
        j= activevars[i]
        hasactive[grep(paste("X\\[, ",j,"\\]",sep=''),topmodel), i]= 1
        hasactive[grep(paste("X\\[, ",j,"\\].spline",sep=''),topmodel), i]= 1
    }
    nactive= rowSums(hasactive)
    #Number of truly inactive variables that were selected
    hasinactive= matrix(0,nrow=length(topmodel),ncol=length(inactivevars))
    for (i in 1:length(inactivevars)) {
        j= inactivevars[i]
        hasinactive[grep(paste("X\\[, ",j,"\\]",sep=''),topmodel), i]= 1
        hasinactive[grep(paste("X\\[, ",j,"\\].spline",sep=''),topmodel), i]= 1
    }
    ninactive= rowSums(hasinactive)
    #Correct model selections
    correct= (nactive==length(activevars)) & (ninactive==0)
    ans= matrix(c(mean(correct), mean(pp[correct]), mean(nactive), mean(ninactive)),nrow=1)
    colnames(ans)= c('correct model','pp when correct','nactive','ninactive')
    return(ans)
}


correctselCoxiMOM= function(simout,p,activevars,uncens) {
    #Return proportion of correct model selections
    # - simout: list where each element is simulation output returned by simScenarioCoxiMOM
    # - p: number of variables
    # - activevars: indexes truly active variables, e.g. activevars=1:2
    # - uncens: if TRUE results are returned for uncensored data, if FALSE results are for censored data
    # Output: a vector with the following elements
    # - Proportion of correct model selections
    # - Average posterior probability of the correct model, conditional on it being selected
    # - Average number of truly active variables that were selected
    # - Average number of truly inactive variables that were selected
    inactivevars= setdiff(1:p, activevars)
    nn= rownames(simout[[1]]$margpp)
    if (uncens) {
        topmodel= sapply(simout, function(z) { sel= strsplit(as.character(z$topmodeluncens[1,1]), split=',')[[1]]; paste(rownames(z$margpp)[as.numeric(sel)],collapse=', ') } )
        pp= sapply(simout, function(z) z$topmodeluncens[1,'pp'])
    } else {
        topmodel= sapply(simout, function(z) { sel= strsplit(as.character(z$topmodel[1,1]), split=',')[[1]]; paste(rownames(z$margpp)[as.numeric(sel)],collapse=', ') } )
        pp= sapply(simout, function(z) z$topmodel[1,'pp'])
    }
    #Number of truly active variables that were selected
    hasactive= matrix(0,nrow=length(topmodel),ncol=length(activevars))
    for (i in 1:length(activevars)) {
        j= activevars[i]
        hasactive[grep(paste("X\\[, ",j,"\\]",sep=''),topmodel), i]= 1
        #hasactive[grep(paste("X\\[, ",j,"\\].s",sep=''),topmodel), i]= 1  #not needed, implied by the earlier condition
    }
    nactive= rowSums(hasactive)
    #Number of truly inactive variables that were selected
    hasinactive= matrix(0,nrow=length(topmodel),ncol=length(inactivevars))
    for (i in 1:length(inactivevars)) {
        j= inactivevars[i]
        hasinactive[grep(paste("X\\[, ",j,"\\]",sep=''),topmodel), i]= 1
        #hasinactive[grep(paste("X\\[, ",j,"\\].s",sep=''),topmodel), i]= 1 #not needed, implied by the earlier condition
    }
    ninactive= rowSums(hasinactive)
    #Correct model selections
    correct= (nactive==length(activevars)) & (ninactive==0)
    ans= matrix(c(mean(correct), mean(pp[correct]), mean(nactive), mean(ninactive)),nrow=1)
    colnames(ans)= c('correct model','pp when correct','nactive','ninactive')
    return(ans)
}


correctselLASSO= function(simout,p,activevars,uncens) {
    #Return proportion of correct model selections
    # - simout: list where each element is simulation output returned by either simScenarioAFTLASSO or simScenarioCoxLASSO
    # - p: number of variables
    # - activevars: indexes truly active variables, e.g. activevars=1:2
    # - uncens: if TRUE results are returned for uncensored data, if FALSE results are for censored data
    # Output: a vector with the following elements
    # - Proportion of correct model selections
    # - Average number of truly active variables that were selected
    # - Average number of truly inactive variables that were selected
    nselgroup= t(sapply(simout, function(z) z$nselgroup))
    inactive= setdiff(1:p, activevars)
    if (length(inactive)>0) {
        ninactive= rowSums((nselgroup[,paste('X',inactive,sep='')] + nselgroup[,paste('S',inactive,sep='')])>0)
    } else {
        ninactive= 0
    }
    nactive= rowSums((nselgroup[,paste('X',activevars,sep='')] + nselgroup[,paste('S',activevars,sep='')])>0)
    correct= (nactive==length(activevars)) & (ninactive==0)
    ans= matrix(c(mean(correct,na.rm=TRUE), NA, mean(nactive,na.rm=TRUE), mean(ninactive,na.rm=TRUE)),nrow=1)
    colnames(ans)= c('correct model','pp when correct','nactive','ninactive')
    return(ans)
}



margppmombf= function(simout,p,activevars) {
    #Return average marginal posterior inclusion prob for each active var
    # - simout: list where each element is simulation output returned by simScenarioWithSpur
    # - p: number of variables
    # - activevars: indexes truly active variables, e.g. activevars=1:2
    margpp= Reduce("+",lapply(simout,function(z) z$margpp)) / length(simout)
    inactive= setdiff(1:p, activevars)
    sel= rownames(margpp) %in% c(paste('X[, ',activevars,']',sep=''),paste('X[, ',activevars,'].spline',sep=''))
    spurlinear= rownames(margpp) %in% paste('X[, ',inactive,']',sep='')
    spurspline= rownames(margpp) %in% paste('X[, ',inactive,'].spline',sep='')
    ans= rbind(margpp[sel,], colMeans(margpp[spurlinear,]), colMeans(margpp[spurspline,]))
    rownames(ans)[5:6]= c('spur.linear','spur.spline')
    return(ans)
}


margppCoxiMOM= function(simout,p,activevars) {
    #Return average marginal posterior inclusion prob for each active var
    # - simout: list where each element is simulation output returned by simScenarioCoxiMOM
    # - p: number of variables
    # - activevars: indexes truly active variables, e.g. activevars=1:2
    margpp= Reduce("+",lapply(simout,function(z) z$margppgroups)) / length(simout)
    inactive= setdiff(1:p, activevars)
    sel= rownames(margpp) %in% c(paste('X',activevars,sep=''),paste('S',activevars,sep=''))
    spurlinear= rownames(margpp) %in% paste('X',inactive,sep='')
    spurspline= rownames(margpp) %in% paste('S',inactive,sep='')
    ans= rbind(margpp[sel,], colMeans(margpp[spurlinear,]), colMeans(margpp[spurspline,]))
    rownames(ans)[5:6]= c('spur.linear','spur.spline')
    return(ans)
}





#####################################################################################################
##
## MODEL-FITTING ROUTINES TO CHECK CODE IN MOMBF
##
#####################################################################################################

#log-likelihood in (alpha, rho=log(tau)) parameterization
loglikAFT <- function(par,y,X,bothterms=FALSE) {
  cens= (y[,2]==0); nuncens= sum(!cens)
  alpha <- matrix(head(par,-1),ncol=1); rho <- tail(par,1); exprho= exp(rho)
  if (ncol(X)>0) {
    var.obs <- -0.5 * sum((exprho*y[!cens,1] - X[!cens,]%*%alpha)^2) - nuncens/2 * log(2*pi/exprho^2)
    var.cens <- sum(pnorm(-exprho*y[cens,1] + X[cens,]%*%alpha,log=T))
  } else {
    var.obs <- -0.5 * sum((exprho*y[!cens,1])^2) - nuncens/2 * log(2*pi/exprho^2)
    var.cens <- sum(pnorm(-exprho*y[cens,1],log=T))
  }
  if (!bothterms) { ans= var.obs + var.cens } else { ans= list(all= var.obs+var.cens, obs=var.obs, cens=var.cens) }
  return(ans)
}


invmillsnorm= function(z) dnorm(z) / pnorm(z)
Dz= function(z) { m= invmillsnorm(-z); return(m^2 - z*m) }

#Gradient and hessian of log-likelihood in (alpha, rho=log(tau)) parameterization
loglikAFTgradhess <- function(par,y,X) {
  cens= (y[,2]==0); nuncens= sum(!cens)
  alpha= head(par,-1); rho= tail(par,1); exprho= exp(rho)
  g= g.obs= g.cens= double(length(par)); H= H.obs= H.cens= matrix(NA,nrow=length(par),ncol=length(par))
  if (ncol(X)>0) {
      res= exprho*y[,1] - X%*%alpha
      r= invmillsnorm(-res[cens])
      D= Dz(res[cens])
      g.obs[-length(g.obs)]= t(X[!cens,]) %*% res[!cens]
      g.cens[-length(g.cens)]= t(X[cens,]) %*% matrix(r,ncol=1)
      H.obs[1:ncol(X),1:ncol(X)]= -(t(X[!cens,]) %*% X[!cens,])
      H.cens[1:ncol(X),1:ncol(X)]= -(t(X[cens,]) %*% (D*X[cens,]))
      H.obs[1:ncol(X),ncol(X)+1]= exprho * t(X[!cens,]) %*% matrix(y[!cens,1],ncol=1)
      H.cens[1:ncol(X),ncol(X)+1]= exprho * t(X[cens,]) %*% matrix(y[cens,1]*D,ncol=1)
      #H[1:ncol(X),ncol(X)+1]= exprho * t(X[!cens,]) %*% matrix(y[!cens,1],ncol=1) + exprho * t(X[cens,]) %*% matrix(y[cens,1]*D,ncol=1)
      #H[1:ncol(X),1:ncol(X)]= H.obs + H.cens
      H[ncol(X)+1,1:ncol(X)]= H[1:ncol(X),ncol(X)+1]
  } else {
      res= exprho*y[,1]
      r= invmillsnorm(-res[cens])
      D= Dz(res[cens])
  }
  nuncens= sum(!cens)
  g.obs[ncol(X)+1]= nuncens - exprho * sum(y[!cens,1]*res[!cens])
  g.cens[ncol(X)+1]= -exprho * sum(y[cens,1]*r)
  H.obs[ncol(X)+1,ncol(X)+1]=  - exprho * sum(y[!cens,1]*res[!cens]) - exprho^2 * sum(y[!cens,1]^2)
  H.cens[ncol(X)+1,ncol(X)+1]= - exprho * sum(y[cens,1]*r) - exprho^2 * sum(y[cens,1]^2 * D)
  #g[ncol(X)+1]= nuncens - exprho * sum(y[!cens,1]*res[!cens]) + exprho * sum(y[cens,1]*r)
  #H[ncol(X)+1,ncol(X)+1]= g[ncol(X)+1] - nuncens - exprho^2 * (sum(y[!cens,1]^2) + sum(y[cens,1]^2 * D))
  g= g.obs + g.cens; H= H.obs + H.cens
  return(list(g=g,g.obs=g.obs,g.cens=g.cens,H=H,H.obs=H.obs,H.cens=H.cens))
}


#pMOM + group Zellner + IG
dpmomgzellIG <- function(par, X, groups= 1:ncol(X), tau, taugroup, a=.01, b=.01, logscale=FALSE)  {
  alpha= head(par,-1); rho= tail(par,1); colX= length(alpha); ans=0
  if (ncol(X)>0) {
      ningroup= table(groups)
      isgroup= groups %in% (as.numeric(names(ningroup)[ningroup>1]))
      ans= sum(dmom(alpha[!isgroup], tau=tau, logscale=TRUE))
      if (any(isgroup)) {
          groupids= unique(groups[isgroup])
          for (j in groupids) {
              Sj= (t(X[,groups==j,drop=FALSE]) %*% X[,groups==j,drop=FALSE]) * (sum(groups==j)/taugroup)
              ans= ans + dmvnorm(par[groups==j], sigma= solve(Sj), log=TRUE)
          }
      }
  }
  ans= ans + mombf:::dinvgamma(exp(-2*rho), a/2, b/2, log=TRUE) + log(2) -2*rho
  if (!logscale) ans= exp(ans)
  return(ans)
}


#Gradient and hessian of log pMOM + log group Zellner + log IG
dpmomgzellIGgradhess <- function(par, X, groups= 1:ncol(X), tau, taugroup, a=.01, b=.01)  {
  alpha= head(par,-1); rho= tail(par,1); colX= length(alpha)
  grad= double(length(par)); hess= matrix(0,nrow=length(par),ncol=length(par))
  if (length(par)>1) {
      ningroup= table(groups)
      isgroup= groups %in% (as.numeric(names(ningroup)[ningroup>1]))
      grad[which(!isgroup)]= 2/par[which(!isgroup)] - par[which(!isgroup)]/tau
      diag(hess)[which(!isgroup)]= -2/par[which(!isgroup)]^2 - 1/tau
      if (any(isgroup)) {
          groupids= unique(groups[isgroup])
          for (j in groupids) {
              Sj= t(X[,groups==j,drop=FALSE]) %*% X[,groups==j,drop=FALSE]
              grad[groups==j]= -(sum(groups==j)/taugroup) * (Sj %*% matrix(par[groups==j],ncol=1))
              hess[groups==j,groups==j]= -(sum(groups==j)/taugroup) * Sj
          }
      }
  }
  grad[length(par)]= a - exp(2*par[length(par)]) * b
  hess[length(par),length(par)]= - exp(2*par[length(par)]) * 2*b
  return(list(grad=grad,H=hess))
}


#Group Zellner + group Zellner + IG
dgzellgzellIG <- function(par, X, groups= 1:ncol(X), tau, taugroup, a=.01, b=.01, logscale=FALSE)  {
  alpha= head(par,-1); rho= tail(par,1); colX= length(alpha); ans=0
  if (ncol(X)>0) {
      ans= 0
      groupids= unique(groups)
      for (j in groupids) {
          if (sum(groups==j)==1) {
              Sj= (t(X[,groups==j,drop=FALSE]) %*% X[,groups==j,drop=FALSE]) * (1/tau)
          } else {
              Sj= (t(X[,groups==j,drop=FALSE]) %*% X[,groups==j,drop=FALSE]) * (sum(groups==j)/taugroup)
          }
          ans= ans + dmvnorm(alpha[groups==j], sigma= solve(Sj), log=TRUE)
      }
  }
  ans= ans + mombf:::dinvgamma(exp(-2*rho), a/2, b/2, log=TRUE) + log(2) -2*rho
  if (!logscale) ans= exp(ans)
  return(ans)
}


#log-posterior under group Zellner + log group Zellner + log IG
fgzellgzellIG= function(par, y, X, groups=1:ncol(X), tau, taugroup, a=.01, b=.01) {
    l= loglikAFT(par, y=y, X=X)
    p= dgzellgzellIG(par, X=X, groups=groups, tau=tau, taugroup=taugroup, a=a, b=b, log=TRUE)
    return(l+p)
}


#Gradient and hessian of log-posterior under group Zellner + log group Zellner + log IG
fgzellgzellIGgradhess= function(par, y, X, groups=1:ncol(X), tau, taugroup, a=.01, b=.01) {
    l= loglikAFTgradhess(par, y=y, X=X)
    p= dgzellgzellIGgradhess(par, X=X, groups=groups, tau=tau, taugroup=taugroup, a=a, b=b)
    H= l$H + p$H
    H[lower.tri(H)]= t(H)[lower.tri(H)]
    ans= list(g=l$g + p$g, H=H)
    return(ans)
}


#Prior gradient and hessian of log group Zellner + log group Zellner + log IG
dgzellgzellIGgradhess <- function(par, X, groups= 1:ncol(X), tau, taugroup, a=.01, b=.01)  {
  alpha= head(par,-1); rho= tail(par,1); colX= length(alpha)
  grad= double(length(par)); hess= matrix(0,nrow=length(par),ncol=length(par))
  if (length(par)>1) {
      ningroup= table(groups)
      groupids= unique(groups)
      for (j in groupids) {
          Sj= (t(X[,groups==j,drop=FALSE]) %*% X[,groups==j,drop=FALSE])
          if (sum(groups==j)==1) {
            grad[groups==j]= -(1/tau) * (Sj %*% matrix(alpha[groups==j],ncol=1))
            hess[groups==j,groups==j]= -(1/tau) * Sj
          } else {
            sel= c(groups==j,FALSE)
            grad[sel]= -(1/taugroup) * (Sj %*% matrix(alpha[groups==j],ncol=1))
            hess[sel,sel]= -(sum(groups==j)/taugroup) * Sj
          }
      }
  }
  grad[length(par)]= a - exp(2*par[length(par)]) * b
  hess[length(par),1:(length(par)-1)]=  hess[1:(length(par)-1),length(par)]= 0
  hess[length(par),length(par)]= - exp(2*par[length(par)]) * 2*b
  return(list(grad=grad,H=hess))
}



### ADAPTATION OF JAVIER'S FUNCTION


marglhood.aft= function(y,X,a=.01,b=.01) {
  source('~/Dropbox/projects/fxrubio/DavidJavier/BVSC/R Codes/g2/VarSelection.R')
  if (class(y) != 'Surv') stop("class(y) must be 'Surv'")
  status <- y[,2]==1; n.index= y[,2]==0; y <- y[,1]; p= ncol(X); n= nrow(X)
  log.marglik.mom = AIC = BIC = log.mod.prior = models= vector()
  # Variable selection
  for(i in 1:p){
    subs = combn(p,i)
    nl = dim(subs)[2]; nv = dim(subs)[1]
    for(j in 1:nl){
      Xs <- X[,subs[,j]]
      models= c(models, paste(subs[,j],collapse=','))
      if(nv==1){
        lmr <- lm(y~Xs-1) # LSE estimator as initial values
        init0 <- c(lmr$coefficients/sigma(lmr),1/sigma(lmr))
        # MLE using NR
        init1<-Newton.Raphson2(y,Xs,init0,100,1e-12,grad.rlog.lik1,Hess.rlog.lik1,status=status)
        al <- -2*rlog.lik1(init1,y,Xs,status=status) +2*(1+1)
        bl <- -2*rlog.lik1(init1,y,Xs,status=status) +log(n)*(1+1)
        # MAP using NR
        MAP.mom <- Newton.Raphson2(y,Xs,init0,100,1e-12,grad.rlog.post.mom1,Hess.rlog.post.mom1,status=status)
        lml1 <- log.marg.lik.lap(y=y,X=Xs,MAP.mom,Hess.rlog.post.mom1,rlog.post.mom1,status=status)
        gam.prior <- log.prior.gamma(p,nv ,a,b)
      }
      if(nv>1){
        lmr <- lm(y~Xs-1)
        init0 <- c(lmr$coefficients/sigma(lmr),1/sigma(lmr))
        init1<-Newton.Raphson2(y,Xs,init0,100,1e-12,grad.rlog.lik,Hess.rlog.lik,status=status)
        al <- -2*rlog.lik(init1,y,Xs,status=status) +2*(dim(Xs)[2]+1)
        bl <- -2*rlog.lik(init1,y,Xs,status=status) +log(n)*(dim(Xs)[2]+1)
        MAP.mom <- Newton.Raphson2(y,Xs,init0,100,1e-12,grad.rlog.post.mom,Hess.rlog.post.mom,status=status)
        lml1 <- log.marg.lik.lap(y=y,X=Xs,MAP.mom,Hess.rlog.post.mom,rlog.post.mom,status=status)
        gam.prior <- log.prior.gamma(p,nv,a,b)
      }
      log.marglik.mom <- c(log.marglik.mom,lml1);
      AIC <- c(AIC,al); BIC <- c(BIC,al);
      log.mod.prior <- c(log.mod.prior,gam.prior)
    }
  }
  # Model posterior probabilities
  var.mom <- exp(log.marglik.mom + log.mod.prior - max(log.marglik.mom))
  post.prob.mom <- var.mom/sum(var.mom)
  ans= data.frame(models,post.prob.mom)
  ans= ans[order(post.prob.mom,decreasing=TRUE),]
  return(ans)
}



