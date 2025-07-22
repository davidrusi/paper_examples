
##########################################################################3
# FUNCTIONS TO PERFORM CROSS-VALIDATION
##########################################################################3


# Return BMA out-of-sample predictions using K-fold cross-validation
# Input
# - y: outcome variable, as passed onto modelSelection & modelSelection_eBayes
# - x: covariates, as passed onto modelSelection & modelSelection_eBayes
# - Z: meta-covariates, as passed onto modelSelection_eBayes
# - priorCoef: prior on the coefficients, as passed onto modelSelection & modelSelection_eBayes
# - priorDelta: only used if eBayes==FALSE. prior on the models, as passed onto modelSelection & modelSelection_eBayes
# - priorVar: prior on the variance parameter
kfoldCV.bma <- function(y, x, Z, priorCoef=momprior(), priorDelta=modelbbprior(), priorVar=igprior(.01, .01), eBayes,  K=10, seed, center=TRUE, scale=TRUE, method='auto', niter=5000, niter.mstep=1000, mc.cores=1, verbose=TRUE) {
  ## K-fold cross-validation for BMA predictions
  if (K > nrow(x)) stop("The number of folds cannot be larger than nrow(x)")
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  f <- function(k,...) {
    sel <- subset==k
    if (eBayes) {
      ms <- modelSelection_eBayes(y=y[!sel], x=x[!sel,,drop=FALSE], Z=Z, niter.mcmc=niter, niter.mstep=niter.mstep, verbose=FALSE)
      hyperpar <- ms$eBayes_hyperpar
    } else {
      ms <- modelSelection(y=y[!sel], x=x[!sel,,drop=FALSE], center=center, scale=scale, niter=niter, priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar, method=method, verbose=FALSE, ...)
    }
    pred <- predict(ms, newdata= x[sel,,drop=FALSE])[,'mean']
    if (verbose) cat('.')
    return(list(pred=pred, hyperpar=hyperpar))
  }
  if (verbose) cat('Running cross-validation')
  if (mc.cores > 1) {
    if ("parallel" %in% loadedNamespaces())  {
      allpred <- parallel::mclapply(1:K, f, mc.preschedule=FALSE)
    } else {
      stop("Did not find mclapply. Please load parallel package")
    }
  } else {
    allpred <- lapply(1:K, f)
 }
  pred= double(nrow(x))
  hyperpar= matrix(NA, nrow=K, ncol=nrow(allpred[[1]]$hyperpar))
  for (k in 1:K) {
    pred[subset==k] <- allpred[[k]]$pred
    if (eBayes) hyperpar[k,]= allpred[[k]]$hyperpar
  }
  if (verbose) cat('.')
  ans= list(pred=pred, hyperpar=hyperpar)
  return(ans)
}





