
###########################################################################################
## 1. DEFINE SIMULATION ROUTINES
###########################################################################################

## MAIN SIMULATION ROUTINES ##
##############################

# Simulate data
# Input
# - simid: simulation id, used to set random number generator
# - theta: true parameter values
# - phi: true residual variance
# - n: sample size of the sample to be simulated
# Output
# - y: outcome
# - x: covariates
sim_data <- function(theta, n, rho, phi=1) {
  p <- length(theta)
  S <- diag(p); S[upper.tri(S)] <- rho; S[lower.tri(S)] <- rho
  x <- rmvnorm(n,sigma=S)
  x <- scale(x)
  y <- x %*% matrix(theta,ncol=1) + rnorm(n,0,sd=sqrt(phi))
  return(list(y=y, x=x))
}

# Simulate regression parameter values
# Let theta be the regression parameter
# - The logit prior probability logit P(theta[,j] != 0) = Z[j,] * w
#
#   where Z[,1] == 1 and Z[j,-1] ~ N(0, S), with S[i,i]= 1 and S[i,j]=rho
#
# - The non-zero entries in theta are uniformly spaced between theta.min and theta.max
#
# Input
# - p: number of covariates
# - w: coefficients used to generate theta
# - rho: correlation between Z's
#
# Output
# - Z: meta-covariates used to compute logit prior probabilities
# - priorpp: prior probability that theta != 0
# - theta: simulated regression parameter values
sim_parameters <- function(theta.min, theta.max, p, w, rho) {
  S <- diag(length(w)-1); S[upper.tri(S)] <- rho; S[lower.tri(S)] <- rho
  Z <- cbind(1,rmvnorm(p,sigma=S))
  priorpp <- 1 / (1 + exp(-Z %*% matrix(w, ncol=1)))
  nonzero <- (runif(length(priorpp)) < priorpp)
  theta <- rep(0, length(nonzero))
  theta[nonzero] <- seq(theta.min, theta.max, length=sum(nonzero))
  ans <- list(theta=theta, priorpp=priorpp, Z=Z)
  return(ans)
}


simBMA <- function(simid, p, n, eBayes, priorCoef, priorDelta, theta.min, theta.max, rho, phi=1, niter=5000, niter.mstep=1000, ...) {
  #Simulate data, run Bayesian model selection and averaging
  #
  # Simulate y ~ N(X theta, phi I), where X[i,] ~ N(0, S) with S[i,i]=1 and S[i,j]= rho
  # 
  # theta is simulated such that logit prior probability logit P(theta[,j] != 0) = Z[j,] * w
  #
  #      where Z[,1] == 1 and Z[j,-1] ~ N(0, S)
  #
  # The Empirical Bayes analysis only observes discretized Zcut = (Z > 0)
  #
  # The non-zero entries in theta are uniformly spaced between theta.min and theta.max
  #
  #
  # Input
  # - simid: simulation id, used to set random number generator
  # - p: number of covariates
  # - n: sample size of the sample to be simulated
  # - eBayes: if 'all' modelSelection_eBayes is called using all columns in Z, if 'intercept' only Z[,1] is used, if 'none' modelSelection is called (so Z is not used). If 'none' then priorDelta must be specified
  # - priorCoef: prior on the regression coefficients
  # - priorDelta: prior on the model space. Ignored unless eBayes == 'none'
  # - rho: covariates are draw from a multivariate Normal centered at 0, unit variances and all pairwise correlations equal to rho
  # - phi: true residual variance
  # - niter: number of MCMC iterations, passed onto modelSelection & modelSelection_eBayes
  # - niter.mstep: number of MCMC iterations in each M-step of the empirical Bayes algorithm
  # - ...: other arguments passed onto modelSelection_eBayes and modelSelection
  #
  # Output
  # - theta: regression parameter values used to simulate y ~ N(X theta, phi I)
  # - summary: a data frame with the posterior mode (most probable model a posteriori) and BMA results
  # - pcorrect: posterior mean of the proportion of correct variable inclusion/exclusions, and posterior probability of the true model (proportion of MCMC visits to the true model)
  # - hyperpar: if eBayes == TRUE, hyperpar contains the estimated hyper-parameters
  require(mvtnorm)
  set.seed(simid)
  if (missing(priorDelta)) priorDelta <- modelbbprior(1,1)
  simpar <- sim_parameters(theta.min=theta.min, theta.max=theta.max, p=p, w=w, rho=rho)
  sim <- sim_data(theta=simpar$theta, n=n, rho=rho, phi=phi)
  y <- sim$y; x <- sim$x
  #Zcut <- (simpar$Z > 0)
  #
  if (eBayes == 'all') {
    ms <- modelSelection_eBayes(y=y, x=x, Z=simpar$Z, niter.mcmc=niter, niter.mstep=niter.mstep, verbose=FALSE, ...)
    hyperpar <- ms$eBayes_hyperpar
  } else if (eBayes == 'intercept') {
    ms <- modelSelection_eBayes(y=y, x=x, Z=simpar$Z[,1,drop=FALSE], niter.mcmc=niter, niter.mstep=niter.mstep, verbose=FALSE, ...)
    hyperpar <- ms$eBayes_hyperpar
  } else if (eBayes == 'none') {
    ms <- modelSelection(y=y, x=x, niter=niter, priorCoef=priorCoef, priorDelta=priorDelta, verbose=FALSE, ...)
    hyperpar <- NULL
  } else {
    stop("Invalid vale of eBayes. It should be 'all', 'intercept' or 'none'")
  }
  b <- coef(ms) #BMA estimates and posterior inclusion probabilities
  b <- b[!(rownames(b) %in% c('intercept','phi')),]
  postsummary <- data.frame(selection=ms$postMode, b=b)
  names(postsummary) <- c('selection', 'estimate', 'ci.low', 'ci.up', 'margpp')
  pcorrect <- colMeans(t(ms$postSample)==(simpar$theta!=0))
  pcorrect <- c(propCorrectVars=mean(pcorrect), propTrueModel=mean(pcorrect==1))
  ans <- list(theta=simpar$theta, summary=postsummary, pcorrect=pcorrect, priorpp=simpar$priorpp, hyperpar=hyperpar)
  cat('.')
  return(ans)
}


simPenLhood <- function(simid, p, n, theta.min, theta.max, rho, phi=1, method) {
  #Simulate data and obtain parameter estimates with penalized likelihood
  # simid: simulation id, used to set random number generator
  # theta: true parameter values
  # phi: true residual variance
  # n: sample size of the sample to be simulated
  # rho: covariates are draw from a multivariate Normal centered at 0, unit variances and all pairwise correlations equal to rho
  # niter: number of MCMC iterations to explore the model space and also number of posterior samples for theta to be obtained
  # method: set method=='SCAD' for SCAD, method=='LASSO' for LASSO
  require(mvtnorm)
  set.seed(simid)
  simpar <- sim_parameters(theta.min=theta.min, theta.max=theta.max, p=p, w=w, rho=rho)
  sim <- sim_data(theta=simpar$theta, n=n, rho=rho, phi=phi)
  y <- sim$y; x <- sim$x
  #
  if (method=='SCAD') {
    require(ncvreg)
    mode <- scadcvfit(y=y,x=x)
  } else if (method=='LASSO') {
    require(lars)
    require(parcor)
    mode <- larscvfit(y=y,x=x)
  } else if (method=='ADALASSO') {
    require(parcor)
    mode <- adalassocvfit(y=y, x=x)
  } else {
    stop("Invalid value of argument method")
  }
  summary <- data.frame(selection=(mode!=0), estimate=mode, margpp=rep(NA,length(mode)))
  pcorrect <- mean((mode!=0)==(simpar$theta!=0))
  pcorrect <- c(propCorrectVars=pcorrect,propTrueModel=(pcorrect==1))
  ans <- list(theta=simpar$theta, summary=summary, pcorrect=pcorrect, priorpp=simpar$priorpp)
  cat('.')
  return(ans)
}


scadcvfit <- function(y,x) {
  #Auxiliary function fitting SCAD via 10-fold cross-validation
  cvscad <- cv.ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
  beta <- cvscad$fit$beta[-1,cvscad$min]
}

larscvfit <- function(y,x) {
  #Auxiliary function fitting LASSO via 10-fold cross-validation
  mylars(X=x,y=y,k=10,use.Gram=TRUE,normalize=TRUE)$coef
}


adalassocvfit <- function(y,x) {
  #Auxiliary function fitting Adaptive LASSO via 10-fold cross-validation
  adalasso(X=x, y=y, k=10, use.Gram=TRUE)$coefficients.adalasso
}



simOracle <- function(simid,theta,n,rho=0,phi) {
  #Simulate data and obtain Oracle predictions using MLE and its 95% CI obtained under the true model
  require(mvtnorm)
  set.seed(simid)
  p <- length(theta)
  S <- diag(p); S[upper.tri(S)] <- rho; S[lower.tri(S)] <- rho
  x <- rmvnorm(n,sigma=S)
  y <- x %*% matrix(theta,ncol=1) + rnorm(n,0,sd=sqrt(phi))
  #
  sel <- theta!=0
  lm1 <- lm(y ~ -1 + x[,theta!=0,drop=FALSE])
  mode <- rep(0,length(theta))
  mode[sel] <- coef(lm1)
  ci <- matrix(0,nrow=length(theta),ncol=2)
  colnames(ci) <- c("X2.5.","X97.5.")
  ci[sel,] <- confint(lm1,level=0.95)
  summary <- data.frame(selection=as.numeric(theta!=0), margpp=as.numeric(theta!=0), mode=mode, mean=mode, ci)
  ans <- list(theta=theta,summary=summary,pcorrect=1)
  cat('.')
  return(ans)
}


summarizeSim <- function(sim, threshold) {
  # Auxiliary function to extract FDR, power and MSE from simulation output
  # For Bayesian methods, the average posterior inclusionp probability for true zeroes and true non-zeroes is also returned
  mse <- (sim$theta - sim$summary$estimate)^2
  if (missing(threshold)) {
      sel <- (sim$summary$estimate != 0) #for frequentist methods, selected model is given by point estimate
  } else {
      sel <- (sim$summary[,'margpp'] > threshold) #for Bayesian methods, PIP is thresholded
  }
  fp <- sum(sel & (sim$theta == 0))
  tp <- sum(sel & (sim$theta != 0))
  fdr <- ifelse(fp + tp == 0, 0, fp / (fp + tp))
  pow <- tp / sum(sim$theta != 0)
  if (!is.null(sim$summary$margpp)) {
    meanpp0 <- mean(sim$summary[sim$theta==0,'margpp'])
    meanpp1 <- mean(sim$summary[sim$theta!=0,'margpp'])
  } else {
    meanpp <- c(NA, NA)
  }
  ans <- c(fdr=fdr, pow=pow, mse=sum(mse), meanpp0=meanpp0, meanpp1=meanpp1)
  return(ans)
}




##########################################################################3
# 2. FUNCTIONS TO PERFORM CROSS-VALIDATION
##########################################################################3


# Return BMA out-of-sample predictions using K-fold cross-validation
# Input
# - y: outcome variable, as passed onto modelSelection & modelSelection_eBayes
# - x: covariates, as passed onto modelSelection & modelSelection_eBayes
# - Z: meta-covariates, as passed onto modelSelection_eBayes
# - priorCoef: prior on the coefficients, as passed onto modelSelection & modelSelection_eBayes
# - priorDelta: only used if eBayes==FALSE. prior on the models, as passed onto modelSelection & modelSelection_eBayes
# - eBayes: if TRUE modelSelection_eBayes is called, else modelSelection is called
# - K: number of folds (set K=length(y) for leave-one-out cross-validation)
# - seed: optionally, you may specify the random number generator seed
# - niter: number of MCMC iterations, passed onto modelSelection & modelSelection_eBayes
# - niter.mstep: number of MCMC iterations in each M-step of the empirical Bayes algorithm
kfoldCV.bma <- function(y, x, Z, priorCoef=momprior(), priorDelta=modelbbprior(), eBayes,  K=10, seed, niter=5000, niter.mstep=1000, mc.cores=1, verbose=TRUE) {
  ## K-fold cross-validation for BMA predictions
  if (K > nrow(x)) stop("The number of folds cannot be larger than nrow(x)")
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  f <- function(k,...) {
    sel <- subset==k
    if (eBayes) {
      ms <- modelSelection_eBayes(y=y[!sel], x=x[!sel,,drop=FALSE], Z=Z, niter.mcmc=niter, niter.mstep=niter.mstep, verbose=FALSE, ...)
      hyperpar <- ms$eBayes_hyperpar
    } else {
      ms <- modelSelection(y=y[!sel], x=x[!sel,,drop=FALSE], niter=niter, priorCoef=priorCoef, priorDelta=priorDelta, verbose=FALSE, ...)
      hyperpar <- NULL
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
  if (eBayes) {
      hyperpar= matrix(NA, nrow=K, ncol=nrow(allpred[[1]]$hyperpar))
  } else {
      hyperpar= matrix(NA, nrow=K, ncol=0)
  }
  for (k in 1:K) {
    pred[subset==k] <- allpred[[k]]$pred
    if (eBayes) hyperpar[k,]= allpred[[k]]$hyperpar
  }
  if (verbose) cat('.')
  ans= list(pred=pred, hyperpar=hyperpar)
  return(ans)
}





