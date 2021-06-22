#FILE OBTAINED FROM https://github.com/SoyeonKimStat/Grouplasso
cv.grplasso <- function(x, y, index, model = LinReg(), nfolds=5, foldid=NULL, nlambda=20,parallel = FALSE) {
  library(grplasso)
  library(parallel)
  library(doParallel)
  if(is.null(foldid)) { foldid = sample(rep(seq(nfolds), length = length(y)))}
  lambda.max <- lambdamax(x, y, index = index, model = model)
  lambda.min = 0.0001 * lambda.max 
  
  lambda=exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda))
  
  ## Fit the solution path on the lambda grid
  se<-function(x) sqrt(var(x)/length(x)) 
  pred.mat <- matrix(NA, ncol=length(lambda),nrow=length(y))
  #from<-proc.time()

  if (parallel) {
    cl = makeCluster(nfolds, outfile="")
    registerDoParallel(cl)
    
    outlist = foreach(i = seq(nfolds),.combine = list,.multicombine = TRUE, .packages = c("grplasso")) %dopar% 
      {
        print(paste(i,"th fold",sep=""))
        which <- foldid == i
        y_sub <- y[!which]
        
        fit <- grplasso(x =x[!which, , drop = FALSE], y = y_sub, index = index, lambda = lambda, model = model)
        predict(fit, x[which,,drop=FALSE])
        
      }
    stopImplicitCluster()
    for(i in 1:nfolds) {
      which <- foldid == i
      pred.mat[which,]<- outlist[[i]]
    }
  } else {     
    for (i in seq(nfolds)) {
      print(paste(i,"th fold",sep=""))
      which <- foldid == i
      y_sub <- y[!which]
      
      fit <- grplasso(x =x[!which, , drop = FALSE], y = y_sub, index = index, lambda = lambda, model = model)
      pred.mat[which,]<- predict(fit, x[which,,drop=FALSE])
      
    }  
  }
  #to <- from-proc.time()
  mse <- sapply(1:length(lambda), function(curi) { 
    sq.error<-(y-pred.mat[,curi])^2
    c(cvm=mean(sq.error),cvse=se(sq.error))
  })
  mse2<- data.frame(t(mse))
  mse2$lambda <- lambda
  
  id<-which.min(mse2$cvm)
  lambda.min <- mse2$lambda[id]
  lambda.1se <- mse2$lambda[which.max(mse2$cvm <= mse2$cvm[id] +mse2$cvse[id])]
  
  test <- grplasso(x =x, y=y, index=index, lambda=c(lambda.1se,lambda.min),model=model)
  list(lambda=lambda, cvm=mse2$cvm, cvse=mse2$cvse, grplasso.fit=test, lambda.min=lambda.min,lambda.1se=lambda.1se,foldid=foldid)  
}


lassopost  <- function(y,x, verbose = FALSE) {
  #Run LASSO + post-selection inference on selected coefficients. lambda set via 10-fold cross-validation
  #Input
  # - y: response variable
  # - x: design matrix
  #
  # Output: matrix with parameter estimates, 95% CIs and P-values for the variables selected by LASSO (lambda set via cross-validation), using the post-selection inference method of Lee et al (2016)
  require(glmnet)
  require(selectiveInference)
  x= x[,apply(x,2,'sd') > 0] #remove intercept
  xstd= scale(x); ystd= scale(y) #important!
  if (is.null(colnames(x))) {
      colnames(xstd)= paste('X',1:ncol(xstd),sep='')
  } else {
      colnames(xstd)= colnames(x)
  }
  sigmahat= estimateSigma(x=xstd, y=ystd, standardize=FALSE)$sigmahat
  cvfit= cv.glmnet(x=xstd, y=ystd, nfolds=10)
  gfit= glmnet(x=xstd, y=ystd, standardize=FALSE, lambda=cvfit$lambda.min)
  b= coef(gfit)[-1]
  lcv= fixedLassoInf(x=xstd,y=ystd,beta=b,lambda=cvfit$lambda.min*nrow(x),sigma=sigmahat,verbose=verbose)
  sel.lasso= (b!=0)
  ci.lassopost= matrix(0,nrow=length(b),ncol=4)
  colnames(ci.lassopost)= c('estimate','ci.low','ci.up','pvalue'); rownames(ci.lassopost)= colnames(xstd)
  ci.lassopost[b != 0,]= cbind(lcv$coef0, lcv$ci, lcv$pv)
  ci.lassopost[b == 0,4]= NA
  return(ci.lassopost)
}


