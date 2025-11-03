###############################################################################
##
## THIS SCRIPT REPRODUCES THE RESULTS FROM SECTION 5.1 "SIMULATION WITH INDEPENDENT ERRORS"
##
###############################################################################

library(mombf)
library(parallel)
library(tidyverse)
#library(devtools)
#devtools::install_github(repo = "skdeshpande91/VCBART/VCBART")
library(VCBART)
library(here)
setwd("~/github/localnulltest")



###############################################################################
## DEFINE FUNCTIONS
###############################################################################

truemean= function(x,z) {
    ans= double(nrow(x))
    group1= (x[,1]==1)
    ans[group1]= ifelse(z[group1] <=0, cos(z[group1]), 1)
    ans[!group1]= ifelse(z[!group1]<=0, cos(z[!group1]), 1/(z[!group1]+1)^2)
    return(ans)
}


#Simulate data where y = truemean(x,z) + e and
# x[,1] is binary with 0.5 observations equal to 1 and 0
# x[,j] ~ N(x[,1], 1) for j=2,...,p
# e ~ N(0, sd.error) iid
#
# Input
# - seed: random number generator seed
# - n: sample size
# - p: number of variables
# - sd.error: standard deviation of the errors
# Output
# - y: simulated outcome values
# - x: simulated covariate values
# - z: simulated z values (a uniform grid in [-3,3])
# - m: true mean of E
simdata= function(seed, n, p, sd.error) {
    if (!missing(seed)) set.seed(seed)
    x= matrix(NA, nrow=n, ncol=p)
    x[,1]= rep(0:1,c(n/2,n/2))
    if (p>1) { for (j in 2:p) x[,j]= x[,1] + rnorm(n) }
    z= rep(seq(-3,3,length=n/2), 2)
    m= truemean(x,z)
    y= truemean(x,z) + rnorm(n, 0, sd=sd.error)
    ans= list(y=y, x=x, z=z, m=m)
    return(ans)
}


#Simulate data as in simdata and compute posterior probabilities for local null tests computed under several basis functions
# Input: same as for simdata. Also:
# - mc.cores: number of cores to use for parallel processing (if it causes issues, set mc.cores=1)
# - uncut: if TRUE, results are obtained for cut and uncut basis. If FALSE, only for cut basis (our recommended methodology)
# Output: matrix with z coordinates and posterior probabilities under a cut 0-degree basis (cut0) and standard cubic B-spline (uncut3)
# Note: results for a cut cubic basis were very similar than for the 0-degree basis, hence not returned
sim_localnulltest= function(seed, n, p, sd.error=0.25, mc.cores, localknots, uncut=TRUE, ...) {
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    if (missing(localknots)) {
      nlocalknots= c(7,9,11)
      fit0= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=0, nlocalknots=nlocalknots, mc.cores = mc.cores, ...) #0-degree, multi-resolution
    } else {
      fit0= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=0, localknots=localknots, mc.cores = mc.cores, ...) #0-degree, multi-resolution
    }
    b0= coef(fit0)
    mse.cut0= sum((predict(fit0) - sim$m)^2)       #MSE in estimating the mean
    if (uncut) {
        #uncut cubic, multi-resolution
        if (missing(localknots)) {
            nlocalknots= c(7,9,11)
            fit.uncut= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=3, usecutbasis=FALSE, nlocalknots=nlocalknots, mc.cores = mc.cores, ...)
        } else {
            fit.uncut= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=3, usecutbasis=FALSE, localknots=localknots, mc.cores = mc.cores, ...)           
        }
        b.uncut= coef(fit.uncut)
        mse.uncut3= sum((predict(fit.uncut) - sim$m)^2)
        margpp= cbind(b0[,c('covariate','z1')], b0[,'margpp'], b.uncut[,'margpp'])
    } else {
        mse.uncut3= NA
        margpp= cbind(b0[,c('covariate','z1')], b0[,'margpp'], NA)
    }
    #Return output
    colnames(margpp)= c('covariate','z','cut0','uncut3')
    ans= list(margpp=margpp, mse=c(cut0=mse.cut0, uncut3=mse.uncut3), pp_localknots= fit0$pp_localknots)
    return(ans)
}



###############################################################################
## GAM-BASED LOCAL NULL TESTS
###############################################################################

# Return simultaneous confidence interval from a gam fit, estimated via simulations
# INPUT
# - fit: object returned by gam
# - parm: name of the smooth term in fit for which we want the confidence intervals
# - level: confidence level
# - nsim: number of simulations (multivariate normal draws)
# - data: data.frame used to fit the object
# OUTPUT: point estimates, lower and upper ends of the confidence interval. This is similar to using predict.gam, except that the standard error is multiplied by a critical level related to the maximum of Gaussians, rather than the standard normal quantile (e.g. 1.96 for level=0.95), resulting in wider intervals
# NOTES
# Function adapted from package gratia, following the description in https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited
#
# REFERENCES
# Marra, G., and Wood, S. N. (2012). Coverage properties of confidence intervals for generalized additive model components. Scandinavian journal of statistics, theory and applications 39, 53–74. doi:10.1111/j.1467-9469.2011.00760.x.
# Nychka, D. (1988). Bayesian confidence intervals for smoothing splines. Journal of the American Statistical Association 83, 1134–1143. doi:10.1080/01621459.1988.10478711.
# Ruppert, D., Wand, M. P., and Carroll, R. J. (2003). Semiparametric regression. Cambridge University Press.
simultaneous_ci= function(fit, parm, level=0.95, nsim=10000, data) {
    #Auxiliary function
    sim_interval <- function(smooth, level, data, se.fit) {
      start <- smooth[["first.para"]]
      end <- smooth[["last.para"]]
      para.seq <- start:end
      Cg <- PredictMat(smooth, data)
      simDev <- Cg %*% t(buDiff[, para.seq])
      absDev <- abs(sweep(simDev, 1L, se.fit, FUN = "/")) #David modified
      masd <- apply(absDev, 2L, max, na.rm=TRUE)
      unname(quantile(masd, probs = level, type = 8))
    }
    V <- gratia:::get_vcov(fit, unconditional = FALSE)  # covariance matrix of fitted object
    buDiff <- mvnfast::rmvn(n = nsim, mu = rep(0, nrow(V)), sigma = V, ncores = 1)     ## simulate un-biased deviations given bayesian covar matrix
    #Select smooth term specified by parm
    S <- gratia:::smooths(fit)
    select <- gratia:::check_user_select_smooths(smooths = S, select = parm, partial_match = FALSE)
    #out <- gratia::smooth_estimates(fit, select = S, n = n, data = data, partial_match = FALSE)
    smooth <- gratia:::get_smooth(fit, parm)
    #Find critical level
    pred= predict(fit, type='iterms', se.fit=TRUE, terms=parm) #David added
    crit <- sim_interval(smooth, level = level, data = data, se.fit= pred$se.fit)
    #Return simultaneous confidence interval
    ans= data.frame(pred$fit, pred$fit - crit * pred$se.fit, pred$fit + crit * pred$se.fit)
    names(ans)= c('estimate','ci.low','ci.up')
    return(ans)
}


sim_localnulltest_gam= function(seed, n, p, sd.error=0.25, k, level=0.95) {
    require(mgcv)
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    df= data.frame(y= sim$y, x= sim$x, z= sim$z)
    names(df)= c("y", paste("x",1:p,sep=""), "z")
    if (missing(k)) {  #use gam's default number of knots
        txt_main= "y ~ s(z)"
        txt_int= paste(" + s(z, by=x",1:ncol(sim$x),")",sep='',collapse="")
    } else {           #use user-specified number of knots k
        txt_main= "y ~ s(z, k=k)"
        txt_int= paste(" + s(z, k=k, by=x",1:ncol(sim$x),")",sep='',collapse="")
    }
    f= as.formula(paste(txt_main, txt_int))
    fit= gam(f, data=df)
    #pred= predict(fit, type='lpmatrix', terms='s(z, by=x[, 1])') #shows that by default, k=12 knots (basis dimension is 9)
    mse= sum((predict(fit) - sim$m)^2)  
    ans= vector("list", p)
    for (i in 1:p) {
        terms= paste("s(z):x",i,sep="")
        ci= simultaneous_ci(fit, parm=terms, level=level, nsim=10^4, data=df)
        ci= data.frame(i, sim$z, ci)
        names(ci)= c('varname','z','estimate','ci.low','ci.up')
        ans[[i]]= ci[sim$x[,1] == 1,] #select 2nd half of observations, to avoid duplicating z's
    }
    ans= do.call(rbind, ans)
    ans= list(ci=ans, mse=mse)
    return(ans)
}



###############################################################################
## FUSED LASSO BASED LOCAL NULL TESTS
###############################################################################

sim_fusedlasso= function(seed, n, p, sd.error=0.25, mc.cores, ...) {
    require(genlasso)
    require(igraph)
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    fit= localnulltest_fusedlasso(y=sim$y, x=sim$x, z=sim$z, localgridsize=100, nbaseknots=20, nlocalknots=7, basedegree=3)
    return(fit)
}


#Return BIC for the best lambda parameter in fused LASSO, for a given value of the other regularization parameter gamma
fusedlasso_bic= function(loggamma, y, X, graph, only_minbic=TRUE, ebic) {
    fit= genlasso::fusedlasso(y=y, X=X, gamma = exp(loggamma), graph=graph, verbose = FALSE)
    #Obtain BIC
    b= coef(fit)$beta
    npar= colSums(b != 0)
    bic= nrow(w) * log(colSums((y - fit$fit)^2) / nrow(w)) + nrow(w) * (log(2 * pi) + 1) + log(nrow(w)) * npar
    if (ebic) bic = bic + npar * log(ncol(X))
    #if (ebic) bic = bic + lchoose(ncol(X), npar)
    bic= data.frame(lambda=fit$lambda, npar=npar, bic=bic) |>
           filter(!is.infinite(lambda))
    b= b[,!is.infinite(fit$lambda)]
    if (only_minbic) {
        ans = bic[which.min(bic$bic),'bic']
    } else {
        sel = which.min(bic$bic)
        ans = list(bic= bic[sel,], beta=b[,sel])
    }
    return(ans)
}


#Estimate local covariate effects by applying fused LASSO to our cut basis
localnulltest_fusedlasso= function(y, x, z, localgridsize=100, localgrid, nbaseknots=20, nlocalknots=7, basedegree=3) {
    require(dplyr)
    cutdegree= 0
    #Check & format input requirements
    check= mombf:::checkargs_localnulltest(y=y, x=x, z=z)
    y= check$y; x= check$x; z= check$z
    # Define local tests
    if (missing(localgrid)) localgrid= mombf:::define_localgrid(localgridsize=localgridsize, z=z)
    # Define knots and corresponding knot-based testing regions
    kk= mombf:::define_knots_localnulltest(z=z, localgrid=localgrid, nbaseknots=nbaseknots, basedegree=basedegree, nlocalknots=nlocalknots, localknots=NULL)
    knots= kk$knots; regionbounds= kk$regionbounds; region=kk$region; regioncoord= kk$regioncoord; testov= kk$testov; testxregion= kk$testxregion; testIntervals= kk$testIntervals
    #Create design matrix & run least-squares regression
    desnew= mombf:::estimationPoints(x=x, regioncoord=regioncoord, regionbounds=regionbounds, testov=testov) #Points at which local effects will be estimated
    des= mombf:::createDesignLocaltest(x=rbind(x,desnew$x), z=rbind(z,desnew$z), y=y, region=c(region,desnew$region), regionbounds=regionbounds, basedegree=basedegree, cutdegree=cutdegree, knots=knots, usecutbasis=TRUE, useSigma=FALSE)
    w= des$w[1:nrow(x),]; wnew= des$w[-1:-nrow(x),]
    #Define edges between neighboring local tests
    sel= (des$w1varname[-1] == des$w1varname[-length(des$w1varname)])
    edges1= ((des$ncolw0+1):(ncol(des$w)-1))[sel]
    edges2= edges1 + 1
    edges= as.vector(rbind(edges1, edges2))
    gr = graph(edges=edges,directed=FALSE)
    #Fused LASSO
    options('warn'=-1)
    f= function(loggamma) fusedlasso_bic(loggamma, y=y, X=w, graph=gr, ebic=TRUE)
    opt= optimize(f, interval=c(log(.001), log(100)), tol=0.01)
    loggamma = opt$minimum
    #system.time(loggamma <- cmna::goldsectmin(f, log(0.001), log(100), tol=.01))
    fit= fusedlasso_bic(loggamma, y=y, X=w, graph=gr, only_minbic=FALSE, ebic=TRUE)
    b= matrix(fit$beta, ncol=1)
    options('warn'=0)
    #Select estimates corresponding to local covariate effects
    b= data.frame(varname=des$w1varname, estimate= b[-1:-des$ncolw0,])
    b$region= sapply(strsplit(des$vargroups[-1:-des$ncolw0], split="\\."), "[", 2)
    b$region= as.numeric(sub("R","",b$region))   
    #Indicate region coordinates and format output
    nregions= length(regionbounds[[1]]) - 1
    region_intervals= matrix(NA, nrow=nregions, ncol=2)
    colnames(region_intervals)= c("zmin","zmax")
    for (i in 1:nregions) region_intervals[i,]= regionbounds[[1]][c(i,i+1)]
    region_intervals= data.frame(region=1:nregions, region_intervals)
    b= merge(b, region_intervals, by='region')
    b= b[order(b$varname, b$region),]
    b= select(b, varname, region, zmin, zmax, estimate) |>
         group_by(varname, region) |>
         summarize(zmin=mean(zmin), zmax=mean(zmax), estimate=mean(estimate)) |>
         as.data.frame()
    return(b)
}



###############################################################################
## P-VALUE BASED LOCAL NULL TESTS
###############################################################################


sim_localnulltest_pvalue= function(seed, n, p, sd.error=0.25, p.adjust.method='BH', ...) {
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    fit= localnulltest_pvalue(y=sim$y, x=sim$x, z=sim$z, p.adjust.method=p.adjust.method, ...)
    return(fit)
}

localnulltest_pvalue= function(y, x, z, p.adjust.method='BH', localgridsize=100, localgrid, nbaseknots=20, nlocalknots=7, basedegree=3) {
    cutdegree= 0
    #Check & format input requirements
    check= mombf:::checkargs_localnulltest(y=y, x=x, z=z)
    y= check$y; x= check$x; z= check$z
    # Define local tests
    if (missing(localgrid)) localgrid= mombf:::define_localgrid(localgrid=localgrid, localgridsize=localgridsize, z=z)
    # Define knots and corresponding knot-based testing regions
    kk= mombf:::define_knots_localnulltest(z=z, localgrid=localgrid, nbaseknots=nbaseknots, basedegree=basedegree, nlocalknots=nlocalknots, localknots=NULL)
    knots= kk$knots; regionbounds= kk$regionbounds; region=kk$region; regioncoord= kk$regioncoord; testov= kk$testov; testxregion= kk$testxregion; testIntervals= kk$testIntervals
    #Create design matrix & run least-squares regression
    desnew= mombf:::estimationPoints(x=x, regioncoord=regioncoord, regionbounds=regionbounds, testov=testov) #Points at which local effects will be estimated
    des= mombf:::createDesignLocaltest(x=rbind(x,desnew$x), z=rbind(z,desnew$z), y=y, region=c(region,desnew$region), regionbounds=regionbounds, basedegree=basedegree, cutdegree=cutdegree, knots=knots, usecutbasis=TRUE, useSigma=FALSE)
    w= des$w[1:nrow(x),]; wnew= des$w[-1:-nrow(x),]
    fit= lm(y ~ -1 + w)
    mysummary= summary(fit)
    #Store pvalues
    pvals= matrix(NA, nrow=ncol(w), ncol=2); colnames(pvals)= c('estimate','pvalue')
    pvals[!mysummary$aliased,]= mysummary$coef[,c(1,4)]
    rownames(pvals)= names(coef(fit))
    pvals= data.frame(varname=des$w1varname, pvals[-1:-des$ncolw0,])
    pvals$region= sapply(strsplit(des$vargroups[-1:-des$ncolw0], split="\\."), "[", 2)
    pvals$region= as.numeric(sub("R","",pvals$region))
    #Indicate region coordinates and format output
    nregions= length(regionbounds[[1]]) - 1
    region_intervals= matrix(NA, nrow=nregions, ncol=2)
    colnames(region_intervals)= c("zmin","zmax")
    for (i in 1:nregions) region_intervals[i,]= regionbounds[[1]][c(i,i+1)]
    region_intervals= data.frame(region=1:nregions, region_intervals)
    pvals= merge(pvals, region_intervals, by='region')
    pvals= pvals[order(pvals$varname, pvals$region),]
    pvals= pvals[,c("varname","region","zmin","zmax","estimate","pvalue")]
    pvals$pvalue.adjusted= p.adjust(pvals$pvalue, method=p.adjust.method)
    return(pvals)
}


###############################################################################
## BART FUNCTIONS
###############################################################################

sim_bart= function(seed, n, p, sd.error=0.25, mc.cores) {
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    bartfit= mc.gbart(x.train= cbind(sim$z,sim$x), y.train=sim$y, ntree=200, mc.cores=mc.cores, seed=1234) #default is 200 trees
    pp.used= colMeans(bartfit$varcount[,-1] > 0) #posterior probability of the covariate being used in >=1 trees
    mean.usage= colMeans(bartfit$varcount[,-1] / 200) #posterior mean of the proportion of trees (out of 200) in which each covariate is used
    ans= list(pp.used= pp.used, mean.usage=mean.usage)
    return(ans)
}


###############################################################################
## VC-BART FUNCTIONS
###############################################################################


#Format VCBART posterior summaries.
#Input: output of summarize_beta
#Ouput: data.frame indicating, for each variable and z coordinate, the posterior mean and 0.95 posterior interval. 
format_vcbart_beta_summary= function(beta, ztest) {
    ans= matrix(nrow=0, ncol=5)
    colnames(ans)= c('variable','z','mean','low','up')
    for (j in 2:dim(beta)[3]) {
      ans= rbind(ans, cbind(j-1, ztest, beta[,,j]))  #posterior mean and 0.95 interval for beta[,j] at all values of z
    }
    ans= data.frame(ans)
    ans$variable= factor(ans$variable)
    return(ans)
}


#Run VCBART on simulated data
sim_vcbart= function(seed, n, p, sd.error=0.25, cutoff=1, mc.cores) {
    if (cutoff > 100) stop("cutoff must be <=100")
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    train= 1:n  #use all data as training
    test= integer(0)
    xtrain= sim$x; xtest= sim$x[test,,drop=FALSE]
    cutpoints= list(seq(min(sim$z), max(sim$z), length=1000))
    nd= 1000 #number of MCMC iterations
    ztest= matrix(seq(min(sim$z), max(sim$z), length=100), ncol=1)
    xtest= matrix(1, nrow=nrow(ztest), ncol=ncol(sim$x))
    #Run VCBART and obtain posterior summaries
    bartfit= VCBART_ind(Y_train=sim$y, X_train=xtrain, Z_cont_train=matrix(sim$z,ncol=1), ni_train=rep(1,length(train)), subj_id_train=1:length(train), X_test=xtest, Z_cont_test=ztest, nd=nd, cutpoints=cutpoints, verbose=FALSE)
    betahat= summarize_beta(bartfit$betahat.train)
    betahat.test= summarize_beta(bartfit$betahat.test)
    betahat.test= format_vcbart_beta_summary(betahat.test, ztest=ztest)
    #MSE for true expectation
    muhat= betahat[,1,1] + rowSums(xtrain * betahat[,1,-1])
    mse= sum((muhat - sim$m[train])^2)  #MSE in estimating the mean
    #Return output
    ans= list(beta.ci=betahat.test, mse=mse)
    return(ans)
}



###############################################################################
## RUN SIMULATIONS
###############################################################################

#Simulations with n=100, Normal shrinkage prior on coefficients
n=100; p=10; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))

#Simulations with n=1000, Normal shrinkage prior on coefficients
n=1000; p=10; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))


#Simulations with n=100, ICAR+ prior on coefficients
n=100; p=10; sd.error=0.25; nsims= 100; priorCoef= icarplusprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef, uncut=FALSE); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_icarplus_n", n,".RData")))

#Simulations with n=1000, ICAR+ prior on coefficients
n=1000; p=10; sd.error=0.25; nsims= 100; priorCoef= icarplusprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef, uncut=FALSE); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_icarplus_n", n,".RData")))


###############################################################################
## RUN SIMULATIONS FOR FUSED LASSO
###############################################################################

#Simulations with n=100
n=100; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_fusedlasso(seed=i, n=n, p=p, sd.error=sd.error); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_fusedlasso_n", n,".RData")))

#Simulations with n=1000
n=1000; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_fusedlasso(seed=i, n=n, p=p, sd.error=sd.error); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_fusedlasso_n", n,".RData")))


###############################################################################
## RUN SIMULATIONS FOR LEAST SQUARES WITH P-VALUE ADJUSTMENT
###############################################################################

#Simulations with n=100
n=100; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest_pvalue(seed=i, n=n, p=p, sd.error=sd.error, p.adjust.method='BH'); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_ols_n", n,".RData")))

#Simulations with n=1000
n=1000; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest_pvalue(seed=i, n=n, p=p, sd.error=sd.error, p.adjust.method='BH'); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_ols_n", n,".RData")))


###############################################################################
## RUN SIMULATIONS FOR GAM
###############################################################################

#Simulations with n=100, default knots (k=12)
n=100; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest_gam(seed=i, n=n, p=p, sd.error=sd.error, level=0.95); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_gam_n", n,".RData")))


#Simulations with n=1000, default knots (k=12)
n=1000; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest_gam(seed=i, n=n, p=p, sd.error=sd.error, level=0.95); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_gam_n", n,".RData")))


#Simulations with n=100, k=24 knots
n=100; p=10; sd.error=0.25; nsims= 100; k=24
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest_gam(seed=i, n=n, p=p, sd.error=sd.error, level=0.95, k=k); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_gam_knots",k,"_n", n,".RData")))


#Simulations with n=1000, k=24 knots
n=1000; p=10; sd.error=0.25; nsims= 100; k=24
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest_gam(seed=i, n=n, p=p, sd.error=sd.error, level=0.95, k=k); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_gam_knots",k,"_n", n,".RData")))



###############################################################################
## RUN SIMULATIONS FOR VC-BART
###############################################################################

#Simulations with n=100
n=100; p=10; sd.error=0.25; nsims= 100; cutoff=1; mc.cores= 4
sim.vcbart= mclapply(1:nsims, function(i) { ans= sim_vcbart(seed=i, n=n, p=p, sd.error=sd.error, cutoff=cutoff); return(ans) }, mc.cores=mc.cores)
save(sim.vcbart, file=here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))

#Simulations with n=1000
n=1000; p=10; sd.error=0.25; nsims= 100; cutoff=1; mc.cores= 4
sim.vcbart= mclapply(1:nsims, function(i) { ans= sim_vcbart(seed=i, n=n, p=p, sd.error=sd.error, cutoff=cutoff); return(ans) }, mc.cores=mc.cores)
save(sim.vcbart, file=here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))





###############################################################################
## SUMMARIZE TYPE I ERROR AND POWER
###############################################################################


# n=100, normalid prior
load("code/simulation_iid_output/sim_n100.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

df= data.frame(id, pow0, pow3) |>
  mutate(region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)

tab= summarize(df, cut0=mean(pow0), uncut3=mean(pow3))

# n=100, ICAR+ prior
load("code/simulation_iid_output/sim_icarplus_n100.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pow0= rowMeans(pp0 > 0.95)

df= data.frame(id, pow0) |>
  mutate(region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)

tab.icarplus= summarize(df, cut0.icarplus=mean(pow0))

# n=100, fused LASSO
load("code/simulation_iid_output/sim_fusedlasso_n100.RData")
bhat= do.call(cbind, lapply(sim, function(z) z$estimate))
pow= rowMeans(abs(bhat) > 0.1, na.rm=TRUE)
df= data.frame(sim[[1]][,c('varname','zmin','zmax')], pow) |>
  mutate(region= factor(paste0('(',zmin,',',zmax,']')), covariate1= (varname=='x1')) |>
  select(-zmin, -zmax) |>
  group_by(covariate1, region)
tab.fusedlasso= summarize(df, fusedlasso=mean(pow, na.rm=TRUE))


# n=100, BH
load("code/simulation_iid_output/sim_ols_n100.RData")
pval= do.call(cbind, lapply(sim, function(z) z$pvalue.adjusted))
pow= rowMeans(pval < .05, na.rm=TRUE)

df= data.frame(sim[[1]][,c('varname','zmin','zmax')], pow) |>
  mutate(region= factor(paste0('(',zmin,',',zmax,']')), covariate1= (varname=='x1')) |>
  select(-zmin, -zmax) |>
  group_by(covariate1, region)

tab.pval= summarize(df, pval=mean(pow, na.rm=TRUE))


# n=100, GAM default knots
load("code/simulation_iid_output/sim_gam_n100.RData")
pow.gam= rowMeans(sapply(sim, function(z) rej= (z$ci$ci.low > 0) | (z$ci$ci.up < 0)))

df= data.frame(sim[[1]]$ci[,c('varname','z')], pow.gam) |>
  mutate(region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3.01)), covariate1= (varname==1)) |>
  group_by(covariate1, region)

tab.gam= summarize(df, pow.gam= mean(pow.gam))


# n=100, GAM k=24 knots
load("code/simulation_iid_output/sim_gam_knots24_n100.RData")
pow.gam= rowMeans(sapply(sim, function(z) rej= (z$ci$ci.low > 0) | (z$ci$ci.up < 0)))

df= data.frame(sim[[1]]$ci[,c('varname','z')], pow.gam) |>
  mutate(region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3.01)), covariate1= (varname==1)) |>
  group_by(covariate1, region)

tab.gam= summarize(df, pow.gam= mean(pow.gam))



# n=100, VC-BART
n=100
load("code/simulation_iid_output/sim_vcbart_n100.RData")
names(sim.vcbart[[1]])

threshold= 0
pow.vcbart= rowMeans(sapply(sim.vcbart, function(z, threshold=threshold) rej= (z$beta.ci[,'low'] > 0) | (z$beta.ci[,'up']<0)))

df= data.frame(sim.vcbart[[1]]$beta.ci[,c('variable','z')], pow.vcbart) |>
  mutate(region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3.01)), covariate1= (variable==1)) |>
  group_by(covariate1, region)

tab.vcbart= summarize(df, pow.vcbart= mean(pow.vcbart))

tab= cbind(tab, tab.icarplus[,'cut0.icarplus'], tab.fusedlasso[,'fusedlasso'], tab.pval[,'pval'], tab.vcbart[,'pow.vcbart'], tab.gam[,'pow.gam'])
cbind(filter(tab, covariate1)[,-1], filter(tab, !covariate1)[,-1:-2])


#n=1000, normalid prior
load("code/simulation_iid_output/sim_n1000.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

df= data.frame(id, pow0, pow3) |>
  mutate(region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)

tab= summarize(df, cut0=mean(pow0), uncut3=mean(pow3))

# n=1000, ICAR+ prior
load("code/simulation_iid_output/sim_icarplus_n1000.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pow0= rowMeans(pp0 > 0.95)

df= data.frame(id, pow0) |>
  mutate(region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)

tab.icarplus= summarize(df, cut0.icarplus=mean(pow0))


# n=1000, fused LASSO
load("code/simulation_iid_output/sim_fusedlasso_n1000.RData")
bhat= do.call(cbind, lapply(sim, function(z) z$estimate))
pow= rowMeans(abs(bhat) > 0.1, na.rm=TRUE)
df= data.frame(sim[[1]][,c('varname','zmin','zmax')], pow) |>
  mutate(region= factor(paste0('(',zmin,',',zmax,']')), covariate1= (varname=='x1')) |>
  select(-zmin, -zmax) |>
  group_by(covariate1, region)
tab.fusedlasso= summarize(df, fusedlasso=mean(pow, na.rm=TRUE))


# n=1000, BH
load("code/simulation_iid_output/sim_ols_n1000.RData")
pval= do.call(cbind, lapply(sim, function(z) z$pvalue.adjusted))
pow= rowMeans(pval < .05, na.rm=TRUE)

df= data.frame(sim[[1]][,c('varname','zmin','zmax')], pow) |>
  mutate(region= factor(paste0('(',zmin,',',zmax,']')), covariate1= (varname=='x1')) |>
  select(-zmin, -zmax) |>
  group_by(covariate1, region)

tab.pval= summarize(df, pval=mean(pow, na.rm=TRUE))


# n=1000, GAM default knots
load("code/simulation_iid_output/sim_gam_n1000.RData")
pow.gam= rowMeans(sapply(sim, function(z) rej= (z$ci$ci.low > 0) | (z$ci$ci.up < 0)))

df= data.frame(sim[[1]]$ci[,c('varname','z')], pow.gam) |>
  mutate(region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3.01)), covariate1= (varname==1)) |>
  group_by(covariate1, region)

tab.gam= summarize(df, pow.gam= mean(pow.gam))


# n=1000, GAM k=24 knots
load("code/simulation_iid_output/sim_gam_knots24_n1000.RData")
pow.gam= rowMeans(sapply(sim, function(z) rej= (z$ci$ci.low > 0) | (z$ci$ci.up < 0)))

df= data.frame(sim[[1]]$ci[,c('varname','z')], pow.gam) |>
  mutate(region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3.01)), covariate1= (varname==1)) |>
  group_by(covariate1, region)

tab.gam= summarize(df, pow.gam= mean(pow.gam))


# n=1000, VC-BART
n=1000
load("code/simulation_iid_output/sim_vcbart_n1000.RData")
names(sim.vcbart[[1]])

threshold= 0
pow.vcbart= rowMeans(sapply(sim.vcbart, function(z, threshold=threshold) rej= (z$beta.ci[,'low'] > 0) | (z$beta.ci[,'up']<0)))

df= data.frame(sim.vcbart[[1]]$beta.ci[,c('variable','z')], pow.vcbart) |>
  mutate(region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3.01)), covariate1= (variable==1)) |>
  group_by(covariate1, region)

tab.vcbart= summarize(df, pow.vcbart= mean(pow.vcbart))

tab= cbind(tab, tab.icarplus[,'cut0.icarplus'], tab.fusedlasso[,'fusedlasso'], tab.pval[,'pval'], tab.vcbart[,'pow.vcbart'], tab.gam[,'pow.gam'])
cbind(filter(tab, covariate1)[,-1], filter(tab, !covariate1)[,-1:-2])


###############################################################################
## MSE
###############################################################################

mse= matrix(NA, nrow=2, ncol=6)
rownames(mse)= c('n=100', 'n=1000')
colnames(mse)= c('cut0', 'uncut3', 'ICAR+','VC-BART','GAM','GAM (k=24)')

n=100
load(here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))
mse["n=100",1:2]= colMeans(do.call(rbind,lapply(sim, "[[", "mse"))) / n
load(here("code", "simulation_iid_output", paste0("sim_icarplus_n", n,".RData")))
mse["n=100",'ICAR+']= colMeans(do.call(rbind,lapply(sim, "[[", "mse")))[1] / n
load(here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))
mse["n=100","VC-BART"]= mean(sapply(sim.vcbart, "[[", "mse")) / n
load("code/simulation_iid_output/sim_gam_n100.RData")
mse["n=100","GAM"]= mean(sapply(sim, "[[", "mse")) / n
load("code/simulation_iid_output/sim_gam_knots24_n100.RData")
mse["n=100","GAM (k=24)"]= mean(sapply(sim, "[[", "mse")) / n

n=1000
load(here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))
mse["n=1000",1:2]= colMeans(do.call(rbind,lapply(sim, "[[", "mse"))) / n
load(here("code", "simulation_iid_output", paste0("sim_icarplus_n", n,".RData")))
mse["n=1000",'ICAR+']= colMeans(do.call(rbind,lapply(sim, "[[", "mse")))[1] / n
load(here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))
mse["n=1000","VC-BART"]= mean(sapply(sim.vcbart, "[[", "mse")) / n
load("code/simulation_iid_output/sim_gam_n1000.RData")
mse["n=1000","GAM"]= mean(sapply(sim, "[[", "mse")) / n
load("code/simulation_iid_output/sim_gam_knots24_n1000.RData")
mse["n=1000","GAM (k=24)"]= mean(sapply(sim, "[[", "mse")) / n


library(xtable)
xtable(sqrt(mse), digits=c(0,3,3,3,3,3))




###############################################################################
## PLOTS FOR n=100
###############################################################################

textsize= 30
load("code/simulation_iid_output/sim_n100.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))

#Average posterior probability
pp0m= rowMeans(pp0)
pp3m= rowMeans(pp3)

#Power for test based on pp > 0.95
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

#Power for VC-BART
load("code/simulation_iid_output/sim_vcbart_n100.RData")

threshold=0
pow.vcbart= rowMeans(sapply(sim.vcbart, function(z, threshold=threshold) rej= (z$beta.ci[,'low'] > 0) | (z$beta.ci[,'up']<0)))
df.vcbart= data.frame(sim.vcbart[[1]]$beta.ci[,c('variable','z')], pow.vcbart)

#Power for GAM
load("code/simulation_iid_output/sim_gam_n100.RData")
pow.gam= rowMeans(sapply(sim, function(z) rej= (z$ci$ci.low > 0) | (z$ci$ci.up < 0)))
df.gam= data.frame(sim[[1]]$ci[,c('varname','z')], pow.gam)




## PLOT POSTERIOR PROBABILITIES ##
##################################

df1= tibble(id, pp0m, pp3m) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pp0m','pp3m'), names_to= "Basis", values_to="pp") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Basis= recode(df1$Basis, pp0m="Degree 0", pp3m="Cubic")

df2= tibble(id, pp0m, pp3m) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pp0m), `Cubic`=mean(pp3m)) |>
    gather(key=Basis, value=pp, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Basis,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(color=Basis, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Posterior probability of covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.75,'cm'))

ggsave("drafts/figs/simiid_pp_n100.pdf")




## PLOT POWER FUNCTION ##
#########################

df1= tibble(id, pow0, pow3) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pow0','pow3'), names_to= "Method", values_to="pow") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Method= recode(df1$Method, pow0="Degree 0", pow3="Cubic")

df1.vcbart= filter(df.vcbart, variable==1) |>
  transform(Method='VC-BART', pow=pow.vcbart, Covariate=paste('x',variable,sep='')) |>
  select(z, Method, pow, Covariate)

df1.gam= filter(df.gam, varname==1) |>
  transform(Method='GAM', pow=pow.gam, Covariate=paste('x',varname,sep='')) |>
  select(z, Method, pow, Covariate)

df1= rbind(df1, df1.vcbart, df1.gam)

df2= tibble(id, pow0, pow3) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pow0), `Cubic`=mean(pow3)) |>
    gather(key=Method, value=pow, 2:3) |>
    transform(Covariate="x2-x10")

df2.vcbart= filter(df.vcbart, variable!=1) |>
  group_by(z) |>
  summarize(pow= mean(pow.vcbart)) |>
  transform(Method='VC-BART', Covariate="x2-x10") |>
  select(z, Method, pow, Covariate)

df2.gam= filter(df.gam, varname!=1) |>
  group_by(z) |>
  summarize(pow= mean(pow.gam)) |>
  transform(Method='GAM', Covariate="x2-x10") |>
  select(z, Method, pow, Covariate)

df2= rbind(df2, df2.vcbart, df2.gam)

df= rbind(df1, df2) |>
    transform(group= paste(Method,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pow)) +
    geom_line(aes(color=Method, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Power for covariate inclusion') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))

ggsave("drafts/figs/simiid_pow_n100.pdf")



###############################################################################
## PLOTS FOR n=1000
###############################################################################

load("code/simulation_iid_output/sim_n1000.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))

#Average posterior probability
pp0m= rowMeans(pp0)
pp3m= rowMeans(pp3)

#Power for test based on pp > 0.95
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

#Power for VC-BART
load("code/simulation_iid_output/sim_vcbart_n1000.RData")

threshold=0
pow.vcbart= rowMeans(sapply(sim.vcbart, function(z, threshold=threshold) rej= (z$beta.ci[,'low'] > 0) | (z$beta.ci[,'up']<0)))
df.vcbart= data.frame(sim.vcbart[[1]]$beta.ci[,c('variable','z')], pow.vcbart)

#Power for GAM
load("code/simulation_iid_output/sim_gam_n1000.RData")
pow.gam= rowMeans(sapply(sim, function(z) rej= (z$ci$ci.low > 0) | (z$ci$ci.up < 0)))
df.gam= data.frame(sim[[1]]$ci[,c('varname','z')], pow.gam)



## PLOT POSTERIOR PROBABILITIES ##
##################################

df1= tibble(id, pp0m, pp3m) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pp0m','pp3m'), names_to= "Basis", values_to="pp") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Basis= recode(df1$Basis, pp0m="Degree 0", pp3m="Cubic")

df2= tibble(id, pp0m, pp3m) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pp0m), `Cubic`=mean(pp3m)) |>
    gather(key=Basis, value=pp, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Basis,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(color=Basis, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Posterior probability of covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.75,'cm'))

ggsave("drafts/figs/simiid_pp_n1000.pdf")




## PLOT POWER FUNCTION ##
#########################

df1= tibble(id, pow0, pow3) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pow0','pow3'), names_to= "Method", values_to="pow") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Method= recode(df1$Method, pow0="Degree 0", pow3="Cubic")

df1.vcbart= filter(df.vcbart, variable==1) |>
  transform(Method='VC-BART', pow=pow.vcbart, Covariate=paste('x',variable,sep='')) |>
  select(z, Method, pow, Covariate)

df1.gam= filter(df.gam, varname==1) |>
  transform(Method='GAM', pow=pow.gam, Covariate=paste('x',varname,sep='')) |>
  select(z, Method, pow, Covariate)

df1= rbind(df1, df1.vcbart, df1.gam)

df2= tibble(id, pow0, pow3) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pow0), `Cubic`=mean(pow3)) |>
    gather(key=Method, value=pow, 2:3) |>
    transform(Covariate="x2-x10")

df2.vcbart= filter(df.vcbart, variable!=1) |>
  group_by(z) |>
  summarize(pow= mean(pow.vcbart)) |>
  transform(Method='VC-BART', Covariate="x2-x10") |>
  select(z, Method, pow, Covariate)

df2.gam= filter(df.gam, varname!=1) |>
  group_by(z) |>
  summarize(pow= mean(pow.gam)) |>
  transform(Method='GAM', Covariate="x2-x10") |>
  select(z, Method, pow, Covariate)

df2= rbind(df2, df2.vcbart, df2.gam)

df= rbind(df1, df2) |>
    transform(group= paste(Method,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pow)) +
    geom_line(aes(color=Method, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Power for covariate inclusion') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))

ggsave("drafts/figs/simiid_pow_n1000.pdf")





###############################################################################
## SIMULATION WITH UNALIGNED KNOTS
###############################################################################

#True group differences now start at z=0.15
truemean= function(x,z) {
    ans= double(nrow(x))
    group1= (x[,1]==1)
    ans[group1]= ifelse(z[group1]  <=0.15,  cos(z[group1]), cos(0.15))
    ans[!group1]= ifelse(z[!group1]<=0.15, cos(z[!group1]), cos(0.15)/(z[!group1] + 0.85)^2)
    return(ans)
}

## Plot mean for both groups ##
###############################

set.seed(1)
n= 2000
x= matrix(rep(0:1,c(n/2,n/2)), ncol=1)
z= rep(seq(-3,3,length=n/2), 2)
m= truemean(x,z)
y= truemean(x,z) + rnorm(n, 0, .25)

nlocalknots= c(7,9,11)
fit0= localnulltest(y, x=x, z=z, cutdegree=0, nlocalknots=nlocalknots, mc.cores = 1)


k1= fit0$regionbounds[[1]][[1]] #knots of the local test basis
k2= fit0$regionbounds[[2]][[1]] #knots of the local test basis
k3= fit0$regionbounds[[3]][[1]] #knots of the local test basis


pdf('~/github/localnulltest/drafts/figs/simiid_unaligned.pdf')
par(mar=c(4.2,4.2,.1,.1))
sel= x[,1]==1
plot(z[sel], m[sel], ylim=range(m), xlab='z', ylab='y', type='l', col='gray', lwd=5, cex.axis=1.5, cex.lab=1.5)
lines(z[!sel], m[!sel], col='gray', lwd=5)
abline(v= 0.15)
abline(h= min(m) + c(-.04,0,.04))
points(k1, rep(min(m)+.04, length(k1)), pch=16, cex=2)
points(k2, rep(min(m), length(k2)), pch=17, cex=2)
points(k3, rep(min(m)-.04, length(k3)), pch=18, cex=2)
legend('topleft',c('Resolution 1','Resolution 2','Resolution 3'), pch=c(16,17,18), cex=1.8)
dev.off()


## Run simulations ##
#####################

#Simulations with n=100, Normal shrinkage prior on coefficients
n=100; p=10; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_unaligned_n", n,".RData")))

#Simulations with n=1000, Normal shrinkage prior on coefficients
n=1000; p=10; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_unaligned_n", n,".RData")))


## Plot posterior probabilities and power function (n=100) ##
#############################################################

textsize= 25
load("code/simulation_iid_output/sim_unaligned_n100.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))

#Average posterior probability
pp0m= rowMeans(pp0)
pp3m= rowMeans(pp3)

#Power for test based on pp > 0.95
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

## PLOT POSTERIOR PROBABILITIES ##
df1= tibble(id, pp0m, pp3m) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pp0m','pp3m'), names_to= "Basis", values_to="pp") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Basis= recode(df1$Basis, pp0m="Degree 0", pp3m="Cubic")

df2= tibble(id, pp0m, pp3m) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pp0m), `Cubic`=mean(pp3m)) |>
    gather(key=Basis, value=pp, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Basis,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(color=Basis, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Posterior probability of covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.75,'cm'))

ggsave("drafts/figs/simiid_unaligned_pp_n100.pdf")


## PLOT POWER FUNCTION 
df1= tibble(id, pow0, pow3) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pow0','pow3'), names_to= "Method", values_to="pow") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Method= recode(df1$Method, pow0="Degree 0", pow3="Cubic")


df2= tibble(id, pow0, pow3) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pow0), `Cubic`=mean(pow3)) |>
    gather(key=Method, value=pow, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Method,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pow)) +
    geom_line(aes(color=Method, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Power for covariate inclusion') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))

ggsave("drafts/figs/simiid_unaligned_pow_n100.pdf")



## Plot posterior probabilities and power function (n=1000) ##
#############################################################

textsize= 30
load("code/simulation_iid_output/sim_unaligned_n1000.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))

#Average posterior probability
pp0m= rowMeans(pp0)
pp3m= rowMeans(pp3)

#Power for test based on pp > 0.95
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

## PLOT POSTERIOR PROBABILITIES ##
df1= tibble(id, pp0m, pp3m) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pp0m','pp3m'), names_to= "Basis", values_to="pp") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Basis= recode(df1$Basis, pp0m="Degree 0", pp3m="Cubic")

df2= tibble(id, pp0m, pp3m) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pp0m), `Cubic`=mean(pp3m)) |>
    gather(key=Basis, value=pp, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Basis,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(color=Basis, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Posterior probability of covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.75,'cm'))

ggsave("drafts/figs/simiid_unaligned_pp_n1000.pdf")


## PLOT POWER FUNCTION 
df1= tibble(id, pow0, pow3) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pow0','pow3'), names_to= "Method", values_to="pow") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Method= recode(df1$Method, pow0="Degree 0", pow3="Cubic")


df2= tibble(id, pow0, pow3) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pow0), `Cubic`=mean(pow3)) |>
    gather(key=Method, value=pow, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Method,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pow)) +
    geom_line(aes(color=Method, lty=Covariate), lwd=3) +
    ylim(0,1) +
    labs(y='Power for covariate inclusion') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))

ggsave("drafts/figs/simiid_unaligned_pow_n1000.pdf")





###############################################################################
## SIMULATION WITH UNALIGNED KNOTS, NON EQUI-SPACED KNOTS
###############################################################################

#True group differences now start at z=0.15
truemean= function(x,z) {
    ans= double(nrow(x))
    group1= (x[,1]==1)
    ans[group1]= ifelse(z[group1]  <=0.15,  cos(z[group1]), cos(0.15))
    ans[!group1]= ifelse(z[!group1]<=0.15, cos(z[!group1]), cos(0.15)/(z[!group1] + 0.85)^2)
    return(ans)
}

#Specify knots
localknots= vector("list",4)
localknots[[1]]= unique(c(seq(-3,0,length=4),seq(0,3,length=4)))  #7 knots, equi-spaced
localknots[[2]]= unique(c(seq(-3,0,length=4),seq(0,3,length=8)))  #11 knots, twice as many knots for z>0
localknots[[3]]= unique(c(seq(-3,0,length=6),seq(0,3,length=6)))  #11 knots, equi-spaced
localknots[[4]]= unique(c(seq(-3,0,length=6),seq(0,3,length=12))) #17 knots, twice as many knots for z>0


## Plot mean for both groups ##
###############################

set.seed(1)
n= 2000
x= matrix(rep(0:1,c(n/2,n/2)), ncol=1)
z= rep(seq(-3,3,length=n/2), 2)
m= truemean(x,z)
y= truemean(x,z) + rnorm(n, 0, .25)

fit0= localnulltest(y, x=x, z=z, cutdegree=0, localknots=localknots, mc.cores = 3)

k1= fit0$regionbounds[[1]][[1]] #knots of the local test basis
k2= fit0$regionbounds[[2]][[1]] #knots of the local test basis
k3= fit0$regionbounds[[3]][[1]] #knots of the local test basis
k4= fit0$regionbounds[[4]][[1]] #knots of the local test basis


pdf('~/github/localnulltest/drafts/figs/simiid_unaligned_nonequispaced.pdf')
par(mar=c(4.2,4.2,.1,.1))
sel= x[,1]==1
plot(z[sel], m[sel], ylim=range(m), xlab='z', ylab='y', type='l', col='gray', lwd=5, cex.axis=1.5, cex.lab=1.5)
lines(z[!sel], m[!sel], col='gray', lwd=5)
abline(v= 0.15)
abline(h= min(m) + c(-.04,0,.04,.08))
points(k1, rep(min(m)+.08, length(k1)), pch=16, cex=2, col='black')
points(k2, rep(min(m)+.04, length(k2)), pch=17, cex=2, col='black')
points(k3, rep(min(m)+.00, length(k3)), pch=16, cex=2, col='gray')
points(k4, rep(min(m)-.04, length(k4)), pch=17, cex=2, col='gray')
legend('topleft',c('Resolution 1','Resolution 2','Resolution 3','Resolution 4'), pch=c(16,17,16,17), col=c('black','black','gray','gray'), cex=1.7)
dev.off()



## Run simulations ##
#####################

#Simulations with n=100, Normal shrinkage prior on coefficients
n=100; p=1; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef, localknots=localknots, uncut=FALSE); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_unaligned_nonequispaced_n", n,".RData")))

#Simulations with n=1000, Normal shrinkage prior on coefficients
n=1000; p=1; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef, localknots=localknots, uncut=FALSE); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_unaligned_nonequispaced_n", n,".RData")))

#Simulations with n=5000, Normal shrinkage prior on coefficients
n=5000; p=1; sd.error=0.25; nsims= 100; priorCoef= normalidprior()
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores=3, priorCoef=priorCoef, localknots=localknots, uncut=FALSE); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_unaligned_nonequispaced_n", n,".RData")))


## Posterior probabilities of each resolution ##
################################################

pp_localknots= matrix(NA, nrow=4, ncol=3)
colnames(pp_localknots)= c('n=100','n=1000','n=5000')
rownames(pp_localknots)= c('7 knots, equi-spaced','11 knots (3 regions for z<0, 7 for z>0)','11 knots, equi-spaced','17 knots (5 regions for z<0, 11 for z>0')

load("code/simulation_iid_output/sim_unaligned_nonequispaced_n100.RData")
pp_localknots[,'n=100']= rowMeans(sapply(sim, '[[', 'pp_localknots'))

load("code/simulation_iid_output/sim_unaligned_nonequispaced_n1000.RData")
pp_localknots[,'n=1000']= rowMeans(sapply(sim, '[[', 'pp_localknots'))

load("code/simulation_iid_output/sim_unaligned_nonequispaced_n5000.RData")
pp_localknots[,'n=5000']= rowMeans(sapply(sim, '[[', 'pp_localknots'))


xtable(pp_localknots)


## Plot posterior probabilities for local tests (n=100) ##
##########################################################

textsize= 25
load("code/simulation_iid_output/sim_unaligned_nonequispaced_n100.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp0m= rowMeans(pp0)
df100= tibble(id, pp0m) |>
         filter(covariate == 1) |>
         select(-covariate) |>
         transform(Covariate="x1",n=100)

load("code/simulation_iid_output/sim_unaligned_nonequispaced_n1000.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp0m= rowMeans(pp0)
df1000= tibble(id, pp0m) |>
          filter(covariate == 1) |>
          select(-covariate) |>
          transform(Covariate="x1",n=1000)

load("code/simulation_iid_output/sim_unaligned_nonequispaced_n5000.RData")
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp0m= rowMeans(pp0)
df5000= tibble(id, pp0m) |>
          filter(covariate == 1) |>
          select(-covariate) |>
          transform(Covariate="x1",n=5000)


df= rbind(df100,df1000,df5000) |>
      group_by(z,n) |>
      summarise(pp= mean(pp0m)) |>
      transform(n= factor(n)) |>
      filter(pp > 0)

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(color=n, lty=n), lwd=3) +
    ylim(0,1) +
    labs(y='') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.75,'cm'))

ggsave("drafts/figs/simiid_unaligned_nonequispaced_pp.pdf")


