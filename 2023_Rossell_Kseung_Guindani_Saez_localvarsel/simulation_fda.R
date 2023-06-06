###############################################################################
##
## THIS SCRIPT REPRODUCES THE RESULTS FROM SECTION 5.2 "FUNCTIONAL DATA SIMULATION"
##
###############################################################################

library(mombf)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(rgen)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(here)
library(ggrepel)
library(ggplot2)
setwd("~/github/localnulltest/code")
source("routines.R")
library(reshape2)
library(gridExtra)
source("lfmm/R/mcmc_fcts.R")
sourceCpp("lfmm/src/mcmc_fcts.cpp")




###############################################################################
## DEFINE FUNCTIONS
###############################################################################

truemean= function(x,z) {
    if (x[1]==1) {
        ans= ifelse(z <=0, cos(z), 1)
    } else {
        ans= ifelse(z<=0, cos(z), 1/(z+1)^2)
    }
    return(ans)
}

#Simulate functional data for nindiv individuals, each function being a grid over npoints
# For each function generate y = truemean(x,z) + e, where
# x[,1] is binary with 0.5 observations equal to 1 and 0
# x[,j] ~ N(x[,1], 1) for j=2,...,p
# e ~ N(0, Sigma), and Sigma is an AR1 covariance, i.e. Sigma[i,i]=1 and Sigma[i,j]= rho^(i-j)
# Input
# - seed: random number seed
# - ninvid: number of individuals (functions) to generate
# - npoints: grid size
# - p: number of covariates
# - rho: AR1 autocorrelation parameter
# Ouput
# - y: outcome
# - x: covariates
# - z: coordinates on the grid
# - function_id: function identifier 
simdata_fda= function(seed, nindiv, npoints, p, rho) {
    set.seed(seed)
    #Generate covariate values for each individual
    xind= matrix(NA, nrow=nindiv, ncol=p)
    xind[,1]= rep(0:1,c(nindiv/2,nindiv/2))
    for (j in 2:p) xind[,j]= xind[,1] + rnorm(nindiv)
    #Error covariance matrix
    Sigma= diag(npoints)
    for (i in 1:nrow(Sigma)) for (j in 1:i) Sigma[i,j]= Sigma[j,i]=rho^(i-j)
    #Time grid
    zseq= seq(-3,3,length=npoints)
    #Generate outcome and covariates
    y= x= z= m= vector("list",nindiv)
    for (i in 1:nindiv) {
        m[[i]]= truemean(xind[i,],zseq)
        e= as.vector(rmvnorm(1, sigma=Sigma)) 
        y[[i]]= m[[i]] + 0.5 * e
        x[[i]]= matrix(rep(xind[i,],nrow(Sigma)), ncol=p, byrow=TRUE)
        z[[i]]= zseq
    }
    y= do.call(c,y)
    x= do.call(rbind,x)
    z= do.call(c,z)
    m= do.call(c,m)
    function_id= rep(1:nindiv, each=npoints)
    ans= list(y=y, x=x, z=z, function_id=function_id, m=m)
    return(ans)
}




#Simulate functional data and run our method with 0 degree cut basis
sim_localnulltest_fda= function(seed, nindiv, npoints, p, rho, single=FALSE, priorCoef=momprior(), priorGroup=groupmomprior()) {
    set.seed(seed)
    sim= simdata_fda(seed=seed, nindiv=nindiv, npoints=npoints, p=p, rho=rho)
    if (single) x= sim$x[,1,drop=FALSE] else x= sim$x
    fit= localnulltest_fda(sim$y, x=x, z=sim$z, function_id=sim$function_id, Sigma='AR/MA', priorCoef=priorCoef, priorGroup=priorGroup)
    b= coef(fit)
    #MSE
    mhat= predict(fit)
    mse= mean((mhat - sim$m)^2)
    #Return output
    margpp= b[,c('covariate','z1','estimate','margpp')]
    colnames(margpp)= c('covariate','z','estimate','cut0')
    ans= list(margpp=margpp, mse=mse)
    return(ans)
}


###############################################################################
## LFMM FUNCTIONS
###############################################################################

#Simulate functional data and run the LFMM method (single covariate case) from Paulon, Sarkar & Mueller 2023
sim_LFMM_single= function(seed, nindiv, npoints, p, rho) {
    require(tidyverse)
    set.seed(seed)
    sim= simdata_fda(seed=seed, nindiv=nindiv, npoints=npoints, p=p, rho=rho)
    #Run MCMC
    X= sim$x[,1,drop=FALSE] + 1 #use covariate 1 only (loc_clust_reff requires index to start at 1) 
    zseq= unique(sim$z)
    time= as.numeric(as.factor(sim$z))
    xgrid= min(time):max(time)
    hypers= list(sigma_s= 0.1, a_sigma= 1, b_sigma= 1, a_phi=10, b_phi=1, a_alpha=1, b_alpha=1)
    fit= loc_clust_reff(y=sim$y, time=time, ind=sim$function_id, X=X, xgrid=xgrid, hypers=hypers, Niter=5000, burnin=500)
    #Posterior probabilities of local effect
    pp= double(length(zseq))
    for (t in 1:length(pp)) {
      z= fit$Z[,t,]
      ngroups= apply(z, 2, function(zz) length(unique(zz)))
      pp[t]= mean(ngroups > 1)
    }
    pp= data.frame(z=zseq, pp=pp)
    #MSE
    m= unique(data.frame(X1=sim$x[,1], z=sim$z, m=sim$m)) |>
      arrange(X1, z)
    mhat= c(rowMeans(fit$post_mean[,1,]), rowMeans(fit$post_mean[,2,]))
    mse= mean((mhat - m$m)^2)
    #Return output
    ans= list(pp=pp, mse=mse)
    return(ans)
}


#Simulate functional data and run the LFMM method (multiple covariates case) from Paulon, Sarkar & Mueller 2023
sim_LFMM_multi= function(seed, nindiv, npoints, p, rho) {
    require(tidyverse)
    set.seed(seed)
    sim= simdata_fda(seed=seed, nindiv=nindiv, npoints=npoints, p=p, rho=rho)
    #Discretize X's
    df= as_tibble(data.frame(sim$x))
    breaks= c(-Inf, 0.5, Inf)  #cut at the mean (roughly 0.5)
    X= mutate(df, across(X2:X10, cut, breaks=breaks)) |>
      mutate(across(X1:X10, as.integer)) |>
      as.matrix()
    X[,1]= X[,1] + 1  #loc_clust_multi_reff requires levels to start at 1
    Xgrid= as.matrix(expand.grid(apply(X, 2, unique, simplify=FALSE)))
    #Hyper-parameters (default values from Paulon, Sarkar & Mueller 2023)
    hypers= list(sigma_s= 0.1, a_sigma= 1, b_sigma= 1, a_phi=10, b_phi=1, a_alpha=1, b_alpha=1)
    #Run MCMC
    time= as.numeric(as.factor(sim$z)); t_grid= min(time):max(time)
    fit= loc_clust_multi_reff(y=sim$y, time=time, ind=sim$function_id, X=X, t_grid=t_grid, Xgrid=Xgrid, hypers=hypers, Niter=5000, burnin=500)
    #MSE in estimating the mean
    mhat= colMeans(fit$post_mean) #predicted mean(y)
    mse= mean((sim$m - mhat)^2)
    #Posterior probability of local effect
    pp= pp_localtest_lfmm(fit, time=time, Xgrid=Xgrid) #check: is it post_mean or post_pred that should be used
    #Return output
    ans= list(pp=pp, mse=mse)
    return(ans)
}


#Compute posterior probabilities for local tests from LFMM model output
#Input
# - fit: output from loc_clust_multi_reff
# - time: time argument given to loc_clust_multi_reff
# - Xgrid: Xgrid argument given to loc_clust_multi_reff
#Output: for each covariate and time, the posterior probability that the covariate has an effect (proportion of MCMC draws where it affected the prediction)
pp_localtest_lfmm = function(fit, time, Xgrid) {
  pp= matrix(NA, nrow=max(time), ncol=ncol(Xgrid))
  colnames(pp)= paste0("X",1:ncol(pp))
  for (t in 1:nrow(pp)) {
    postpred= fit$post_pred[,,t]
    for (j in 1:ncol(Xgrid)) {
      #For each row in postpred, obtain the mean for each value of Xj
      df= data.frame(cbind(Xgrid[,j], t(postpred))) |> 
        group_by(X1) |> 
        summarize_all(mean) |> 
        select(-1)
      #Proportion of iterations in which the group means differ
      pp[t,j]= summarize_all(df, sd) |> 
        as.matrix() |>
        mean()
    }
  }
  #Format output
  pp= data.frame(time= zseq, pp) |>
    pivot_longer(cols=-1, names_to="covariate", values_to="pp") |>
    mutate(covariate= sub("X", "", covariate))
  return(pp)
}

###############################################################################
## PINI-VANTINI FUNCTIONS
###############################################################################

library(fdatest)

sim_piva= function(seed, nindiv, npoints, p, rho, single=FALSE, pval_thresh=0.05, priorCoef=momprior(), priorGroup=groupmomprior()) {
  set.seed(seed)
  sim= simdata_fda(seed=seed, nindiv=nindiv, npoints=npoints, p=p, rho=rho)
  zseq= unique(sim$z)
  if (single) x= sim$x[,1,drop=FALSE] else x= sim$x
  reshaped_y=matrix(sim$y, nrow=nindiv, ncol=npoints, byrow=T)
  x1=sim$x[seq(1, nrow(sim$x), npoints), ][,1]
  ITP_fit <- ITP2bspline(data1 = reshaped_y[x1==0,], data2 = reshaped_y[x1==1,])
  mse= mean((ITP_fit$data.eval- sim$m)^2)
  margpp=ITP_fit$corrected.pval # corrected p-values
  gammahat=ifelse(ITP_fit$corrected.pval < pval_thresh, 1, 0)
  #Return output
  ans= list(margpp=margpp, mse=mse, gammahat=gammahat, z=zseq)
  return(ans)
}

###############################################################################
## RUN SIMULATIONS FOR OUR METHOD
###############################################################################

# Using all covariates

nindiv= 50; npoints= 100; p=10; rho=0.99; nsims=100
pb= txtProgressBar()
sim= vector("list", nsims)
for (i in 1:nsims) {
    sim[[i]]= sim_localnulltest_fda(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho, priorCoef=momprior(), priorGroup=groupmomprior())
    setTxtProgressBar(pb, i/nsims)
}
#save(sim, file=paste0("~/github/localnulltest/code/simulation_fda_output/sim_pmom_nindiv_",nindiv,".RData"))


nindiv= 100; npoints= 100; p=10; rho=0.99; nsims=100
pb= txtProgressBar()
sim= vector("list", nsims)
for (i in 1:nsims) {
    sim[[i]]= sim_localnulltest_fda(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho, priorCoef=momprior(), priorGroup=groupmomprior())
    setTxtProgressBar(pb, i/nsims)
}
#save(sim, file=paste0("~/github/localnulltest/code/simulation_fda_output/sim_pmom_nindiv_",nindiv,".RData"))


# Only using X1

nindiv= 50; npoints= 100; p=10; rho=0.99; nsims=100
pb= txtProgressBar()
sim= vector("list", nsims)
for (i in 1:nsims) {
    sim[[i]]= sim_localnulltest_fda(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho, single=TRUE)
    setTxtProgressBar(pb, i/nsims)
}
save(sim, file=paste0("~/github/localnulltest/code/simulation_fda_output/sim_pmom_single_nindiv_",nindiv,".RData"))


nindiv= 100; npoints= 100; p=10; rho=0.99; nsims=100
pb= txtProgressBar()
sim= vector("list", nsims)
for (i in 1:nsims) {
    sim[[i]]= sim_localnulltest_fda(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho, single=TRUE)
    setTxtProgressBar(pb, i/nsims)
}
save(sim, file=paste0("~/github/localnulltest/code/simulation_fda_output/sim_pmom_single_nindiv_",nindiv,".RData"))



###############################################################################
## RUN SIMULATIONS FOR LFMM (ONLY USING TRULY ACTIVE COVARIATE X1)
###############################################################################

nindiv= 50; npoints= 100; p=10; rho=0.99; nsims=100
sim= vector("list", nsims)
for (i in 1:length(sim)) sim[[i]]= try(sim_LFMM_single(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho))
sel= (sapply(sim, 'class') == 'list'); sim= sim[sel] #remove the few simulations in which code crashed
save(sim, file=here("code", "simulation_fda_output", paste0("sim_lfmm_single_nindiv", nindiv,".RData")))

###############################################################################
## RUN SIMULATIONS FOR PINI-VANTINI (ONLY USING TRULY ACTIVE COVARIATE X1)
##############################################################################

nindiv= 50; npoints= 100; p=10; rho=0.99; nsims=100
pb= txtProgressBar()
sim= vector("list", nsims)
for (i in 1:nsims) {
  sim[[i]]= sim_piva(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho, priorCoef=normalidprior(), priorGroup=normalidprior())
  setTxtProgressBar(pb, i/nsims)
}
save(sim, file=paste0("~/github/localnulltest/code/simulation_fda_output/sim_vantini_nindiv_",nindiv,".RData"))


nindiv= 100; npoints= 100; p=10; rho=0.99; nsims=100
pb= txtProgressBar()
sim= vector("list", nsims)
for (i in 1:nsims) {
  sim[[i]]= sim_piva(seed=i, nindiv=nindiv, npoints=npoints, p=p, rho=rho, priorCoef=normalidprior(), priorGroup=normalidprior())
  setTxtProgressBar(pb, i/nsims)
}
save(sim, file=paste0("~/github/localnulltest/code/simulation_fda_output/sim_vantini_nindiv_",nindiv,".RData"))



###############################################################################
## PLOTS
###############################################################################

textsize= 22
#textsize= 35

## POSTERIOR PROBABILITIES ##
#############################

# Using all 10 covariates
nindiv= 50
load(paste0("simulation_fda_output/sim_pmom_nindiv_",nindiv,".RData"))
pp= sapply(sim, function(z) z$margpp$cut0)
df= data.frame(sim[[1]]$margpp[,1:2], pp=rowMeans(pp))
txt= filter(df, z==max(z)) |> mutate(covariate= paste0("X", covariate))
ggplot(df, aes(z, pp)) +
    geom_line(aes(group=covariate)) +
    geom_text_repel(aes(x=z, y=pp, label=covariate), data=txt, size=textsize/4, max.overlaps=20) +
    labs(y='Posterior probability of local covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(5,'cm'))
ggsave(paste0("../drafts/figs/simfda_pp_nindiv",nindiv,".pdf"))


## Using only X1 ##
###################

nindiv= 50
load(paste0("simulation_fda_output/sim_pmom_single_nindiv_",nindiv,".RData"))
pp= sapply(sim, function(z) z$margpp$cut0)
df= data.frame(sim[[1]]$margpp[,1:2], pp=rowMeans(pp))

#LFMM method
load(paste0("simulation_fda_output/sim_lfmm_single_nindiv", nindiv,".RData"))
z.lfmm= sim[[1]]$pp$z
pp.lfmm= rowMeans(sapply(sim, function(zz) zz$pp$pp))
df.lfmm= data.frame(z=z.lfmm, pp=pp.lfmm)

txt= rbind(df[which.max(df$z),c('z','pp')], df.lfmm[which.max(df$z),c('z','pp')]) |>
  mutate(pp= pp+.02)
txt$method= c('0-degree cut basis','LFMM')

ggplot() +
    geom_line(aes(z, pp), data=df) +
    geom_line(aes(z, pp), data=df.lfmm, color='gray') +
    geom_text_repel(aes(x=z, y=pp, label=method), data=txt, size=textsize/4) +
    labs(y='Posterior probability of local covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(5,'cm'))
ggsave(paste0("../drafts/figs/simfda_pp_single_nindiv",nindiv,".pdf"))





###############################################################################
## SUMMARIZE TYPE I ERROR AND POWER
###############################################################################


## OUR METHOD, USING ONLY COVARIATE 1 ##
########################################

nindiv= 50
load(paste0("simulation_fda_output/sim_pmom_single_nindiv_",nindiv,".RData"))
pp= sapply(sim, function(z) z$margpp$cut0)
pow= data.frame(sim[[1]]$margpp[,1:2], pow=rowMeans(pp > 0.95))
df= mutate(pow, region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3))) |>
  group_by(region)
tab= summarize(df, cut0=mean(pow))


## LFMM, USING ONLY COVARITE 1 ##
#################################

load(paste0("simulation_fda_output/sim_lfmm_single_nindiv", nindiv,".RData"))
z.lfmm= sim[[1]]$pp$z
pow.lfmm= rowMeans(sapply(sim, function(zz) zz$pp$pp > 0.95))
pow.lfmm= data.frame(z=z.lfmm, pow=pow.lfmm)

df.lfmm= mutate(pow.lfmm, region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3))) |>
  group_by(region)
tab= summarize(df.lfmm, lfmm=mean(pow))


## PINI & VANTINI, USING ONLY COVARIATE 1 ##
############################################

nindiv= 50
load(paste0("simulation_fda_output/sim_vantini_nindiv_",nindiv,".RData"))
z.piva= zseq; #sim[[1]]$pp$z
pow.piva= rowMeans(sapply(sim, function(zz) zz$margpp < 0.05))
pow.piva= data.frame(z=z.piva, pow=pow.piva)
df.piva= mutate(pow.piva, region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3))) |>
  group_by(region)
tab= summarize(df.piva, piva=mean(pow))

nindiv= 100
load(paste0("simulation_fda_output/sim_vantini_nindiv_",nindiv,".RData"))
z.piva= zseq; #sim[[1]]$pp$z
pow.piva= rowMeans(sapply(sim, function(zz) zz$margpp < 0.05))
pow.piva= data.frame(z=z.piva, pow=pow.piva)
df.piva= mutate(pow.piva, region= cut(z, breaks=c(-3.01,-2,-1,0,1,2,3))) |>
  group_by(region)
tab= summarize(df.piva, piva=mean(pow))


## OUR METHOD, USING ALL COVARIATES ##
######################################

nindiv= 50
load(paste0("simulation_fda_output/sim_pmom_nindiv_",nindiv,".RData"))
pp= sapply(sim, function(z) z$margpp$cut0)
pow= data.frame(sim[[1]]$margpp[,1:2], pow=rowMeans(pp > 0.95))
df= mutate(pow, region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)
tab= summarize(df, cut0=mean(pow))


nindiv= 100
load(paste0("simulation_fda_output/sim_pmom_nindiv_",nindiv,".RData"))
pp= sapply(sim, function(z) z$margpp$cut0)
pow= data.frame(sim[[1]]$margpp[,1:2], pow=rowMeans(pp > 0.95))
df= mutate(pow, region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)
tab= summarize(df, cut0=mean(pow))



