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
source(here('code/routines.R'))
setwd(here('drafts/figs'))



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
    for (j in 2:p) x[,j]= x[,1] + rnorm(n)
    z= rep(seq(-3,3,length=n/2), 2)
    m= truemean(x,z)
    y= truemean(x,z) + rnorm(n, 0, sd=sd.error)
    ans= list(y=y, x=x, z=z, m=m)
    return(ans)
}


#Simulate data as in simdata and compute posterior probabilities for local null tests computed under several basis functions
# Input: same as for simdata
# Output: matrix with z coordinates and posterior probabilities under a cut 0-degree basis (cut0) and standard cubic B-spline (uncut3)
# Note: results for a cut cubic basis were very similar than for the 0-degree basis, hence not returned
sim_localnulltest= function(seed, n, p, sd.error=0.25, mc.cores) {
    sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)
    nlocalknots= c(7,9,11)
    fit0= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=0, nlocalknots=nlocalknots, mc.cores = mc.cores) #0-degree, multi-resolution
    b0= coef(fit0)
    fit.uncut= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=3, usecutbasis=FALSE, nlocalknots=nlocalknots, mc.cores = mc.cores) #uncut cubic, multi-resolution
    b.uncut= coef(fit.uncut)
    #MSE in estimating the mean
    mse.cut0= sum((predict(fit0) - sim$m)^2)  
    mse.uncut3= sum((predict(fit.uncut) - sim$m)^2)
    #Return output
    margpp= cbind(b0[,c('covariate','z1')], b0[,'margpp'], b.uncut[,'margpp'])
    colnames(margpp)= c('covariate','z','cut0','uncut3')
    ans= list(margpp=margpp, mse=c(cut0=mse.cut0, uncut3=mse.uncut3))
    return(ans)
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
    check= checkargs_localnulltest(y=y, x=x, z=z)
    y= check$y; x= check$x; z= check$z
    # Define local tests
    if (missing(localgrid)) localgrid= define_localgrid(localgrid=localgrid, localgridsize=localgridsize, z=z)
    # Define knots and corresponding knot-based testing regions
    kk= define_knots_localnulltest(z=z, localgrid=localgrid, nbaseknots=nbaseknots, basedegree=basedegree, nlocalknots=nlocalknots)
    knots= kk$knots; regionbounds= kk$regionbounds; region=kk$region; regioncoord= kk$regioncoord; testov= kk$testov; testxregion= kk$testxregion; testIntervals= kk$testIntervals
    #Create design matrix & run least-squares regression
    desnew= estimationPoints(x=x, regioncoord=regioncoord, regionbounds=regionbounds, testov=testov) #Points at which local effects will be estimated
    des= createDesignLocaltest(x=rbind(x,desnew$x), z=rbind(z,desnew$z), y=y, region=c(region,desnew$region), regionbounds=regionbounds, basedegree=basedegree, cutdegree=cutdegree, knots=knots, usecutbasis=TRUE, useSigma=FALSE)
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


#Modification of Sameer Deshpande's get_beta_support at https://github.com/skdeshpande91/VCBART
# - Corrects a bug that caused it to crash when R=1 (dimension of z)
# - Requires a single MCMC run, rather reporting the average across 2 MCMC runs
# Output: exceed_cutoff_probs are marginal posterior inclusion probabilities, support are the selected covariates
my_get_beta_support= function (chain1, burn, max_cutoff=1) {
    R= dim(chain1$var_counts_samples)[1]
    p= dim(chain1$var_counts_samples)[2]
    exceed_cutoff_probs <- array(dim = c(R, p, max_cutoff))
    support <- list()
    for (cutoff in 1:max_cutoff) {
        tmp_support <- list()
        exceed_cutoff_probs[, , cutoff]= apply(chain1$var_counts_samples[,,-(1:burn),drop=FALSE] >= cutoff, MARGIN = c(1, 2), FUN = mean)
        for (k in 1:p) tmp_support[[k]] <- which(exceed_cutoff_probs[, k, cutoff] >= 0.5)
        support[[cutoff]] <- tmp_support
    }
    return(list(exceed_cutoff_probs = exceed_cutoff_probs, support = support))
}

#Modification of Sameer Deshpande's get_beta_support at https://github.com/skdeshpande91/VCBART
# - Corrects a bug that caused it to crash when R=1 (dimension of z)
# - Requires a single MCMC run, rather reporting the average across 2 MCMC runs
# Input: chain1 is the output from VCBART, burn the number of burn-in iterations
# Output: posterior mean, sd, and 0.95 interval for each observation's mean
my_summarize_beta= function (chain1, burn) {
    N_train <- dim(chain1$beta_train_samples)[1]
    N_test <- dim(chain1$beta_test_samples)[1]
    p <- dim(chain1$beta_train_samples)[2]
    beta_summary_train <- array(dim = c(N_train, 4, p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
    beta_summary_test <- array(dim = c(N_test, 4, p), dimnames = list(c(), c("MEAN", "SD", "L95", "U95"), c()))
    for (k in 1:p) {
        beta_summary_train[, "MEAN", k] <- apply(chain1$beta_train_samples[, k, ], FUN=mean, MARGIN=1, na.rm=TRUE)
        beta_summary_train[, "SD", k] <- apply(chain1$beta_train_samples[, k, ], FUN = sd, MARGIN = 1, na.rm = TRUE)
        beta_summary_train[, "L95", k] <- apply(chain1$beta_train_samples[, k, ], FUN = quantile, MARGIN = 1, probs = 0.025)
        beta_summary_train[, "U95", k] <- apply(chain1$beta_train_samples[, k, ], FUN = quantile, MARGIN = 1, probs = 0.975)
        beta_summary_test[, "MEAN", k] <- apply(chain1$beta_test_samples[, k, ], FUN = mean, MARGIN = 1, na.rm = TRUE)
        beta_summary_test[, "SD", k] <- apply(chain1$beta_test_samples[, k, ], FUN = sd, MARGIN = 1, na.rm = TRUE)
        beta_summary_test[, "L95", k] <- apply(chain1$beta_test_samples[, k, ], FUN = quantile, MARGIN = 1, probs = 0.025)
        beta_summary_test[, "U95", k] <- apply(chain1$beta_test_samples[, k, ], FUN = quantile, MARGIN = 1, probs = 0.975)
    }
    return(list(train = beta_summary_train, test = beta_summary_test))
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
    bartfit= VCBART(Y_train=sim$y, X_train=xtrain, Z_train=matrix(sim$z,ncol=1), n_train=length(train), X_test=xtest, Z_test=ztest, n_test=nrow(xtest), error_structure="ind", split_probs_type="adaptive", nd=nd, cutpoints=cutpoints, verbose=FALSE)
    beta_summary= my_summarize_beta(bartfit, burn=500)
    betahat= beta_summary$train[,"MEAN",]
    betahat.test= beta_summary$test[,"MEAN",]
    #MSE for true expectation
    muhat= betahat[,1] + rowSums(betahat[,-1] * xtrain)
    mse= sum((muhat - sim$m[train])^2)  #MSE in estimating the mean
    #Posterior inclusion probability of z for each covariate
    beta_support= my_get_beta_support(bartfit, burn=500, max_cutoff=100)
    margpp= beta_support$exceed_cutoff_probs[1,,]  #select first row since z is univariate
    #margpp= beta_support$exceed_cutoff_probs[1,,cutoff]  #select first row since z is univariate
    #Return output
    ans= list(margpp=margpp, mse=mse, df=df)
    return(ans)
}



#Create test set consisting of a range of z values, where each covariate is set to 0 or 1 for each z-value and the remaining covariates are set to their mean
createTestset= function(x, z) {
    zseq= rep(seq(min(z), max(z), length=100), 2)
    m= colMeans(x)
    xseq= rep(0:1, each=100)
    ztest= matrix(rep(zseq, ncol(x)), ncol=1)
    xtest= lapply(1:ncol(x), function(i) matrix(NA, nrow=length(zseq), ncol=ncol(x)))
    covariate= rep(1:ncol(x), each=length(zseq))
    for (i in 1:ncol(x)) {
        xtest[[i]][,i]= xseq
        xtest[[i]][,-i]= matrix(rep(m[-i], length(zseq)), nrow=length(zseq), ncol=ncol(x)-1, byrow=TRUE)
    }
    xtest= do.call(rbind, xtest)
    ans= list(xtest=xtest, ztest=ztest, covariate=covariate)
    return(ans)
}



###############################################################################
## RUN SIMULATIONS
###############################################################################

#Simulations with n=100
n=100; p=10; sd.error=0.25; nsims= 100; mc.cores=3
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores = mc.cores); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))

#Simulations with n=1000
n=1000; p=10; sd.error=0.25; nsims= 100
sim= lapply(1:nsims, function(i) { ans= sim_localnulltest(seed=i, n=n, p=p, sd.error=sd.error, mc.cores = mc.cores); cat("."); return(ans) })
save(sim, file=here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))


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



n=100
load(here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))
names(sim.vcbart[[1]])
max_cutoff= 50
truly_active= 2


#Proportion of simulations in which a variable j for various cutoffs, ranging from 1 to max_cutoff
# Variable inclusion is based on the marginal posterior probability of being in >=cutoff trees. When said probability is > 0.95, the variable in considered included
# - j: index of the variable
# - max_cutoff: maximum cutoff to consider
prop_rejected_vcbart= function(j, max_cutoff=50) {
    ans= rowMeans(sapply(sim.vcbart, function(zz) zz$margpp[j+1, 1:max_cutoff] > 0.95))  #use j+1 since 1st variable is the intercept
    return(ans)
}

#Power
p= 10; max_cutoff= 50
prop_rejected_vcbart(2, max_cutoff=max_cutoff)

#Type I error
prop_rej= matrix(NA, nrow=p-1, ncol=max_cutoff)
for (j in 2:p) prop_rej[j-1,]= prop_rejected_vcbart(j, max_cutoff=max_cutoff)
colMeans(prop_rej)



pow.vcbart= sapply(sim.vcbart, function(zz) zz$margpp[2, 1:max_cutoff] > 0.95)  #number of rejected null hypothesis for various cutoffs (post prob > 0.95 of being in cutoff trees)


typeI.vcbart= sapply(sim.vcbart, function(zz) zz$margpp[-2, 1:max_cutoff] > 0.95)  #number of rejected null hypothesis for various cutoffs (post prob > 0.95 of being in cutoff trees)
rowMeans(pow.vcbart)



###############################################################################
## SUMMARIZE TYPE I ERROR AND POWER
###############################################################################


# n=100

load(here("code/simulation_iid_output/sim_n100.RData"))
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

df= data.frame(id, pow0, pow3) |>
  mutate(region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)

tab= summarize(df, cut0=mean(pow0), uncut3=mean(pow3))


load(here("code/simulation_iid_output/sim_ols_n100.RData"))
pval= do.call(cbind, lapply(sim, function(z) z$pvalue.adjusted))
pow= rowMeans(pval < .05, na.rm=TRUE)

df= data.frame(sim[[1]][,c('varname','zmin','zmax')], pow) |>
  mutate(region= factor(paste0('(',zmin,',',zmax,']')), covariate1= (varname=='x1')) |>
  select(-zmin, -zmax) |>
  group_by(covariate1, region)

tab.pval= summarize(df, pval=mean(pow, na.rm=TRUE))

tab= cbind(tab, tab.pval[,'pval'])
cbind(filter(tab, covariate1)[,-1], filter(tab, !covariate1)[,-1:-2])


# n=1000

load(here("code/simulation_iid_output/sim_n1000.RData"))
id= sim[[1]]$margpp[,c('covariate','z')]
pp0= do.call(cbind, lapply(sim, function(z) z$margpp[,'cut0']))
pp3= do.call(cbind, lapply(sim, function(z) z$margpp[,'uncut3']))
pow0= rowMeans(pp0 > 0.95)
pow3= rowMeans(pp3 > 0.95)

df= data.frame(id, pow0, pow3) |>
  mutate(region= cut(z, breaks=c(-3,-2,-1,0,1,2,3)), covariate1= (covariate==1)) |>
  group_by(covariate1, region)

tab= summarize(df, cut0=mean(pow0), uncut3=mean(pow3))


load(here("code/simulation_iid_output/sim_ols_n1000.RData"))
pval= do.call(cbind, lapply(sim, function(z) z$pvalue.adjusted))
pow= rowMeans(pval < .05, na.rm=TRUE)

df= data.frame(sim[[1]][,c('varname','zmin','zmax')], pow) |>
  mutate(region= factor(paste0('(',zmin,',',zmax,']')), covariate1= (varname=='x1')) |>
  select(-zmin, -zmax) |>
  group_by(covariate1, region)

tab.pval= summarize(df, pval=mean(pow, na.rm=TRUE))

tab= cbind(tab, tab.pval[,'pval'])
cbind(filter(tab, covariate1)[,-1], filter(tab, !covariate1)[,-1:-2])


###############################################################################
## MSE
###############################################################################

mse= matrix(NA, nrow=2, ncol=3)
rownames(mse)= c('n=100', 'n=1000')
colnames(mse)= c('cut0', 'uncut3', 'VC-BART')

n=100
load(here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))
load(here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))
mse["n=100",1:2]= colMeans(do.call(rbind,lapply(sim, "[[", "mse"))) / n
mse["n=100","VC-BART"]= mean(sapply(sim.vcbart, "[[", "mse")) / n

n=1000
load(here("code", "simulation_iid_output", paste0("sim_n", n,".RData")))
load(here("code", "simulation_iid_output", paste0("sim_vcbart_n", n,".RData")))
mse["n=1000",1:2]= colMeans(do.call(rbind,lapply(sim, "[[", "mse"))) / n
mse["n=1000","VC-BART"]= mean(sapply(sim.vcbart, "[[", "mse")) / n


library(xtable)
xtable(sqrt(mse), digits=c(0,3,3,3))




###############################################################################
## PLOTS FOR n=100
###############################################################################

textsize= 18
#textsize= 35
load(here("code/simulation_iid_output/sim_n100.RData"))
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
    geom_line(aes(lty=Basis, color=Covariate), lwd=3) +
    scale_colour_grey() +
    ylim(0,1) +
    labs(y='Posterior probability of covariate effect') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.75,'cm'))

ggsave("simiid_pp_n100.pdf")




## PLOT POWER FUNCTION ##
#########################

df1= tibble(id, pow0, pow3) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pow0','pow3'), names_to= "Basis", values_to="pp") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Basis= recode(df1$Basis, pow0="Degree 0", pow3="Cubic")

df2= tibble(id, pow0, pow3) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pow0), `Cubic`=mean(pow3)) |>
    gather(key=Basis, value=pp, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Basis,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(lty=Basis, color=Covariate), lwd=3) +
    scale_colour_grey() +
    ylim(0,1) +
    labs(y='Power for covariate inclusion') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))


ggsave("simiid_pow_n100.pdf")



###############################################################################
## PLOTS FOR n=1000
###############################################################################

load(here("code/simulation_iid_output/sim_n1000.RData"))
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
    geom_line(aes(lty=Basis, color=Covariate), lwd=3) +
    scale_colour_grey() +
    ylim(0,1) +
    labs(y='Posterior probability of local null test') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))

ggsave("simiid_pp_n1000.pdf")


## PLOT POWER FUNCTION ##
#########################

df1= tibble(id, pow0, pow3) |>
    filter(covariate == 1) |>
    pivot_longer(cols=c('pow0','pow3'), names_to= "Basis", values_to="pp") |>
    select(-covariate) |>
    transform(Covariate="x1")
df1$Basis= recode(df1$Basis, pow0="Degree 0", pow3="Cubic")

df2= tibble(id, pow0, pow3) |>
    filter(covariate != 1) |>
    group_by(z) |>
    summarise(`Degree 0`=mean(pow0), `Cubic`=mean(pow3)) |>
    gather(key=Basis, value=pp, 2:3) |>
    transform(Covariate="x2-x10")

df= rbind(df1, df2) |>
    transform(group= paste(Basis,Covariate,sep=", "))

ggplot(df, aes(x=z, y=pp)) +
    geom_line(aes(lty=Basis, color=Covariate), lwd=3) +
    scale_colour_grey() +
    ylim(0,1) +
    labs(y='Power for covariate inclusion') +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.17,.8), legend.key.width=unit(1.7,'cm'))


ggsave("simiid_pow_n1000.pdf")




