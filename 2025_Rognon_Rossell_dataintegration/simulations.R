###########################################################################################
##
##  R CODE TO REPRODUCE SIMULATION STUDY
##
###########################################################################################

library(modelSelection)
#library(BMS)
library(parallel)
library(mvtnorm)
library(lars)
library(parcor)
#Install parcor from https://cran.r-project.org/src/contrib/Archive/parcor/
#It requires first installing package ppls, from https://cran.r-project.org/src/contrib/Archive/ppls/
library(tidyverse)
source("routines.R")


###########################################################################################
## 1. RUN SIMULATIONS
###########################################################################################

# Define simulation settings
# n: sample size
# p: total number of variables. All have coefficient=0 except the last five with coefficients=0.6, 1.2, 1.8, 2.4, 3.0
# rho: correlation between predictors
# theta.min: smallest non-zero regression coefficient
# theta.min: largest non-zero regression coefficient
#
# Further simulation parameters
# w: hyper-parameter value (shared across all settings)
# nsim: number of simulations to perform in each simulation setting

simpars <- matrix(NA,nrow=3,ncol=5); colnames(simpars) <- c('n','p','rho','theta.min','theta.max')
simpars[1,]  <- c(200,  100, 0.5, 1/3, 2/3)
simpars[2,]  <- c(200,  200, 0.5, 1/3, 2/3)
simpars[3,]  <- c(100,  200, 0.5, 1/3, 2/3)
logit= function(x) log(x/(1-x))
#w= c(logit(.05), 2, 0); scenario= 'scenario1'
#w= c(logit(.05), 1, 0); scenario= 'scenario2'
#w= c(logit(.05), 0, 0); scenario= 'scenario3'
#w= c(logit(.05), 1.5, 0); scenario= 'scenario4'
w= c(logit(.05), 0.75, 0); scenario= 'scenario5'
nsim <- 100


# Run simulations
mc.cores <- 1  #if package parallel is not available, set mc.cores <- 1
for (i in 1:nrow(simpars)) {
  n <- simpars[i,'n']; p <- simpars[i,'p']; rho <- simpars[i,'rho']
  theta.min <- simpars[i,'theta.min']; theta.max <- simpars[i,'theta.max']
  fname= paste('../data/lmsim_',scenario,'_p',p,'_n',n,'_rho',rho,'.RData',sep='')
  # Bayesian methods with pMOM prior on regression coefficient
  sim.eb= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='all', theta.min=theta.min, theta.max=theta.max, rho=rho, niter=5000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.eb0= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='intercept', theta.min=theta.min, theta.max=theta.max, rho=rho, niter=5000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.bb= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='none', priorModel=modelbbprior(), theta.min=theta.min, theta.max=theta.max, rho=rho, niter=5000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  # Bayesian methods with Zellner's prior on regression coefficient
  sim.ebz= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='all', priorCoef=zellnerprior(), theta.min=theta.min, theta.max=theta.max, rho=rho, niter=5000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.ebz0= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='intercept', priorCoef=zellnerprior(), theta.min=theta.min, theta.max=theta.max, rho=rho, niter=5000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.bbz= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='none', priorCoef=zellnerprior(), priorModel=modelbbprior(), theta.min=theta.min, theta.max=theta.max, rho=rho, niter=5000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  # Penalized likelihood methods
  sim.lasso <- mclapply(1:nsim, simPenLhood, p=p, n=n, theta.min=theta.min, theta.max=theta.max, rho=rho, method='LASSO', mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.alasso <- mclapply(1:nsim, simPenLhood, p=p, n=n, theta.min=theta.min, theta.max=theta.max, rho=rho, method='ADALASSO', mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.scad <- mclapply(1:nsim, simPenLhood, p=p, n=n, theta.min=theta.min, theta.max=theta.max, rho=rho, method='SCAD', mc.cores=mc.cores, mc.preschedule=FALSE)
  save(sim.eb, sim.eb0, sim.bb, sim.ebz, sim.ebz0, sim.bbz, sim.lasso, sim.alasso, sim.scad, file=fname, compress='bzip2')
  cat("\n")
}



# Run simulations assessing sensitivity to g-prior's dispersion parameter

zp1 <- zellnerprior(0.1)
zp2 <- zellnerprior(10)
for (i in 1:nrow(simpars)) {
  cat("Simulation scenario ", i, "\n")
  n <- simpars[i,'n']; p <- simpars[i,'p']; rho <- simpars[i,'rho']
  theta.min <- simpars[i,'theta.min']; theta.max <- simpars[i,'theta.max']
  fname= paste('../data/lmsim_',scenario,'_p',p,'_n',n,'_rho_sensitivity',rho,'.RData',sep='')
  sim.ebz1= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='all', priorCoef=zp1, theta.min=theta.min, theta.max=theta.max, rho=rho, niter=1000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  sim.ebz2= mclapply(1:nsim, simBMA, p=p, n=n, eBayes='all', priorCoef=zp1, theta.min=theta.min, theta.max=theta.max, rho=rho, niter=1000, niter.mstep=1000, mc.cores=mc.cores, mc.preschedule=FALSE)
  save(sim.ebz1, sim.ebz2, file=fname, compress='bzip2')
  cat("\n")
}




###########################################################################################
## 2. BUILD FIGURES FROM SIMULATION OUTPUT
###########################################################################################

library(xtable)
library(ggrepel)

# Create tables summarizing all simulation scenarios

simpars= matrix(NA,nrow=3,ncol=5); colnames(simpars)= c('n','p','rho','theta.min','theta.max')
simpars[1,] = c(200,  100, 0.5, 1/3, 2/3)
simpars[2,] = c(200,  200, 0.5, 1/3, 2/3)
simpars[3,] = c(100,  200, 0.5, 1/3, 2/3)
logit= function(x) log(x/(1-x))
wmat= matrix(NA, nrow=5, ncol=3)
wmat[1,]= c(logit(.05), 2, 0) 
wmat[2,]= c(logit(.05), 1, 0) 
wmat[3,]= c(logit(.05), 0, 0) 
wmat[4,]= c(logit(.05), 1.5, 0)
wmat[5,]= c(logit(.05), 0.75, 0)
scenarios= c('scenario1','scenario2','scenario3','scenario4','scenario5')

simoutput= taboutput= vector("list", length(scenarios) * nrow(simpars))
for (j in 1:length(scenarios)) {
  scenario= scenarios[j]
  w= wmat[j,]
  for (i in 1:nrow(simpars)) {
    n= simpars[i,'n']; p= simpars[i,'p']; rho= simpars[i,'rho']
    fname= paste('../data/lmsim_',scenario,'_p',p,'_n',n,'_rho',rho,'.RData',sep='')
    load(fname)
    #
    simoutput.eb  = data.frame(t(sapply(sim.eb, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes')
    #
    simoutput.eb0  = data.frame(t(sapply(sim.eb0, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes (Intercept)')
    #
    simoutput.bb  = data.frame(t(sapply(sim.bb, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Beta-Binomial')
    #
    simoutput.ebz  = data.frame(t(sapply(sim.ebz, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes (Zellner)')
    # 
    simoutput.ebz0  = data.frame(t(sapply(sim.ebz0, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes (Zellner, Intercept)')
    # 
    simoutput.bbz  = data.frame(t(sapply(sim.bbz, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Beta-Binomial (Zellner)')
    # 
    simoutput.lasso  = data.frame(t(sapply(sim.lasso, summarizeSim))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='LASSO')
    # 
    simoutput.alasso  = data.frame(t(sapply(sim.alasso, summarizeSim))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='ALASSO')
    # 
    simoutput.scad  = data.frame(t(sapply(sim.scad, summarizeSim))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='SCAD')
    #
    idx= (j-1) * nrow(simpars) + i
    simoutput[[idx]]= rbind(simoutput.eb, simoutput.eb0, simoutput.bb, simoutput.ebz, simoutput.ebz0, simoutput.bbz, simoutput.lasso, simoutput.alasso, simoutput.scad)
    simoutput[[idx]]= data.frame(scenario=scenario, n=n, p=p, simoutput[[idx]])
    #
    taboutput[[idx]]= group_by(simoutput[[idx]], metric, method) |>
      summarize(mean=mean(value)) |>
      pivot_wider(id_cols=c('method'), names_from = metric, values_from = mean)
    taboutput[[idx]]= data.frame(scenario=scenario, n=n, p=p, taboutput[[idx]])
  }
}
taboutput= do.call(rbind, taboutput)  #put results from all settings into a single table


# Display table with MSE, Power and FDR for each scenario in latex format
method.levels= c("Empirical Bayes", "Empirical Bayes (Intercept)", "Beta-Binomial", "Empirical Bayes (Zellner)", "Empirical Bayes (Zellner, Intercept)", "Beta-Binomial (Zellner)", "LASSO", "ALASSO", "SCAD")
    
for (i in 1:length(scenarios)) {
    # Create the table
    tabsel= filter(taboutput, scenario == scenarios[i], method %in% method.levels) |>
      select(method, n, p, mse, pow, fdr) |>
      transform(np= paste("(n=",n,", p=",p,")",sep="")) |>
      transform(method = factor(method, levels=method.levels)) |>
      transform(np = factor(np, levels=c("(n=100, p=200)","(n=200, p=200)","(n=200, p=100)"))) |>
      arrange(np, method) |>
      select(-np)
    # Output the table in latex style
    xres= xtable(tabsel, digits= c(0,0,0,0,2,2,2))
    cat("SCENARIO", i, "\n\n")
    print(xres, include.rownames = FALSE)
    cat("-----------------------------------------------------------------------------\n\n")
}


# Create figures
# To avoid cluttering, the intercept-only eBayes results are not plotted
# They can be added to the plots by uncommenting the lines featuring #to see also intercept-only eBayes

## Define colors and line types for each method
style_mapping = data.frame(rbind(c("Empirical Bayes", "black", "solid"), c("Beta-Binomial", "black", "dashed"), c("LASSO", "gray50", "solid"), c("ALASSO", "gray50", "dashed"), c("SCAD", "gray50", "dotted")))
#style_mapping = data.frame(rbind(c("Empirical Bayes", "black", "solid"), c("Empirical Bayes (Intercept)", "black", "dashed"), c("Beta-Binomial", "black", "dotted"), c("LASSO", "gray50", "solid"), c("ALASSO", "gray50", "dashed"), c("SCAD", "gray50", "dotted"))) #to see also intercept-only eBayes
names(style_mapping)= c("method","color","linetype")

methods= c('Beta-Binomial', 'Empirical Bayes', 'LASSO', 'ALASSO', 'SCAD') 
#methods= c('Beta-Binomial', 'Empirical Bayes', 'Empirical Bayes (Intercept)', 'LASSO', 'ALASSO', 'SCAD') #to see also intercept-only eBayes

cex.axis= 40; cex.title= 40; cex.text= 13; point.size= 6; lwd= 3

## Scenario 1

scenario.sel = "scenario1"
tab2plot= transform(taboutput, np= paste("(n=",n,", p=",p,")",sep="")) |>
    transform(method = factor(method), np = factor(np, levels=c("(n=100, p=200)","(n=200, p=200)","(n=200, p=100)"))) |>
    filter(scenario==scenario.sel, method %in% methods) |>
    merge(style_mapping, by='method')

ggplot(tab2plot, aes(x=np, y=mse, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd=lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    scale_y_log10() +
    labs(y="MSE", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_MSE.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_MSE0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)



ggplot(tab2plot, aes(x=np, y=fdr, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd=lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    labs(y="FDR", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_FDR.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_FDR0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)


ggplot(tab2plot, aes(x=np, y=pow, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd= lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    labs(y="Power", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_power.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_power0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)


## Scenario 2

scenario.sel = "scenario2"
tab2plot= transform(taboutput, np= paste("(n=",n,", p=",p,")",sep="")) |>
    transform(method = factor(method), np = factor(np, levels=c("(n=100, p=200)","(n=200, p=200)","(n=200, p=100)"))) |>
    filter(scenario==scenario.sel, method %in% methods) |>
    merge(style_mapping, by='method')

ggplot(tab2plot, aes(x=np, y=mse, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd=lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    scale_y_log10() +
    labs(y="MSE", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_MSE.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_MSE0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)


ggplot(tab2plot, aes(x=np, y=fdr, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd=lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    labs(y="FDR", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_FDR.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_FDR0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)


ggplot(tab2plot, aes(x=np, y=pow, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd= lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    labs(y="Power", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_power.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_power0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)



## Scenario 3

scenario.sel = "scenario3"
tab2plot= transform(taboutput, np= paste("(n=",n,", p=",p,")",sep="")) |>
    transform(method = factor(method), np = factor(np, levels=c("(n=100, p=200)","(n=200, p=200)","(n=200, p=100)"))) |>
    filter(scenario==scenario.sel, method %in% methods) |>
    merge(style_mapping, by='method')

ggplot(tab2plot, aes(x=np, y=mse, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd=lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    scale_y_log10() +
    labs(y="MSE", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_MSE.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_MSE0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)


ggplot(tab2plot, aes(x=np, y=fdr, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd=lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    labs(y="FDR", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_FDR.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_FDR0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)


ggplot(tab2plot, aes(x=np, y=pow, group=method, color=color)) +
    geom_point(aes(shape = method), size = point.size) +
    geom_line(aes(linetype = linetype), lwd= lwd) +
    scale_color_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    scale_linetype_identity(guide = "legend", labels = style_mapping$method, breaks = style_mapping$color) +
    geom_text_repel(aes(label=method), direction= "y", hjust= 1, position= position_nudge(x=-0.5), x=0.95, size= cex.text, show.legend=FALSE, data= filter(tab2plot, np == "(n=100, p=200)")) +
    labs(y="Power", x="") +
    theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size=cex.title), legend.position='none')

fname = paste("../paper/figs/",scenario.sel, "_power.pdf", sep="")
#fname = paste("../paper/figs/",scenario.sel, "_power0.pdf", sep="") #to see also intercept-only eBayes
ggsave(fname)




## Table for the sensitivity analysis wrt Zellner's prior dispersion g

simpars= matrix(NA,nrow=3,ncol=5); colnames(simpars)= c('n','p','rho','theta.min','theta.max')
simpars[1,] = c(200,  100, 0.5, 1/3, 2/3)
simpars[2,] = c(200,  200, 0.5, 1/3, 2/3)
simpars[3,] = c(100,  200, 0.5, 1/3, 2/3)
logit= function(x) log(x/(1-x))
scenarios= c('scenario5')
wmat= matrix(NA, nrow=5, ncol=3)
wmat[1,]= c(logit(.05), 0.75, 0)

simoutput= taboutput= vector("list", length(scenarios) * nrow(simpars))
for (j in 1:length(scenarios)) {
  scenario= scenarios[j]
  w= wmat[j,]
  for (i in 1:nrow(simpars)) {
    n= simpars[i,'n']; p= simpars[i,'p']; rho= simpars[i,'rho']
    fname= paste('../data/lmsim_',scenario,'_p',p,'_n',n,'_rho',rho,'.RData',sep='')
    fname_sens= paste('../data/lmsim_',scenario,'_p',p,'_n',n,'_rho_sensitivity',rho,'.RData',sep='')
    load(fname)
    load(fname_sens)
    #
    simoutput.ebz  = data.frame(t(sapply(sim.ebz, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes (Zellner, g = 1)')
    #
    simoutput.ebz1  = data.frame(t(sapply(sim.ebz1, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes (Zellner, g = 0.1)')
    #
    simoutput.ebz2  = data.frame(t(sapply(sim.ebz2, summarizeSim, threshold=0.95))) |>
      pivot_longer(cols=c('fdr','pow','mse','meanpp0','meanpp1'), names_to='metric') |>
      transform(method='Empirical Bayes (Zellner, g = 10)')
    #
    idx= (j-1) * nrow(simpars) + i
    simoutput[[idx]]= rbind(simoutput.ebz, simoutput.ebz1, simoutput.ebz2)
    simoutput[[idx]]= data.frame(scenario=scenario, n=n, p=p, simoutput[[idx]])
    #
    taboutput[[idx]]= group_by(simoutput[[idx]], metric, method) |>
      summarize(mean=mean(value)) |>
      pivot_wider(id_cols=c('method'), names_from = metric, values_from = mean)
    taboutput[[idx]]= data.frame(scenario=scenario, n=n, p=p, taboutput[[idx]])
  }
}
taboutput= do.call(rbind, taboutput)  #put results from all settings into a single table


# Display table with MSE, Power and FDR for each scenario in latex format
method.levels= c("Empirical Bayes (Zellner, g = 1)", "Empirical Bayes (Zellner, g = 0.1)", "Empirical Bayes (Zellner, g = 10)")
    
for (i in 1:length(scenarios)) {
    # Create the table
    tabsel= filter(taboutput, scenario == scenarios[i], method %in% method.levels) |>
      select(method, n, p, mse, pow, fdr) |>
      transform(np= paste("(n=",n,", p=",p,")",sep="")) |>
      transform(method = factor(method, levels=method.levels)) |>
      transform(np = factor(np, levels=c("(n=100, p=200)","(n=200, p=200)","(n=200, p=100)"))) |>
      arrange(np, method) |>
      select(-np)
    # Output the table in latex style
    xres= xtable(tabsel, digits= c(0,0,0,0,3,3,3))
    cat("SCENARIO", i, "\n\n")
    print(xres, include.rownames = FALSE)
    cat("-----------------------------------------------------------------------------\n\n")
}





###########################################################################################
## 3. RUN TIMES
###########################################################################################

simpars <- matrix(NA, nrow=4, ncol=5); colnames(simpars) <- c('n','p','rho','theta.min','theta.max')
simpars[1,]  <- c(500,  100, 0.5, 1/3, 2/3)
simpars[2,]  <- c(500,  500, 0.5, 1/3, 2/3)
simpars[3,]  <- c(500,  1000, 0.5, 1/3, 2/3)
simpars[4,]  <- c(500,  2000, 0.5, 1/3, 2/3)
logit= function(x) log(x/(1-x))
w= c(logit(.01), 0.5, 0)
nsim <- 5


# Run simulations
time_mat <- matrix(NA, nrow=nsim, ncol=4)
colnames(time_mat) <- c("bb (niter=5,000)", "eb (niter= 5,000; niter.mstep=1,000 iter)", "eb (niter= 5,000; niter.mstep=100)", "eb (niter= 1,000; niter.mstep=100)")
runtime <- lapply(1:nrow(simpars), function(i) time_mat)
dif_mat <- matrix(NA, nrow=nsim, ncol=2)
colnames(dif_mat) <- c("niter.mstep 1,000 to 100", "niter 5,000 to 1,000")
ppdif <- lapply(1:nrow(simpars), function(i) dif_mat)
nonzero <- matrix(NA, nrow=nrow(simpars), ncol=nsim)
rownames(nonzero) <- paste('simpars', 1:nrow(nonzero), sep='')

for (i in 1:nrow(simpars)) {
  cat("Simulation", i)
  n <- simpars[i,'n']; p <- simpars[i,'p']; rho <- simpars[i,'rho']
  theta.min <- simpars[i,'theta.min']; theta.max <- simpars[i,'theta.max']
  for (j in 1:nsim) {
    set.seed(j)
    simpar <- sim_parameters(theta.min=theta.min, theta.max=theta.max, p=p, w=w, rho=rho)
    sim <- sim_data(theta=simpar$theta, n=n, rho=rho, phi=1)
    df <- data.frame(y=sim$y, sim$x)
    nonzero[i,j] <- sum(simpar$theta != 0)
    # Beta-Binomial (5,000 iter)
    runtime[[i]][j,1] <- system.time(ms.bb <- modelSelection(y ~ ., data=df, priorCoef=zellnerprior(), niter=5000, verbose=FALSE))['elapsed']
    # eBayes (5,000 iter, 1,000 M-step iter)
    runtime[[i]][j,2] <- system.time(ms.eb1 <- modelSelection_eBayes(y=sim$y, x=sim$x, data=df, Z=simpar$Z, niter.mcmc=5000, niter.mstep=1000, verbose=FALSE))['elapsed']
    # eBayes (5,000 iter, 100 M-step iter)
    runtime[[i]][j,3] <- system.time(ms.eb2 <- modelSelection_eBayes(y=sim$y, x=sim$x, data=df, Z=simpar$Z, niter.mcmc=5000, niter.mstep=100, verbose=FALSE))['elapsed']
    # eBayes (1,000 iter, 100 M-step iter)
    runtime[[i]][j,4] <- system.time(ms.eb3 <- modelSelection_eBayes(y=sim$y, x=sim$x, data=df, Z=simpar$Z, niter.mcmc=1000, niter.mstep=100, verbose=FALSE))['elapsed']
    # MAE in estimated posterior marginal inclusion probabilities
    exclude <- c(1, p+2) #rows corresponding to intercept and error variance (always = 1)
    pp <- matrix(NA, nrow=p, ncol=3)
    pp[,1] <- coef(ms.eb1)[-exclude, 'margpp']
    pp[,2] <- coef(ms.eb2)[-exclude, 'margpp']
    pp[,3] <- coef(ms.eb3)[-exclude, 'margpp']
    ppdif[[i]][j,] <- colMeans(abs(pp[,1:2] - pp[,2:3]))
    cat(".")
  }
  cat("\n")
}



runtime_mean <- t(sapply(runtime, colMeans))
ppdif_mean <- t(sapply(ppdif, colMeans))
runtime_mean <- cbind(simpars[,c('p','n')], nonzero= rowMeans(nonzero), runtime_mean, ppdif_mean)
xtable::xtable(t(runtime_mean))





###########################################################################################
## 4. SENSITIVITY TO HYPER-PARAMETER INITIALIZATION
###########################################################################################

simpars <- matrix(NA, nrow=2, ncol=5); colnames(simpars) <- c('n','p','rho','theta.min','theta.max')
simpars[1,]  <- c(500,  100, 0.5, 1/3, 2/3)
simpars[2,]  <- c(500,  500, 0.5, 1/3, 2/3)
logit= function(x) log(x/(1-x))
w= c(logit(.01), 0.5, 0)
nsim <- 5


# Run simulations
w_hat <- as.data.frame(matrix(NA, nrow= nrow(simpars) * nsim * 3, ncol=6))
names(w_hat) <- c("p", "simulation", "method", "w0", "w1", "w2")
w_hat$p <- rep(simpars[,'p'], each=nsim * 3)
w_hat$simulation <- rep(rep(1:nsim, each= 3), nrow(simpars))
w_hat$method <- rep(c('Zero init, M=1000', 'Beta-Binomial init, M=1000', 'Zero init, M=100'), nrow(simpars) * nsim)

for (i in 1:nrow(simpars)) {
  cat("Simulation", i)
  n <- simpars[i,'n']; p <- simpars[i,'p']; rho <- simpars[i,'rho']
  theta.min <- simpars[i,'theta.min']; theta.max <- simpars[i,'theta.max']
  for (j in 1:nsim) {
    set.seed(j)
    simpar <- sim_parameters(theta.min=theta.min, theta.max=theta.max, p=p, w=w, rho=rho)
    sim <- sim_data(theta=simpar$theta, n=n, rho=rho, phi=1)
    df <- data.frame(y=sim$y, sim$x)
    # eBayes (init hyper-parameters at w=0)
    ms.eb0 <- modelSelection_eBayes(y=sim$y, x=sim$x, data=df, Z=simpar$Z, wini= c(0,0,0), niter.mcmc=5000, niter.mstep=1000, verbose=FALSE)
    # eBayes (default init from Beta-Binomial)
    ms.eb1 <- modelSelection_eBayes(y=sim$y, x=sim$x, data=df, Z=simpar$Z, niter.mcmc=5000, niter.mstep=1000, verbose=FALSE)
    # eBayes (init hyper-parameters at w=0, less iterations)
    ms.eb2 <- modelSelection_eBayes(y=sim$y, x=sim$x, data=df, Z=simpar$Z, wini= c(0,0,0), niter.mcmc=5000, niter.mstep=100, verbose=FALSE)
    # Save hyper-parameter estimates
    what <- t(cbind(ms.eb0$eBayes_hyperpar, ms.eb1$eBayes_hyperpar, ms.eb2$eBayes_hyperpar))
    w_hat[(w_hat$p == p) & (w_hat$simulation == j), c('w0','w1','w2')] <- what
    cat(".")
  }
  cat("\n")
}


# Function to obtain mean absolute difference in w_hat between two methods
mad_what <- function(w_hat, method1, method2) {
  # Compute absolute differences in w_hat
  w0 <- dplyr::filter(w_hat, method == method1) |>
    dplyr::select(-simulation, -method)
  w1 <- dplyr::filter(w_hat, method == method2) |>
    dplyr::select(-simulation, -method)
  wdif <- abs(dplyr::select(w1, -p) - dplyr::select(w0, -p))
  # Add (p, simulation) identifiers
  wdif$p <- w0$p
  wdif$simulation <- w0$simulation
  # Compute mean for each p
  ans <- dplyr::group_by(wdif, p) |>
    dplyr::summarize(across(w0:w2, mean))
  return(ans)
}


mad_what(w_hat, method1="Beta-Binomial init, M=1000", method2="Zero init, M=1000")
mad_what(w_hat, method1="Beta-Binomial init, M=1000", method2="Zero init, M=100")

