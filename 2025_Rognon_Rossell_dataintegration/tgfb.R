#devtools::install_github("https://github.com/davidrusi/mombf")
library(tidyverse)
library(modelSelection)
library(hgu133plus2.db)
library(xtable)
source("routines.R")
priorCoef= momprior()
#priorCoef= zellnerprior()  #uncomment this line for Zellner prior results

# Import data

x= read.table("~/github/dataintegration_test/data/tgfb.txt", header=TRUE)
shortlist= as.character(read.table("~/github/dataintegration_test/data/mouse_shortlist.txt", header=TRUE)[,1])
y= x[,1]
x= x[,-1]
inshort= (colnames(x) %in% shortlist)
Z= matrix(inshort, ncol=1)

# Bayesian variable selection with eBayes and Beta-Binomial prior

system.time(ms.eb <- modelSelection_eBayes(y=y, x=x, Z=Z, priorCoef=priorCoef, niter.mcmc= 5000, niter.mstep= 1000))
system.time(ms.bb <- modelSelection(y=y, x=x, niter= 5000, priorCoef=priorCoef, priorModel= modelbbprior()))
b.eb= coef(ms.eb) #Get BMA estimates and PIP
b.bb= coef(ms.bb) #Get BMA estimates and PIP

# Leave-one-out cross-validation

msloo.eb= kfoldCV.bma(y=y, x=x, Z=Z, eBayes=TRUE,  priorCoef=priorCoef, K=length(y), verbose=TRUE, niter=5000, niter.mstep=1000)
msloo.bb= kfoldCV.bma(y=y, x=x, eBayes=FALSE, priorCoef=priorCoef, priorModel=modelbbprior(), K=length(y), verbose=TRUE)


# Report results

## BMA results. Genes with post prob > 0.5

b= round(cbind(b.eb, b.bb), 3)
b= b[c(-1,-nrow(b)),] #exclude intercept and error variance
sel= (b[,4] > 0.5) | (b[,8] > 0.5)  #covariates selected by median probability model under either analysis
b= round(b[sel,], 2)
b= b[order(b[,4],decreasing=TRUE),]
rownames(b) %in% shortlist #which of the selected covariates where in the short list defining Z

probe_ids= sub("X", "", rownames(b))
annotations= select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
tab= cbind(annotations[,c(-1,-3)], b)
rownames(tab)= NULL
tab

cor(msloo.eb$pred, y)^2
cor(msloo.bb$pred, y)^2 


# Prior inclusion probability for genes out/in the short list

priorprob= 1/(1 + exp(-ms.eb$Z %*% ms.eb$eBayes_hyperpar))
data.frame(inshort, priorprob) |> group_by(inshort) |> summarize(mean(priorprob))


# Confidence interval for hyper-parameters
hist(msloo.eb$hyperpar[,2])
quantile(msloo.eb$hyperpar[,2], prob=c(0.025,0.975))


# Compare posterior inclusion probabilities eBayes vs. Beta-Binomial

par(mar=c(4.5,4.5,.1,.1))
col= ifelse(inshort, 'black', 'gray'); pch= ifelse(inshort, 15, 2)
plot(ms.bb$margpp, ms.eb$margpp, ylab="Posterior inclusion probability (EBayes)", xlab="Posterior inclusion probability (Beta-Binomial)", cex.lab= 1.5, cex.axis=1.5, col=col, pch=pch, xlim=c(0,1), ylim=c(0,1))
abline(0,1)
legend('topleft', c('In mouse study','Not in mouse study'), col=c('black','gray'), pch=c(15,2), cex=1.5)


# Exploratory analysis
r= as.vector(cor(y, x))
assoc= data.frame(inshort, margpp=ms.bb$margpp, r=r) |>
  transform(inshort= case_match(as.character(inshort), "TRUE" ~ "Yes", "FALSE" ~ "No"))


## Marginal correlations with TGFB vs. found in mouse study yes/no

ggplot(assoc, aes(x=inshort, y=r)) +
    geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
    geom_point() +
    labs(x= "Found in mouse study", y= "Correlation with TGFB") +
    theme(axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16))









