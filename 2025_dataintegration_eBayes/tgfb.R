#devtools::install_github("https://github.com/davidrusi/mombf")
library(tidyverse)
library(mombf)
library(hgu133plus2.db)
library(xtable)
source("routines.R")

# Import data (first you should unzip tgfb.zip)

x= read.table("~/github/dataintegration_test/data/tgfb.txt", header=TRUE)
shortlist= as.character(read.table("~/github/dataintegration_test/data/mouse_shortlist.txt", header=TRUE)[,1])
inshort= (colnames(x) %in% shortlist)
Z= matrix(inshort, ncol=1)


# Exploratory analysis. Marginal correlations with TGFB vs. found in mouse study yes/no
r= as.vector(cor(y, x))
assoc= data.frame(inshort, r=r) |>
  transform(inshort= case_match(as.character(inshort), "TRUE" ~ "Yes", "FALSE" ~ "No"))

ggplot(assoc, aes(x=inshort, y=r)) +
    geom_violin(draw_quantiles=c(0.25,0.5,0.75)) +
    geom_point() +
    labs(x= "Found in mouse study", y= "Correlation with TGFB") +
    theme(axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16))


# Bayesian variable selection with eBayes and Beta-Binomial prior

ms.eb= modelSelection_eBayes(y=y, x=x, Z=Z, niter.mcmc= 5000, niter.mstep= 1000)
ms.bb= modelSelection(y=y, x=x, niter= 5000, priorDelta= modelbbprior())
b.eb= coef(ms.eb) #Get BMA estimates and PIP
b.bb= coef(ms.bb) #Get BMA estimates and PIP

# Leave-one-out cross-validation

msloo.eb= kfoldCV.bma(y=y, x=x, Z=Z, eBayes=TRUE,  K=length(y), verbose=TRUE, niter=5000, niter.mstep=1000)
msloo.bb= kfoldCV.bma(y=y, x=x, eBayes=FALSE, priorDelta=modelbbprior(), K=length(y), verbose=TRUE)


# Report results

## BMA. Genes with posterior inclusion probability > 0.5

b= round(cbind(b.eb, b.bb), 3)
b= b[c(-1,-nrow(b)),] #exclude intercept and error variance
sel= (b[,4] > 0.5) | (b[,8] > 0.5)  #covariates selected by median probability model under either analysis
b= b[sel,]
b= b[order(b[,4],decreasing=TRUE),]
rownames(b) %in% colnames(x_short) #which of the selected covariates where in the short list defining Z

probe_ids= sub("X", "", rownames(b))
annotations= select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
cbind(annotations, b)


## Leave-one-out cross-validation. Squared correlation between observations and predictions

cor(msloo.eb$pred, y)^2  # Empirical Bayes
cor(msloo.bb, y)^2       # Beta-Binomial


# Prior inclusion probability for genes out/in the short list

priorprob= 1/(1 + exp(-ms.eb$Z %*% ms.eb$eBayes_hyperpar))
data.frame(inshort, priorprob) |> group_by(inshort) |> summarize(mean(priorprob))


# Compare posterior inclusion probabilities eBayes vs. Beta-Binomial

par(mar=c(4.5,4.5,.1,.1))
col= ifelse(inshort, 'black', 'gray'); pch= ifelse(inshort, 15, 2)
plot(ms.bb$margpp, ms.eb$margpp, ylab="Posterior inclusion probability (EBayes)", xlab="Posterior inclusion probability (Beta-Binomial)", cex.lab= 1.5, cex.axis=1.5, col=col, pch=pch, xlim=c(0,1), ylim=c(0,1))
abline(0,1)
legend('topleft', c('In mouse study','Not in mouse study'), col=c('black','gray'), pch=c(15,2), cex=1.5)






