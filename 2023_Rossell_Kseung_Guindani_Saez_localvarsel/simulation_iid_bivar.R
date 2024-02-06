###############################################################################
##
## THIS SCRIPT REPRODUCES THE EXAMPLE WITH BIVARIATE COORDINATES z
##
###############################################################################

library(mombf)
library(parallel)
library(tidyverse)


###############################################################################
## DEFINE FUNCTIONS
###############################################################################

truemean= function(x,z) {
    ans= double(nrow(x))
    group1= (x[,1]==1)
    sel= (z[,1] <=0) | (z[,2] <=0)
    ans[group1]= ifelse(z[group1,1] <=0, cos(z[group1,1]), 1) + ifelse(z[group1,2] <=0, cos(z[group1,2]), 1)
    ans[!group1]= ifelse(z[!group1,1]<=0, cos(z[!group1,1]), 1/(z[!group1,1]+1)^2) + ifelse(z[!group1,2]<=0, cos(z[!group1,2]), 1/(z[!group1,2]+1)^2)
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
    zseq= seq(-3,3,length=round(sqrt(n/2)))
    z= expand.grid(zseq, zseq)
    x= matrix(NA, nrow=2*nrow(z), ncol=p)
    x[,1]= rep(0:1,each=nrow(z))
    for (j in 2:p) x[,j]= x[,1] + rnorm(nrow(x))
    z= rbind(z,z)
    m= truemean(x,z)
    y= m + rnorm(length(m), 0, sd=sd.error)
    ans= list(y=y, x=x, z=z, m=m)
    return(ans)
}


seed=1; n=1000; p=2; sd.error=0.25
sim= simdata(seed=seed, n=n, p=p, sd.error=sd.error)

fit0= localnulltest(sim$y, x=sim$x, z=sim$z, cutdegree=0, nlocalknots=8, niter=1000, verbose=TRUE)
b= coef(fit0)



#Plot true mean group differences
zseq= seq(-3,3,length=50)
z= expand.grid(zseq, zseq)
x= matrix(NA, nrow=2*nrow(z), ncol=1)
x[,1]= rep(0:1,each=nrow(z))
z= rbind(z,z)
m= truemean(x, z)
mat= data.frame(m, x1=x[,1], z1=z[,1], z2=z[,2])
mat$m= mat$m + rnorm(nrow(mat), sd=0.0001)
m0= filter(mat, x1==0) |> rename(m0=m) |> select(z1, z2, m0)
m1= filter(mat, x1==1) |> rename(m1=m) |> select(m1)
mdif= cbind(m0, select(m1, m1)) |> mutate(mdif= m1 - m0)

ggplot(mdif, aes(z1, z2)) +
    geom_tile(aes(fill=mdif)) +
    scale_fill_gradient2(low = "white", mid = "gray", high = "black",  midpoint = .02) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"), legend.title=element_text(size=20), legend.text=element_text(size=20), legend.position = c(0.8, 0.2)) +
    guides(fill=guide_legend(title="Group differences"))
ggsave("bspline_fit_cubic_bivar.pdf")


#Plot posterior probabilities
b2= mutate(b, reject= ifelse(margpp > 0.95, 'Yes', 'No')) |> filter(covariate==1)
ggplot(b2, aes(z1, z2, color=reject, pch=reject)) +
    geom_point(show.legend= FALSE) +
    geom_label(aes(label=round(margpp,3)), size=7, label.size=1, show.legend=FALSE) +
    scale_colour_grey(start= 0.8, end=0.2) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
ggsave("bspline_cubic_pp_bivar.pdf")


