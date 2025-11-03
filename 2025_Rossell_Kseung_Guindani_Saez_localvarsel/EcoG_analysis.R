library(mombf)
library(tidyverse)
setwd("~/github/localnulltest/code")

load("~/github/localnulltest/data/ofc3_prewindowed.RData")

timings_w <- seq(-900, 1900, 50)

# n patients and electrodes by patient
n_pts <- max(gamble_data$patient)
n_elecs <- t(unlist(lapply(buttonpressw, function(x) dim(x)[3])))

### fda local null test 

## set patient number, time window 
## and electrode of interest
## (5,29) from -450 to 1400

patient_num=5 
elec_num=29 

pt_gamble_data <- gamble_data %>% filter(patient==patient_num)
pt_gamble_data$trial <- 1:nrow(pt_gamble_data)
pt_ecog <- buttonpressw[[patient_num]]
pt_elec_dat <- pt_ecog[,,elec_num]

npoints=length(timings_w)
ntrials=max(pt_gamble_data$trial)

## centering the trials to a common 
## baseline intercept

intercepts <- rep(NA, ntrials)
center_elec_dat <- matrix(NA, dim(pt_elec_dat)[1], dim(pt_elec_dat)[2])
for(i in 1:ntrials)
{
  intercepts[i] <- lm(pt_elec_dat[i,] ~ timings_w)$coefficients[1]
  center_elec_dat[i,] <- pt_elec_dat[i, ]-intercepts[i] 
}


## subsetting the data to
## specified time window

subset_timings <- timings_w[timings_w>-450 & timings_w<1500]
subset_timings_index <- which(timings_w>-450 & timings_w<1500)
center_elec_dat <- center_elec_dat[, subset_timings_index]
timings_w <- subset_timings
npoints <- length(timings_w)

## reshaping data to fit 
## localnulltest function args

y=c(t(center_elec_dat))
x=as.numeric(rep(pt_gamble_data[, "win.ind"], each=npoints))
id=as.numeric(rep(pt_gamble_data$trial, each=npoints))
z=rep(timings_w, times=ntrials)
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))

##########################################
## local null test for single covariates 
##########################################

## win/loss

set.seed(1234)
fitwin=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                  function_id=df$id, Sigma='AR/MA', data=df,
                 priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitwin$covareffects$margpp

# rpe reward prediction error (RPE) 
# (difference between the obtained reward 
# and the expected value of the gamble),

x <- as.numeric(rep(pt_gamble_data[, "rpe"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(1234)
fitrpe=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                       function_id=df$id, Sigma='AR/MA', data=df,
                       priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitrpe$covareffects$margpp

# regret 
# the amount of extra money that would 
# have been won for the non-chosen option 

x <- as.numeric(rep(pt_gamble_data[, "regret"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(1234)
fitregret=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                       function_id=df$id, Sigma='AR/MA', data=df,
                       priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitregret$covareffects$margpp

## gamble.ind
## choice of gamble (e.g., whether subject chose to gamble, w)

x <- as.numeric(rep(pt_gamble_data[, "gamble.ind"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(123)
fitgamble=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                         function_id=df$id, Sigma='AR/MA', data=df,
                         priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitgamble$covareffects$margpp

## previous.risk 

x <- as.numeric(rep(pt_gamble_data[, "previous.risk"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(1234)
fitprevrisk=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                            function_id=df$id, Sigma='AR/MA', data=df,
                            priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitprevrisk$covareffects$margpp

## previous.win.ind

x <- as.numeric(rep(pt_gamble_data[, "previous.win.ind"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(1234)
fitprevrisk=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                              function_id=df$id, Sigma='AR/MA', data=df,
                              priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitprevrisk$covareffects$margpp

## previous.rpe

x <- as.numeric(rep(pt_gamble_data[, "previous.rpe"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(1234)
fitprevrpe=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                              function_id=df$id, Sigma='AR/MA', data=df,
                              priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitprevrpe$covareffects$margpp


## previous.regret

x <- as.numeric(rep(pt_gamble_data[, "previous.regret"], each=npoints))
df <- data.frame(y=matrix(y), x=matrix(as.numeric(x)), id=as.matrix(id), z=as.matrix(z))
set.seed(1234)
fitprevregret=localnulltest_fda(y=df$y, x=df$x, z=df$z,
                             function_id=df$id, Sigma='AR/MA', data=df,
                             priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
#fitprevregret$covareffects$margpp



##########################################
##  Multivariable analysis
##########################################

# outcome & choice related variables 


x_win=as.numeric(rep(pt_gamble_data[, "win.ind"], each=npoints))
x_rpe <- as.numeric(rep(pt_gamble_data[, "rpe"], each=npoints))
x_regret <- as.numeric(rep(pt_gamble_data[, "regret"], each=npoints))
x_gamble <- as.numeric(rep(pt_gamble_data[, "gamble.ind"], each=npoints))

# past trial characteristics

x_previous.risk <- as.numeric(rep(pt_gamble_data[, "previous.risk"], each=npoints)) 
x_previous.win.ind <- as.numeric(rep(pt_gamble_data[, "previous.win.ind"], each=npoints)) 
x_previous.rpe <- as.numeric(rep(pt_gamble_data[, "previous.rpe"], each=npoints)) 
x_previous.regret <- as.numeric(rep(pt_gamble_data[, "previous.regret"], each=npoints)) 

# building dataframe for design matrix

x=cbind(x_win, x_rpe, x_regret, x_gamble, x_previous.risk, 
        x_previous.win.ind, x_previous.rpe, x_previous.regret)

df <- data.frame(y=matrix(y), id=as.matrix(id), z=as.matrix(z))

df <- mutate(df, cov=x) 

# multivariate localnull test 

set.seed(1234)
fitmulti=localnulltest_fda(y=df$y, x=df$cov, z=df$z,
                    function_id=df$id, Sigma='AR/MA', data=df,
                    priorCoef=momprior(), priorGroup=groupmomprior(), mc.cores=1)
b=coef(fitmulti)
#save(fitmulti,b, file="Ecog_results.RData")

##########################################
## Figures 
##########################################

##########################################
## In the main manuscript 
##########################################

library(ggplot2)


## Figure 4 
## EcoG Data as a function of Win
## Differences mean loss
## A figure with mean and 80%CI

# wins: did the gamble resulted in a win or a loss 
# (covariate of interest)

ga_win <- pt_gamble_data[, "win.ind"]==1 
no_win <- !ga_win

ga_win_meanhfa <- apply(center_elec_dat[ga_win, ], 2, mean)
no_win_meanhfa <- apply(center_elec_dat[no_win, ], 2, mean)
ga_win_meanhfa_upper <- ga_win_meanhfa + 1.282 * apply(center_elec_dat[ga_win, ], 2, sd)/sqrt(ntrials)
ga_win_meanhfa_lower <- ga_win_meanhfa - 1.282 * apply(center_elec_dat[ga_win, ], 2, sd)/sqrt(ntrials)
no_win_meanhfa_upper <- no_win_meanhfa + 1.282 * apply(center_elec_dat[no_win, ], 2, sd)/sqrt(ntrials)
no_win_meanhfa_lower <- no_win_meanhfa - 1.282 * apply(center_elec_dat[no_win, ], 2, sd)/sqrt(ntrials)

textsize= 20  
ggplot()+
  geom_line(data = data.frame(x = timings_w, y = no_win_meanhfa), aes(x, y), 
            color = "black", linetype="dashed", size = 0.75) +
  geom_line(data = data.frame(x = timings_w, y = ga_win_meanhfa), aes(x, y), 
            color = "black", linetype="solid", size = 0.75) +
  geom_ribbon(aes(x=timings_w, y=ga_win_meanhfa, ymin = ga_win_meanhfa_lower, ymax = ga_win_meanhfa_upper), 
              fill = "azure4", alpha = 0.3) +
  geom_ribbon(aes(x=timings_w, y=no_win_meanhfa, ymin = no_win_meanhfa_lower, ymax = no_win_meanhfa_upper), 
              fill = "gray68", alpha = 0.3) + #gray68 #lightgray #azure4
  labs(title = "", x = "Time", y = "HFA") +
  guides(shape=FALSE)+
  theme_minimal() +
  theme(axis.text=element_text(size=textsize-8), 
        axis.title=element_text(size=textsize), 
        legend.title=element_text(size=textsize), 
        legend.text=element_text(size=textsize), 
        legend.position=c(.83,.2), legend.key.width=unit(3,'cm'), 
        strip.text.x = element_text(size=textsize))
#ggsave("../drafts/figs/Ecog_trials1.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)

## Plot the marginal posterior probabilities
## for all covariates

b=coef(fitmulti)
b=data.frame(covariate=1:ncol(df$cov), 
            covariaten=c('Win', 'RPE', 'Regret', 'Gamble', 
                           'Prev. Risk', 'Prev. Win', 
                           'Prev. RPE', 'Prev. Regret')) |>
  right_join(b) |>
  select(-covariate) |>
  rename(covariate= covariaten)


textsize= 20
mycols= rep(c('black','darkgrey'), ncol(x)/2 + 1)[1:ncol(x)]
ggplot(b, aes(z1, margpp)) +
  geom_line() +
  #geom_line(aes(group=covariate, lty=covariate)) +
  labs(x='Time', y='Posterior probability of a covariate effect') +
  facet_wrap( ~ covariate) +
  ylim(0,1) + 
  scale_colour_grey() +
  theme_bw() +
  theme(axis.text=element_text(size=textsize-8), 
        axis.title=element_text(size=textsize), 
        legend.title=element_text(size=textsize), 
        legend.text=element_text(size=textsize), 
        legend.position=c(.83,.2), legend.key.width=unit(3,'cm'), 
         strip.text.x = element_text(size=textsize)
  )
#ggsave("../drafts/figs/Ecog_pp_localtests.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)

##########################################
## In the Supplementary Material
##########################################

### Time-varying effects in fitwin

b=coef(fitwin)
b=data.frame(covariate=1:ncol(df$cov), 
             covariaten=c('Win')) |>
  right_join(b) |>
  select(-covariate) |>
  rename(covariate= covariaten)
bsel= filter(b, covariate %in% c('Win'))
ggplot(bsel, aes(z1, estimate)) +
  geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
  labs(x='Time', y='Estimated covariate effect') +
  scale_colour_grey() +
  theme_bw() 
#ggsave("../drafts/figs/Ecog_estimate_fitwin.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)

### Time-varying effects in fitrpe

b=coef(fitrpe)
b=data.frame(covariate=1:ncol(df$cov), 
             covariaten=c('RPE')) |>
  right_join(b) |>
  select(-covariate) |>
  rename(covariate= covariaten)
bsel= filter(b, covariate %in% c('RPE'))
ggplot(bsel, aes(z1, estimate)) +
  geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
  labs(x='Time', y='Estimated covariate effect') +
  scale_colour_grey() +
  theme_bw() 
#ggsave("../drafts/figs/Ecog_estimate_fitrpe.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)

### Time-varying effects in fitregret

b=coef(fitregret)
b=data.frame(covariate=1:ncol(df$cov), 
             covariaten=c('Regret')) |>
  right_join(b) |>
  select(-covariate) |>
  rename(covariate= covariaten)
bsel= filter(b, covariate %in% c('Regret'))
ggplot(bsel, aes(z1, estimate)) +
  geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
  labs(x='Time', y='Estimated covariate effect') +
  scale_colour_grey() +
  theme_bw() 
#ggsave("../drafts/figs/Ecog_estimate_fitregret.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)

### Time-varying effects in fitgamble

b=coef(fitgamble)
b=data.frame(covariate=1:ncol(df$cov), 
             covariaten=c('Gamble')) |>
  right_join(b) |>
  select(-covariate) |>
  rename(covariate= covariaten)
bsel= filter(b, covariate %in% c('Gamble'))
ggplot(bsel, aes(z1, estimate)) +
  geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
  labs(x='Time', y='Estimated covariate effect') +
  scale_colour_grey() +
  theme_bw() 
#ggsave("../drafts/figs/Ecog_estimate_fitgamble.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)



### Time-varying effects  in fitmulti

b=coef(fitmulti)
b=data.frame(covariate=1:ncol(df$cov), 
             covariaten=c('Win', 'RPE', 'Regret', 'Gamble', 
                          'Prev. Risk', 'Prev. Win', 
                          'Prev. RPE', 'Prev. Regret')) |>
  right_join(b) |>
  select(-covariate) |>
  rename(covariate= covariaten)
bsel= filter(b, covariate %in% c('Win', 'RPE'))
ggplot(bsel, aes(z1, estimate)) +
  geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
  labs(x='Time', y='Estimated covariate effect') +
  scale_colour_grey() +
  theme_bw() 
#ggsave("../drafts/figs/Ecog_estimate_fitmulti.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)


##########################################
## EcoG Data as a function of RPE
## RPE positive (values >0)
## RPE negative (values <=0)
## A figure with mean and 80%CI

rpe_pos <-  pt_gamble_data[,"rpe"]>0
rpe_neg <-  pt_gamble_data[,"rpe"]<=0

rpe_pos_meanhfa <- apply(center_elec_dat[rpe_pos, ], 2, mean)
rpe_neg_meanhfa <- apply(center_elec_dat[rpe_neg, ], 2, mean)
rpe_pos_meanhfa_upper <- rpe_pos_meanhfa + 1.282 * apply(center_elec_dat[rpe_pos, ], 2, sd)/sqrt(ntrials)
rpe_pos_meanhfa_lower <- rpe_pos_meanhfa - 1.282 * apply(center_elec_dat[rpe_pos, ], 2, sd)/sqrt(ntrials)
rpe_neg_meanhfa_upper <- rpe_neg_meanhfa + 1.282 * apply(center_elec_dat[rpe_neg, ], 2, sd)/sqrt(ntrials)
rpe_neg_meanhfa_lower <- rpe_neg_meanhfa - 1.282 * apply(center_elec_dat[rpe_neg, ], 2, sd)/sqrt(ntrials)

textsize= 20  
ggplot()+
  geom_line(data = data.frame(x = timings_w, y = rpe_neg_meanhfa), aes(x, y), 
            color = "black", linetype="dashed", size = 0.75) +
  geom_line(data = data.frame(x = timings_w, y = rpe_pos_meanhfa), aes(x, y), 
            color = "black", linetype="solid", size = 0.75) +
  geom_ribbon(aes(x=timings_w, y=rpe_pos_meanhfa, ymin = rpe_pos_meanhfa_lower, ymax = rpe_pos_meanhfa_upper), 
              fill = "azure4", alpha = 0.3) +
  geom_ribbon(aes(x=timings_w, y=rpe_neg_meanhfa, ymin = rpe_neg_meanhfa_lower, ymax = rpe_neg_meanhfa_upper), 
              fill = "gray68", alpha = 0.3) + #gray68 #lightgray #azure4
  labs(title = "", x = "Time", y = "HFA") +
  guides(shape=FALSE)+
  theme_minimal() +
  theme(axis.text=element_text(size=textsize-8), 
        axis.title=element_text(size=textsize), 
        legend.title=element_text(size=textsize), 
        legend.text=element_text(size=textsize), 
        legend.position=c(.83,.2), legend.key.width=unit(3,'cm'), 
        strip.text.x = element_text(size=textsize))
#ggsave("../drafts/figs/Ecog_rpe_s1.pdf", 
#       width = 6, height = 6, units = "in", dpi = 300)
