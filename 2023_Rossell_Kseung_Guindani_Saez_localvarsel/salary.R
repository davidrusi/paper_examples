#These data were obtained from the Current Population Survey at https://www.census.gov/programs-surveys/cps.html
#
#To use these data you should add a reference to:
#
#Flood, Sarah and King, Miriam and Rodgers, Renae and Ruggles, Steven and Warren, J Robert (2020).
#Integrated Public Use Microdata Series, Current Population Survey: Version 7.0 [dataset]
#Minneapolis, MN: IPUMS, 10.18128/D030.V7.0
#
# See https://cps.ipums.org/cps/citation.shtml for the rules regarding citation and use of the CPS database

library(mombf)
library(parallel)
library(tidyverse)
setwd("~/github/localnulltest/code")
#source('routines.R')
load("../data/salary.RData")

#SELECT CASES. AGED 18-65, SINGLE RACE, NON-MILITARY EMPLOYED 35-40H/WEEK
data= filter(salary, incomewage>=500, age>=18, age<=65, employment=="At work", !(occ %in% c("unemployed/neverworked","military")), hoursworked>=35, hoursworked<=40, wkstat=="Full-time") |>
  mutate(occ= factor(occ), logincome= log10(incomewage), black= ifelse(race=='Black',1,0), college= ifelse(edu=='CG',1,0), classworker=factor(classworker)) |>
  mutate(female_black= female * black) |>
  select(c(logincome, age, female, hispanic, black, college, occ, classworker, female_black))


y= data$logincome
z= data$age
x.adjust= model.matrix(~ occ, data=data)
x= data[,c('female','hispanic','black','college','classworker')]
x= mutate(x, government= as.numeric(classworker=='Government employee'), selfemployed= as.numeric(classworker=='Self-employed')) |>
  select(-classworker)

## Analysis without interactions ## 
###################################

fit= localnulltest(y, x=x, z=z, x.adjust=x.adjust, nlocalknots=c(5,10,15,20), verbose=TRUE, niter=2000, mc.cores=1)
b= coef(fit)
save(fit, b, file="salary_results.RData")

b= data.frame(covariate=1:ncol(x), covariaten= names(x)) |>
     right_join(b) |>
     select(-covariate) |>
     rename(covariate= covariaten)

textsize= 35
mycols= rep(c('black','darkgrey'), ncol(x)/2 + 1)[1:ncol(x)]
ggplot(b, aes(z1, margpp)) +
    geom_line() +
    #geom_line(aes(group=covariate, lty=covariate)) +
    labs(x='Age (years)', y='Posterior probability of a covariate effect') +
    facet_wrap( ~ covariate) +
    ylim(0,1) + 
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.83,.2), legend.key.width=unit(5,'cm'), strip.text.x = element_text(size=textsize))
ggsave("../drafts/figs/salary_pp_localtests.pdf")

textsize= 35
bsel= filter(b, covariate %in% c('black','female','hispanic','college'))

ggplot(bsel, aes(z1, estimate)) +
    geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
    labs(x='Age (years)', y='Estimated covariate effect (log10)') +
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.83,.7), legend.key.width=unit(5,'cm'))

ggsave("../drafts/figs/salary_estimated_effects.pdf")



## SET 20 KNOTS ##

fit2= localnulltest(y, x=x, z=z, x.adjust=x.adjust, nbaseknots=20, nlocalknots=20, verbose=TRUE)  #results in low power
b2= coef(fit2)

b2= data.frame(covariate=1:ncol(x), covariaten= names(x)) |>
     right_join(b2) |>
     select(-covariate) |>
     rename(covariate= covariaten)

textsize= 35
mycols= rep(c('black','darkgrey'), ncol(x)/2 + 1)[1:ncol(x)]
ggplot(b2, aes(z1, margpp)) +
    geom_line() +
    #geom_line(aes(group=covariate, lty=covariate)) +
    labs(x='Age (years)', y='Posterior probability of a covariate effect') +
    facet_wrap( ~ covariate) +
    ylim(0,1) + 
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.83,.2), legend.key.width=unit(5,'cm'), strip.text.x = element_text(size=textsize))
ggsave("../drafts/figs/salary_pp_localtests_20knots.pdf")

textsize= 35
bsel= filter(b2, covariate %in% c('black','female','hispanic','college'))
ggplot(bsel, aes(z1, estimate)) +
    geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
    labs(x='Age (years)', y='Estimated covariate effect (log10)') +
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.83,.7), legend.key.width=unit(5,'cm'))
ggsave("../drafts/figs/salary_estimated_effects_20knots.pdf")



## SET 30 KNOTS ##

fit2= localnulltest(y, x=x, z=z, x.adjust=x.adjust, nbaseknots=30, nlocalknots=30, verbose=TRUE)  #results in low power
b2= coef(fit2)

b2= data.frame(covariate=1:ncol(x), covariaten= names(x)) |>
     right_join(b2) |>
     select(-covariate) |>
     rename(covariate= covariaten)

textsize= 35
mycols= rep(c('black','darkgrey'), ncol(x)/2 + 1)[1:ncol(x)]
ggplot(b2, aes(z1, margpp)) +
    geom_line() +
    #geom_line(aes(group=covariate, lty=covariate)) +
    labs(x='Age (years)', y='Posterior probability of a covariate effect') +
    facet_wrap( ~ covariate) +
    ylim(0,1) + 
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.83,.2), legend.key.width=unit(5,'cm'), strip.text.x = element_text(size=textsize))
ggsave("../drafts/figs/salary_pp_localtests_30knots.pdf")

textsize= 35
bsel= filter(b2, covariate %in% c('black','female','hispanic','college'))
ggplot(bsel, aes(z1, estimate)) +
    geom_line(aes(group=covariate, col=covariate, lty=covariate), lwd=2) +
    labs(x='Age (years)', y='Estimated covariate effect (log10)') +
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.position=c(.83,.7), legend.key.width=unit(5,'cm'))
ggsave("../drafts/figs/salary_estimated_effects_30knots.pdf")



## ANALYSIS WITH INTERACTIONS ##
################################

#Add interaction female:black
x= data[,c('female','hispanic','black','college','classworker','female_black')]
x= mutate(x, government= as.numeric(classworker=='Government employee'), selfemployed= as.numeric(classworker=='Self-employed')) |>
    select(-classworker)

fit= localnulltest(y, x=x, z=z, x.adjust=x.adjust, nlocalknots=c(5,10,15,20), verbose=TRUE, niter=2000, mc.cores=1)
b= coef(fit)

save(fit, b, file="salary_results_femaleblack.RData")

b= data.frame(covariate=1:ncol(x), covariaten= names(x)) |>
     right_join(b) |>
     select(-covariate) |>
     rename(covariate= covariaten) |>
     filter(covariate %in% c('female','black','female_black'))

textsize= 35
mycols= rep(c('black','darkgrey'), ncol(x)/2 + 1)[1:ncol(x)]
ggplot(b, aes(z1, margpp)) +
    geom_line() +
    #geom_line(aes(group=covariate, lty=covariate)) +
    labs(x='Age (years)', y='Posterior probability of a covariate effect') +
    facet_wrap( ~ covariate) +
    ylim(0,1) + 
    scale_colour_grey() +
    theme_bw() +
    theme(axis.text=element_text(size=textsize), axis.title=element_text(size=textsize), legend.title=element_text(size=textsize), legend.text=element_text(size=textsize), legend.key.width=unit(5,'cm'), strip.text.x = element_text(size=textsize))

ggsave("../drafts/figs/salary_pp_localtests_femaleblack.pdf")
