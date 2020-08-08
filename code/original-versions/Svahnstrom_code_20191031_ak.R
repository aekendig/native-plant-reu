#### REU 2019 data analysis ####
#Authors: Vida Svahnstrom
#Goal: 

#Clear existing data 
rm(list=ls())
#set working directory 
#setwd("~/Desktop/REU/R Stuff")
#load library
library("tidyverse")
library(brms)

#import data set with biomass results
dat1 <- read.csv("./data/Svahnstrom_greenhouse_biomass_20190912 - Sheet1.csv")
#rename columns
#colnames(dat1) <- c("treatment","density","species",
#                    "replicate","mv_biomass","nat_biomass", "notes"
#) #I don't think this is necessary because these are the column names already 
#import data set with infection results
dat2 <- read.csv("./data/Infection quantification data - Sheet1.csv")
#rename columns
colnames(dat2) <- c("treatment","density","species",
                    "replicate","mv_leaves_inf","mv_leaves_1","mv_leaves_2",
                    "mv_leaves_3", "nat_leaves_inf","nat_leaves","Notes"
)

#Check if the relationship between native biomass and density varies by treatment or sp
ggplot(dat1,aes(x=as.factor(density), y=nat_biomass, color=treatment)) +
  geom_point(alpha=.1) +
  facet_wrap(~species) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 

# density as an integer
ggplot(dat1,aes(x=density, y=nat_biomass, color=treatment)) +
  geom_point(alpha=.1) +
  facet_wrap(~species) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 

#check this statistically
m1 <- lm(nat_biomass ~ density*treatment*species, data=dat1)
summary(m1)

# estimate average max biomass
dat1 %>%
  filter(density == 0) %>%
  summarise(mean_bio = mean(nat_biomass),
            sd_bio = sd(nat_biomass))

# estimate competitive effect
a <- 1
dat1 %>%
  select(density) %>%
  unique() %>%
  mutate(nat_bio = 1.8/(1 + a*density)) %>%
  ggplot(aes(x = density, y = nat_bio)) +
  geom_line()

#### work on priors
# non-linear model
m1n <- brm(data = dat1, family = gaussian,
           bf(nat_biomass ~ b0/(1 + alpha * density), b0 ~ treatment * species, alpha ~ treatment * species, nl = T),
           prior <- c(prior(normal(1.8, 10), nlpar = "b0"),
                      prior(gamma(1, 10), nlpar = "alpha")),
           iter = 6000, warmup = 1000, chains = 1, cores = 1,
           control = list(adapt_delta = 0.99999))

# check out model
summary(m1n) # all Rhat values are 1
plot(m1n) # check that lines are generally overlapping
pp_check(m1n, nsamples = 50) # check that samples line up with mean
dat1_comp <- dat1 %>%
  filter(!is.na(nat_biomass)) %>%
  mutate(predicted = predict(m1n)[, 1])
cor.test(dat1_comp$nat_biomass, dat1_comp$predicted) # high correlation between predicted and observed values

# create a dataset to look at predicted values
dat1_pred <- tibble(density = rep(c(0:100), 6), 
                    species = rep(c("C", "S", "V"), each = 202),
                    treatment = rep(rep(c("W", "F"), each = 101), 3))


# ad predicted values to data
dat1_pred <- dat1_pred %>%
  cbind(fitted(m1n, newdata = dat1_pred, nsamples = 100)) %>%
  rename(pred = Estimate,
         lower = Q2.5,
         upper = Q97.5) %>%
  as_tibble()






#NEW PLOT I made with Amy
ggplot(dat3,aes(mv_inf_prop, y=mv_biomass, color=treatment)) +
  geom_point(alpha=.5) +
  facet_wrap(~species) 

#TO DO:
#make lm for above with only fungus data
#do one with y=mv_inf_prop and predictor as density*treatment*species

#Merge dataframes 
dat3 <- merge(dat1, dat2, all = T)

#Check if infection of natives effects native biomass
ggplot(dat3,aes(x=nat_leaves_inf, y=nat_biomass)) +
  geom_point(alpha=.1) +
  facet_wrap(~species) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 
ma <- lm(nat_biomass ~ nat_leaves_inf*species, data=dat3)
summary(ma)


dat3 <- dat3 %>% 
  mutate(sum_leaves_mv = rowSums(cbind(mv_leaves_1, mv_leaves_2, mv_leaves_3), na.rm = T),
         av_leaves_mv = sum_leaves_mv/3,
         av_leaves_mv_pot = av_leaves_mv*density,
         mv_inf_prop = mv_leaves_inf/av_leaves_mv_pot,
         mv_inf_factor = as.factor(ifelse(mv_leaves_inf == 0, "no", "yes")),
         nat_inf_factor = as.factor(ifelse(nat_leaves_inf == 0, "no","yes")))

#Check if infection of natives (as a factor) effects native biomass                                 
ggplot(dat3,aes(nat_inf_factor, y=nat_biomass)) +
  geom_point(alpha=.1) +
  facet_wrap(~species) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 

ma <- lm(nat_biomass ~ nat_inf_factor*species, data=dat3)
summary(ma)

#Check if proportion of Mv leaves infection effects native biomass
ggplot(dat3,aes(x=mv_inf_prop, y=nat_biomass)) +
  geom_point(alpha=.5) +
  facet_wrap(~species, scales="free") 
   

m2 <- lm(nat_biomass ~ mv_inf_prop*species, data=dat3)
summary(m2)

#Check if Mv infection (as a factor) effects native biomass
ggplot(dat3,aes(x=mv_inf_factor, y=nat_biomass)) +
  geom_point(alpha=.1) +
  facet_wrap(~species) +
  stat_summary(geom = "point", fun.y = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 

m3 <- lm(nat_biomass ~ as.factor(mv_inf_factor)*species, data=dat3)
summary(m3)
