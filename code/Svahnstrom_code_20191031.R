#### REU 2019 data analysis ####
#Authors: Vida Svahnstrom
#Goal: 

#Clear existing data 
rm(list=ls())
#set working directory 
setwd("~/Desktop/REU/R Stuff")
#load library
library("tidyverse")

#import data set with biomass results
dat1 <- read.csv("Svahnstrom_greenhouse_biomass_20190912 - Sheet1.csv")
#rename columns
colnames(dat1) <- c("treatment","density","species",
                    "replicate","mv_biomass","nat_biomass", "notes"
)
#import data set with infection results
dat2 <- read.csv("Infection quantification data - Sheet1.csv")
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
#check this statistically
m1 <- lm(nat_biomass ~ density*treatment*species, data=dat1)
summary(m1)

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
