#### info ####

# Goal: How do Microstegium and Bipolaris interact to affect native plant biomass?

# Authors: Amy Kendig and Vida Svahnstrom

# Previous version: Svahnstrom_code_20190424.R


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)
library(grid)
library(gridExtra)

# import data
dat <- read_csv("./data/biomass_20190912.csv")
# previous file name: Svahnstrom_greenhouse_biomass_20190912 - Sheet1


#### edit data ####

# notes
unique(dat$notes)
filter(dat, !is.na(notes))

# min values
min(dat$nat_biomass, na.rm = T)
min(dat$mv_biomass, na.rm = T)

# add columns
dat1 <- dat %>%
  mutate(fungus = recode(treatment, "F" = 1, "W" = 0),
         log_nat_biomass = log(nat_biomass),
         log_mv_biomass = log(mv_biomass))


#### visualize ####

# native biomass by density
ggplot(dat1, aes(x = density, y = nat_biomass, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 

# log-transformed native biomass by density
ggplot(dat1, aes(x = density, y = log_nat_biomass, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2) 

# mv biomass by density
ggplot(dat1, aes(x = density, y = mv_biomass, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2)

# native biomass by mv biomass
filter(dat1, !is.na(mv_biomass)) %>%
  ggplot(aes(x = mv_biomass, y = nat_biomass, color = treatment)) +
  geom_point() +
  facet_wrap(~ species)


#### native biomass by density ####

# all species together

# intercept
filter(dat1, species == "C" & density == 0 & fungus == 0) %>%
  summarise(b0 = mean(nat_biomass, na.rm = T),
            log_b0 = mean(log_nat_biomass, na.rm = T))

# distributions
x <- seq(0, 10, length.out = 100)
y <- dexp(x, 2)
y <- dgamma(x, shape = 1.75, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# remove NA's from data
dat2 <- filter(dat1, !is.na(nat_biomass))

# model
nat_mod_1 <- brm(data = dat2, family = gaussian,
                 bf(nat_biomass ~ b0/(1 + alpha * density), 
                    b0 ~ fungus * species, 
                    alpha ~ fungus * species, 
                    nl = T),
                 prior <- c(prior(gamma(2, 1), nlpar = "b0", coef = "Intercept"),
                            prior(normal(0, 1), nlpar = "b0"),
                            prior(gamma(1, 1), nlpar = "alpha", coef = "Intercept"),
                            prior(normal(0, 1), nlpar = "alpha")),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)
plot(nat_mod_1)
# this causes divergent transitions becuase combinations of the intercept and treatments lead to negative estimates
# placing a lower bound on the intercept applies it to the other coefficients as well, which is not ideal

# log-transform biomass
nat_mod_2 <- brm(data = dat2, family = gaussian,
                 bf(log_nat_biomass ~ b0/(1 + alpha * density), 
                    b0 ~ fungus * species, 
                    alpha ~ fungus * species, 
                    nl = T),
                 prior <- c(prior(gamma(0.7, 1), nlpar = "b0", coef = "Intercept"),
                            prior(normal(0, 1), nlpar = "b0"),
                            prior(gamma(1, 1), nlpar = "alpha", coef = "Intercept"),
                            prior(normal(0, 1), nlpar = "alpha")),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)
plot(nat_mod_2)
# same thing happened

# fit each species and treatment separately because you can only have positive effects of fungus/water or species when they are all in the same model


#### native biomass by density - separate models ####

# simulation to check fit
nat_sim_fun <- function(dat, mod){
  
  # predicted data
  pdat <- tibble(density = seq(0, 100, length.out = 300)) %>%
    mutate(nat_biomass = fitted(mod, newdata = .)[, "Estimate"])
  
  # figure
  ggplot(dat, aes(x = density, y = nat_biomass)) +
    stat_summary(geom = "point", fun = "mean", size = 2) +
    stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
    geom_line(data = pdat, color = "pink")
}

# separate data
cfdat <- filter(dat2, species == "C" & fungus == 1)
cwdat <- filter(dat2, species == "C" & fungus == 0)
sfdat <- filter(dat2, species == "S" & fungus == 1)
swdat <- filter(dat2, species == "S" & fungus == 0)
vfdat <- filter(dat2, species == "V" & fungus == 1)
vwdat <- filter(dat2, species == "V" & fungus == 0)

# overall b0
filter(dat2, density == 0) %>%
  summarise(b0 = mean(nat_biomass))
# 1.75

# c fungus nat model
nat_mod_cf <- brm(data = cfdat, family = gaussian,
                  bf(nat_biomass ~ b0/(1 + alpha * density), 
                     b0 + alpha ~ 1, 
                     nl = T),
                  prior <- c(prior(gamma(1.75, 1), nlpar = "b0", lb = 0),
                             prior(gamma(1, 1), nlpar = "alpha", lb = 0)),
                  iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(nat_mod_cf)                 
plot(nat_mod_cf)
nat_sim_fun(cfdat, nat_mod_cf)

# c water nat model
nat_mod_cw <- update(nat_mod_cf, newdata = cwdat)
summary(nat_mod_cw)                 
plot(nat_mod_cw)
nat_sim_fun(cwdat, nat_mod_cw)

# s fungicide nat model
nat_mod_sf <- update(nat_mod_cf, newdata = sfdat)
summary(nat_mod_sf)                 
plot(nat_mod_sf)
nat_sim_fun(sfdat, nat_mod_sf)

# s water nat model
nat_mod_sw <- update(nat_mod_cf, newdata = swdat)
summary(nat_mod_sw)                 
plot(nat_mod_sw)
nat_sim_fun(swdat, nat_mod_sw)

# v fungicide nat model
nat_mod_vf <- update(nat_mod_cf, newdata = vfdat)
summary(nat_mod_vf)                 
plot(nat_mod_vf)
nat_sim_fun(vfdat, nat_mod_vf)

# v water nat model
nat_mod_vw <- update(nat_mod_cf, newdata = vwdat)
summary(nat_mod_vw)                 
plot(nat_mod_vw)
nat_sim_fun(vwdat, nat_mod_vw)


#### mv biomass by density ####

# remove missing data
# reorder species
dat3 <- filter(dat1, !is.na(mv_biomass)) %>%
  mutate(species = fct_relevel(species, "V", "S", "C"))

# model
mv_mod <- brm(data = dat3, family = gaussian,
              mv_biomass ~ density*fungus*species + I(density^2)*fungus*species,
              prior <- c(prior(normal(0, 1), class = "b")),
              iter = 6000, warmup = 1000, chains = 3, cores = 2)
summary(mv_mod)                 
plot(mv_mod)

# simulate data
mv_sim_dat <- tibble(density = seq(0, 100, length.out = 300)) %>%
  merge(tibble(fungus = rep(c(0, 1), 3),
               species = rep(c("C", "S", "V"), each = 2)),
        all = T) %>%
  as_tibble() %>%
  mutate(mv_biomass = fitted(mv_mod, newdata = .)[, "Estimate"],
         mv_biomass_lower = fitted(mv_mod, newdata = .)[, "Q2.5"],
         mv_biomass_upper = fitted(mv_mod, newdata = .)[, "Q97.5"],
         Treatment = recode(fungus, "1" = "pathogen inoculation", "0" = "control (water)")) 

# plot model
ggplot(dat3, aes(x = density, y = mv_biomass)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  geom_line(data = mv_sim_dat, color = "pink") +
  facet_grid(as.factor(fungus) ~ species)


#### save models ####

save(nat_mod_cf, file = "./output/c_biomass_fungus_density_model.rda")
save(nat_mod_cw, file = "./output/c_biomass_water_density_model.rda")
save(nat_mod_sf, file = "./output/s_biomass_fungus_density_model.rda")
save(nat_mod_sw, file = "./output/s_biomass_water_density_model.rda")
save(nat_mod_vf, file = "./output/v_biomass_fungus_density_model.rda")
save(nat_mod_vw, file = "./output/v_biomass_water_density_model.rda")
save(mv_mod, file = "./output/mv_biomass_density_model.rda")


#### native biomass figure ####

# model fits
nat_fit_fun <- function(mod){
  
  dat <- tibble(density = seq(0, 100, length.out = 300)) %>%
    mutate(nat_biomass = fitted(mod, newdata = .)[, "Estimate"],
           nat_biomass_lower = fitted(mod, newdata = .)[, "Q2.5"],
           nat_biomass_upper = fitted(mod, newdata = .)[, "Q97.5"])
  
  return(dat)
}
nat_fit_cf <- nat_fit_fun(nat_mod_cf) %>%
  mutate(treatment = "F")
nat_fit_cw <- nat_fit_fun(nat_mod_cw) %>%
  mutate(treatment = "W")
nat_fit_sf <- nat_fit_fun(nat_mod_sf) %>%
  mutate(treatment = "F")
nat_fit_sw <- nat_fit_fun(nat_mod_sw) %>%
  mutate(treatment = "W")
nat_fit_vf <- nat_fit_fun(nat_mod_vf) %>%
  mutate(treatment = "F")
nat_fit_vw <- nat_fit_fun(nat_mod_vw) %>%
  mutate(treatment = "W")

# combine fits by species
nat_fit_c <- full_join(nat_fit_cf, nat_fit_cw) %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))
nat_fit_s <- full_join(nat_fit_sf, nat_fit_sw) %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))
nat_fit_v <- full_join(nat_fit_vf, nat_fit_vw) %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))

# coefficients
nat_coef_cf <- fixef(nat_mod_cf) %>%
  as_tibble() %>%
  mutate(Parameter = c("b0", "alpha"),
         treatment = "F")
nat_coef_cw <- fixef(nat_mod_cw) %>%
  as_tibble() %>%
  mutate(Parameter = c("b0", "alpha"),
         treatment = "W")
nat_coef_sf <- fixef(nat_mod_sf) %>%
  as_tibble() %>%
  mutate(Parameter = c("b0", "alpha"),
         treatment = "F")
nat_coef_sw <- fixef(nat_mod_sw) %>%
  as_tibble() %>%
  mutate(Parameter = c("b0", "alpha"),
         treatment = "W")
nat_coef_vf <- fixef(nat_mod_vf) %>%
  as_tibble() %>%
  mutate(Parameter = c("b0", "alpha"),
         treatment = "F")
nat_coef_vw <- fixef(nat_mod_vw) %>%
  as_tibble() %>%
  mutate(Parameter = c("b0", "alpha"),
         treatment = "W")

# combine coefficients by species
nat_coef_c <- full_join(nat_coef_cf, nat_coef_cw) %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"),
         Disease = recode(Parameter, alpha = "Indirect", b0 = "Direct"))
nat_coef_s <- full_join(nat_coef_sf, nat_coef_sw) %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"),
         Disease = recode(Parameter, alpha = "Indirect", b0 = "Direct"))
nat_coef_v <- full_join(nat_coef_vf, nat_coef_vw) %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"),
         Disease = recode(Parameter, alpha = "Indirect", b0 = "Direct"))

# combine data by species
cdat <- filter(dat2, species == "C") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))
sdat <- filter(dat2, species == "S") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))
vdat <- filter(dat2, species == "V") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))

# default theme
theme_def <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.title = element_text(size = 12, face = "italic", hjust = 0.5))

# colors
col_pal = c("#018571", "#a6611a")

# parameter plots
par_plot_c <- ggplot(nat_coef_c, aes(x = Disease, y = Estimate, color = Treatment)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), position = position_dodge(0.2), width = 0.1) +
  geom_point(size = 2, position = position_dodge(0.2))+
  scale_color_manual(values = col_pal) +
  theme_def +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10))

par_plot_s <- par_plot_c %+%
  nat_coef_s

par_plot_v <- par_plot_c %+%
  nat_coef_v

# density plots
dens_plot_c <- ggplot(cdat, aes(x = density, y = nat_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = nat_fit_c, aes(ymin = nat_biomass_lower, ymax = nat_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = nat_fit_c, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  annotation_custom(ggplotGrob(par_plot_c), xmin = 20, xmax = 95, ymin = 1, ymax = 4.1) +
  ggtitle("Panicum clandestinum") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def
  
dens_plot_s <- ggplot(sdat, aes(x = density, y = nat_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = nat_fit_s, aes(ymin = nat_biomass_lower, ymax = nat_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = nat_fit_s, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  annotation_custom(ggplotGrob(par_plot_s), xmin = 20, xmax = 95, ymin = 0.32, ymax = 1.33) +
  ggtitle("Eragrostis spectabilis") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_y_continuous(breaks = c(0, 1)) +
  theme_def 

dens_plot_v <- ggplot(vdat, aes(x = density, y = nat_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = nat_fit_v, aes(ymin = nat_biomass_lower, ymax = nat_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = nat_fit_v, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  annotation_custom(ggplotGrob(par_plot_v), xmin = 20, xmax = 95, ymin = 0.56, ymax = 2.3) +
  ggtitle("Elymus virginicus") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_y_continuous(breaks = c(0, 1, 2)) +
  theme_def
  
# combine plots
dens_plot_nat <- plot_grid(dens_plot_v + theme(legend.position = "none"),
                           dens_plot_s + theme(legend.position = "none"),
                           dens_plot_c + theme(legend.position = "none"),
                           nrow = 1,
                           labels = LETTERS[1:3],
                           label_size = 12)

# axes
y_plot_nat <- textGrob("Biomass (g)", gp = gpar(fontsize = 12), rot = 90)
x_plot_nat <- textGrob(expression(paste(italic(Microstegium), " density", sep = "")), gp = gpar(fontsize = 12))

# combine
dens_plot_nat2 <- grid.arrange(arrangeGrob(dens_plot_nat, bottom = x_plot_nat, left = y_plot_nat))

# legend
leg_nat <- get_legend(dens_plot_v)

# save plot
tiff("./output/native_biomass_density_figure.tiff", width = 7.5, height = 3, units = "in", res = 300)
grid.arrange(arrangeGrob(dens_plot_nat2, bottom = leg_nat, padding = unit(1, "line")))
dev.off()


#### native biomass table ####

# combine and select
nat_coef <- nat_coef_c %>%
  mutate(Species = "Panicum clandestinum") %>%
  full_join(nat_coef_s %>%
              mutate(Species = "Eragrostis spectabilis")) %>%
  full_join(nat_coef_v %>%
              mutate(Species = "Elymus virginicus")) %>%
  arrange(Species, Treatment, Parameter) %>%
  select(Species, Treatment, Parameter, Estimate:Q97.5) %>%
  mutate(Estimate = round(Estimate, 2),
         Est.Error = round(Est.Error, 2),
         Q2.5 = round(Q2.5, 2),
         Q97.5 = round(Q97.5, 2))

# save
write_csv(nat_coef, "./output/native_biomass_density_parameters.csv")

#### Microstegium biomass figure ####

# separate data
mvdatc <- filter(dat3, species == "C") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))
mvdats <- filter(dat3, species == "S") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))
mvdatv <- filter(dat3, species == "V") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"))

mv_sim_datc <- filter(mv_sim_dat, species == "C")
mv_sim_dats <- filter(mv_sim_dat, species == "S")
mv_sim_datv <- filter(mv_sim_dat, species == "V")

# summarise data to specify y limits
mvdatc_sum <- mvdatc %>%
  group_by(density, Treatment) %>%
  summarise(mv_biomass_lower = mean_cl_boot(mv_biomass)[, 2],
            mv_biomass_upper = mean_cl_boot(mv_biomass)[, 3],
            mv_biomass = mean_cl_boot(mv_biomass)[, 1])

mvdatv_sum <- mvdatv %>%
  group_by(density, Treatment) %>%
  summarise(mv_biomass_lower = mean_cl_boot(mv_biomass)[, 2],
            mv_biomass_upper = mean_cl_boot(mv_biomass)[, 3],
            mv_biomass = mean_cl_boot(mv_biomass)[, 1])

# mv plots
mv_plot_c <- ggplot(mvdatc_sum, aes(x = density, y = mv_biomass, 
                                    ymin = mv_biomass_lower, ymax = mv_biomass_upper,
                                    color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = mv_sim_datc, alpha = 0.3, color = NA) +
  geom_line(data = mv_sim_datc, aes(linetype = Treatment)) +
  geom_errorbar(width = 0.1, alpha = 0.5) +
  geom_point(size = 2) +
  ggtitle("Panicum clandestinum") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 12))

mv_plot_s <- ggplot(mvdats, aes(x = density, y = mv_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = mv_sim_dats, aes(ymin = mv_biomass_lower, ymax = mv_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = mv_sim_dats, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  ggtitle("Eragrostis spectabilis") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 12))

mv_plot_v <- ggplot(mvdatv_sum, aes(x = density, y = mv_biomass, 
                                    ymin = mv_biomass_lower, ymax = mv_biomass_upper,
                                    color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = mv_sim_datv, alpha = 0.3, color = NA) +
  geom_line(data = mv_sim_datv, aes(linetype = Treatment)) +
  geom_errorbar(width = 0.1, alpha = 0.5) +
  geom_point(size = 2) +
  ggtitle("Elymus virginicus") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 12))

# combine plots
dens_plot_mv <- plot_grid(mv_plot_v + theme(legend.position = "none"),
                           mv_plot_s + theme(legend.position = "none"),
                           mv_plot_c + theme(legend.position = "none"),
                           nrow = 1,
                           labels = LETTERS[1:3],
                           label_size = 12)

# axes
y_plot_mv <- textGrob(expression(paste(italic(Microstegium), " biomass (g)", sep = "")), gp = gpar(fontsize = 12), rot = 90)
x_plot_mv <- x_plot_nat

# combine
dens_plot_mv2 <- grid.arrange(arrangeGrob(dens_plot_mv, bottom = x_plot_mv, left = y_plot_mv))

# save plot
tiff("./output/mv_biomass_density_figure.tiff", width = 7.5, height = 3, units = "in", res = 300)
grid.arrange(arrangeGrob(dens_plot_mv2, bottom = leg_nat, padding = unit(1, "line")))
dev.off()


#### Microstegium biomass table ####

# extract coefficients
mv_coef <- fixef(mv_mod) %>%
  as_tibble(rownames = "Coefficient")

# save
write_csv(mv_coef, "./output/mv_biomass_density_coefficients.csv")
