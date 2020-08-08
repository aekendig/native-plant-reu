#### info ####

# Goal: How do Microstegium and Bipolaris interact to affect native plant biomass?

# Authors: Amy Kendig and Vida Svahnstrom


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse) # version 1.3.0
library(brms) # version 2.13.0
library(cowplot) # version 1.0.0
library(grid) # version 4.0.1
library(gridExtra) # version 2.3
library(tidybayes) # version 2.0.3
library(extrafont) # version 0.17

# import data
dat <- read_csv("./data/biomass_20190912.csv")


#### edit data ####

# notes
unique(dat$notes)
filter(dat, !is.na(notes))

# min values
min(dat$nat_biomass, na.rm = T)
min(dat$mv_biomass, na.rm = T)

# add columns
# remove pot that was accidentally inoculated
dat1 <- dat %>%
  mutate(fungus = recode(treatment, "F" = 1, "W" = 0),
         fungusF = recode(treatment, "F" = "inoculation", "W" = "control"),
         log_nat_biomass = log(nat_biomass),
         log_mv_biomass = log(mv_biomass)) %>%
  filter(!(treatment == "W" & density == 100 & species == "C" & replicate == 1))


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
y <- dexp(x, 0.5)
y <- dgamma(x, shape = 2, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# remove NA's from data
dat2 <- filter(dat1, !is.na(nat_biomass))

# model
# nat_mod_1 <- brm(data = dat2, family = gaussian,
#                  bf(nat_biomass ~ b0/(1 + alpha * density), 
#                     b0 ~ 0 + fungusF * species, 
#                     alpha ~ 0 + fungusF * species, 
#                     nl = T),
#                  prior <- c(prior(gamma(2, 1), nlpar = "b0", lb = 0),
#                             prior(exponential(0.5), nlpar = "alpha", lb = 0),
#                             prior(cauchy(0, 1), class = sigma)),
#                  iter = 6000, warmup = 1000, chains = 1, cores = 1)
# summary(nat_mod_1)
# plot(nat_mod_1)
# prior_summary(nat_mod_1)
# this isn't quite right - the effects of species identity are constrained to be positive

# combine species and fungus treatments into one variable
dat3 <- dat2 %>%
  mutate(spfungus = paste(species, fungusF, sep = ""),
         Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "mock inoculation (control)"))

# model
nat_mod_2 <- brm(data = dat3, family = gaussian,
                 bf(nat_biomass ~ b0/(1 + alpha * density), 
                    b0 ~ 0 + spfungus, 
                    alpha ~ 0 + spfungus, 
                    nl = T),
                 prior <- c(prior(gamma(2, 1), nlpar = "b0", lb = 0),
                            prior(exponential(0.5), nlpar = "alpha", lb = 0),
                            prior(cauchy(0, 1), class = sigma)),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)
prior_summary(nat_mod_2)
summary(nat_mod_2)
plot(nat_mod_2)

# increase the chains and cores
nat_mod_3 <- update(nat_mod_2, chains = 3, cores = 1)
summary(nat_mod_3)
plot(nat_mod_3)
pp_check(nat_mod_3, nsamples = 100)

# simulate data
nat_sim_dat <- tibble(density = seq(0, 100, length.out = 300)) %>%
  merge(tibble(fungusF = rep(c("control", "inoculation"), 3),
               species = rep(c("C", "S", "V"), each = 2)),
        all = T) %>%
  as_tibble() %>%
  mutate(spfungus = paste(species, fungusF, sep = "")) %>%
  mutate(nat_biomass = fitted(nat_mod_3, newdata = .)[, "Estimate"],
         nat_biomass_lower = fitted(nat_mod_3, newdata = .)[, "Q2.5"],
         nat_biomass_upper = fitted(nat_mod_3, newdata = .)[, "Q97.5"],
         Treatment = recode(fungusF, "inoculation" = "pathogen inoculation", "control" = "mock inoculation (control)")) 

# plot model
ggplot(dat3, aes(x = density, y = nat_biomass)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  geom_line(data = nat_sim_dat, color = "pink") +
  facet_grid(fungusF ~ species)


#### mv biomass by density ####

# remove missing data
# reorder species
dat4 <- filter(dat1, !is.na(mv_biomass)) %>%
  mutate(species = fct_relevel(species, "C", "V", "S"))

# intercept prior
filter(dat4, species == "C" & fungus == 0 & density == 2) %>%
  summarise(mean = mean(mv_biomass),
            sd = sd(mv_biomass))

# model
mv_mod <- brm(data = dat4, family = gaussian,
              mv_biomass ~ density*fungus*species + I(density^2)*fungus*species,
              prior <- c(prior(normal(0, 10), class = "b"),
                         prior(normal(2, 10), class = "Intercept"),
                         prior(cauchy(0, 1), class = sigma)),
              iter = 6000, warmup = 1000, chains = 3, cores = 1, 
              control = list(max_treedepth = 15))
prior_summary(mv_mod)
summary(mv_mod)                 
plot(mv_mod)
pp_check(mv_mod, nsamples = 100)

# simulate data
mv_sim_dat <- tibble(density = seq(0, 100, length.out = 300)) %>%
  merge(tibble(fungus = rep(c(0, 1), 3),
               species = rep(c("C", "S", "V"), each = 2)),
        all = T) %>%
  as_tibble() %>%
  mutate(mv_biomass = fitted(mv_mod, newdata = .)[, "Estimate"],
         mv_biomass_lower = fitted(mv_mod, newdata = .)[, "Q2.5"],
         mv_biomass_upper = fitted(mv_mod, newdata = .)[, "Q97.5"],
         Treatment = recode(fungus, "1" = "pathogen inoculation", "0" = "mock inoculation (control)")) 

# plot model
ggplot(dat4, aes(x = density, y = mv_biomass)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  geom_line(data = mv_sim_dat, color = "pink") +
  facet_grid(as.factor(fungus) ~ species)


#### native biomass figure ####

# coefficients
nat_coef <- fixef(nat_mod_3) %>%
  as_tibble() %>%
  mutate(Parameter = rep(c("b0", "alpha"), each = 6),
         Treatment = rep(c("mock inoculation (control)", "pathogen inoculation"), 6),
         species = rep(rep(c("C", "S", "V"), each = 2), 2)) %>%
  mutate(Disease = recode(Parameter, "b0" = "Direct", "alpha" = "Indirect"))

# split by species
nat_coef_c <- filter(nat_coef, species == "C")
nat_coef_v <- filter(nat_coef, species == "V")
nat_coef_s <- filter(nat_coef, species == "S")

# split data by species
datc <- filter(dat3, species == "C")
datv <- filter(dat3, species == "V")
dats <- filter(dat3, species == "S")

# split fit by species
nat_sim_datc <- filter(nat_sim_dat, species == "C")
nat_sim_datv <- filter(nat_sim_dat, species == "V")
nat_sim_dats <- filter(nat_sim_dat, species == "S")

# default theme
theme_def <- theme_bw(base_family = "Arial") +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.title = element_text(size = 10, face = "italic", hjust = 0.5))

# colors
col_pal = c("#55A48B", "#C0A76D")

# parameter plots
par_plot_c <- ggplot(nat_coef_c, aes(x = Disease, y = Estimate, color = Treatment)) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), position = position_dodge(0.5), width = 0) +
  geom_point(size = 2, position = position_dodge(0.5))+
  scale_color_manual(values = col_pal) +
  scale_y_continuous(limits = c(0, 8)) +
  theme_def +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8))

par_plot_s <- par_plot_c %+%
  nat_coef_s

par_plot_v <- par_plot_c %+%
  nat_coef_v

# density plots
dens_plot_c <- ggplot(datc, aes(x = density, y = nat_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = nat_sim_datc, aes(ymin = nat_biomass_lower, ymax = nat_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = nat_sim_datc, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  annotation_custom(ggplotGrob(par_plot_c), xmin = 10, xmax = 105, ymin = 1, ymax = 4.2) +
  ggtitle("Dichanthelium") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_y_continuous(limits = c(0, 4.2)) +
  theme_def
  
dens_plot_s <- ggplot(dats, aes(x = density, y = nat_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = nat_sim_dats, aes(ymin = nat_biomass_lower, ymax = nat_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = nat_sim_dats, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  annotation_custom(ggplotGrob(par_plot_s), xmin = 10, xmax = 105, ymin = 1, ymax = 4.2) +
  ggtitle("Eragrostis") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_y_continuous(limits = c(0, 4.2)) +
  theme_def 

dens_plot_v <- ggplot(datv, aes(x = density, y = nat_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = nat_sim_datv, aes(ymin = nat_biomass_lower, ymax = nat_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = nat_sim_datv, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  annotation_custom(ggplotGrob(par_plot_v), xmin = 10, xmax = 105, ymin = 1, ymax = 4.2) +
  ggtitle("Elymus") +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  scale_y_continuous(limits = c(0, 4.2)) +
  theme_def
  
# combine plots
dens_plot_nat <- plot_grid(dens_plot_c + theme(legend.position = "none"),
                           dens_plot_v + theme(legend.position = "none"),
                           dens_plot_s + theme(legend.position = "none"),
                           nrow = 1,
                           labels = LETTERS[1:3],
                           label_size = 10)

# axes
y_plot_nat <- textGrob("Biomass (g)", gp = gpar(fontsize = 10), rot = 90)
x_plot_nat <- textGrob(expression(paste(italic(Microstegium), " density", sep = "")), gp = gpar(fontsize = 10))

# combine
dens_plot_nat2 <- grid.arrange(arrangeGrob(dens_plot_nat, bottom = x_plot_nat, left = y_plot_nat))

# legend
leg_nat <- get_legend(dens_plot_v)

# save plot
tiff("./output/Fig5.tiff", width = 5.2, height = 2.5, units = "in", res = 300)
grid.arrange(arrangeGrob(dens_plot_nat2, bottom = leg_nat, padding = unit(1, "line")))
dev.off()


#### native biomass table ####

# extract posterior samples
nat_post <- posterior_samples(nat_mod_3)

# calculate differences
nat_diff <- nat_post %>%
  transmute(dfv = b_b0_spfungusVinoculation - b_b0_spfungusVcontrol,
            dfs = b_b0_spfungusSinoculation - b_b0_spfungusScontrol,
            dfc = b_b0_spfungusCinoculation - b_b0_spfungusCcontrol,
            ifv = b_alpha_spfungusVinoculation - b_alpha_spfungusVcontrol,
            ifs = b_alpha_spfungusSinoculation - b_alpha_spfungusScontrol,
            ifc = b_alpha_spfungusCinoculation - b_alpha_spfungusCcontrol,
            dvc = b_b0_spfungusVcontrol - b_b0_spfungusCcontrol,
            dsc = b_b0_spfungusScontrol - b_b0_spfungusCcontrol,
            dev = b_b0_spfungusScontrol - b_b0_spfungusVcontrol,
            ivc = b_alpha_spfungusVcontrol - b_alpha_spfungusCcontrol,
            isc = b_alpha_spfungusScontrol - b_alpha_spfungusCcontrol,
            iev = b_alpha_spfungusScontrol - b_alpha_spfungusVcontrol) %>%
  gather() %>%
  mutate(Disease = case_when(substr(key, 1, 1) == "d" ~ "Direct",
                             TRUE ~ "Indirect"),
         Comparison = case_when(substr(key, 2, 2) == "f" ~ "Pathogen inoculation - mock inoculation",
                                substr(key, 2, 2) == "s" ~ "Eragrostis - Dichanthelium",
                                substr(key, 2, 2) == "v" ~ "Elymus - Dichanthelium",
                                substr(key, 2, 2) == "e" ~ "Eragrostis - Elymus"),
         Species = case_when(substr(key, 3, 3) == "v"  ~ "Elymus virginicus",
                             substr(key, 3, 3) == "s"  ~ "Eragrostis spectabilis",
                             substr(key, 3, 3) == "c"  ~ "Dichanthelium clandestinum")) %>%
  select(-key)
head(nat_diff)

# Mean values
nat_table <- nat_diff %>%
  group_by(Comparison, Species, Disease) %>%
  mean_hdi() %>%
  ungroup() %>%
  rename("Estimate" = value, "Q2.5" = .lower, "Q97.5" = .upper)

# save
write_csv(nat_table, "./output/native_biomass_density_table.csv")


#### mv biomass figure ####

# separate data
mvdatc <- filter(dat4, species == "C") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "mock inoculation (control)"))
mvdats <- filter(dat4, species == "S") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "mock inoculation (control)"))
mvdatv <- filter(dat4, species == "V") %>%
  mutate(Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "mock inoculation (control)"))

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
  ggtitle(expression(paste("Native: ", italic(Dichanthelium), sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def +
  theme(plot.title = element_text(size = 9.5, hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 14))

mv_plot_s <- ggplot(mvdats, aes(x = density, y = mv_biomass, color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = mv_sim_dats, aes(ymin = mv_biomass_lower, ymax = mv_biomass_upper), alpha = 0.3, color = NA) +
  geom_line(data = mv_sim_dats, aes(linetype = Treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  ggtitle(expression(paste("Native: ", italic(Eragrostis), sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def +
  theme(plot.title = element_text(size = 9.5, hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 14))

mv_plot_v <- ggplot(mvdatv_sum, aes(x = density, y = mv_biomass, 
                                    ymin = mv_biomass_lower, ymax = mv_biomass_upper,
                                    color = Treatment, fill = Treatment)) + 
  geom_ribbon(data = mv_sim_datv, alpha = 0.3, color = NA) +
  geom_line(data = mv_sim_datv, aes(linetype = Treatment)) +
  geom_errorbar(width = 0.1, alpha = 0.5) +
  geom_point(size = 2) +
  ggtitle(expression(paste("Native: ", italic(Elymus), sep = ""))) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  theme_def +
  theme(plot.title = element_text(size = 9.5, hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 14))

# combine plots
dens_plot_mv <- plot_grid(mv_plot_c + theme(legend.position = "none"),
                          mv_plot_v + theme(legend.position = "none"),
                          mv_plot_s + theme(legend.position = "none"),
                          nrow = 1,
                          labels = LETTERS[1:3],
                          label_size = 10)

# axes
y_plot_mv <- textGrob(expression(paste(italic(Microstegium), " biomass (g)", sep = "")), gp = gpar(fontsize = 10), rot = 90)
x_plot_mv <- x_plot_nat

# combine
dens_plot_mv2 <- grid.arrange(arrangeGrob(dens_plot_mv, bottom = x_plot_mv, left = y_plot_mv))

# legend
leg_mv <- get_legend(mv_plot_v)

# save plot
tiff("./output/Fig4.tiff", width = 5.2, height = 2.5, units = "in", res = 300)
grid.arrange(arrangeGrob(dens_plot_mv2, bottom = leg_mv, padding = unit(1, "line")))
dev.off()


#### mv biomass table ####

# extract coefficients
mv_coef <- fixef(mv_mod) %>%
  as_tibble(rownames = "Coefficient")

# save
write_csv(mv_coef, "./output/mv_biomass_density_coefficients.csv")


#### save models ####

save(nat_mod_3, file = "./output/nat_biomass_density_model.rda")
save(mv_mod, file = "./output/mv_biomass_density_model.rda")