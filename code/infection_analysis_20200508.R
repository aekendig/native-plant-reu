#### info ####

# Goal: How does Microstegium density affect infection?

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
dat <- read_csv("./data/infection_20190809.csv")
# previous file name: Infection quantification data - Sheet1


#### edit data ####

# notes
unique(dat$notes)

# add columns
dat1 <- dat %>%
  mutate(fungus = recode(treatment, "F" = 1, "W" = 0),
         Treatment = recode(treatment, "F" = "pathogen inoculation", "W" = "control (water)"),
         Species = recode(species, "C" = "Panicum\nclandestinum", "V" = "Elymus\nvirginicus", "S" = "Eragrostis\nspectabilis"),
         mv_leaves_avg = rowMeans(cbind(mv_leaves_stem_1, mv_leaves_stem_2, mv_leaves_stem_3), na.rm = T),
         mv_leaves_est = round(mv_leaves_avg) * density,
         mv_prop_infec = mv_leaves_infec / mv_leaves_est,
         mv_leaves_healthy = mv_leaves_est - mv_leaves_infec,
         native_leaves = replace_na(native_leaves, 1),
         native_prop_infec = native_leaves_infec / native_leaves,
         native_leaves_healthy = native_leaves - native_leaves_infec)


#### visualize ####

# mv average leaves
ggplot(dat1, aes(x = density, y = mv_leaves_avg, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2)

# mv est leaves
ggplot(dat1, aes(x = density, y = mv_leaves_est, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2)

# mv infection
ggplot(dat1, aes(x = density, y = mv_prop_infec, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2)

# accidentally inoculated pot
filter(dat1, mv_prop_infec > 0 & treatment == "W")
# W 100 C 1 - repeat biomass analysis without this one

# nat infection
ggplot(dat1, aes(x = density, y = native_prop_infec, color = treatment)) +
  facet_wrap(~ species) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width= .2)

# nat and mv infection
ggplot(dat1, aes(x = mv_prop_infec, y = native_prop_infec, color = treatment)) +
  geom_point()


#### mv leaf model ####

# remove missing data
# make long by replicate
leafdat <- filter(dat1, !is.na(mv_leaves_stem_1)) %>%
  select(c(treatment:replicate, fungus, Treatment, mv_leaves_stem_1:mv_leaves_stem_3)) %>%
  gather(key = stem, value = leaves, -c(treatment:Treatment)) %>%
  mutate(species = fct_relevel(species, "V", "S", "C"),
         stem = str_extract(stem, "[[:digit:]]+"),
         pot = paste(treatment, density, species, replicate, sep = ""),
         density_scaled = (density - mean(density)) / sd(density)) %>%
  filter(!is.na(leaves))

# check stems
unique(leafdat$stem)

# distribution
hist(leafdat$leaves)
mean(leafdat$leaves)
var(leafdat$leaves)

# models attempted
# Poisson brms - lots of divergent transitions and really large estimates in the fitted 95% interval
# Poisson glmer - very small predicted values
# Beverton-Holt with a single species-fungus treatment combination with Gaussian response - underestimated values because it's going towards 0 at high density
# Gaussian linear model - better fit, misses non-linearity

# model (didn't end up using this)
mv_mod_leaves <- brm(data = leafdat, family = gaussian,
                    leaves ~ density*fungus*species + I(density^2)*fungus*species,
                    prior <- c(prior(normal(0, 1), class = "b")),
                    iter = 6000, warmup = 1000, chains = 1, cores = 1)
prior_summary(mv_mod_leaves)              
summary(mv_mod_leaves)                 
plot(mv_mod_leaves)

# simulate data
leafdat_pred <- tibble(density = seq(0, 100, length.out = 300)) %>%
  merge(tibble(fungus = rep(c(0, 1), 3),
               species = rep(c("C", "S", "V"), each = 2)),
        all = T) %>%
  as_tibble() %>%
  mutate(leaves = fitted(mv_mod_leaves, newdata = ., re_formula = NA)[, "Estimate"],
         leaves_lower = fitted(mv_mod_leaves, newdata = ., re_formula = NA)[, "Q2.5"],
         leaves_upper = fitted(mv_mod_leaves, newdata = ., re_formula = NA)[, "Q97.5"],
         Treatment = recode(fungus, "1" = "pathogen inoculation", "0" = "control (water)"))

# plot model
ggplot(leafdat, aes(x = density, y = leaves, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = leafdat_pred, aes(ymin = leaves_lower, ymax = leaves_upper), alpha = 0.5, color = NA) +
  geom_line(data = leafdat_pred) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ species)
# similar trends to biomass figure - not sure it's helpful. mostly shows Microstegium size response to total biomass and density


#### mv infection model ####

# remove missing data
# reorder species
mvdat <- filter(dat1, !is.na(mv_leaves_est) & treatment == "F") %>%
  mutate(species = fct_relevel(species, "V", "S", "C"))

# model
mv_mod_infec <- brm(data = mvdat, family = binomial,
                    mv_leaves_infec | trials(mv_leaves_est) ~ density*species,
                    prior <- c(prior(normal(0, 10), class = "b"),
                               prior(normal(0, 10), class = "Intercept")),
                    iter = 6000, warmup = 1000, chains = 3, cores = 2)
prior_summary(mv_mod_infec)              
summary(mv_mod_infec)                 
plot(mv_mod_infec)
pp_check(mv_mod_infec, nsamples = 100)

# simulate data
mv_sim_dat <- mvdat %>%
  mutate(mv_leaves_infec = fitted(mv_mod_infec, newdata = ., type = "response")[, "Estimate"],
         mv_leaves_infec_lower = fitted(mv_mod_infec, newdata = ., type = "response")[, "Q2.5"],
         mv_leaves_infec_upper = fitted(mv_mod_infec, newdata = ., type = "response")[, "Q97.5"],
         mv_prop_infec = mv_leaves_infec / mv_leaves_est,
         mv_prop_infec_lower = mv_leaves_infec_lower / mv_leaves_est,
         mv_prop_infec_upper = mv_leaves_infec_upper / mv_leaves_est) 

# plot model
ggplot(mvdat, aes(x = density, y = mv_prop_infec)) +
  geom_ribbon(data = mv_sim_dat, aes(ymin = mv_prop_infec_lower, ymax = mv_prop_infec_upper), alpha = 0.5, color = NA) +
  geom_line(data = mv_sim_dat) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ species)
# need to smooth prediction line more, but otherwise looks okay


#### native infection figure ####

# colors
col_pal = c("#80cdc1", "#dfc27d")

# number of plants
dat1 %>% 
  filter(treatment == "F") %>%
  group_by(species) %>%
  count()

# counts per species
dat1 %>%
  filter(native_leaves_infec > 0) %>%
  group_by(species, density, treatment) %>%
  summarise(n = n(),
            native_prop_infec = mean(native_prop_infec))

# format data
natdat <- filter(dat1, treatment == "F") %>%
  mutate(native_leaves_healthy = case_when(native_leaves_infec == 0 ~ NA_real_,
                                           TRUE ~ native_leaves_healthy),
         native_leaves_infec = case_when(native_leaves_infec == 0 ~ NA_real_,
                                         TRUE ~ native_leaves_infec)) %>%
  select(species, Species, density, replicate, native_leaves_infec, native_leaves_healthy) %>%
  gather(key = infection, value = native_leaves, -c(species:replicate)) %>%
  mutate(densityf = as.factor(density),
         Lesions = recode(infection, native_leaves_infec = "yes", native_leaves_healthy = "no"))

# separate by species
natdatc <- filter(natdat, species == "C")
natdatv <- filter(natdat, species == "V")
natdats <- filter(natdat, species == "S")

# default theme
theme_def <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = c(0.75, 0.73), 
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.title = element_text(size = 12, face = "italic", hjust = 0.5))

# figures
nat_plot_c <- ggplot(natdatc, aes(x = densityf, y = native_leaves, fill = Lesions)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = col_pal) +
  ggtitle("Panicum") +
  theme_def

nat_plot_v <- ggplot(natdatv, aes(x = densityf, y = native_leaves, fill = Lesions)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = col_pal) +
  ggtitle("Elymus") +
  theme_def

nat_plot_s <- ggplot(natdats, aes(x = densityf, y = native_leaves, fill = Lesions)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = col_pal) +
  ggtitle("Eragrostis") +
  theme_def

# combine plots
nat_plot_comb <- plot_grid(nat_plot_v + theme(legend.position = "none"),
                           nat_plot_s,
                           nat_plot_c + theme(legend.position = "none"),
                           nrow = 1,
                           labels = LETTERS[1:3],
                           label_size = 12)

# axes
y_plot_nat <- textGrob("Leaves per plant with lesions", gp = gpar(fontsize = 12), rot = 90)
x_plot_nat <- textGrob(expression(paste(italic(Microstegium), " density", sep = "")), gp = gpar(fontsize = 12))

# save plot
tiff("./output/Fig2.tiff", width = 7.5, height = 2.5, units = "in", res = 300)
grid.arrange(arrangeGrob(nat_plot_comb, bottom = x_plot_nat, left = y_plot_nat))
dev.off()


#### mv infection figure ####

# proportion by density
mv_prop <- mvdat %>%
  group_by(species, density) %>%
  summarise(mv_prop_infec = mean(mv_prop_infec) %>% round(2))

# make data long
mvdatl <- mvdat %>%
  select(species, density, mv_leaves_infec, mv_leaves_healthy) %>%
  gather(key = infection, value = mv_leaves, -c(species:density)) %>%
  mutate(densityf = as.factor(density),
         Lesions = recode(infection, mv_leaves_infec = "yes", mv_leaves_healthy = "no")) %>%
  left_join(mv_prop)

# separate by species
mvdatc <- filter(mvdatl, species == "C")
mvdatv <- filter(mvdatl, species == "V")
mvdats <- filter(mvdatl, species == "S")

# figures
mv_plot_c <- ggplot(mvdatc, aes(x = densityf, y = mv_leaves, fill = Lesions)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  geom_text(y = 840, aes(label = mv_prop_infec), check_overlap = T, size = 2.5) +
  scale_fill_manual(values = col_pal) +
  ggtitle(expression(paste("Native: ", italic(Panicum), sep = ""))) +
  theme_def +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 840))

mv_plot_v <- ggplot(mvdatv, aes(x = densityf, y = mv_leaves, fill = Lesions)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  geom_text(y = 840, aes(label = mv_prop_infec), check_overlap = T, size = 2.5) +
  scale_fill_manual(values = col_pal) +
  ggtitle(expression(paste("Native: ", italic(Elymus), sep = ""))) +
  theme_def +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 840))

mv_plot_s <- ggplot(mvdats, aes(x = densityf, y = mv_leaves, fill = Lesions)) +
  stat_summary(geom = "bar", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  geom_text(y = 840, aes(label = mv_prop_infec), check_overlap = T, size = 2.5) +
  scale_fill_manual(values = col_pal) +
  ggtitle(expression(paste("Native: ", italic(Eragrostis), sep = ""))) +
  theme_def +
  theme(legend.position = c(0.2, 0.5),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 840))

# combine plots
mv_plot_comb <- plot_grid(mv_plot_v,
                           mv_plot_s,
                           mv_plot_c,
                           nrow = 1,
                           labels = LETTERS[1:3],
                           label_size = 12)

# axes
y_plot_mv <- textGrob(expression(paste(italic(Microstegium), " leaves per pot", sep = "")), gp = gpar(fontsize = 12), rot = 90)
x_plot_mv <- x_plot_nat

# save plot
tiff("./output/Fig4.tiff", width = 7.5, height = 2.5, units = "in", res = 300)
grid.arrange(arrangeGrob(mv_plot_comb, bottom = x_plot_mv, left = y_plot_mv))
dev.off()


#### mv infection table ####

# extract coefficients
mv_coef <- fixef(mv_mod_infec) %>%
  as_tibble(rownames = "Coefficient")

# save
write_csv(mv_coef, "./output/mv_infection_density_coefficients.csv")


#### mv leaf figure ####

# colors
col_pal2 = c("#018571", "#a6611a")

# separate by species
ldatc <- filter(leafdat, species == "C")
ldatv <- filter(leafdat, species == "V")
ldats <- filter(leafdat, species == "S")

# figures
leaf_plot_c <- ggplot(ldatc, aes(x = density, y = leaves, color = Treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  scale_color_manual(values = col_pal2) +
  ggtitle(expression(paste("Native: ", italic(Panicum), sep = ""))) +
  theme_def +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))

leaf_plot_v <- ggplot(ldatv, aes(x = density, y = leaves, color = Treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  scale_color_manual(values = col_pal2) +
  ggtitle(expression(paste("Native: ", italic(Elymus), sep = ""))) +
  theme_def +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))

leaf_plot_s <- ggplot(ldats, aes(x = density, y = leaves, color = Treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, alpha = 0.5) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  scale_color_manual(values = col_pal2) +
  ggtitle(expression(paste("Native: ", italic(Eragrostis), sep = ""))) +
  theme_def +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(size = 12, hjust = 0.5))

# combine plots
leaf_plot_comb <- plot_grid(leaf_plot_v,
                          leaf_plot_s + theme(legend.position = "none"),
                          leaf_plot_c,
                          nrow = 1,
                          labels = LETTERS[1:3],
                          label_size = 12)

# axes
y_plot_leaf <- textGrob(expression(paste(italic(Microstegium), " leaves per tiller", sep = "")), gp = gpar(fontsize = 12), rot = 90)
x_plot_leaf <- x_plot_nat

# combine
leaf_plot_comb2 <- grid.arrange(arrangeGrob(leaf_plot_comb, bottom = x_plot_leaf, left = y_plot_leaf))

# legend
leg_leaf <- get_legend(leaf_plot_s)

# save plot
tiff("./output/S1Fig.tiff", width = 7.5, height = 3, units = "in", res = 300)
grid.arrange(arrangeGrob(leaf_plot_comb2, bottom = leg_leaf, padding = unit(1, "line")))
dev.off()
