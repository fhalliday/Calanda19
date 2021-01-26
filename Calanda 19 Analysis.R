
# packages----
require(tidyverse)
require(piecewiseSEM)
require(nlme)
require(interactions)
require(ggbeeswarm)
require(lme4)
require(V.PhyloMaker)
require(car)
require(ggeffects)
require(MuMIn)
require(cowplot)

# ggplot theme----
theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) 
theme_set(theme_figs)

update_geom_defaults('point', c(size = 1.5))
update_geom_defaults('errorbar', c(size = 0.5))

pd <- position_dodge(0.3)



# data----
calanda19.com <- read.csv("Disease/Community_disease_load.csv") %>% 
      mutate(plot = str_sub(SubplotID, 1, 4)) %>% 
  
  # add plot level data
  left_join(
    read.csv("Plots/Plot_elevation.csv") %>% 
              select(-Notes)) %>%
  
  # add plant diversity data
  left_join(
    read.csv("Vegetation/Calanda19_diversity.csv") %>% 
      select(SubplotID = PlotID, plant.richness:hill.simpson)
  ) %>% 
  
  # add phylogenetic diversity
  left_join(
    read.csv("Phylogeny/Calanda_19_mpd_no_trees.csv") %>% 
      select(SubplotID = PlotID, mpd.obs, mpd.obs.z)
  ) %>% 
  
  # add traits
  left_join(
    read.csv("Traits/Calanda_2019_community_traits_no_trees.csv") %>% 
      select(SubplotID = PlotID, Chlorophyll:RaoQ)
  ) %>% 
  
  # add data-loggers
  left_join(
    read.csv("TMS/Calanda_summarized_tms_2019.csv") %>% 
      select(site = Site, mean.min.soil.t:mean.moisture)
  ) %>% 
  
  # some transformations
  mutate(
    l.richness = car::logit(plant.richness), 
    l.simpson = car::logit(hill.simpson),
    Chlo.z = scale(Chlorophyll),
         CN.z = scale(CN),
         f1.z = scale(f1),
         H.z = scale(Height),
         Lon.z = scale(Leaf_lifespan),
         N.z = scale(Leaf_N),
         P.z = scale(Leaf_P),
         Amax.z = scale(Photosynthetic_rate),
         SLA.z = scale(SLA),
         SM.z = scale(Seed_mass),
         elevation_centered = scale(elevation, scale = F),
         soil_surf_t_scaled = scale(mean.soil.surface.t),
         soil_m_scaled = scale(mean.moisture),
         richness_centered = scale(plant.richness, scale =F),
         l.richness_centered = scale(l.richness, scale = F),
         mpd.obs_scaled = scale(mpd.obs),
         sqrt.Disease.l = sqrt(Disease.l)) %>% 
  
  # some additional variables
  mutate(
    Meadow = meadow,
    Site = site,
    PlotID = plot
  )

summary(calanda19.com)
names(calanda19.com)
# correlations ----
# look at data
calanda19.com %>% 
  dplyr::select(elevation, mean.soil.t, mean.soil.surface.t, mean.air.t, mean.moisture, Disease.l, plant.richness, hill.shannon, f1, Seed_mass, Height, mpd.obs.z, FRic) %>% 
  pairs()

calanda19.com %>% 
  dplyr::select(elevation, mean.soil.t, mean.soil.surface.t, mean.air.t, mean.moisture, Disease.l, plant.richness, hill.shannon, f1, Seed_mass, Height, mpd.obs.z, FRic) %>% 
  cor()

# shannon and richness are colinear r = .75
# functional richness and plant richness are positively correlated, but not colinear r = 0.42
# temp and elevation are colinear, so probably best not to include temp in a model, but note that increasing elevation was strongly associated with reduced temperatures.
# moisture is almost entirely independent of elevation, so probably not useful here either.

# check model assumptions ----
# check that component models fit assumptions
com.dis.ri <- lme(Disease.l ~ plant.richness * elevation + f1*elevation + mpd.obs.z*elevation, random = ~1|Meadow/Site/PlotID,  data = calanda19.com)
# check residuals
resids.fig <- function(mod, df) {
  residdf <- dplyr::mutate(df, resids = residuals(mod, type = 'normalized'),
                           fits = fitted(mod))
  fig2 <-ggplot(residdf, aes(x = fits, y = resids)) + geom_point() +
    labs(x = 'Fitted values', y = '')
  
  fig3 <- ggplot(residdf) + stat_qq(aes(sample = resids)) +
    labs(x = 'Theoretical Quantiles', y = 'Sample Quantiles')
  
  # qqline plot = FALSE, according to James should work
  
  fig4 <- ggplot(residdf, aes(x = resids)) + geom_histogram(aes(y=..density..), colour = 'grey50') +
    labs(x = 'Residuals', y = 'Frequency') + scale_y_continuous(expand = c(0, 0)) +
    stat_function(fun = dnorm, color = "red", args = list(mean = mean(residdf$resids),
                                                          sd = sd(residdf$resids)))
  grid::grid.draw(rbind(ggplotGrob(fig2), ggplotGrob(fig3), ggplotGrob(fig4), size = 'first'))
  
  return(summary(mod))
}
variables.fig <- function(df, mod, variable){
  df %>% 
    mutate(resids = residuals(mod, type = 'normalized')) %>% 
    ggplot(aes_string(x=variable, y="resids", color=variable)) +
    xlab(variable) +
    ylab("residuals") +
    # geom_boxplot() +
    geom_point(position=position_jitter(h=0, w=0.4))
}

resids.fig(com.dis.ri, calanda19.com)
# appears to violate assumptions of heteroscedasticity and normality

variables.fig(calanda19.com, com.dis.ri, "elevation")
variables.fig(calanda19.com, com.dis.ri, "plant.richness")
variables.fig(calanda19.com, com.dis.ri, "mpd.obs.z")
variables.fig(calanda19.com, com.dis.ri, "f1")

# square-root transform the response
com.dis.ri2 <- lme(sqrt(Disease.l) ~ plant.richness * elevation + 
                     f1*elevation + 
                     mpd.obs.z*elevation, 
                   random = ~1|Meadow/Site/PlotID,  data = calanda19.com)
# check residuals
resids.fig(com.dis.ri2, calanda19.com)
# that normalized the residuals somewhat, but still appears to violate assumptions of heterscedasticity
# still appears to violate assumptions of heteroscedasticity

variables.fig(calanda19.com, com.dis.ri2, "plant.richness")
variables.fig(calanda19.com, com.dis.ri2, "elevation")
variables.fig(calanda19.com, com.dis.ri2, "mpd.obs.z")
variables.fig(calanda19.com, com.dis.ri2, "f1")

# model variance separately among sites
com.dis.ri3 <- lme(sqrt.Disease.l ~ plant.richness * elevation + 
                     f1.z*elevation + 
                     mpd.obs.z*elevation,
                   random = ~1|Meadow/Site/PlotID,
                   weights = varIdent(form = ~ 1 | Site),
                   data = calanda19.com,
                   control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
resids.fig(com.dis.ri3, calanda19.com)
# much much better
variables.fig(calanda19.com, com.dis.ri3, "plant.richness")
variables.fig(calanda19.com, com.dis.ri3, "elevation")
variables.fig(calanda19.com, com.dis.ri3, "mpd.obs.z") # this one still isn't perfect -- possibly due to outliers?
variables.fig(calanda19.com, com.dis.ri3, "f1")

# plant richness model
pr.ri <- lme(plant.richness ~ elevation, 
             random = ~1|Meadow/Site/PlotID,
             weights = varIdent(form = ~ 1 | Site),
             data = calanda19.com,
             control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

resids.fig(pr.ri, calanda19.com)
# looks ok...

variables.fig(calanda19.com, pr.ri, "elevation")
# looks pretty good
variables.fig(calanda19.com, pr.ri, "mean.moisture")

# pace of life model
f1.ri <- lme(f1 ~ elevation, 
             random = ~1|Meadow/Site/PlotID,
             weights = varIdent(form = ~ 1 | Site),
             data = calanda19.com,
             control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

resids.fig(f1.ri, calanda19.com)
# seems fine

# phylogenetic diversity model
mpd.ri <- lme(mpd.obs.z ~ elevation, 
             random = ~1|Meadow/Site/PlotID,
             weights = varIdent(form = ~ 1 | Site),
             data = calanda19.com,
             control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

resids.fig(mpd.ri, calanda19.com)
# this one seems a bit odd, possibly due to the presence of outliers?

# 1. does elevation modify community structure?----

calanda19.fmod <- lme(f1 ~ elevation, 
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.fmod)
# marginal effect of elevation on functional traits
MuMIn::r.squaredGLMM(calanda19.fmod)

calanda19.mpdmod <- lme(mpd.obs.z ~ elevation, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.mpdmod)
# significant effect of elevation on mean-pairwise phylogenetic diversity
MuMIn::r.squaredGLMM(calanda19.mpdmod)

calanda19.rmod <- lme(plant.richness ~ elevation, 
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.rmod)
# significant but very weak effect of elevation on richness
MuMIn::r.squaredGLMM(calanda19.rmod)

r.e <- ggeffects::ggpredict(calanda19.rmod, "elevation") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = elevation, y = plant.richness), position = position_jitter(), shape = 1)+
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "\nHost richness\n", x = "Elevation", title = "")

# pdf("Figures/richness_elevation.pdf", height = 4, width = 6)
# r.e
# dev.off()


mpd.e <- ggeffects::ggpredict(calanda19.mpdmod, "elevation") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = elevation, y = mpd.obs.z), position = position_jitter(), shape = 1) + 
  geom_hline(yintercept = 0, lty = 2, color = "grey50") +
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Host phylogenetic \ndiversity (mpd.obs.z)", x = "Elevation", title = "")

# pdf("Figures/mpdz_elevation.pdf", height = 4, width = 6)
# mpd.e
# dev.off()

f.e <- ggeffects::ggpredict(calanda19.fmod, "elevation") %>% plot(rawdata = F) +
  geom_point(data = calanda19.com, aes(x = elevation, y = f1), position = position_jitter(), shape = 1) + 
  theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Host community \npace of life", x = "Elevation", title = "")

# pdf("Figures/f1_elevation.pdf", height = 4, width = 6)
# f.e
# dev.off()

# pdf("Figures/moderators_elevation.pdf", height = 7, width = 4)
# plot_grid(r.e, mpd.e, f.e, ncol = 1)
# dev.off()

# 2. does elevation modify how community structure alters disease risk? ----
calanda19.dmod <- lme(sqrt.Disease.l ~ elevation + plant.richness + mpd.obs.z + f1 +
              plant.richness*elevation + mpd.obs.z*elevation + f1*elevation,
            random = ~1|Meadow/Site/PlotID,
            weights = varIdent(form = ~ 1 | Site),
            data = calanda19.com,
            control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.dmod)
# r squared
MuMIn::r.squaredGLMM(calanda19.dmod)

# although a spatial
# autocorrelation term could be used here instead of the nested blocks, we did not GPS each subplot, so we stick with this analysis instead of something spatially explicit

# Yes-- elevation modifies the relationship between host community functional traits and disease
# show this effect

f1_modplot <- emmeans::emtrends(calanda19.dmod, ~elevation, var = "f1", at = list(elevation = seq(600,1800, 1))) %>% 
  data.frame() %>% 
  mutate(in.data = case_when(
    elevation < min(calanda19.com$elevation) ~ "No",
    elevation > max(calanda19.com$elevation) ~ "No",
    TRUE ~ "Yes"
  )) %>% 
  ggplot(aes(x = elevation, y = f1.trend)) +
  geom_ribbon(fill = "grey", aes(ymin = lower.CL, ymax = upper.CL)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.com, aes(y = f1), sides = "b") +
  labs(y = "Effect of community pace-of-life on disease", x = "Elevation")

# pdf("Figures/F1_moderation.pdf", height = 4, width = 5)
# f1_modplot
# dev.off()
# 

f1_eleplot <- calanda19.com %>% 
  mutate(elevation_grouped = case_when(
    elevation < 1000 ~ "low elevation (< 1000m)",
    elevation >= 1000 & elevation <= 1500 ~ "mid elevation (1000m - 1500m)",
    elevation > 1500 ~ "high elevation (> 1500m)"
  )) %>% 
  ggplot(aes(x = f1, y = sqrt.Disease.l,  group = fct_reorder(elevation_grouped, elevation))) + 
  facet_wrap(~fct_reorder(elevation_grouped, elevation), ncol = 1, scales = "free_y")+
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", color = "black") +
  theme(legend.position = c(.4, .9)) +
  labs(y = "Community disease load (square-root transformed)", x = "Host community pace-of-life")

# pdf("Figures/F1_disease_elevation_raw_data.pdf", height = 7, width = 4)
# f1_eleplot
# dev.off()
# 

# finally to make it super clear, drop mid-elevation
calanda19.com %>% 
  mutate(elevation_grouped = case_when(
    elevation <= 1000 ~ "low elevation",
    elevation > 1000 & elevation <= 1500 ~ "mid elevation",
    elevation > 1500 ~ "high elevation"
  )) %>% 
  filter(elevation_grouped != "mid elevation") %>% 
  ggplot(aes(x = f1.z, y = sqrt.Disease.l, color = fct_reorder(elevation_grouped, elevation), group = fct_reorder(elevation_grouped, elevation))) + 
  # facet_wrap(~fct_reorder(elevation_grouped, elevation), ncol = 1, scales = "free_y")+
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", show.legend = F) +
  labs(y= "Community disease", x = "Community pace of life") +
  colorblindr::scale_colour_OkabeIto() +
  theme(legend.position = c(.5, .85))
# but this is pretty misleading since it ignores 1/3 of the data

# What's happening with species richness?

richness_eleplot <- calanda19.com %>% 
  mutate(elevation_grouped = case_when(
    elevation < 1000 ~ "low elevation (< 1000m)",
    elevation >= 1000 & elevation <= 1500 ~ "mid elevation (1000m - 1500m)",
    elevation > 1500 ~ "high elevation (> 1500m)"
  )) %>% 
  ggplot(aes(x = plant.richness, y = sqrt.Disease.l,  group = fct_reorder(elevation_grouped, elevation))) + 
  facet_wrap(~fct_reorder(elevation_grouped, elevation), ncol = 1, scales = "free_y")+
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", color = "black") +
  theme(legend.position = c(.4, .9)) +
  labs(y = "Community disease load (square-root transformed)", x = "Host richness")

# pdf("Figures/richness_disease_elevation_raw_data.pdf", height = 7, width = 4)
# richness_eleplot
# dev.off()

r.d <- calanda19.com %>% 
  ggplot(aes(x = plant.richness, y = sqrt.Disease.l)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", color = "black") +
  theme(legend.position = c(.4, .9)) +
  labs(y = "Community disease load (square-root transformed)", x = "Host richness")

# pdf("Figures/richness_disease_raw_data.pdf", height = 4, width = 6)
# r.d
# dev.off()


# and the main effect of elevation?
d.e <- calanda19.com %>% 
  ggplot(aes(x = elevation, y = sqrt.Disease.l)) + 
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", color = "black") +
  theme(legend.position = c(.4, .9)) +
  labs(y = "Community disease load (square-root transformed)", x = "Elevation")
# 
# pdf("Figures/disease_elevation_raw_data.pdf", height = 4, width = 6)
# d.e
# dev.off()


# Do we trust this model? leave one out cross validation
# this takes a long time to run, so I'm commenting it out 
# xval <- data.frame()
# for (i in 1:nrow(calanda19.com)){
#   train.data  <- calanda19.com[-i, ]
#   test.data <- calanda19.com[i, ]
#   model <- update(calanda19.dmod, data = train.data)
#   xval <- rbind(xval, data.frame(obs = test.data$sqrt.Disease.l, pred = predict(model, test.data)))
# 
#   }
# 
# # calculate the RMSE
# xval %>%
#   mutate(error = obs - pred) %>%
#   summarize(rmse = sqrt(mean(error^2)))
# 
# # compare that with the RMSE of the full model
# sqrt(mean(calanda19.dmod$residuals^2))

# So for this model I could report the pseudo R2 (marginal = .238, conditional = .460), RMSE = 0.295, and LOOCV RMSE = 0.315

# 2. can this effect be attributed to changing distributions of host composition with elevation?


# 3. What is the relative of each path? ----
model.c19 <- psem(
 
  lme(sqrt.Disease.l ~ elevation + plant.richness + mpd.obs.z + f1 + 
        f1*elevation,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(f1 ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(mpd.obs.z ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  # residual correlations among variables
  # mpd.obs.z is independent of species richness, so no need to include that here
  f1 %~~% plant.richness,
  f1 %~~% mpd.obs.z)
# note that warning about NAs refers to data that we won't ever use, so i think it's ok.

model.c19.fit <- summary(model.c19, .progressBar = T)
model.c19.fit

# visualize the standardized direct, indirect, effects, 
# I realize that this code is very ugly. That's because I was figuring out what I wanted as I wrote the code.
c19.effect.summary <- data.frame(model.c19.fit$coefficients) %>% 
  filter(Response == "sqrt.Disease.l", Predictor != "elevation", Predictor != "elevation:f1",
         # exclude mpd.obs.z since that path is not supported by the model
         Predictor != "mpd.obs.z") %>% 
  mutate(Effect = "Biotic drivers") %>% 
  group_by(Effect) %>% 
  summarize(net_abs_effect = sum(abs(Std.Estimate)),
            net_effect = sum(Std.Estimate)) %>% 
  rbind(
    data.frame(model.c19.fit$coefficients) %>% 
      filter(Response == "sqrt.Disease.l", Predictor == "elevation") %>% 
      mutate(Effect = "Elevation (direct)") %>% 
      group_by(Effect) %>% 
      summarize(net_abs_effect = sum(abs(Std.Estimate)),
                net_effect = sum(Std.Estimate))
    
  ) %>% 
  rbind(
    data.frame(model.c19.fit$coefficients) %>% 
      filter(Response == "sqrt.Disease.l", Predictor == "elevation:f1") %>% 
      mutate(Effect = "Elevation (interactive)") %>% 
      group_by(Effect) %>% 
      summarize(net_abs_effect = sum(abs(Std.Estimate)),
                net_effect = sum(Std.Estimate))
  ) %>% 
  rbind(
    data.frame(model.c19.fit$coefficients) %>% 
      filter(Response == "sqrt.Disease.l", Predictor != "elevation", Predictor != "elevation:f1") %>% 
      mutate(moderator = Predictor) %>% 
      rbind(
        data.frame(model.c19.fit$coefficients) %>% 
          filter(Response != "sqrt.Disease.l", Predictor == "elevation") %>% 
          mutate(moderator = Response)
      ) %>% 
      mutate(Effect = "Elevation (indirect)") %>% 
      group_by(Effect, moderator) %>% 
      summarize(Std.Estimate = prod(Std.Estimate)) %>% 
      # drop f1 and phylogenetic diversity paths, since these aren't supported by the model
      # no significant effect of elevation on F1, no significant effect of phylogenetic diversity on disease
      filter(moderator == "plant.richness") %>% 
      group_by(Effect) %>% 
      summarize(net_abs_effect = sum(abs(Std.Estimate)),
                net_effect = sum(Std.Estimate))
  ) %>% 
  gather(
    key = "Type", value = "Estimate", net_abs_effect:net_effect
  ) %>% 
  rbind(
    data.frame(model.c19.fit$coefficients) %>% 
      filter(Response == "sqrt.Disease.l", Predictor != "elevation", Predictor != "elevation:f1") %>% 
      transmute(Effect = "Biotic drivers", Type = Predictor, Estimate = Std.Estimate)) %>% 
  rbind(
    data.frame(model.c19.fit$coefficients) %>% 
      filter(Response == "sqrt.Disease.l", Predictor != "elevation", Predictor != "elevation:f1") %>% 
      mutate(Type = Predictor) %>% 
      rbind(
        data.frame(model.c19.fit$coefficients) %>% 
          filter(Response != "sqrt.Disease.l", Predictor == "elevation") %>% 
          mutate(Type = Response)
      ) %>% 
      mutate(Effect = "Elevation (indirect)") %>% 
      group_by(Effect, Type) %>% 
      summarize(Estimate = prod(Std.Estimate)) %>% ungroup()
  )


c19.effect.summary %>% 
  mutate(Driver = ifelse(Effect == "Biotic drivers", "Biotic drivers", "Elevation")) %>% 
  ggplot(aes(x = fct_relevel(Effect, "Elevation (indirect)"), y = Estimate, fill = fct_relevel(Type, "net_abs_effect", "net_effect", "f1", "mpd.obs.z", "plant.richness"))) +
  facet_grid(~Driver, scales = "free_x") + 
  geom_col(position = position_dodge())


c19.effect.summary %>% 
  mutate(Driver = ifelse(Effect == "Biotic drivers", "Biotic drivers", "Elevation")) %>% 
  ggplot(aes(x = fct_relevel(Type, "net_abs_effect", "net_effect", "f1", "mpd.obs.z", "plant.richness"), y = Estimate)) +
  facet_grid(~Driver, scales = "free_x") + 
  geom_col(position = position_dodge())


# this is even more absurd code because I just made it on the fly, but it shows the effects that I'm interested in
std.effecs.fig <- c19.effect.summary %>%  
  filter(Effect != "Biotic drivers", Type == "net_abs_effect"| Type == "net_effect") %>% 
  group_by(Type) %>% 
  summarize(Estimate = sum(Estimate)) %>% 
  mutate(Type = case_when(
    Type == "net_abs_effect" ~ "Net |effect|",
    Type == "net_effect" ~ "Net effect",
    Type == "plant.richness" ~ "Rich",
    Type == "mpd.obs.z" ~ "Phylo",
    Type == "f1" ~ "Pace"
  ),
  Driver = "Elevation") %>% 
  rbind(
    c19.effect.summary %>%  
      filter(Effect == "Elevation (direct)", Type == "net_effect") %>% 
      transmute(Type = case_when(
        Type == "net_effect" ~ "Direct"
      ),
      Estimate = Estimate,
      Driver = "Elevation")
  ) %>% 
  rbind(
    c19.effect.summary %>% 
      filter(Effect == "Biotic drivers") %>% 
      mutate(Type = case_when(
        Type == "net_abs_effect" ~ "Net |effect|",
        Type == "net_effect" ~ "Net effect",
        Type == "plant.richness" ~ "Rich",
        Type == "mpd.obs.z" ~ "Phylo",
        Type == "f1" ~ "Pace"
      ),
             Driver = "Community structure") %>% 
      select(-Effect) 
  ) %>% 
  rbind(
    c19.effect.summary %>% 
      filter(Effect != "Biotic drivers", Type != "net_abs_effect", Type != "net_effect") %>% 
      transmute(Type = case_when(
        Type == "net_abs_effect" ~ "Net |effect|",
        Type == "net_effect" ~ "Net effect",
        Type == "plant.richness" ~ "Rich",
        Type == "mpd.obs.z" ~ "Phylo",
        Type == "f1" ~ "Pace"
      ),
                Estimate = Estimate, Driver = "Elevation")
    
  ) %>% 
  rbind(
    c19.effect.summary %>% 
      filter(Effect == "Elevation (interactive)", Type != "net_abs_effect") %>% 
      transmute(Type = case_when(
        Type == "net_effect" ~ "Interactive"
      ),
                Estimate = Estimate, Driver = "Elevation")
  ) %>% 
  mutate(effect = case_when(
    Type == "Net |effect|" ~  "Net |effect|",
    Type == "Net effect" ~ "Net effect",
    TRUE ~ "Path coefficent"),
    # change type Net |effect| to just "Net"
    Type = case_when(
      Type == "Net |effect|" ~ "Net",
      TRUE ~ as.character(Type)
    )
    ) %>% 
  filter(Type != "Net effect") %>% 
  # annotate the bars that aren't actually supported by the model
  mutate(n.s = c("", "", "", "", "n.s.", "", "n.s.", "n.s.","","" )) %>% 
  ggplot(aes(x = fct_relevel(Type, "Net"), y = Estimate, fill = effect)) + 
  facet_wrap(~Driver, scales = "free",strip.position = "bottom") +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = c("black", "grey")) +
  # annotate the bars that aren't actually supported by the model
  geom_text(aes(label = n.s), vjust = 1) +
  labs(y = "Standardized estimate", x = "") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-.62, 1)) +
  geom_text(data = data.frame(
    Driver = c("Community structure", "Elevation"),
    effect = c("Path coefficent", "Path coefficent"),
    lbl = c("A)", "B)")),
    aes(x = .6, y = 1, label = lbl), vjust = -.5) 

# pdf("Figures/SEM_effects_22Dec20.pdf", height = 5, width = 10)
# std.effecs.fig
# dev.off()

# this figure shows indirect vs direct standardized path coefficients
# paths that were not supported in the model are annotated with "n.s."
# The net effect (red) is the sum of the absolute value of all paths, excluding those that were not supported in the model.
# A) paths associated with host community structure. Black bars are the direct effects of host community structure on parasite community load
# B) paths associated with elevation, "Direct" is the direct effect of elevation on community parasite load, "Interactive" is the interaction between Elevation and host pace of life, the pace, phylo, and rich bars are the product of indirect pathways involving these three mediators.



# Not currently included in the manuscript----
# Would a functional diversity metric improve this model?
# I'm not sure what the rationale for this would be. More functionally diverse assemblages have lower density of any particular type of host?
model.c19_2 <- psem(
  
  lme(sqrt.Disease.l ~ elevation + plant.richness + mpd.obs.z + f1 + FRic +
        f1*elevation,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(f1 ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(FRic ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(mpd.obs.z ~ elevation, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  # residual correlations among variables
  f1 %~~% plant.richness,
  f1 %~~% mpd.obs.z,
  FRic %~~% plant.richness,
  FRic %~~% mpd.obs.z,
  FRic %~~% f1)

summary(model.c19_2, .progressBar = T)

# Incorporating functional diversity does not improve this model, though
# elevation is a significant predictor of functional diversity -- with
# functional diversity increasing with elevation. and functional diversity is
# correlated with community-level pace-of-life


# Is it possible that individual traits are driving the host pace-of-life response?----
MuMIn::AICc(calanda19.dmod)
# So the model with pace-of-life had AICc 207.99, marginal R2 = .238, conditional R2 = .460, RMSE = 0.295, and LOOCV RMSE = 0.315

# I'll start with a base model, to make this a bit simpler.
base.model <- lme(sqrt.Disease.l ~ plant.richness*elevation + mpd.obs.z*elevation,
                  random = ~1|Meadow/Site/PlotID,
                  weights = varIdent(form = ~ 1 | Site),
                  data = calanda19.com,
                  control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

# First model replaces f1 with Chlorophyll, Leaf_lifespan, Leaf_N, Leaf_P, SLA
traitmod <- update(base.model, .~. + Chlo.z*elevation + 
                     Lon.z*elevation +
                     N.z*elevation +
                     P.z*elevation +
                     SLA.z*elevation)
car::Anova(traitmod)
MuMIn::AICc(traitmod)
MuMIn::r.squaredGLMM(traitmod)

# LOOCV function
LOOCV <- function(m) {
  for (i in 1:nrow(calanda19.com)){
    train.data  <- calanda19.com[-i, ]
    test.data <- calanda19.com[i, ]
    model <- update(m, data = train.data)
    xval <- rbind(xval, data.frame(obs = test.data$sqrt.Disease.l, pred = predict(model, test.data)))
    
  }
  
  # calculate the RMSE
  LOOCV.rmse <- xval %>%
    mutate(error = obs - pred) %>%
    summarize(LOOCV.rmse = sqrt(mean(error^2)))
  
  # compare that with the RMSE of the full model
  rmse <- sqrt(mean(m$residuals^2))
  
  cbind(rmse = rmse, LOOCV.rmse = LOOCV.rmse)
}

# commented out becasue this takes a long time to run
# LOOCV(traitmod)

# AIC 319.42, marginal R2 = .252, conditional R2 = .432, RMSE = 0.293, and LOOCV RMSE = 0.319

# model with Chlorophyll content only
Chlomod <- update(base.model, .~. + Chlo.z * elevation)
car::Anova(Chlomod)
MuMIn::AICc(Chlomod)
MuMIn::r.squaredGLMM(Chlomod)
# LOOCV(Chlomod)
#AICc 220.28, marginal R2 = .167, conditional R2 = .393, RMSE = 0.296, and LOOCV RMSE = 0.318

# Model with leaf longevity only
Lonmod <- update(base.model, .~. + Lon.z * elevation)
car::Anova(Lonmod)
MuMIn::AICc(Lonmod)
MuMIn::r.squaredGLMM(Lonmod)
# LOOCV(Lonmod)
# AICc 217.93, marginal R2 = .183, conditional R2 = .434, RMSE = 0.292, and LOOCV RMSE = 0.315

# Model with leaf nitrogen only
Nmod <- update(base.model, .~. + N.z * elevation)
car::Anova(Nmod)
MuMIn::AICc(Nmod)
MuMIn::r.squaredGLMM(Nmod)
# LOOCV(Nmod)
# AICc 221.72, marginal R2 = .166, conditional R2 = .393, RMSE = 0.294, and LOOCV RMSE = 0.317

# Model with leaf phosphorus only
Pmod <- update(base.model, .~. + P.z * elevation)
car::Anova(Pmod)
MuMIn::AICc(Pmod)
MuMIn::r.squaredGLMM(Pmod)
# LOOCV(Pmod)
# AICc 220.09, marginal R2 = .176, conditional R2 = .399, RMSE = 0.298, and LOOCV RMSE = 0.318

# Model with SLA only
SLAmod <- update(base.model, .~. + SLA.z * elevation)
car::Anova(SLAmod)
MuMIn::AICc(SLAmod)
MuMIn::r.squaredGLMM(SLAmod)
# LOOCV(SLAmod)
# AICc 212.41, marginal R2 = .251, conditional R2 = .391, RMSE = 0.298, and LOOCV RMSE = 0.316

# So of the traits associated with pace-of-life, only one, SLA, interacted with elevation to influence community disease load.

# this is interesting. But remember the trait-continuum was calculated at the
# host species level, and leaf chlorophyll content loaded even more strongly
# onto that trait than SLA. but when we scale that up to the community level, we
# do not see the same patterns.
calanda19.com %>% 
  gather(trait, trait.value, Chlorophyll, Leaf_lifespan, Leaf_N, Leaf_P, SLA) %>% 
  ggplot(aes(x = f1, y = trait.value)) +
  facet_wrap(~trait, scales = "free_y") + 
  geom_point(shape = 1) +
  geom_smooth(method = lm, formula = y~x, color = "black", se = F) +
  labs(y = "Community weighted mean trait value", x = "Community weigthed mean pace of life")

# here's the same relationship at the species level
calanda.traits <- read.csv("Traits/Calanda_2019_plant_traits_no_trees.csv")

calanda.traits %>% 
  gather(trait, trait.value, Chlorophyll, Leaf_lifespan, Leaf_N, Leaf_P, SLA) %>% 
  ggplot(aes(x = f1, y = trait.value)) +
  facet_wrap(~trait, scales = "free_y") + 
  geom_point(shape = 1) +
  geom_smooth(method = lm, formula = y~x, color = "black", se = F) +
  labs(y = "Community weighted mean trait value", x = "Community weigthed mean pace of life")

# kind-of points to the value of omitting single traits analyses, since some
# traits were missing for many taxa, making the CWM value of those traits a poor
# representation of the community as a whole. because we used several traits to
# do our factor analysis, we were able to average over those missing traits at
# the species level, which gives us a much stronger estimate of CWM at the
# community level. note that our community-level SLA and F1 estimates are likely
# so similar because there's almost no missing data on SLA at the species level, 
# and SLA had a very high factor loading onto pace of life.
