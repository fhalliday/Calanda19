
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

# data----
calanda19.com <- read.csv("Disease/Community_disease_load.csv", sep = ",", strip.white = T) %>% 
      mutate(plot = str_sub(SubplotID, 1, 4)) %>% 
  
  # add plot level data
  left_join(
    read.csv("Plots/Plot_elevation.csv", sep = ",", strip.white = T) %>% 
              select(-Notes)) %>%
  
  # add plant diversity data
  left_join(
    read.csv("Vegetation/Calanda19_diversity.csv", sep = ",", strip.white = T) %>% 
      select(SubplotID = PlotID, plant.richness:hill.simpson)
  ) %>% 
  
  # add phylogenetic diversity
  left_join(
    read.csv("Phylogeny/Calanda_19_mpd_no_trees.csv", sep = ",", strip.white = T) %>% 
      select(SubplotID = PlotID, mpd.obs, mpd.obs.z)
  ) %>% 
  
  # add traits
  left_join(
    read.csv("Traits/Calanda_2019_community_traits_no_trees.csv", sep = ",", strip.white = T) %>% 
      select(SubplotID = PlotID, Chlorophyll:RaoQ)
  ) %>% 
  
  # add data-loggers
  left_join(
    read.csv("TMS/Calanda_summarized_tms_2019.csv", sep = ",", strip.white = T) %>% 
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
    # rescale variables
    soil_surf_t_lowsd_scaled = mean.soil.surface.t-(mean(mean.soil.surface.t)-sd(mean.soil.surface.t)),
    soil_surf_t_highsd_scaled = mean.soil.surface.t-(mean(mean.soil.surface.t)+sd(mean.soil.surface.t)),
    soil_surf_t_centered = scale(mean.soil.surface.t, scale = F),
    soil_t_centered = scale(mean.soil.t, scale = F),
    air_t_centered = scale(mean.air.t, scale = F),
    moisture_centered = scale(mean.moisture, scale = F),
    elevation_centered = scale(elevation, scale = F),
    sqrt.Disease.l = sqrt(Disease.l),
    f1_scaled = scale(f1, scale = F),
    richness_scaled = scale(plant.richness, scale = F)) %>% 
  
  # some additional variables
  mutate(
    Meadow = meadow,
    Site = site,
    PlotID = plot
  )

summary(calanda19.com)
names(calanda19.com)

calanda19.site <-  read.csv("TMS/Calanda_summarized_tms_2019_plotIDs.csv") %>% 
  select(PlotID = plot, Site, mean.min.soil.t:mean.moisture) %>% 
  # add elevation
  left_join(
    read.csv("Plots/Plot_elevation.csv") %>% 
      select(PlotID = plot, Meadow = meadow, elevation)) %>% 
  # standardization
  mutate(soil_surf_t_lowsd_scaled = mean.soil.surface.t-(mean(mean.soil.surface.t)-sd(mean.soil.surface.t)),
         soil_surf_t_highsd_scaled = mean.soil.surface.t-(mean(mean.soil.surface.t)+sd(mean.soil.surface.t)),
         soil_surf_t_centered = scale(mean.soil.surface.t, scale = F)) %>% 
  # to this, I need to add site-level measures of other variables so that I can use it in the SEM
  left_join(calanda19.com %>% 
  group_by(Site) %>% 
  summarize(sqrt.Disease.l = mean(sqrt.Disease.l),
            plant.richness = mean(plant.richness),
            f1_scaled = mean(f1_scaled)))


# check model assumptions ----

# Analysis for eLife revision -- replacing Elevation with Temperature in the models----
# Major changes to this version
# 1. I dropped phylogenetic diversity following reviewer comments
# 2. I replaced elevation with soil moisture and temperature in mixed models
# 3. I incorporated temperature into the SEM, though I left out soil moisture, since it was unrelated to any other variable in the analysis.

# how does elevation influence local microclimate?
calanda19.stmod <- lm(mean.soil.t ~ elevation, data = calanda19.site)
anova(calanda19.stmod)
summary(calanda19.stmod)

calanda19.sstmod <- lm(mean.soil.surface.t ~ elevation, data = calanda19.site)
anova(calanda19.sstmod)
summary(calanda19.sstmod)

calanda19.atmod <- lm(mean.air.t ~ elevation, data = calanda19.site)
anova(calanda19.atmod)
summary(calanda19.atmod)

calanda19.mmod <- lm(mean.moisture ~ elevation, data = calanda19.site)
anova(calanda19.mmod)
summary(calanda19.mmod)

# 1. do abiotic conditions modify community structure?----
# traits
calanda19.fmod.t <- lme(f1 ~ mean.soil.surface.t + mean.moisture, 
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.fmod.t)
# highly significant effect of temp on functional traits
# that's interesting -- very different from the elevation mod
MuMIn::r.squaredGLMM(calanda19.fmod.t)
summary(calanda19.fmod.t)

# soil t
calanda19.fmod.t2 <- lme(f1 ~ mean.soil.t + mean.moisture, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.fmod.t2)
MuMIn::r.squaredGLMM(calanda19.fmod.t2)
summary(calanda19.fmod.t2)

# air t
calanda19.fmod.t3 <- lme(f1 ~ mean.air.t + mean.moisture, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.fmod.t3)
MuMIn::r.squaredGLMM(calanda19.fmod.t3)
summary(calanda19.fmod.t3)

# elevation
calanda19.fmod.e <- lme(f1 ~ elevation + mean.moisture, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.fmod.e)
MuMIn::r.squaredGLMM(calanda19.fmod.e)
summary(calanda19.fmod.e)

# richness
calanda19.rmod.t <- lme(plant.richness ~ mean.soil.surface.t + mean.moisture, 
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.rmod.t)
# significant effect of t on plant richness too
MuMIn::r.squaredGLMM(calanda19.rmod.t)
summary(calanda19.rmod.t)

# soil temp
calanda19.rmod.t2 <- lme(plant.richness ~ mean.soil.t + mean.moisture, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
summary(calanda19.rmod.t2)
car::Anova(calanda19.rmod.t2)
MuMIn::r.squaredGLMM(calanda19.rmod.t2)

# air temp
calanda19.rmod.t3 <- lme(plant.richness ~ mean.air.t + mean.moisture, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
summary(calanda19.rmod.t3)
car::Anova(calanda19.rmod.t3)
MuMIn::r.squaredGLMM(calanda19.rmod.t3)

# elevation
calanda19.rmod.e <- lme(plant.richness ~ elevation + mean.moisture, 
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
summary(calanda19.rmod.e)
car::Anova(calanda19.rmod.e)
MuMIn::r.squaredGLMM(calanda19.rmod.e)

# 2. does temp modify how community structure alters disease risk? ----
calanda19.dmod.t <- lme(sqrt.Disease.l ~ soil_surf_t_centered + richness_scaled + f1_scaled +
                          richness_scaled*soil_surf_t_centered + f1_scaled*soil_surf_t_centered+
                          richness_scaled*moisture_centered + f1_scaled*moisture_centered,
                      random = ~1|Meadow/Site/PlotID,
                      weights = varIdent(form = ~ 1 | Site),
                      data = calanda19.com,
                      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.dmod.t)
# r squared
MuMIn::r.squaredGLMM(calanda19.dmod.t)
summary(calanda19.dmod.t)

# commented out because it takes a long time to run
# LOOCV(calanda19.dmod.t)

# pseudo R2 (marginal = .228, conditional = .497), RMSE = 0.292, and LOOCV RMSE = 0.311

# soil temp
calanda19.dmod.t2 <- lme(sqrt.Disease.l ~ soil_t_centered + richness_scaled + f1_scaled +
                           richness_scaled*soil_t_centered + f1_scaled*soil_t_centered+
                           richness_scaled*moisture_centered + f1_scaled*moisture_centered,
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.dmod.t2)
summary(calanda19.dmod.t2)

# air temp
calanda19.dmod.t3 <- lme(sqrt.Disease.l ~ air_t_centered + richness_scaled + f1_scaled +
                           richness_scaled*air_t_centered + f1_scaled*air_t_centered+
                           richness_scaled*moisture_centered + f1_scaled*moisture_centered,
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.dmod.t3)
summary(calanda19.dmod.t3)

# elevation
calanda19.dmod.e <- lme(sqrt.Disease.l ~ elevation_centered + richness_scaled + f1_scaled +
                          richness_scaled*elevation_centered + f1_scaled*elevation_centered+
                          richness_scaled*moisture_centered + f1_scaled*moisture_centered,
                        random = ~1|Meadow/Site/PlotID,
                        weights = varIdent(form = ~ 1 | Site),
                        data = calanda19.com,
                        control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))
car::Anova(calanda19.dmod.e)
summary(calanda19.dmod.e)

# Using soil-surface temp, test whether it is possible that individual traits are driving the host pace-of-life response----
MuMIn::AICc(calanda19.dmod.t)
# So the model with pace-of-life had AICc 151.1249, marginal R2 = .227, conditional R2 = .497, RMSE = 0.292, and LOOCV RMSE = 0.311

# I'll start with a base model, to make this a bit simpler.
base.model <- lme(sqrt.Disease.l ~ soil_surf_t_centered + richness_scaled +
                    richness_scaled * soil_surf_t_centered +
                    richness_scaled * moisture_centered,
                  random = ~1|Meadow/Site/PlotID,
                  weights = varIdent(form = ~ 1 | Site),
                  data = calanda19.com,
                  control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7))

# First model replaces f1 with Chlorophyll, Leaf_lifespan, Leaf_N, Leaf_P, SLA
traitmod <- update(base.model, .~. + Chlo.z*soil_surf_t_centered + Chlo.z*moisture_centered + 
                     Lon.z*soil_surf_t_centered + Lon.z*moisture_centered +
                     N.z*soil_surf_t_centered + N.z*moisture_centered +
                     P.z*soil_surf_t_centered + P.z*moisture_centered +
                     SLA.z*soil_surf_t_centered + SLA.z*moisture_centered)
car::Anova(traitmod)
MuMIn::AICc(traitmod)
MuMIn::r.squaredGLMM(traitmod)

# commented out becasue this takes a long time to run
# LOOCV(traitmod)

# AIC 220.11, marginal R2 = .257, conditional R2 = 0.294, RMSE = 0.294, and LOOCV RMSE = 0.320

# model with Chlorophyll content only
Chlomod <- update(base.model, .~. + Chlo.z * soil_surf_t_centered + Chlo.z * moisture_centered)
car::Anova(Chlomod)
MuMIn::AICc(Chlomod)
MuMIn::r.squaredGLMM(Chlomod)
# LOOCV(Chlomod)
#AICc 164.57, marginal R2 = .163, conditional R2 = 0.422, RMSE = 0..2941, and LOOCV RMSE = 0.3135

# Model with leaf longevity only
Lonmod <- update(base.model, .~. + Lon.z * soil_surf_t_centered + Lon.z * moisture_centered)
car::Anova(Lonmod)
MuMIn::AICc(Lonmod)
MuMIn::r.squaredGLMM(Lonmod)
# LOOCV(Lonmod)
# AICc 161.70, marginal R2 = .181, conditional R2 = .441, RMSE = 0.2883, and LOOCV RMSE = 0.3092

# Model with leaf nitrogen only
Nmod <- update(base.model, .~. + N.z * soil_surf_t_centered + N.z * moisture_centered)
car::Anova(Nmod)
MuMIn::AICc(Nmod)
MuMIn::r.squaredGLMM(Nmod)
# LOOCV(Nmod)
# AICc 167.70, marginal R2 = .161, conditional R2 = .402, RMSE = 0.2916, and LOOCV RMSE = 0.3127

# Model with leaf phosphorus only
Pmod <- update(base.model, .~. + P.z * soil_surf_t_centered + P.z * moisture_centered)
car::Anova(Pmod)
MuMIn::AICc(Pmod)
MuMIn::r.squaredGLMM(Pmod)
# LOOCV(Pmod)
# AICc 165.38, marginal R2 = 0.168, conditional R2 = 0.433, RMSE = 0.2958, and LOOCV RMSE = 0.3136

# Model with SLA only
SLAmod <- update(base.model, .~. + SLA.z * soil_surf_t_centered + SLA.z * moisture_centered)
car::Anova(SLAmod)
MuMIn::AICc(SLAmod)
MuMIn::r.squaredGLMM(SLAmod)
# LOOCV(SLAmod)
# AICc 159.30, marginal R2 = 0.229, conditional R2 = 0.452, RMSE = 0.294, and LOOCV RMSE = 0.311


# 3. What is the relative contribution of each path? ----
# here, I'm omitting soil moisture, since it is unrelated to the other variables in the data
# and I'm using mean-centered temperature pace of life
model.c19.t <- psem(
  
  lme(sqrt.Disease.l ~ soil_surf_t_centered + plant.richness + f1_scaled + 
        f1_scaled*soil_surf_t_centered,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ soil_surf_t_centered, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(f1_scaled ~ soil_surf_t_centered, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),

  # fit as a linear model at the site level since there is only a single temp measurement per site
  lm(soil_surf_t_centered ~ elevation, 
      data = calanda19.site),
  
  # residual correlations among variables
  f1_scaled %~~% plant.richness)
# note that warning about NAs refers to data that we won't ever use, so i think it's ok.

model.c19.t.fit <- summary(model.c19.t, .progressBar = T)
model.c19.t.fit


# explore parameter estimates for pace of life at varying temperatures
model.c19.low.t <- psem(
  
  lme(sqrt.Disease.l ~ soil_surf_t_lowsd_scaled + plant.richness + f1_scaled + 
        f1_scaled*soil_surf_t_lowsd_scaled,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ soil_surf_t_lowsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(f1_scaled ~ soil_surf_t_lowsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lm(soil_surf_t_lowsd_scaled ~ elevation, 
     data = calanda19.site),
  
  # residual correlations among variables
  f1_scaled %~~% plant.richness)
# note that warning about NAs refers to data that we won't ever use, so i think it's ok.

model.c19.low.t.fit <- summary(model.c19.low.t, .progressBar = T)
model.c19.low.t.fit

# explore parameter estimates for pace of life at varying temperatures
model.c19.high.t <- psem(
  
  lme(sqrt.Disease.l ~ soil_surf_t_highsd_scaled + plant.richness + f1_scaled + 
        f1_scaled*soil_surf_t_highsd_scaled,
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(plant.richness ~ soil_surf_t_highsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lme(f1_scaled ~ soil_surf_t_highsd_scaled, 
      random = ~1|Meadow/Site/PlotID,
      weights = varIdent(form = ~ 1 | Site),
      data = calanda19.com,
      control=list(maxIter=1000, msMaxIter = 1000, tolerance = 1e-7)),
  
  lm(soil_surf_t_highsd_scaled ~ elevation, 
     data = calanda19.site),
  
  # residual correlations among variables
  f1_scaled %~~% plant.richness)
# note that warning about NAs refers to data that we won't ever use, so i think it's ok.

model.c19.high.t.fit <- summary(model.c19.high.t, .progressBar = T)
model.c19.high.t.fit


# remake figure 3

f1_modplot_t <- emmeans::emtrends(calanda19.dmod.t, ~soil_surf_t_centered, var = "f1_scaled", at = list(soil_surf_t_centered = seq(min(calanda19.com$soil_surf_t_centered), max(calanda19.com$soil_surf_t_centered), .1))) %>% 
  data.frame() %>% 
  mutate(in.data = case_when(
    soil_surf_t_centered < min(calanda19.com$soil_surf_t_centered) ~ "No",
    soil_surf_t_centered > max(calanda19.com$soil_surf_t_centered) ~ "No",
    TRUE ~ "Yes"
  ),
  mean_soil_surf_t = soil_surf_t_centered + mean(calanda19.com$mean.soil.surface.t)) %>% 
  ggplot(aes(x = mean_soil_surf_t, y = f1_scaled.trend)) +
  geom_ribbon(fill = "grey", aes(ymin = lower.CL, ymax = upper.CL)) +
  geom_line() + 
  geom_hline(yintercept = 0, lty =2, alpha = .6) +
  geom_rug(data = calanda19.com, aes(x = mean.soil.surface.t, y = f1), sides = "b") +
  labs(y = "Effect of community pace-of-life on disease", x = "Soil-surface Temperature")

# pdf("Figures/F1_moderation_revision.pdf", height = 4, width = 5)
# f1_modplot_t
# dev.off()
# 


# Remake fig s3 with multiple env conditions

# function to construct the plot
calanda.pl <- function(mod){
  
  yv = case_when(
    attributes(mod$terms)$variables[[2]] == "plant.richness" ~ "Host richness",
    attributes(mod$terms)$variables[[2]] == "f1" ~ "\nHost community \npace of life"
  )
  
  xv = case_when(
    attributes(mod$terms)$variables[[3]] == "mean.soil.surface.t" ~ "Soil-surface Temperature",
    attributes(mod$terms)$variables[[3]] == "mean.soil.t" ~ "Soil Temperature",
    attributes(mod$terms)$variables[[3]] == "mean.air.t" ~ "Air Temperature",
    attributes(mod$terms)$variables[[3]] == "elevation" ~ "Elevation"
  )
  
  ggeffects::ggpredict(mod, attributes(mod$terms)$variables[[3]]) %>% plot(rawdata = T) +
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
    labs(y = yv, x = xv, title = "")
}

r.t1 <- calanda.pl(calanda19.rmod.t)

# host richness and soil temp is not significant, so I will just plot the raw data
# r.t2 <- calanda.pl(calanda19.rmod.t2)
# 
r.t2 <- calanda19.com %>% 
  ggplot(aes(x = mean.soil.t, y = plant.richness)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Host richness", x = "Soil Temperature", title = "")


r.t3 <- calanda.pl(calanda19.rmod.t3)
r.e <- calanda.pl(calanda19.rmod.e)

f.t1 <- calanda.pl(calanda19.fmod.t)
f.t2 <- calanda.pl(calanda19.fmod.t2)
f.t3 <- calanda.pl(calanda19.fmod.t3)
f.e <- calanda.pl(calanda19.fmod.e)

# two moisture plots that don't have fits because there was no relationship
r.m <- calanda19.com %>% 
  ggplot(aes(x = mean.moisture, y = plant.richness)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Host richness", x = "Soil Moisture", title = "")

f.m <- calanda19.com %>% 
  ggplot(aes(x = mean.moisture, y = f1)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "\nHost communiy \npace of life", x = "Soil Moisture", title = "")

#  
# pdf("Figures/fig_s3_8_apr_21.pdf", height = 10, width = 6)
# plot_grid(
#   r.t1, f.t1,
#   r.t2, f.t2,
#   r.t3, f.t3,
#   r.e, f.e,
#   r.m, f.m,
#   ncol = 2
#   )
# dev.off()

# temp variables raw correlations with disease
sm.d <- calanda19.com %>% 
  ggplot(aes(x = mean.moisture, y = Disease.l)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  scale_y_sqrt() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Community parasite load", x = "Soil Moisture", title = "")

t1.d <- calanda19.com %>% 
  ggplot(aes(x = mean.soil.surface.t, y = Disease.l)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  scale_y_sqrt() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Community parasite load", x = "Soil surface temperature", title = "")

t2.d <- calanda19.com %>% 
  ggplot(aes(x = mean.soil.t, y = Disease.l)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  scale_y_sqrt() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Community parasite load", x = "Soil temperature", title = "")

t3.d <- calanda19.com %>% 
  ggplot(aes(x = mean.air.t, y = Disease.l)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  scale_y_sqrt() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Community parasite load", x = "Air temperature", title = "")

e.d <- calanda19.com %>% 
  ggplot(aes(x = elevation, y = Disease.l)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  scale_y_sqrt() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(y = "Community parasite load", x = "Elevation", title = "")

plot_grid(
  t1.d,
  t2.d,
  t3.d,
  # e.d,
  sm.d
  )

# abiotic conditions vs elevation plots
sm.e <- calanda19.site %>% 
  ggplot(aes(y = mean.moisture, x = elevation)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  scale_y_sqrt() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(x = "Elevation", y = "Soil Moisture", title = "")

t1.e <- calanda19.site %>% 
  ggplot(aes(y = mean.soil.surface.t, x = elevation)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  geom_smooth(method = "lm", color = "black") + 
  scale_y_sqrt() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(x = "Elevation", y = "Soil surface temperature", title = "")

t2.e <- calanda19.site %>% 
  ggplot(aes(y = mean.soil.t, x = elevation)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  geom_smooth(method = "lm", color = "black") + 
  scale_y_sqrt() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(x = "Elevation", y = "Soil temperature", title = "")

t3.e <- calanda19.site %>% 
  ggplot(aes(y = mean.air.t, x = elevation)) + 
  geom_point(alpha = .35, position = position_jitter(width = .2)) +
  theme_classic() +
  scale_color_brewer() + 
  geom_smooth(method = "lm", color = "black") + 
  scale_y_sqrt() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) +
  labs(x = "Elevation", y = "Air temperature", title = "")

plot_grid(
  t1.e,
  t2.e,
  t3.e,
  sm.e
)
