# Summarize Calanda 2019 TMS data

c19_temp <- read.csv("TMS/TMS_data_clean_calibrated.csv") %>% 
  # summarize as mean daily avg, min, and max
  mutate(Site = substr(site,1,2)) %>% 
  # first calculate daily min, max, avg
  group_by(Site, Date) %>% 
  summarize(min.soil.t = min(soil_temp),
            max.soil.t = max(soil_temp),
            avg.soil.t = mean(soil_temp),
            min.soil.surface.t = min(soil_surface_temp),
            max.soil.surface.t = max(soil_surface_temp),
            avg.soil.surface.t = mean(soil_surface_temp),
            min.air.t = min(air_temp),
            max.air.t = max(air_temp),
            avg.air.t = mean(air_temp),
            min.moisture = min(volumetric_soil_moisture),
            max.moisture = max(volumetric_soil_moisture),
            avg.moisture = mean(volumetric_soil_moisture)) %>% 
  ungroup() %>% 
  # now average across all days that they were in the field
  group_by(Site) %>% 
  summarize(mean.min.soil.t = mean(min.soil.t),
            mean.max.soil.t = mean(max.soil.t),
            mean.soil.t = mean(avg.soil.t),
            mean.min.soil.surface.t = mean(min.soil.surface.t),
            mean.max.soil.surface.t = mean(max.soil.surface.t),
            mean.soil.surface.t = mean(avg.soil.surface.t),
            mean.min.air.t = mean(min.air.t),
            mean.max.air.t = mean(max.air.t),
            mean.air.t = mean(avg.air.t),
            mean.min.moisture = mean(min.moisture),
            mean.max.moisture = mean(max.moisture),
            mean.moisture = mean(avg.moisture))

# write.csv(c19_temp, "TMS/Calanda_summarized_tms_2019.csv", row.names = F)
