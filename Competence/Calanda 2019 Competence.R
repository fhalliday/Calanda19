# Calanda_2019_community_competence
calanda.spec.competence <- read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Community_disease_load_species_level_fixed16122019.csv") %>% 
  # remove temp and moisture columns
  dplyr::select(-mean.min.soil.t,
                -mean.max.soil.t,
                -mean.soil.t,
                -mean.min.soil.surface.t,
                -mean.max.soil.surface.t,
                -mean.soil.surface.t,
                -mean.min.air.t,
                -mean.max.air.t,
                -mean.air.t,
                -mean.min.moisture,
                -mean.max.moisture,
                -mean.moisture) %>% 
  group_by(SPP) %>% 
  summarize(
    Disease = mean(Disease),
    rel.cover = mean(rel.cover),
    perc.plots = length(unique(SubplotID)) / 220) %>% 
  mutate(
    # diluters and amplifiers are present in at least 5% of plots
    # scored based on how diseased they get as a species, on average
    disease2 = case_when(
      perc.plots >= 0.05 ~ Disease,
      TRUE ~ NA_real_
    ),
    competence.rank = cume_dist(disease2),
    # name them "amplifiers" if they are in the top 20% of diseased and diluters if they are in the bottom 20%
    # but note that this definition doesn't really work for diluters who tend to be rare species
    competence.grouped = factor(case_when(
      competence.rank >= 0.8 ~ "amplifier",
      competence.rank <= 0.2 ~ "diluter",
      perc.plots < 0.05 ~ "NA",
      TRUE ~ "intermediate competence")))

calanda.competence <- read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Community_disease_load_species_level_fixed16122019.csv") %>% 
  left_join(calanda.spec.competence %>% 
              dplyr::select(SPP, competence.rank, competence.grouped)) %>% 
  # in each SubplotID 
  group_by(SubplotID, competence.grouped) %>% 
  summarize(
    rel.abund = mean(rel.cover)) %>% 
  ungroup() %>% 
  filter(competence.grouped != "NA") %>% 
  mutate(competence.grouped = case_when(
    competence.grouped == "amplifier" ~ "amplifier_rel_abund",
    competence.grouped == "diluter" ~ "diluter_rel_abund",
    competence.grouped == "intermediate competence" ~ "intermediate_rel_abund")) %>% 
  spread(competence.grouped, c(rel.abund)) %>% 
  mutate(amplifier_rel_abund = replace_na(amplifier_rel_abund, 0),
         diluter_rel_abund = replace_na(diluter_rel_abund, 0)) %>% 
  # calculate absolute abundance
  left_join(
    calanda19.spec %>% 
      left_join(calanda.spec.competence %>% 
                  dplyr::select(SPP, competence.rank, competence.grouped)) %>% 
      # in each SubplotID 
      group_by(SubplotID, competence.grouped) %>% 
      summarize(
        abs.abund = mean(abs.cover)) %>% 
      ungroup() %>% 
      filter(competence.grouped != "NA") %>% 
      mutate(competence.grouped = case_when(
        competence.grouped == "amplifier" ~ "amplifier_abs_abund",
        competence.grouped == "diluter" ~ "diluter_abs_abund",
        competence.grouped == "intermediate competence" ~ "intermediate_abs_abund")) %>% 
      spread(competence.grouped, c(abs.abund)) %>% 
      mutate(amplifier_abs_abund = replace_na(amplifier_abs_abund, 0),
             diluter_abs_abund = replace_na(diluter_abs_abund, 0)) 
  )

write.csv(calanda.competence, row.names = F, "/Users/fletcherhalliday/Desktop/Laine lab postdoc/Calanda data/Calanda 2019 competence/Calanda_2019_community_competence.csv")


# Some species level analyses----
calanda19.spec <- calanda19.spec <- read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Community_disease_load_species_level_fixed16122019.csv") %>% 
  # remove temp and moisture columns
  dplyr::select(-mean.min.soil.t,
                -mean.max.soil.t,
                -mean.soil.t,
                -mean.min.soil.surface.t,
                -mean.max.soil.surface.t,
                -mean.soil.surface.t,
                -mean.min.air.t,
                -mean.max.air.t,
                -mean.air.t,
                -mean.min.moisture,
                -mean.max.moisture,
                -mean.moisture) %>% 
  # add calibrated temp and soil moisture data
  left_join(
    read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/TMS_data_clean_calibrated.csv") %>% 
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
                mean.moisture = mean(avg.moisture)))

Calanda.traits <- read.csv("/Users/fletcherhalliday/Desktop/Laine lab postdoc/Calanda data/TRY data/Calanda_2019_plant_traits_no_trees.csv")

# who are the most and least competent hosts?
calanda.spec.competence %>% filter(competence.grouped == "amplifier" | competence.grouped == "diluter") %>% 
  arrange(competence.grouped)

# how diseased are the most competent hosts and does that change with elevation?
calanda19.spec %>% 
  left_join(calanda.spec.competence %>% 
              dplyr::select(SPP, competence.grouped), by = c("SPP" = "SPP")) %>% 
  # filter the species-level dataset to just the amplifiers
  filter(competence.grouped == "amplifier") %>% 
  mutate(Meadow = factor(str_sub(Site, 1,1), levels = c("I", "A", "N", "O", "U"))) %>% 
  ggplot(aes(x = SPP, y = asinh.Disease, color = SPP)) + 
  facet_wrap(~Meadow, scales = "free") +
  geom_boxplot() + 
  labs(x = "Host species", y = "Mean leaf area damaged (asinh-transformed)\n\n")+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 1),
        legend.position = c(1,.1))

# density of amplifiers?
calanda19.spec %>% 
  left_join(calanda.spec.competence %>% 
              dplyr::select(SPP, competence.grouped), by = c("SPP" = "SPP")) %>% 
  # filter the species-level dataset to just the amplifiers
  filter(competence.grouped == "amplifier") %>% 
  mutate(Meadow = factor(str_sub(Site, 1,1), levels = c("I", "A", "N", "O", "U"))) %>% 
  ggplot(aes(x = SPP, y = abs.cover, color = SPP)) + 
  facet_wrap(~Meadow, scales = "free") +
  geom_boxplot() + 
  labs(x = "Host species", y = "Host density\n\n")+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 1),
        legend.position = c(1,.1))

# What traits are associated with being an amplifier vs a diluter?----
calanda.trait.by.competence <- calanda19.spec %>% 
  left_join(calanda.spec.competence %>% 
              dplyr::select(SPP, competence.rank, competence.grouped)) %>% 
  
  # need to edit the dataset so that the species names match up properly with the traits database
  mutate(Species = SPP) %>% 
  mutate(Species = str_replace(Species, "\\.", "_")) %>% 
  # removed unknown Dicots from analysis because I can not assign them to any meaningful taxa
  filter(Species != "Dicot_sp",
         Species != "Dicot_sp.1.N1.5.2",
         Species != "Dicot_sp.2",
         Species != "Dicot_sp.2.N1.5.2",
         Species != "Dicot_sp.5",
         Species != "Dicot_sp.N4.1.1",
         Species != "Dicot_sp.U2.3.3") %>% 
  
  mutate(Species = case_when(
    Species == "Agrostis_capillaris...A..schraderiana" ~ "Agrostis_sp",
    Species == "Agrostis_schraderiana...A..stolonifera" ~ "Agrostis_sp",
    Species == "Alchemilla_sp.1" ~ "Alchemilla_sp",
    Species == "Alchemilla_sp.2" ~ "Alchemilla_sp",
    # I just assumed that the arabis and arabidopsis combo unknow was arabis, because that made things easier.
    # need to check with mikko that it's an ok thing to do
    Species == "Arabis_sp...Arabidopsis.sp" ~ "Arabis_sp",
    # simplified this by assuming that this was always a bromus, 
    # because there were not any flowering Koeleria observed on Calanda
    Species == "Bromus_sp...Koeleria.sp" ~ "Bromus_sp",
    Species == "Campanula_schleuzerii...C..rotundifolia" ~ "Campanula_sp_",
    Species == "Campanula_sp" ~ "Campanula_sp_",
    Species == "Carex_sp.1" ~ "Carex_sp",
    Species == "Carex_sp.2" ~ "Carex_sp",
    Species == "Carex_sp.3" ~ "Carex_sp",
    Species == "Carex_sp.4" ~ "Carex_sp",
    Species == "Carex_sp.U2.3.1" ~ "Carex_sp",
    Species == "Carex_sp.U2.8.2" ~ "Carex_sp",
    # just called this a carlina -- need to check with Mikko that this is ok
    Species == "Carlina_sp...Cirsium.sp" ~ "Carlina_sp",
    # just called this a crepis
    Species == "Crepis_sp...Leontodon.sp" ~ "Crepis_sp",
    # just called this a daucus
    Species == "Daucus_carota...Carum.carvi" ~ "Daucus_carota",
    Species == "Hieracium_pilosella...Hieracium.hoppeanum" ~ "Hieracium_sp",
    Species == "Poa_sp.A1.6.1" ~ "Poa_sp",
    Species == "Poaceae_sp.1" ~ "Poaceae_sp",
    Species == "Poaceae_sp.2" ~ "Poaceae_sp",
    Species == "Potentilla_sp" ~ "Potentilla_sp_",
    Species == "Potentilla_sp.I7.1.1" ~ "Potentilla_sp_",
    # just called this Taraxacum
    Species == "Taraxacum_sp...Leontodon.sp" ~ "Taraxacum_sp",
    Species == "Vaccinium_vitis.idaea" ~ "Vaccinium_vitis-idaea",
    TRUE ~ as.character(Species)
  )) %>% 
  # sum across the taxa that are now combined, so I only have one observation of each species in each plot
  group_by(SubplotID, Site, Species, elevation, plant.richness, hill.shannon, hill.simpson, mean.min.soil.t, mean.max.soil.t, mean.soil.t, mean.min.soil.surface.t, mean.max.soil.surface.t, mean.soil.surface.t, mean.min.air.t, mean.max.air.t, mean.air.t, mean.min.moisture, mean.max.moisture, mean.moisture, competence.rank, competence.grouped, tot.cover) %>% 
  dplyr::select(-SPP, -asinh.Disease.l, -asinh.Disease, -rel.cover, -PL:-focal.plantago) %>% 
  summarize_all(sum) %>% ungroup() %>% 
  # replace the _ with " " to match with the names in TRY
  mutate(Species = trimws(gsub("_"," ",Species)),
         asinh.Disease.l = asinh(Disease.l),
         asinh.Disease = asinh(Disease),
         rel.cover = abs.cover/tot.cover,
         is.amplifier = ifelse(competence.grouped == "amplifier", "Amplifier", "Not amplifier")) %>% 
  # add traits
  left_join(Calanda.traits)

calanda.trait.by.competence %>% 
  gather(trait, trait.value, Chlorophyll:f1) %>% 
  ggplot(aes(x = is.amplifier, y = trait.value)) +
  facet_wrap(~trait, scales = "free_y") + 
  # geom_point(alpha = .5, shape = 1, position = position_quasirandom()) + 
  geom_boxplot()

# at the species-level, amplifiers don't appear to be much different than non-amplifying species



# do trait distributions of the most competent hosts change along the elevational gradient?----


