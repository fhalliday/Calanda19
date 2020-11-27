require(tidyverse)
require(hillR)

# first look at disease data ----
tmp <- read.csv("Disease/Calanda FULL Community Disease survey CLEAN no blanks plotID and speciesID fixed 5.11.2019.csv")
summary(tmp)

# trim white space
tmp <- tmp %>% mutate(
  method = trimws(method, which = "both"),
  SubplotID = trimws(SubplotID, which = "both"),
  SPP = trimws(SPP, which = "both"),
  Phen = trimws(Phen, which = "both"),
  Dam = trimws(Dam, which = "both"))

# look at 90% leaf spot
tmp %>% filter(tot.damage>100)
# 50% chew and 80% mildew. My guess is that the "80% mildew" is for the unchewed portion of the leaf.
# so I vote that we change p mildew to 40%

# there's an instance where leaf spot is unusually high. what is the deal with that?
filter(tmp, Lspot>80)
# this is in subplot A1.6.1 on Potentilla sp

filter(tmp, SubplotID == "A1.6.1", SPP == "Potentilla sp")
filter(tmp, SubplotID == "A1.6.1")
# no obvious cause of this, but 90% seems super unlikely, especially for potentilla

filter(tmp, SPP == "Potentilla sp", Lspot>0)

# I suggest we drop this plant from the dataset

CD.clean <- tmp %>% 
  # create column that is the species id and subplotID combined
  mutate(spp.plot.id = paste(SubplotID, SPP)) %>% 
  # remove the plant with unusually high leafspot damage
  filter(spp.plot.id != "A1.6.1 Potentilla sp") %>%
  # change powdery mildew on  subplot I5.5.1, plant 4, leaf 1, from 80 to 40
  mutate(Pmildew = replace(Pmildew, SubplotID == "I5.5.1" & PlantID == 4 & Leaf == 1, 40)) %>% 
  # Plot I5.8.3 was surveyed using both the "closest" and "touching" method.
  # Only keep the survey that used "touching" for those plots
  filter(!(SubplotID == "I5.8.3" & method == "closest"))

# vegetation
Calanda_19_cover <- read.csv('/Users/fletcherhalliday/Desktop/Laine lab postdoc/Calanda data/Calanda phylogeny/Calanda 19 species list/Calanda_19_Vegetation_data_errors_checked_fh_27_4_2020.csv', stringsAsFactors = F)

# use the same species names that were used to calculate phylogenetic diversity 
c19_spp_plot <- Calanda_19_cover %>% 
  dplyr::select(-Focal.plant..a..No.focal.plant..b., 
                -(Cow.feces:Survey.done.by)) %>% 
  gather(Species, cover, Achillea.millefolium:Viola.sp) %>% 
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
    Species == "Alchemilla_sp.1" ~ "Alchemilla_sp_1",
    Species == "Alchemilla_sp.2" ~ "Alchemilla_sp_2",
    # I just assumed that the arabis and arabidopsis combo unknow was arabis, because that made things easier.
    # need to check with mikko that it's an ok thing to do
    Species == "Arabis_sp...Arabidopsis.sp" ~ "Arabis_sp",
    # simplified this by assuming that this was always a bromus, 
    # because there were not any flowering Koeleria observed on Calanda
    Species == "Bromus_sp...Koeleria.sp" ~ "Bromus_sp",
    Species == "Campanula_schleuzerii...C..rotundifolia" ~ "Campanula_sp_",
    Species == "Campanula_sp" ~ "Campanula_sp_",
    Species == "Carex_sp.1" ~ "Carex_sp_1",
    Species == "Carex_sp.2" ~ "Carex_sp_2",
    Species == "Carex_sp.3" ~ "Carex_sp_3",
    Species == "Carex_sp.4" ~ "Carex_sp_4",
    Species == "Carex_sp.U2.3.1" ~ "Carex_sp_U2.3.1",
    Species == "Carex_sp.U2.8.2" ~ "Carex_sp",
    # just called this a carlina -- need to check with Mikko that this is ok
    Species == "Carlina_sp...Cirsium.sp" ~ "Carlina_sp",
    # just called this a crepis
    Species == "Crepis_sp...Leontodon.sp" ~ "Leontodon_sp",
    # just called this a daucus
    Species == "Daucus_carota...Carum.carvi" ~ "Daucus_carota",
    Species == "Hieracium_pilosella...Hieracium.hoppeanum" ~ "Hieracium_sp",
    Species == "Poa_sp.A1.6.1" ~ "Poa_sp",
    Species == "Poaceae_sp.1" ~ "Poaceae_sp_1",
    Species == "Poaceae_sp.2" ~ "Poaceae_sp",
    Species == "Potentilla_sp" ~ "Potentilla_sp_",
    Species == "Potentilla_sp.I7.1.1" ~ "Potentilla_sp_",
    # just called this Taraxacum
    Species == "Taraxacum_sp...Leontodon.sp" ~ "Taraxacum_sp",
    Species == "Vaccinium_vitis.idaea" ~ "Vaccinium_vitis-idaea",
    TRUE ~ as.character(Species)
  )) %>% 
  # sum across the taxa that are now combined, so I only have one observation of each species in each plot
  group_by(PlotID, Species) %>% 
  summarize(cover = sum(cover)) %>% 
  # replace the _ with " " to match with the names in TRY
  mutate(Species = trimws(gsub("_"," ",Species)))



# create a plot-level disease dataset----
# first, take the vegetation data and put it in long form
# CV.long <- read.csv("Vegetation/Calanda_19_Vegetation_data_errors_checked_fh_27_4_2020.csv") %>% 
#   # long form
#   gather(SPP, Cover, Achillea.millefolium:Viola.sp) %>% 
#   # create plot -species id column so that we can match everything up nicely
#   mutate(
#     SubplotID = trimws(PlotID, which = "both"),
#     SPP = trimws(SPP, which = "both"),
#     spp.plot.id = paste(SubplotID, SPP)) %>% 
#   filter(SPP != "Cow feces",
#          SPP != "Bare ground",
#          SPP != "Rock",
#          SPP != "Bryophytes",
#          SPP != "Litter") 

# calculate total cover to relativise
tot_cover <- c19_spp_plot %>% 
  group_by(PlotID) %>% summarize(tot_cover = sum(cover))

# There are a few taxa where classification was difficult, so we will aggregate them to genus so that we can be confident in our estimates:
unique(c19_spp_plot$Species)
unique(CD.clean$SPP)

# All crepis and Leontodon species can be combined into "Leontodon sp"
# Bromus sp / Koeleria sp in the disease data needs to be changed to "Bromus sp"
# Hieracium pilosella / Hieracium hoppeanum will be combined to 
# Poaceae sp 1 becomes Poaceae sp
# combine the festuca species in this analysis since festuca were difficult to distinguish in the field

CD.cover <- CD.clean %>% 
  mutate(SPP = case_when(
    SPP == "Bromus sp / Koeleria sp" ~ "Bromus sp",
    SPP == "Leontodon autumnalis" ~ "Leontodon sp",
    # SPP == "Carex sp 1" ~ "Carex sp",
    SPP == "Crepis sp / Leontodon sp" ~ "Leontodon sp",
    # SPP == "Carex sp 2" ~ "Carex sp",
    SPP == "Hieracium pilosella / Hieracium hoppeanum" ~ "Hieracium piloselloides",
    # SPP == "Alchemilla sp 1" ~ "Alchemilla sp",
    # SPP == "Alchemilla sp 2" ~ "Alchemilla sp",
    # SPP == "Carex sp 4" ~ "Carex sp",
    SPP == "Vaccinium vitis-idae" ~ "Vaccinium vitis-idaea",
    SPP == "Erica carnea" ~ "Erica herbacea",
    SPP == "Polygonum viviparum" ~ "Persicaria vivipara",
    # SPP == "Carex sp 3" ~ "Carex sp",
    # SPP == "Poaceae sp 1" ~ "Poaceae sp",
    # SPP == "Carex sp U2.3.1" ~ "Carex sp",
    SPP == "Festuca ovina" ~ "Festuca sp",
    SPP == "Festuca rubra" ~ "Festuca sp",
    TRUE ~ as.character(SPP)
  )) %>%
  mutate(Spp_plotID = paste(SPP,SubplotID)) %>% 
  left_join(c19_spp_plot %>% 
              mutate(Species = case_when(
                Species == "Crepis sp" ~ "Leontodon sp",
                Species == "Festuca ovina" ~ "Festuca sp",
                Species == "Festuca rubra" ~ "Festuca sp",
                TRUE ~ as.character(Species)
              )) %>% 
              left_join(tot_cover) %>% 
              mutate(rel_cover = cover/tot_cover,
                     Spp_plotID = paste(Species, PlotID)) %>% 
              select(Spp_plotID, abs_cover = cover, tot_cover, rel_cover)) %>% 
  na.omit()

# something to consider is that very rare species will contribute nothing to community load if they were not observed during the vegetation survey and therefore have a cover of zero

SP.community.load <- CD.cover %>% 
  mutate(Chew.l = Chew * rel_cover, # this is the proportion of leaves with any chewing damage, etc...
    mScrape.l = mScrape * rel_cover,
    iScrape.l = iScrape * rel_cover,
    Mine.l = Mine * rel_cover,
    Gall.l = Gall * rel_cover,
    Window.l = Window * rel_cover,
    Lspot.l = Lspot * rel_cover,
    Rust.l = Rust * rel_cover,
    Pmildew.l = Pmildew * rel_cover,
    Dmildew.l = Dmildew * rel_cover,
    Blight.l = Blight * rel_cover,
    wet.l = wet * rel_cover,
    rhiz.like.l = rhiz.like * rel_cover,
    chlorosis.l = clorosis * rel_cover,
    mildew_like.l = white.mildew.like.thing * rel_cover,
    thrips.l = thrips * rel_cover,
    tent.l = tent * rel_cover,
    scler.like.l = scler.like * rel_cover,
    Phomopsis.l = Phomopsis * rel_cover,
    lf.curl.l = lf.curl * rel_cover,
    necrotic.lesion.l = necrotic.lesion * rel_cover,
    Chlorotic.spots.l = Chlorotic.spots * rel_cover,
    Choking.l = Choking * rel_cover,
    slime_mold_like.l = slime_mold_like * rel_cover, 
    Disease.l = (Lspot + Rust + Pmildew + Dmildew + Blight.l + wet.l + rhiz.like + clorosis + white.mildew.like.thing + scler.like + lf.curl + necrotic.lesion + Chlorotic.spots + Choking + slime_mold_like) * rel_cover,
    Herbivory.l = (Chew + mScrape + iScrape + Mine + Gall + Window + thrips + tent) * rel_cover,
    tot.damage.l = Disease.l + Herbivory.l) %>% 
  select(SubplotID, SPP, Chew.l:tot.damage.l) %>% 
  group_by(SubplotID, SPP) %>% 
  summarize_all(mean) %>%
  ungroup() %>% 
  select(-SPP) %>% 
  group_by(SubplotID) %>% 
  summarize_all(sum)
SP.community.load  

summary(SP.community.load$Disease.l)
hist(SP.community.load$Disease.l)
hist(asinh(SP.community.load$Disease.l))

SP.community.load

# write.csv(SP.community.load, "Disease/Community_disease_load.csv", row.names = F)

# # combine data into a single dataframe for analysis
# Com.load.all <- SP.community.load %>% 
#   mutate(Site = substr(SubplotID, 1,2)) %>% 
#   # add ism data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/ISM_locations_PlotID_fixed.csv", sep = ";") %>% 
#       mutate(SubplotID = paste(substr(meadow,1,1),site,".",plot, ".", subplot, sep = "")) %>% 
#       select(SubplotID, PL, density.type, smallplot, focal.plantago)) %>% 
#   # add site data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/for analysis do not change/Site data 25 June 2019.csv") %>% 
#       mutate(Site = paste(substr(meadow,1,1), site, sep = "")) %>% 
#       select(Site, elevation)
#   ) %>% 
#   # add species richness data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Vegetation_data_plotID_plantID_fixed_with_comma.csv") %>% 
#       mutate(SubplotID = PlotID) %>% 
#       select(SubplotID, Achillea.millefolium:Viola.sp) %>% 
#       # calculate diversity indices
#       transmute(
#         SubplotID = SubplotID,
#         plant.richness = hillR::hill_taxa(.[,-1], q = 0),
#         hill.shannon = hillR::hill_taxa(.[,-1], q = 1),
#         hill.simpson = hillR::hill_taxa(.[,-1], q = 2))) %>% 
#   # transform disease load to normalize it ish
#   mutate(asinh.Disease.l = asinh(Disease.l)) %>% 
#   # add temp and soil moisture data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/TMS_data_clean.csv") %>% 
#       # summarize as mean daily avg, min, and max
#       mutate(Site = substr(site,1,2)) %>% 
#       # first calculate daily min, max, avg
#       group_by(Site, Date) %>% 
#       summarize(min.soil.t = min(soil_temp),
#                 max.soil.t = max(soil_temp),
#                 avg.soil.t = mean(soil_temp),
#                 min.soil.surface.t = min(soil_surface_temp),
#                 max.soil.surface.t = max(soil_surface_temp),
#                 avg.soil.surface.t = mean(soil_surface_temp),
#                 min.air.t = min(air_temp),
#                 max.air.t = max(air_temp),
#                 avg.air.t = mean(air_temp),
#                 min.moisture = min(moisture),
#                 max.moisture = max(moisture),
#                 avg.moisture = mean(moisture)) %>% 
#       ungroup() %>% 
#       # now average across all days that they were in the field
#       group_by(Site) %>% 
#       summarize(mean.min.soil.t = mean(min.soil.t),
#                 mean.max.soil.t = mean(max.soil.t),
#                 mean.soil.t = mean(avg.soil.t),
#                 mean.min.soil.surface.t = mean(min.soil.surface.t),
#                 mean.max.soil.surface.t = mean(max.soil.surface.t),
#                 mean.soil.surface.t = mean(avg.soil.surface.t),
#                 mean.min.air.t = mean(min.air.t),
#                 mean.max.air.t = mean(max.air.t),
#                 mean.air.t = mean(avg.air.t),
#                 mean.min.moisture = mean(min.moisture),
#                 mean.max.moisture = mean(max.moisture),
#                 mean.moisture = mean(avg.moisture)))

# write.csv(Com.load.all, "/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Community_disease_load.csv", row.names = F)


# Make disease dataset at the individual species level----
# Will calculate this both as load (weighted by rel abund) and as %leaf area damaged
SP.community.load.spp <- CD.cover %>% 
  mutate(Chew.l = Chew * rel_cover, # this is the proportion of leaves with any chewing damage, etc...
         mScrape.l = mScrape * rel_cover,
         iScrape.l = iScrape * rel_cover,
         Mine.l = Mine * rel_cover,
         Gall.l = Gall * rel_cover,
         Window.l = Window * rel_cover,
         Lspot.l = Lspot * rel_cover,
         Rust.l = Rust * rel_cover,
         Pmildew.l = Pmildew * rel_cover,
         Dmildew.l = Dmildew * rel_cover,
         Blight.l = Blight * rel_cover,
         wet.l = wet * rel_cover,
         rhiz.like.l = rhiz.like * rel_cover,
         chlorosis.l = clorosis * rel_cover,
         mildew_like.l = white.mildew.like.thing * rel_cover,
         thrips.l = thrips * rel_cover,
         tent.l = tent * rel_cover,
         scler.like.l = scler.like * rel_cover,
         Phomopsis.l = Phomopsis * rel_cover,
         lf.curl.l = lf.curl * rel_cover,
         necrotic.lesion.l = necrotic.lesion * rel_cover,
         Chlorotic.spots.l = Chlorotic.spots * rel_cover,
         Choking.l = Choking * rel_cover,
         slime_mold_like.l = slime_mold_like * rel_cover, 
         Disease.l = (Lspot + Rust + Pmildew + Dmildew + Blight.l + wet.l + rhiz.like + clorosis + white.mildew.like.thing + scler.like + lf.curl + necrotic.lesion + Chlorotic.spots + Choking + slime_mold_like) * rel_cover,
         Herbivory.l = (Chew + mScrape + iScrape + Mine + Gall + Window + thrips + tent) * rel_cover,
         tot.damage.l = Disease.l + Herbivory.l,
         Disease = (Lspot + Rust + Pmildew + Dmildew + Blight.l + wet.l + rhiz.like + clorosis + white.mildew.like.thing + scler.like + lf.curl + necrotic.lesion + Chlorotic.spots + Choking + slime_mold_like),
         Herbivory = (Chew + mScrape + iScrape + Mine + Gall + Window + thrips + tent),
         tot.damage = Disease + Herbivory) %>% 
  select(SubplotID, SPP, Chew:slime_mold_like, Chew.l:tot.damage.l, Disease, Herbivory, tot.damage, abs_cover, tot_cover, rel_cover) %>% 
  group_by(SubplotID, SPP) %>% 
  summarize_all(mean) %>%
  ungroup()

# write.csv(SP.community.load.spp, "Disease/Community_disease_load_species_level.csv")
  
# # combine with other data
# Com.load.all.spp <- SP.community.load.spp %>% 
#   mutate(Site = substr(SubplotID, 1,2)) %>% 
#   # add ism data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/ISM_locations_PlotID_fixed.csv", sep = ";") %>% 
#       mutate(SubplotID = paste(substr(meadow,1,1),site,".",plot, ".", subplot, sep = "")) %>% 
#       select(SubplotID, PL, density.type, smallplot, focal.plantago)) %>% 
#   # add site data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/for analysis do not change/Site data 25 June 2019.csv") %>% 
#       mutate(Site = paste(substr(meadow,1,1), site, sep = "")) %>% 
#       select(Site, elevation)
#   ) %>% 
#   # add species richness data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Vegetation_data_plotID_plantID_fixed_with_comma.csv") %>% 
#       mutate(SubplotID = PlotID) %>% 
#       select(SubplotID, Achillea.millefolium:Viola.sp) %>% 
#       # calculate diversity indices
#       transmute(
#         SubplotID = SubplotID,
#         plant.richness = hillR::hill_taxa(.[,-1], q = 0),
#         hill.shannon = hillR::hill_taxa(.[,-1], q = 1),
#         hill.simpson = hillR::hill_taxa(.[,-1], q = 2))) %>% 
#   # transform disease load to normalize it ish
#   mutate(asinh.Disease.l = asinh(Disease.l),
#          asinh.Disease = asinh(Disease)) %>% 
#   # add temp and soil moisture data
#   left_join(
#     read.csv("/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/TMS_data_clean.csv") %>% 
#       # summarize as mean daily avg, min, and max
#       mutate(Site = substr(site,1,2)) %>% 
#       # first calculate daily min, max, avg
#       group_by(Site, Date) %>% 
#       summarize(min.soil.t = min(soil_temp),
#                 max.soil.t = max(soil_temp),
#                 avg.soil.t = mean(soil_temp),
#                 min.soil.surface.t = min(soil_surface_temp),
#                 max.soil.surface.t = max(soil_surface_temp),
#                 avg.soil.surface.t = mean(soil_surface_temp),
#                 min.air.t = min(air_temp),
#                 max.air.t = max(air_temp),
#                 avg.air.t = mean(air_temp),
#                 min.moisture = min(moisture),
#                 max.moisture = max(moisture),
#                 avg.moisture = mean(moisture)) %>% 
#       ungroup() %>% 
#       # now average across all days that they were in the field
#       group_by(Site) %>% 
#       summarize(mean.min.soil.t = mean(min.soil.t),
#                 mean.max.soil.t = mean(max.soil.t),
#                 mean.soil.t = mean(avg.soil.t),
#                 mean.min.soil.surface.t = mean(min.soil.surface.t),
#                 mean.max.soil.surface.t = mean(max.soil.surface.t),
#                 mean.soil.surface.t = mean(avg.soil.surface.t),
#                 mean.min.air.t = mean(min.air.t),
#                 mean.max.air.t = mean(max.air.t),
#                 mean.air.t = mean(avg.air.t),
#                 mean.min.moisture = mean(min.moisture),
#                 mean.max.moisture = mean(max.moisture),
#                 mean.moisture = mean(avg.moisture)))
# View(Com.load.all.spp)

# write.csv(Com.load.all.spp, "/Users/fletcherhalliday/Dropbox/Laine Group Shared files/Fletcher/Data/Clean data/Community_disease_load_species_level.csv", row.names = F)

