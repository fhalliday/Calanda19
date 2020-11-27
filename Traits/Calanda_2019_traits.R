# Generate calanda 2019 functional composition measurements
# libraries----
require(tidyverse)
require(vroom)
require(umx)
require(ggbeeswarm)

# Hieracium sp
# Hieracium lactucella
# Pulsatilla sp



#Step one get the TRY trait data----
# I have to use a local file here, because the TRY datasets are too large for github storage
try.species <- read.delim("Traits/TryAccSpecies.txt") %>% 
  # create a column for genus
  # this will allow me to aggregate at the genus level for unknown taxa
  mutate(Accgenus = gsub(" .*$", "", AccSpeciesName))

Calanda_19_cover <- read.csv('Vegetation/Calanda_19_Vegetation_data_errors_checked_fh_27_4_2020.csv', stringsAsFactors = F)

# use the same species names that were used to calculate phylogenetic diversity 
# with the exception that Hieraceum becomes "pilosella"
c19_spp_plot <- Calanda_19_cover %>% 
  dplyr::select(-(Cow.feces:Survey.done.by)) %>% 
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
    Species == "Hieracium_pilosella...Hieracium.hoppeanum" ~ "Pilosella_sp",
    Species == "Poa_sp.A1.6.1" ~ "Poa_sp",
    Species == "Poaceae_sp.1" ~ "Poaceae_sp",
    Species == "Poaceae_sp.2" ~ "Poaceae_sp",
    Species == "Potentilla_sp" ~ "Potentilla_sp_",
    Species == "Potentilla_sp.I7.1.1" ~ "Potentilla_sp_",
    # just called this Taraxacum
    Species == "Taraxacum_sp...Leontodon.sp" ~ "Taraxacum_sp",
    Species == "Vaccinium_vitis.idaea" ~ "Vaccinium_vitis-idaea",
    # the try database has more data for the Pilosella name alterantes of Hieraceum, so I'm going to change them here to acommodate that
    Species == "Hieracium_lactucella" ~ "Pilosella_lactucella",
    Species == "Hieracium_piloselloides" ~ "Pilosella_officinarum",
    Species == "Hieracium sp" ~ "Pilosella_sp",
    TRUE ~ as.character(Species)
  )) %>% 
  # sum across the taxa that are now combined, so I only have one observation of each species in each plot
  group_by(PlotID, Species) %>% 
  summarize(cover = sum(cover)) %>% 
  # replace the _ with " " to match with the names in TRY
  mutate(Species = trimws(gsub("_"," ",Species)))

# I will probably end up dropping trees from this (and the phylogeny)
# which of these taxa are trees and how much cover do trees tend to have in plots?
unique(c19_spp_plot$Species)

c19_spp_plot %>% 
  filter(Species == "Juniperus communis"| Species == "Larix decidua"| Species == "Pinus sp") %>% 
  group_by(PlotID) %>% 
  summarize(cover= sum(cover)) %>% 
  summary()
# trees never accounted for more than 7% cover in a plot

c19_spp_all <- c19_spp_plot %>% 
  group_by(Species) %>% summarize(n()) %>% 
  transmute(Species = Species,
            # flag genus-level unknowns
            is.unknown = ifelse(trimws(gsub(".* ", "", Species))=="sp", 1,0),
            # create genus column for unknown species
            genus = gsub(" .*$", "", Species))

c19_spp <- c19_spp_all %>% 
  filter(is.unknown == 0)

c19_genus <- c19_spp_all %>% 
  filter(is.unknown == 1)

# species in our named species dataset that have matches in TRY
try.c <- try.species %>% filter(AccSpeciesName %in% c19_spp$Species)
dim(try.c)

# here are the species in our named species dataset that do not have any matches in TRY
try.not.c <- c19_spp %>% filter(!Species %in% try.species$AccSpeciesName)
try.not.c

# these will get added to the genus-level dataset
c19_genus2 <- rbind(c19_genus, try.not.c)

# create a list of every species observed on calanda from these genera
all_calanda_species <- read.csv("Vegetation/All_Calanda_Species.csv", sep = ",") %>% 
  # not sure why the formatting for this name is all messy, but I can fix it here
  mutate(Species = case_when(
    Species == "Vaccinium\xcavitis-idaea\xcasubsp.\xcaminus" ~ "Vaccinium vitis-idaea subsp. minus",
    TRUE ~ as.character(Species)
  )) %>% 
  # note For one taxon, Rosa sp. We did not have any observations that were categorized to the species level, so we will also include every Rosa species listed on infoflora.ch as candidates for functional trait analysis.
  rbind(
    read.csv("Vegetation/All_CH_Rosa_Species.csv")
  )

all.c <- all_calanda_species %>% filter(Genus %in% c19_genus2$genus)

# species in the full calanda dataset that have matches in TRY
try.c2 <- try.species %>% filter(AccSpeciesName %in% all.c$Species)

# are there any genera left for which we do not have any representative taxa in TRY?
all.c %>% filter(!Genus %in% try.c2$Accgenus)
# no -- so now we have representatives for most species, and when we couldn't find a representative for a species,
# we took the genus-level mean for all species that have been observed on Calanda.
# 

try.accessions <- rbind(try.c, try.c2)

# create a dataframe that I can use to combine the trait data with the species data
c1 <- c19_spp %>% left_join(try.c, by = c("Species" = "AccSpeciesName")) %>% drop_na() %>% 
  transmute(Species = Species, name.in.try = Species, AccSpeciesID = AccSpeciesID)

c2 <- try.c2 %>% mutate(genus = gsub(" .*$", "", AccSpeciesName)) %>% 
  left_join(c19_genus2) %>% 
  transmute(Species = Species, name.in.try = AccSpeciesName, AccSpeciesID = AccSpeciesID)

calanda.try.lookup <- rbind(c1, c2) %>% 
  mutate(
    Species = case_when(
      # need to fix the Anemone pulstatilla which is used as a member of the "Puslatilla sp" group
      name.in.try == "Anemone pulsatilla" ~ "Pulsatilla sp",
      TRUE ~ as.character(Species)
    )
  ) %>% 
  # add in the pilosella species estimates for hieraceum sp (they are synonyms)
  rbind(
    c2 %>% filter(Species == "Pilosella sp") %>% 
      mutate(Species = "Hieracium sp")
  )


# now take the full try databse and select just the taxa that I want from try.accessions
# Note that this database is too big to store on github, so it's stored elsewhere on my computer.
TRYdata <- vroom("Traits/9136.txt",
      col_select = c(AccSpeciesName, AccSpeciesID, TraitID, TraitName, DataName, StdValue, UnitName)) %>%
  filter(AccSpeciesName %in% try.accessions$AccSpeciesName)
#
# write.csv(TRYdata, "Traits/Calanda_TRY_data.csv", row.names = F)

TRYdata <- read.csv("Traits/Calanda_TRY_data.csv")

# step 2 select traits of interest and combine trait data with dataset----
TRY.traits <- TRYdata %>% full_join(calanda.try.lookup) %>% 
  # select just the traits that I want to work with to start
  filter(
    # SLA
    TraitID == 3115| TraitID == 3116| TraitID == 3117|
    # CN
      TraitID == 146 | 
    # Leaf chlorophyll
      TraitID == 164 |
    # leaf lifespan (longevity; 12)
      TraitID == 12 |
    # leaf n per dry mass (14)
      TraitID == 14 |
    # leaf phosphorus per dry mass (15)
      TraitID == 15 & DataName == "Leaf phosphorus content per dry mass (Pmass)" |
    # plant height (vegetative; 3106)
      TraitID == 3106 |
    # seed dry mass (26)
      TraitID == 26 |
    # Photosynthetic rate (40)
      TraitID == 40) %>% 
  # rename traits to allow me to combine the three different measurements of SLA
  mutate(Trait = case_when(
    TraitID == 3115 ~ "SLA",
    TraitID == 3116 ~ "SLA",
    TraitID == 3117 ~ "SLA",
    TraitID == 146 ~ "CN",
    TraitID == 164 ~ "Chlorophyll",
    TraitID == 12 ~"Leaf_lifespan",
    TraitID == 14 ~ "Leaf_N",
    TraitID == 15 ~ "Leaf_P",
    TraitID == 3106 ~ "Height",
    TraitID == 26 ~ "Seed_mass",
    TraitID == 40 ~ "Photosynthetic_rate")) %>% 
  # Create averages that match up with Calanda data
  group_by(Species, Trait) %>% 
  summarize(StdValue = mean(as.numeric(StdValue), na.rm = T)) %>% 
  spread(Trait, StdValue) %>% 
  # drop trees (which never accounted for more than 7% of cover in a plot)
  filter(Species != "Juniperus communis",
         Species != "Larix decidua",
         Species != "Pinus sp")

TRY.traits  
# is there a common axis of traits?
# using FIML to deal with missing data
TRY.traits.fa <- TRY.traits %>% 
  ungroup() %>% 
  # get rid of seed mass since and height, since they are not leaf traits
  # get rid of CN because it's collinear with Leaf N and has poorer coverage in the data
  select(-Species, -Seed_mass, - Height, -CN) %>% 
  data.frame()

umxEFA(TRY.traits.fa, 
       factors = 1, 
       scores = 'Regression',
       minManifests = 2,
       rotation = "varimax") %>% 
  summary(., refModels = mxRefModels(., run = TRUE))

# suggests that photosynthetic rate is not a good factor, so I'll drop that
TRY.traits.fa2 <- TRY.traits %>% 
  ungroup() %>% 
  select(-Species, -Seed_mass, -CN, -Height, -Photosynthetic_rate) %>% 
  data.frame()

umxEFA(TRY.traits.fa2, 
       factors = 1, 
       scores = 'Regression',
       minManifests = 2,
       rotation = "varimax") %>% 
  summary(., refModels = mxRefModels(., run = TRUE))


TRY.fa <- umxEFA(TRY.traits.fa2, 
                 factors = 1, 
                 scores = 'Regression',
                 minManifests = 2,
                 rotation = "varimax",
                 return = "loadings")

umxEFA(TRY.traits.fa2, 
       factors = 1, 
       scores = 'Regression',
       minManifests = 2,
       rotation = "varimax") %>% 
  summary(., refModels = mxRefModels(., run = TRUE))

TRY.fa$F1

Calanda.traits <- cbind(TRY.traits, f1 = (TRY.fa$F1))
# write.csv(Calanda.traits, "Traits/Calanda_2019_plant_traits_no_trees.csv", row.names = F)

# visualize calanda traits at species level
pace.plot <- c19_spp_plot %>% left_join(Calanda.traits) %>% 
  filter(cover > 0) %>% 
  group_by(Species) %>% 
  summarize(mean_cover = mean(cover),
            low_cover = mean(cover) - sd(cover),
            up_cover = mean(cover) + sd(cover),
            f1 = mean(f1)) %>% 
  ungroup() %>% 
  # mutate(Species = str_replace(Species, "\\s+", "\n")) %>% 
  ggplot(aes(x = f1, y = mean_cover, ymin = low_cover, ymax = up_cover, color = Species)) + 
  # geom_errorbar(width = 0, color = "gray") +
  geom_point(size = 3, alpha = .3) +
  geom_text(aes(label = Species), check_overlap = TRUE, hjust=-0.05, size = 6) +
  # geom_text(aes(label=ifelse(abs(f1) > .8, as.character(Species),'')),hjust=0,vjust=0, check_overlap = TRUE) +
  scale_y_continuous(trans = "pseudo_log", limits = c(0, 40)) +
  theme(legend.position = "non") +
  scale_x_continuous(limits = c(-2,4)) +
  labs(y = "Host absolute abundance (% cover)", x = "Host pace-of-life")
pace.plot

# pdf("Figures/host_pace_of_life.pdf", height = 7, width = 14)
# pace.plot
# dev.off()


Calanda.traits %>% filter(f1 > 2)
# calculate community-weighted mean traits
C19_plot_traits <- c19_spp_plot %>% left_join(Calanda.traits) %>% 
  # change shape of dataset to make this easier
  gather(Trait, Trait_value, Chlorophyll:f1) %>% 
  # multiply each trait by the cover of the species in the plot (divided by 100)
  mutate(abs.trait = Trait_value * (cover/100)) %>% 
  # sum across all plants in each plot for each trait 
  group_by(PlotID, Trait) %>% 
  summarize(total.cover = sum(cover, na.rm = T)/100,
            abs.trait = sum(abs.trait, na.rm = T)) %>% 
  mutate(rel.trait = abs.trait/total.cover) %>% 
  select(-total.cover, -abs.trait) %>% 
  spread(Trait, rel.trait)

# Calculate community functional diversity (Rao Q?)
# Functional diversity calculation can not handle multiple taxa with the same values, so those taxa need to be combined.
Calanda.traits.m <- Calanda.traits %>% 
  ungroup() %>% 
  # filter(
  #   # drop redundant taxa (which will be combined in the veg dataset)
  #   Species != "Agrostis capillaris", #combine with Agrostis sp
  #   Species != "Agrostis schraderiana", #combine with Agrostis sp
  #   Species != "Arabis ciliata", #combine with Arabis sp
  #   Species != "Hieracium sp", #combine with Pilosella sp
  #   Species != "Ranunculus tuberosis", #combine with Ranunculus sp
  #   Species != "Leontodon hispidus") %>% #combine with Leontodon sp
  column_to_rownames(var = "Species") %>% 
  # select height, seed mass, and Pace-of-life, to represent leaf (pace of life), h and s. 
  select(-f1) %>% 
  as.matrix()

Calanda.abund <- c19_spp_plot %>% ungroup() %>% 
  # drop trees and taxa that could not be identified to genus (generally quite rare)
  filter(Species != "Juniperus communis", 
         Species != "Larix decidua",
         Species != "Pinus sp",
         Species != "Asteraceae sp",
         Species != "Orchidaceae sp",
         Species != "Poaceae sp"
         ) %>% 
  # mutate(Species2 = case_when(
  #   Species == "Agrostis capillaris" ~ "Agrostis sp",
  #   Species == "Agrostis schraderiana" ~ "Agrostis sp",
  #   Species == "Arabis ciliata" ~ "Arabis sp",
  #   Species == "Hieracium sp" ~ "Pilosella sp",
  #   Species == "Ranunculus tuberosis" ~ "Ranunculus sp",
  #   Species == "Leontodon hispidus" ~ "Leontodon sp",
  #   TRUE ~ as.character(Species))) %>% 
  # group_by(PlotID, Species2) %>% 
  # summarize(cover = sum(cover)) %>% 
  # select(PlotID, Species = Species2, cover) %>% 
  # ungroup() %>% 
  spread(Species, cover) %>% 
  column_to_rownames(var = "PlotID") %>% 
as.matrix()


Calanda.FD <- FD::dbFD(Calanda.traits.m,
          Calanda.abund, corr = "none")

Calanda.fd <- data.frame(
  FRic = Calanda.FD$FRic,
  FEve = Calanda.FD$FEve,
  FDiv = Calanda.FD$FDiv,
  RaoQ = Calanda.FD$RaoQ
) %>% 
  rownames_to_column(var = "PlotID")

Calanda.fd

C19_plot_traits <- C19_plot_traits %>% 
  left_join(Calanda.fd)

# write.csv(C19_plot_traits, "Traits/Calanda_2019_community_traits_no_trees.csv", row.names = F)

# visualize all traits with elevation
C19_plot_traits %>% ungroup() %>% 
  mutate(meadow = substr(PlotID, 1,1),
         site = substr(PlotID, 1, 2),
         elevation = case_when(
           meadow == "I" ~ 600,
           meadow == "A" ~ 1200,
           meadow == "N" ~ 1400,
           meadow == "O" ~ 1600,
           meadow == "U" ~ 1700)) %>% 
  mutate(meadow = fct_reorder(meadow,elevation)) %>% 
  gather(Trait, tv, Chlorophyll:SLA) %>% 
  ggplot(aes(x=meadow, y = tv)) +
  facet_wrap(~Trait, scales = "free_y") + 
  geom_point(position = position_quasirandom(dodge.width = 0.5), shape = 1, alpha = .5) + 
  stat_summary(geom = "point", fun.y = mean, position = position_dodge(width = 0.5), shape = "-", size = 10)

# correlations among leaf traits, height, and seed mass
C19_plot_traits %>% ungroup() %>% select(f1, Height, Seed_mass ) %>% pairs()
# some correlation but nothing particularly strong
C19_plot_traits %>% ungroup() %>% select(f1, Height, Seed_mass ) %>% cor()
# weakly positively correlated

# Visualize functional diversity with elevation
C19_plot_traits %>% ungroup() %>% 
  mutate(meadow = substr(PlotID, 1,1),
         site = substr(PlotID, 1, 2),
         elevation = case_when(
           meadow == "I" ~ 600,
           meadow == "A" ~ 1000,
           meadow == "N" ~ 1400,
           meadow == "O" ~ 1600,
           meadow == "U" ~ 1700)) %>%
  ggplot(aes(x= elevation, y = FRic)) + 
  geom_point(position = position_jitter(), shape = 1) +
  geom_smooth(method = "lm", formula = 'y~ poly(x,2)')

C19_plot_traits %>% ungroup() %>% 
  mutate(meadow = substr(PlotID, 1,1),
         site = substr(PlotID, 1, 2),
         elevation = case_when(
           meadow == "I" ~ 600,
           meadow == "A" ~ 1000,
           meadow == "N" ~ 1400,
           meadow == "O" ~ 1600,
           meadow == "U" ~ 1700)) %>%
  ggplot(aes(x= elevation, y = FEve)) + 
  geom_point(position = position_jitter(), shape = 1) +
  geom_smooth(method = "lm", formula = 'y~ poly(x,2)')

C19_plot_traits %>% ungroup() %>% 
  mutate(meadow = substr(PlotID, 1,1),
         site = substr(PlotID, 1, 2),
         elevation = case_when(
           meadow == "I" ~ 600,
           meadow == "A" ~ 1000,
           meadow == "N" ~ 1400,
           meadow == "O" ~ 1600,
           meadow == "U" ~ 1700)) %>%
  ggplot(aes(x= elevation, y = FDiv)) + 
  geom_point(position = position_jitter(), shape = 1) +
  geom_smooth(method = "lm", formula = 'y~ poly(x,2)')

C19_plot_traits %>% ungroup() %>% 
  mutate(meadow = substr(PlotID, 1,1),
         site = substr(PlotID, 1, 2),
         elevation = case_when(
           meadow == "I" ~ 600,
           meadow == "A" ~ 1000,
           meadow == "N" ~ 1400,
           meadow == "O" ~ 1600,
           meadow == "U" ~ 1700)) %>%
  ggplot(aes(x= elevation, y = RaoQ)) + 
  geom_point(position = position_jitter(), shape = 1) +
  geom_smooth(method = "lm", formula = 'y~ poly(x,2)')
