# note -- i think i need to omit trees from this database otherwise they might dramatically alter my estimates of phylo div

# Calanda 2019 Phylogeny
# libraries----
require(tidyverse)
require(V.PhyloMaker)
require(ape)
require(picante)
require(geiger)

# build phylogeny----
calanda_19_spp <- read.csv("Vegetation/Calanda_19_spp_list_for_phylogeny.csv")

# generate phylogeny
calanda_19_tree <- phylo.maker(sp.list = calanda_19_spp, tree = GBOTB.extended, nodes = nodes.info.1)    

# plot the phylogeny
plot.phylo(calanda_19_tree$scenario.3, cex = .4, type = "p")

# pdf("Figures/Calanda_19_Phylogeny.pdf", height = 7, width = 7)
# plot.phylo(calanda_19_tree$scenario.3, cex = .4, type = "r")
# dev.off()

# save the phylogeny
# save(calanda_19_tree, file = "Phylogeny/Calanda_19_phylogeny_5_may_2020.RData")

# Calculate phylogenetic diversity----
# branch tips -- need to make sure species names are consistent between this and the cover data
calanda_19_tree$scenario.3$tip.label

Calanda_19_cover <- read.csv('Vegetation/Calanda_19_Vegetation_data_errors_checked_fh_27_4_2020.csv', stringsAsFactors = F)
# make names match the tree species. Note that we need to merge some taxa in order to do that

# species that need to be renamed
Calanda_19_cover %>% 
  dplyr::select(-(Cow.feces:Survey.done.by)) %>% 
  # column_to_rownames(var = "PlotID") %>% 
  gather(Species, cover, Achillea.millefolium:Viola.sp) %>% 
  mutate(Species = str_replace(Species, "\\.", "_")) %>% 
  filter(Species %in% calanda_19_tree$scenario.3$tip.label == FALSE) %>% 
  group_by(Species) %>% 
  summarize(cover = sum(cover))

c19_mat <- Calanda_19_cover %>% 
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
    # Calling this leontodon because crepis seems to be quite rare on calanda
    Species == "Crepis_sp...Leontodon.sp" ~ "Leontodon_sp",
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
  # combine the taxa that are now combined, so I only have one observation of each species in each plot
  group_by(PlotID, Species) %>% 
  summarize(cover = sum(cover)) %>% 
  # drop trees (which never accounted for more than 7% of cover in a plot)
  filter(Species != "Juniperus communis",
         Species != "Larix decidua",
         Species != "Pinus sp") %>% 
  # make wide form
  spread(Species, cover) %>% 
  column_to_rownames(var = "PlotID")

phydist <- cophenetic(calanda_19_tree$scenario.3)

mpd.C19 <-ses.mpd(c19_mat, phydist,
                    null.model='taxa.labels',
                    abundance.weighted=T,
                    runs=1000)

mpd.C19$PlotID = rownames(mpd.C19)

# write.csv(mpd.C19,"Phylogeny/Calanda_19_mpd_no_trees.csv", row.names=F)

plot(mpd.C19$mpd.obs~mpd.C19$ntaxa) #mpd and richness are weakly correlated.
plot(mpd.C19$mpd.obs.z~mpd.C19$ntaxa) #these are essentially uncorrelateed (which is good)

###mpd.obs = mean pairwise distance within a plot -- this is a measure of phylogenetic diversity.
###mpd.obs.z = mean pairwise distance in a plot that is not due to chance alone. in other words, to what degree is a plot more or less diverse than random, given the number of species and their relative abundance.

# plot mpd.obs.z against meadow (which is a rough equivalent of elevation)
mpd.C19 %>% mutate(meadow = substr(PlotID, 1,1),
                   site = substr(PlotID, 1, 2),
                   elevation = case_when(
                     meadow == "I" ~ 600,
                     meadow == "A" ~ 1000,
                     meadow == "N" ~ 1400,
                     meadow == "O" ~ 1600,
                     meadow == "U" ~ 1700
                   )) %>% 
  ggplot(aes(x = elevation, y = mpd.obs.z, group = site)) + 
  # geom_point(position = position_jitter(), shape = 1, alpha = .5) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, lty = 2)
# on average, plots apear to be phylogenetically clustered, with more sites showing phylogenetic clustering as elevation increases
# I would have expected something hump-shaped here if there was mixing of high- and low- adapted species at intermediate elevations

# weak increase in phylogenetic diversity with elevation, probably attributable to increase in taxonomic diversity
mpd.C19 %>% mutate(meadow = substr(PlotID, 1,1),
                   site = substr(PlotID, 1, 2),
                   elevation = case_when(
                     meadow == "I" ~ 600,
                     meadow == "A" ~ 1000,
                     meadow == "N" ~ 1400,
                     meadow == "O" ~ 1600,
                     meadow == "U" ~ 1700
                   )) %>% 
  ggplot(aes(x= elevation, y = mpd.obs)) + 
  geom_point(position = position_jitter(), shape = 1) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2)) + 
  geom_hline(yintercept = 0, lty = 2)

