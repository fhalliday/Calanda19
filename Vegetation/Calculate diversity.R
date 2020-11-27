# Create community-level vegetation dataset
require(tidyverse)
require(hillR)

c19_veg <- read.csv("Vegetation/Calanda_19_Vegetation_data_errors_checked_fh_27_4_2020.csv") %>% 
  select(-Cow.feces:-Survey.done.by) %>% 
  transmute(PlotID = PlotID,
            plant.richness = hillR::hill_taxa(.[,-1], q = 0),
            hill.shannon = hillR::hill_taxa(.[,-1], q = 1),
            hill.simpson = hillR::hill_taxa(.[,-1], q = 2))
c19_veg

write.csv(c19_veg, "Vegetation/Calanda19_diversity.csv", row.names = F)
