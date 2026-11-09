# Fig. 1A: sampling locations of the 1,229 wild individuals analysed (ANN 589, ARG 297, PET 144, FAL 199).
# One point per population. Only populations with at least one individual in the four sample lists used for
# the analyses (ann/arg/pet/fal.sample) are plotted, so the map matches the analysed data:
#   - ANN includes the H. annuus subsp. texanus populations (ann_TX), which are part of the 589 ANN individuals;
#   - populations not used in the analyses (ann_AZ, H. niveus subsp. canescens, and five unused PET/FAL
#     populations) are not shown.
# The dashed outlines and the "Local / Species-wide / Inter-specific sweeps" labels were added in Illustrator.
#
# Run from this folder (00_sampling/):  Rscript five_species_distribution.R
# Inputs in data/: States_shapefile/ (US state borders), sunflower_4_sps_pop_info.csv (individual, population,
# taxon, latitude, longitude; from Todesco et al. 2020), and ann/arg/pet/fal.sample (individuals analysed).
library(sf)
library(dplyr)
library(ggplot2)

# US state borders (lower 48)
states <- st_read('data/States_shapefile/States_shapefile.shp', quiet = TRUE) %>%
  filter(!State_Code %in% c('AK', 'HI'))

# individual -> population, with coordinates
pop_info <- read.csv('data/sunflower_4_sps_pop_info.csv')

taxa <- c(ann = 'H. annuus', arg = 'H. argophyllus',
          fal = 'H. petiolaris subsp. fallax', pet = 'H. petiolaris subsp. petiolaris')
pops <- bind_rows(lapply(names(taxa), function(sp) {
  ind <- read.table(paste0('data/', sp, '.sample'))$V1
  pop_info %>% filter(ind_ID %in% ind) %>%
    group_by(Population.ID, Latitude, Longitude) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(taxon = taxa[[sp]])
}))

# check: individuals and populations per taxon (expect 589/297/199/144 individuals; 58/30/21/15 populations)
print(pops %>% group_by(taxon) %>% summarise(individuals = sum(n), populations = n()))

pops$taxon <- factor(pops$taxon, levels = taxa)
labels <- c(expression(italic('H. annuus')), expression(italic('H. argophyllus')),
            expression(italic('H. petiolaris') ~ 'subsp.' ~ italic('fallax')),
            expression(italic('H. petiolaris') ~ 'subsp.' ~ italic('petiolaris')))

p <- ggplot() +
  geom_sf(data = states, fill = NA, colour = 'grey70', linewidth = 0.25) +
  geom_point(aes(x = Longitude, y = Latitude, fill = taxon), data = pops,
             shape = 21, size = 3, colour = 'black', stroke = 0.2) +
  scale_fill_manual(values = c('#af231c', '#fac03d', '#2775b6', '#6EC8BD'), labels = labels) +
  coord_sf(xlim = c(-125, -67), ylim = c(23, 50), expand = TRUE) +
  scale_x_continuous(labels = ~ .x) +
  scale_y_continuous(labels = ~ .x) +
  labs(x = 'Longitude', y = 'Latitude', fill = NULL) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = 'inside', legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_blank(), legend.key = element_blank())

ggsave('fourspecies.pdf', p, width = 7, height = 4.4)
ggsave('fourspecies.png', p, width = 7, height = 4.4, dpi = 300)
