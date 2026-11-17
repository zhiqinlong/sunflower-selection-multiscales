# GEA candidate genes and haploblock intervals for each taxon (run locally after gea2loci.py).
#
# 1. Collinearity: for each taxon, environmental variables with a GEA result are correlated across populations
#    (Spearman); caret::findCorrelation(cutoff = 0.7) removes, from each pair with |rho| > 0.7, the variable with the
#    larger mean absolute correlation (deterministic). Retained: 10 (pet), 12 (fal), 13 (arg), 13 (ann) variables.
#    Variable set: the GEA variables that also have a column in environmental_factors.csv = 21 climate + 14 soil
#    variables + elevation (36). Sodium is not included (GEA file "sodium", column "NA"); MAR and PAS have no column.
# 2. GEA loci of the retained variables (|BF| > 20, from gea2loci.py) -> <sp>_gea_loci.bed; genes whose body +/- 10 kb
#    overlaps a GEA locus -> <sp>_gea_gene.bed (input for enrichment_gea_haploblocks.R and wavelet_gea_gwas.R).
# 3. Haploblocks of Todesco et al. (2020) per taxon -> <sp>.inversion.bed (FAL uses the H. petiolaris haploblocks),
#    and chromosome lengths -> genome.txt (for bedtools shuffle).
library(caret)
library(dplyr)
library(stringr)

dir <- '/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/gea_loci'
gea_dir <- file.path(dir, 'gea_selection')          # one BED per variable and taxon, from gea2loci.py
evs <- read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/Marco_nature/environmental_factors.csv')
gene_file <- '/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_genes_updown10kb.bed'
bedtools <- "/Users/zhiqinlong/Desktop/computer/work/Bio_tools/bedtools2/bin/bedtools"

taxa <- list(pet = c(gea = 'petpet',      taxon = 'petiolaris subsp. petiolaris'),
             fal = c(gea = 'petfal',      taxon = 'fal'),
             arg = c(gea = 'argophyllus', taxon = 'arg'),
             ann = c(gea = 'annuus',      taxon = 'ann'))

for (sp in names(taxa)) {
  ev_sp <- taxa[[sp]][['gea']]
  # variables with a GEA result for this taxon (file names <variable>[_e]_<taxon>.bed)
  EV_names <- list.files(gea_dir, pattern = paste0(ev_sp, '.bed'))
  EV_names <- sapply(strsplit(EV_names, paste0('_', ev_sp, "|_e")), function(x) x[1])
  EV_names <- setdiff(EV_names, c("latitude", "longitude"))

  sp_evs <- evs %>%
    filter(str_detect(Taxon, taxa[[sp]][['taxon']])) %>%
    select(any_of(EV_names)) %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  corrmatrix <- cor(sp_evs, method = "spearman")
  selected_evs <- setdiff(colnames(corrmatrix), findCorrelation(corrmatrix, cutoff = 0.7, names = TRUE))
  cat(sp, ':', ncol(sp_evs), 'candidate variables,', length(selected_evs), 'retained:', sort(selected_evs), '\n')

  # GEA loci of the retained variables, and genes (+/- 10 kb) overlapping them
  gea_loci <- bind_rows(lapply(selected_evs, function(ev) {
    file <- list.files(gea_dir, pattern = paste0('^', ev, ".*", ev_sp), full.names = TRUE)
    read.csv(file, sep = '\t', header = FALSE) %>% mutate(EV = ev)
  }))
  adap_file <- file.path(dir, paste0(sp, '_gea_loci.bed'))
  write.table(distinct(gea_loci), file = adap_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  system(paste(bedtools, 'intersect -a', adap_file, '-b', gene_file, '-wa -wb >', file.path(dir, paste0(sp, '_gea_gene.bed'))))
}

# Haploblocks (Todesco et al. 2020); positions in the table are in Mb
inversions <- read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/Marco_nature/all_inversions_location.csv')
inv_taxon <- c(ann = 'annuus', arg = 'arg', pet = 'pet', fal = 'pet')
for (sp in names(inv_taxon)) {
  inversions %>%
    filter(sps == inv_taxon[[sp]]) %>%
    mutate(START = as.integer(START * 1e6), END = as.integer(END * 1e6),
           chr = sprintf('Ha412HOChr%02d', CHR)) %>%
    select(chr, START, END) %>%
    distinct() %>%
    write.table(file.path(dir, paste0(sp, '.inversion.bed')), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

# Chromosome lengths (17 chromosomes, no unplaced contigs) for bedtools shuffle
read.csv(file.path(dir, 'Ha412HOv2.0-20181130.fasta.fai'), header = FALSE, sep = '\t') %>%
  filter(!grepl('Ha412HOChr00c', V1)) %>%
  select(1, 2) %>%
  write.table(file.path(dir, 'genome.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
