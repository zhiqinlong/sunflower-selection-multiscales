# Number of gene-duplication events for each H. annuus HA412-HOv2 gene (used for Fig. 4B and Figs. S15-S17).
#
# OrthoFinder counts duplication events per orthogroup (gene family) across H. annuus, H. argophyllus and
# H. petiolaris subsp. neglectus. Each HA412-HOv2 gene is given the count of the orthogroup it belongs to, so the
# value describes the duplication history of the gene's family, not of the gene alone. Genes from all four study taxa
# were identified on the HA412-HOv2 annotation, so this one table serves all four taxa.
#
# Same steps as the "ann_dup" table in ../transcriptome/pleiotropy_level_plot.R. Output: gene_duplications.csv
library(tidyverse)

dir <- '/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/duplication_level'   # OrthoFinder output
gene_bed <- '/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_genes.bed'
gtf <- '/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/old_results/2023Aug11_FINAL_INTERPRO.AED_.sorted.interpro.new.gtf'

# 1. HA412-HOv2 genes in each orthogroup (Orthogroups.txt: "OG0000001: gene gene gene ...")
og_ann <- read.table(file.path(dir, 'Orthogroups.txt'), sep = ':', header = FALSE) %>%
  separate_rows(V2, sep = ' ') %>%
  mutate(V2 = trimws(V2)) %>%
  filter(grepl('Ha412HOChr', V2)) %>%
  distinct() %>%
  transmute(protein_id = V2, Orthogroup = V1)

# 2. The protein IDs (Ha412HOChr1g...) differ from the gene IDs in the gene BED (Ha412HOChr01g...), so match them
#    through genomic position: protein ID -> position (GTF) -> gene ID (gene BED).
protein_pos <- read.csv(gtf, sep = '\t', header = FALSE) %>%
  filter(V3 == 'transcript') %>%
  separate(V9, into = c('protein_id', 'ts'), sep = '; transcript_') %>%
  mutate(protein_id = gsub('gene_id mRNA:', '', protein_id)) %>%
  transmute(protein_id, chr = V1, s = V4, e = V5)
genes_pos <- read.csv(gene_bed, sep = '\t', header = FALSE) %>%
  mutate(gene_id = gsub('mRNA:|ncRNA:', '', V4)) %>%
  filter(!grepl('Ha412HOChr00c', gene_id)) %>%      # unplaced contigs
  transmute(gene_id, chr = V1, s = V2, e = V3)

# 3. Duplication events per orthogroup (events with >= 50% support), assigned to each gene.
gene_dup <- read.csv(file.path(dir, 'Duplications_per_Orthogroup.tsv'), sep = '\t') %>%
  rename(n_dup = Duplications..50..support.) %>%
  merge(og_ann, by = 'Orthogroup', all = TRUE) %>%
  merge(protein_pos, by = 'protein_id', all = TRUE) %>%
  select(-protein_id) %>%
  left_join(genes_pos, by = c('chr', 's', 'e')) %>%
  filter(!is.na(n_dup) & !is.na(chr)) %>%
  select(gene_id, chr, s, e, Orthogroup, n_dup)

write.csv(gene_dup, 'gene_duplications.csv', row.names = FALSE)

# summary
cat('genes with a duplication count:', nrow(gene_dup), '\n')
cat('orthogroups:', n_distinct(gene_dup$Orthogroup), '\n')
cat('genes in orthogroups with >= 1 duplication:', sum(gene_dup$n_dup > 0),
    sprintf('(%.1f%%)', 100 * mean(gene_dup$n_dup > 0)), '\n')
print(summary(gene_dup$n_dup))
