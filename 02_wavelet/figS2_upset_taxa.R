# Fig. S2: UpSet plot of locally adaptive (wavelet) genes shared among the four taxa.
# Input: top0001/<sp>_top0001_adap_gene.bed from 03_wavelet_gea_gwas_haploblocks/wavelet_candidates.R
# (genes ANN 4,719, ARG 1,727, PET 3,445, FAL 3,158). Colours were changed in Illustrator.
library(UpSetR)
library(tidyverse)
gene_dir <- '/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/top0001'
sps <- c('ann', 'arg', 'pet', 'fal')
data_upset <- lapply(sps, function(sp)
  unique(read.csv(file.path(gene_dir, paste0(sp, '_top0001_adap_gene.bed')), sep = '\t', header = FALSE)$V4))
names(data_upset) <- toupper(sps)
print(lengths(data_upset))

pdf(file.path(dirname(gene_dir), 'wavelet_upset.pdf'), width = 8, height = 5)
upset(fromList(data_upset),
      order.by = "freq",
      sets = toupper(sps),
      mainbar.y.label = "Gene number",
      sets.x.label = "Gene numbers",
      set_size.show = TRUE,
      point.size = 2,
      line.size = 1,
      mb.ratio = c(0.7, 0.3),
      text.scale = c(1.5, 1.25, 1.25, 1.5, 1.5, 1.5),
      main.bar.color = "violetred4",
      sets.bar.color = "turquoise4",
      matrix.color = "slateblue4",
      shade.color = "wheat4")
dev.off()
