# Wavelet genes that are also GWAS candidates (and GEA genes): Table S1 and Figs. S5, S6 (run locally).
#
#  1. Table S1 (gwas_wavelet_gea.csv): for each taxon, genes found by the wavelet analysis and by the GWAS of
#     Todesco et al. (2020) (Supplementary Data "GWA candidates <sp>"), with their GEA variables if any. One row per
#     distinct combination of GWAS trait, GEA variable and wavelet SNP/scale.
#  2. For each group of such genes (genes within 150 kb are drawn together):
#     A. wavelet P values of the six scales in a 1-Mb window around the genes, with haploblocks shaded (blue) and the
#        genes marked by dotted lines; the title gives the scale of the strongest wavelet SNP;
#     B. gene models with the GWAS traits;
#     C. for genes that are also GEA candidates: allele frequencies of the strongest wavelet SNP per population on a
#        map of the associated environmental variable (map_plot_function.R; rasters from environment_rasters.R).
#  Figs. S5 and S6 were assembled from these panels in Illustrator (trait colours in panel B were added there).
library(readxl)
library(data.table)
library(tidyverse)
library(ggpubr)
library(gggenes)

base     <- '/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis'
out_dir  <- file.path(base, 'gea_gwas_wavelet')
sps      <- c('pet', 'fal', 'arg', 'ann')

######## 1. Table S1 ########
gwa_gea_wavelet <- list()
for (sp in sps) {
  gwas <- read_excel(file.path(out_dir, '41586_2020_2467_MOESM4_ESM.xlsx'), sheet = paste0("GWA candidates ", sp)) %>%
    fill(1, .direction = "down") %>%
    select(1:12) %>%
    rename(gene_id = `Ha412-HO gene ID`) %>% mutate(GWAS = 1)
  gea <- read.csv(file.path(base, 'wavelet/gea_loci', paste0(sp, '_gea_gene.bed')), sep = '\t', header = FALSE) %>%
    select(c(2, 3, 4, 6, 7, 8)) %>%
    rename(GEA_locus_start = V2, GEA_locus_end = V3, gene_id = V8, gene_start = V6, gene_end = V7, EV = V4) %>%
    mutate(GEA = 1)
  wavelet <- read.csv(file.path(base, 'wavelet/separate_scale/top0001', paste0(sp, '_top0001_adap_gene.bed')), sep = '\t', header = FALSE) %>%
    select(c(1, 4, 5, 6, 7, 8)) %>%
    rename(Chr = V1, gene_id = V4, Wavelet_locus_start = V5, Wavelet_locus_end = V6, Wavelet_Pvalue = V7, Wavelet_scale = V8) %>%
    mutate(Wavelet = 1)
  gwa_gea_wavelet[[sp]] <- list(gwas, gea, wavelet) %>%
    reduce(full_join, by = "gene_id") %>%
    mutate(across(c(GEA, Wavelet, GWAS), ~replace_na(.x, 0))) %>%
    mutate(type = rowSums(across(c(GEA, Wavelet, GWAS)))) %>%
    filter(type >= 2)
}
gwa_gea_wavelet_df <- bind_rows(gwa_gea_wavelet, .id = "species") %>%
  mutate(species = toupper(species)) %>%
  filter(Wavelet == 1 & GWAS == 1) %>%
  distinct()
# Table S1: 13 columns, one row per distinct combination (the full join repeats rows that differ only in columns
# not shown, e.g. GEA locus coordinates)
table_s1 <- gwa_gea_wavelet_df %>%
  transmute(Species = species, GWAS_Trait = Trait, Gene_id = gene_id,
            `Closest Arabidopsis ortholog` = `Closest Arabidopsis hortolog`, `Blast score`,
            `Arabidopsis gene description`, `Arabidopsis gene symbol`, GEA_EV = EV, Chromosome = Chr,
            Wavelet_locus_start, Wavelet_locus_end, Wavelet_Pvalue, Wavelet_scale) %>%
  distinct()
write.csv(table_s1, file.path(out_dir, 'gwas_wavelet_gea.csv'), quote = TRUE, row.names = FALSE, fileEncoding = "UTF-8")

######## 2. Figures S5 / S6 ########
gff3 <- read.csv("/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_.sorted.interpro.gff3",
                 sep = "\t", header = FALSE, comment.char = "#") %>%
  mutate(gene_id = str_extract(V9, "Ha412HOChr\\d+g\\d+"))
inversions <- read.csv(file.path(base, 'data/Marco_nature/all_inversions_location.csv'))
pops <- read.csv(file.path(base, 'wavelet/sunflower_4_sps_pop_info.csv'))
scale_index <- list(ann = c(4, 6, 7, 8, 9, 10), pet = c(8, 9, 10, 11, 12, 13),
                    fal = c(6, 7, 8, 9, 10, 12), arg = c(6, 8, 9, 10, 11, 12))
extend <- 150000      # genes closer than this are drawn in one figure
zoom_width <- 1e6     # window on each side of the genes
wavelet_out <- file.path(base, 'wavelet/separate_scale/wavelet_output')
snp_freq_dir <- file.path(base, 'wavelet/separate_scale/evs_snp_freq/snp_var_pvalue_maf0.1')
source(file.path(base, 'wavelet/separate_scale/evs_snp_freq/map_plot_function.R'))   # = map_plot_function.R in this folder
load(file.path(base, 'wavelet/separate_scale/all_dataframe_envs_data.Rdata'))         # 'rasters', from environment_rasters.R

for (sp in sps) {
  tobe_visualized <- gwa_gea_wavelet_df %>%
    filter(species == toupper(sp), !is.na(Wavelet_Pvalue)) %>%
    group_by(gene_id) %>%
    slice_min(Wavelet_Pvalue, n = 1, with_ties = TRUE) %>%   # strongest wavelet SNP per gene
    ungroup()
  if (nrow(tobe_visualized) == 0) next

  # wavelet P values of the six scales for this taxon
  load(file.path(wavelet_out, paste0(sp, ".pvalue.maf0.1.Rdata")))   # jl_out
  load(file.path(wavelet_out, paste0(sp, ".rm_snp.Rdata")))          # rm_snp
  snp <- fread(file.path(wavelet_out, paste0(sp, ".maf0.1.snp.pos")), header = FALSE) %>% select(chr = V1, pos = V4)
  if (length(rm_snp) > 0) snp <- snp[-rm_snp, ]
  samp <- read.csv(file.path(base, 'wavelet', paste0(sp, ".sample")), sep = "\t", header = FALSE) %>% transmute(ind_ID = V1)
  xy <- inner_join(pops, samp, by = "ind_ID") %>% distinct(Population.ID, Latitude, Longitude)
  D <- as.matrix(dist(xy[, c("Latitude", "Longitude")])) * 111.3195
  mn <- quantile(D[upper.tri(D)], .01, na.rm = TRUE); if (mn == 0) mn <- 15
  mx <- quantile(D[upper.tri(D)], .99, na.rm = TRUE) / 2
  scales <- exp(seq(log(mn), log(mx), length.out = 15))[scale_index[[sp]]]
  plot_data <- bind_rows(replicate(6, snp, simplify = FALSE)) %>%
    mutate(pvalue = unlist(lapply(jl_out[1:6], function(x) x[[3]])),
           scale = rep(round(scales, 2), each = nrow(snp)),
           logp = -log10(pmax(pvalue, .Machine$double.xmin)),
           pos = as.numeric(as.character(pos)))
  cutoff <- quantile(plot_data$pvalue, .001, na.rm = TRUE)
  inv_sp <- ifelse(sp == "ann", "annuus", ifelse(sp == "fal", "pet", sp))
  inv <- inversions %>% filter(sps == inv_sp) %>%
    mutate(START = START * 1e6, END = END * 1e6, chr = sprintf("Ha412HOChr%02d", CHR))

  # group nearby genes
  gene_group <- gff3 %>%
    filter(V3 == "gene", gene_id %in% tobe_visualized$gene_id) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(V1, V4) %>%
    group_by(V1) %>%
    mutate(previous_end = lag(cummax(V5)),
           new_group = if_else(is.na(previous_end) | V4 - previous_end > extend, 1L, 0L),
           group = cumsum(new_group)) %>%
    group_by(V1, group) %>%
    mutate(group_id = cur_group_id()) %>%
    ungroup()

  for (g in unique(gene_group$group_id)) {
    genes <- gene_group %>% filter(group_id == g) %>% pull(gene_id)
    sub_gff3 <- gff3 %>% filter(gene_id %in% genes) %>% select(V1, V3, V4, V5, V7, gene_id)
    chr <- unique(sub_gff3$V1)
    gene_start <- min(sub_gff3$V4); gene_end <- max(sub_gff3$V5)
    tag <- paste0(sp, chr, gene_end)

    # B. gene models with GWAS traits
    gene_df <- sub_gff3 %>%
      filter(V3 == "gene") %>%
      left_join(tobe_visualized %>% select(gene_id, Trait), by = "gene_id") %>%
      mutate(forward = V7 == "+", y = factor(gene_id, levels = rev(genes))) %>%
      distinct()
    gene_p <- ggplot(gene_df, aes(xmin = V4, xmax = V5, y = as.factor(y), forward = forward, fill = 'grey')) +
      geom_gene_arrow() +
      geom_gene_label(aes(label = gene_id), size = 3) +
      geom_text(aes(x = (V4 + V5) / 2, y = as.factor(y), label = Trait), vjust = -1.5, size = 3, inherit.aes = FALSE) +
      scale_y_discrete(expand = expansion(add = c(.5, 1))) +
      scale_x_continuous(labels = function(x) paste0(round(x / 1e6, 2), " Mb")) +
      labs(x = "Genomic position", y = NULL) +
      theme_genes() + theme(legend.position = "none")

    # A. zoomed wavelet P values
    manh_plot <- sub_gff3 %>% filter(V3 == "gene") %>% left_join(tobe_visualized, by = "gene_id")
    best <- manh_plot %>% filter(!is.na(Wavelet_Pvalue)) %>% slice_min(Wavelet_Pvalue, n = 1, with_ties = FALSE)
    sm <- scales[which.min(abs(scales - best$Wavelet_scale))]      # scale of the strongest wavelet SNP
    win_start <- max(0, gene_start - zoom_width); win_end <- gene_end + zoom_width
    d <- plot_data %>% filter(chr == !!chr, between(pos, win_start, win_end)) %>% mutate(sig = pvalue <= cutoff)
    iv <- inv %>% filter(chr == !!chr, END >= win_start, START <= win_end)
    manh_p <- ggplot() +
      geom_rect(data = iv, aes(xmin = START, xmax = END), ymin = -Inf, ymax = Inf, fill = "darkblue", alpha = .1, colour = NA) +
      geom_point(data = filter(d, !sig), aes(pos, logp), colour = "grey65", size = .6, alpha = .5) +
      geom_point(data = filter(d, sig), aes(pos, logp), colour = "#D33662", size = 1.2) +
      facet_wrap(~scale, scales = 'free_y', strip.position = 'right', ncol = 1) +
      geom_hline(yintercept = -log10(cutoff), colour = "#3366CC", linetype = "longdash") +
      geom_vline(xintercept = c(gene_start, gene_end), linetype = "dotted") +
      scale_x_continuous(labels = function(x) paste0(round(x / 1e6, 2), " Mb")) +
      coord_cartesian(xlim = c(win_start, win_end)) +
      labs(title = paste(toupper(sp), chr, paste0("Scale = ", round(sm, 2), " km")),
           x = "Genomic position", y = expression(-log[10](italic(P)))) +
      theme_classic() + theme(panel.border = element_rect(fill = NA))
    ggsave(ggarrange(manh_p, gene_p, nrow = 2),
           file = file.path(out_dir, 'wavelet_gea_gwas/wavelet_gwas', paste0(tag, 'manh_gene.pdf')), width = 5, height = 9)

    # C. allele-frequency maps for genes that are also GEA candidates
    map_rows <- manh_plot %>% filter(!is.na(EV), !is.na(Wavelet_locus_end)) %>% distinct(gene_id, EV, .keep_all = TRUE)
    if (nrow(map_rows) == 0) next
    bim <- read.csv(file.path(snp_freq_dir, sp, paste0(chr, '.A.bim')), sep = '\t', header = FALSE)
    freq <- fread(file.path(snp_freq_dir, sp, paste0(chr, '.freq')), header = FALSE) %>% as.data.frame()
    colnames(freq) <- c('Population_ID', bim$V2)
    for (i in seq_len(nrow(map_rows))) {
      snp_id <- paste(chr, map_rows$Wavelet_locus_end[i], sep = ":")
      plot_freq <- freq[, c('Population_ID', snp_id)] %>%
        rename(major = 2) %>% mutate(minor = 1 - major) %>% select(Population_ID, major, minor)
      ev <- map_rows$EV[i]
      map_p <- map_plot(sp, plot_freq, ev, map_rows$Wavelet_scale[i], rasters[[sp]][[ev]], map_rows$gene_id[i])
      ggsave(map_p, file = file.path(out_dir, 'wavelet_gea_gwas', paste0(tag, map_rows$gene_id[i], ev, '.map.png')), width = 6, height = 6)
    }
  }
}
