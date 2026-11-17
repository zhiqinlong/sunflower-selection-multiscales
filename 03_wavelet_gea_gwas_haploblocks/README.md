# Locally adaptive genes and their links to environment, traits and haploblocks

Turns the wavelet results (`02_wavelet/`) into per-scale candidate genes, and compares them with the
genotype–environment associations (GEA), GWAS candidates and large haploblocks of Todesco et al. (2020).
Run the scripts in this order:

| Step | Script | Where | What it does | Output |
|---|---|---|---|---|
| 1 | `wavelet_candidates.R` | local | Lowest 0.1% of GV P values over the six scales = outliers; genes ± 10 kb around them; Manhattan plots | `top0001/<sp>_top0001_adap_loci.bed`, `<sp>_top0001_adap_gene.bed`; **Fig. S1** |
| 2 | `gea2loci.py` | mrfox | GEA outlier SNPs (\|BF\| > 20, BayPass without inversion covariates), one BED per variable | `gea_selection/<variable>_<taxon>.bed` |
| 3 | `gea_genes.R` | local | Removes collinear variables (Spearman \|ρ\| > 0.7, `caret::findCorrelation`); genes ± 10 kb near GEA SNPs of the retained variables; haploblock BEDs; chromosome lengths | `<sp>_gea_gene.bed`, `<sp>.inversion.bed`, `genome.txt` |
| 4 | `enrichment_gea_haploblocks.R` | local | Wavelet genes vs GEA genes (1,000 random gene sets) and vs haploblocks (1,000 shuffled haploblock sets), all scales and per scale | **Figs. S3, S4**; `FigS3_S4_enrichment_results.csv` |
| 5 | `environment_rasters.R` | local, once | Map layers of each variable: 1961–1990 climate normals; soil values interpolated (IDW) between populations | `all_dataframe_envs_data.Rdata` |
| 6 | `wavelet_gea_gwas.R` | local | Wavelet genes that are also GWAS candidates (with GEA variables); zoomed Manhattan plots, gene models, allele-frequency maps (uses `map_plot_function.R`) | **Table S1** (`gwas_wavelet_gea.csv`); **Figs. S5, S6** |

`data/evs_types.csv` lists the environmental variables and whether each is a climate or soil variable.

## Notes
- **Scales.** The six analysed scales per taxon are scale indices (of 15) ann 4, 6, 7, 8, 9, 10; pet 8–13;
  fal 6, 7, 8, 9, 10, 12; arg 6, 8, 9, 10, 11, 12. They match `02_wavelet/snp_var_pvalue.R` and Fig. 1B.
- **Haploblocks.** FAL uses the *H. petiolaris* haploblocks (Todesco et al. report one set for the species).
- **Variables.** Candidates are the GEA variables with a column in `environmental_factors.csv`: 21 climate variables,
  14 soil variables and elevation (36). Sodium is not included (its GEA file is named `sodium`, its column `NA`);
  MAR and PAS have no column. After the collinearity filter, 10 (PET), 12 (FAL), 13 (ARG) and 13 (ANN) variables
  remain; elevation is one of the 13 for ANN.
- **Randomisations** use fixed seeds (`set.seed`, `bedtools shuffle -seed`), so Figs. S3 and S4 are reproducible.
