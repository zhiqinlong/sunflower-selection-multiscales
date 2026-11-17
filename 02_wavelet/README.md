# Wavelet analysis of spatial genetic variation (local adaptation)

Uses WaveletSpatialGenetic (Lasky et al. 2024) to measure how genetic differentiation changes with geographic scale
in each taxon, and to test every SNP for unusually strong differentiation at a given scale.

| Step | File | Where | What it does | Output |
|---|---|---|---|---|
| 1 | `wavelet_scans.sh` | mrfox | All cluster steps, each explained in the script (run with ParaFly) | |
| | `change_raw2freq.R` | (steps 3, 6) | Individual genotypes → allele frequency per population | `<chr>.freq` |
| | `GenomicWavDissim2.R` | (step 4) | Genome-wide wavelet genetic dissimilarity (WGD) at 15 scales; null from 100 permutations of locations | `<sp>.GenomicWavDissim.Rdata` |
| | `snp_var_pvalue.R` | (step 8) | Scale-specific genetic variance (GV) of every SNP at the six analysed scales; χ² fit → P values | `<sp>.pvalue.maf0.1.Rdata` |
| 2 | `fig1b_wgd_scales.R` | local | Mean WGD across scales vs null | **Fig. 1B** |
| 3 | `figS2_upset_taxa.R` | local | Wavelet genes shared among taxa (after `03_.../wavelet_candidates.R`) | **Fig. S2** |

`WaveletSpatialGenetic/` contains functions by Jesse Lasky (https://github.com/jesserlasky/WaveletSpatialGenetic).
`JL_OutflankFuns5_Fst.R` and `Likelihood functions for OutFLANK_JRL.R` were changed slightly (`na.rm = TRUE`
added, extra `print()` calls). Please cite the original work.

## Settings
- **SNPs:** WGD uses LD-pruned SNPs (MAF > 0.1, PLINK `--indep-pairwise 50 5 0.2`); GV uses all SNPs with MAF > 0.1.
- **Scales:** 15 scales evenly spaced on a log axis from the 1st percentile to half the 99th percentile of pairwise
  distances between populations (Euclidean distance on latitude/longitude × 111.3 km).
- **Six analysed scales** (indices of the 15): ann 4, 6, 7, 8, 9, 10; pet 8–13; fal 6, 7, 8, 9, 10, 12; arg 6, 8, 9, 10,
  11, 12. These are the orange circles in Fig. 1B.
- **Masking:** at each scale, locations whose nearest neighbour is farther than half the scale are not used.
- **P values:** χ² distribution fitted by maximum likelihood after trimming 2.5% of GV values at each tail; outliers
  (lowest 0.1% of P over the six scales) and their genes are called in `03_wavelet_gea_gwas_haploblocks/wavelet_candidates.R`.
