#!/bin/bash
# Wavelet analysis of spatial genetic variation in each taxon (ann, arg, fal, pet), with WaveletSpatialGenetic
# (Lasky et al. 2024; functions in WaveletSpatialGenetic/). Run on mrfox from this folder, with R activated.
# Each step prints one command per chromosome (or per taxon); run them with ParaFly and wait for a step to finish:
#
#   bash wavelet_scans.sh <step> > com_<step>
#   ParaFly -c com_<step> -CPU 17          (step 8: -CPU 4, it needs more memory)
#
# Step 7 does its work directly.
#
#   A. Genome-wide wavelet genetic dissimilarity (WGD) across scales, Fig. 1B    steps 1-4  (LD-pruned SNPs)
#   B. Scale-specific genetic variance (GV) and P value for every SNP             steps 5-8  (all SNPs, MAF > 0.1)
# Input: phased, imputed VCFs from 01_filtering_phasing_polarization. Population coordinates:
# sunflower_4_sps_pop_info.csv (= 00_sampling/data/), in this folder.
step=$1
plink=/SundanceScratch/zhiqin/software/miniconda3/envs/plink/bin/plink
vcf() { echo /TEMP_SCRATCH/zhiqin/second_version/beagle/beagle_output/$1/$1.$2.vqsr90.pass.highmappability.biallelic.fbeagle.mis.2.beagled.annotated.vcf.gz; }  # <sp> <chr>
taxa="ann arg fal pet"
chroms=$(echo Ha412HOChr{01..17})
for sp in $taxa; do mkdir -p genome_wide_wavvar/$sp snp_var_pvalue_maf0.1/$sp; done
pruned=genome_wide_wavvar     # A: LD-pruned SNPs
allsnp=snp_var_pvalue_maf0.1  # B: all SNPs with MAF > 0.1

case $step in

# ============ A. Genome-wide WGD (Fig. 1B) ======================================================================
# 1. LD pruning: MAF > 0.1, PLINK --indep-pairwise 50 5 0.2.
1) for sp in $taxa; do for chr in $chroms; do
     echo "$plink --vcf $(vcf $sp $chr) --maf 0.1 --indep-pairwise 50 5 0.2 --allow-extra-chr --allow-no-sex --make-bed --double-id --out $pruned/$sp/$chr"
   done; done ;;

# 2. Genotypes of the pruned SNPs as allele counts per individual (--recode A).
2) for sp in $taxa; do for chr in $chroms; do
     echo "$plink --vcf $(vcf $sp $chr) --extract $pruned/$sp/$chr.prune.in --recode A --allow-extra-chr --allow-no-sex --make-bed --double-id --out $pruned/$sp/$chr.prunein"
   done; done ;;

# 3. Allele frequency per population (change_raw2freq.R; frequency of the major allele of the taxon).
3) for sp in $taxa; do for chr in $chroms; do
     echo "Rscript change_raw2freq.R -s $pruned/$sp/$chr.prunein.fam -g $pruned/$sp/$chr.prunein.raw -o $pruned/$sp/$chr.freq"
   done; done ;;

# 4. Genome-wide WGD at 15 log-spaced scales (1st percentile to half the 99th percentile of pairwise distances
#    between populations), with a null from 100 permutations of population locations (GenomicWavDissim2.R).
#    Output: <sp>.GenomicWavDissim.Rdata (plotted by fig1b_wgd_scales.R).
4) for sp in $taxa; do
     echo "Rscript GenomicWavDissim2.R -s $pruned/$sp/Ha412HOChr01.freq -p $pruned/$sp/ -m 15 -o $sp.GenomicWavDissim.Rdata"
   done ;;

# ============ B. Scale-specific GV and P values ================================================================
# 5. All SNPs with MAF > 0.1, as allele counts.
5) for sp in $taxa; do for chr in $chroms; do
     echo "$plink --vcf $(vcf $sp $chr) --maf 0.1 --allow-extra-chr --allow-no-sex --make-bed --double-id --recode A --out $allsnp/$sp/$chr.A"
   done; done ;;

# 6. Allele frequency per population.
6) for sp in $taxa; do for chr in $chroms; do
     echo "Rscript change_raw2freq.R -s $allsnp/$sp/$chr.A.fam -g $allsnp/$sp/$chr.A.raw -o $allsnp/$sp/$chr.freq"
   done; done ;;

# 7. SNP positions in the same order as the frequency columns: <sp>.maf0.1.snp.pos
7) for sp in $taxa; do
     for chr in $chroms; do cat $allsnp/$sp/$chr.A.bim; done > $allsnp/$sp/$sp.maf0.1.snp.pos
   done ;;

# 8. GV of every SNP at the six analysed scales of the taxon, and its P value from a chi-square distribution fitted
#    by maximum likelihood (snp_var_pvalue.R). Output: <sp>.pvalue.maf0.1.Rdata, <sp>.wavvar.maf0.1.Rdata,
#    <sp>.rm_snp.Rdata (monomorphic SNPs removed). Outliers and genes: 03_wavelet_gea_gwas_haploblocks/wavelet_candidates.R
8) for sp in $taxa; do
     echo "Rscript snp_var_pvalue.R -s $allsnp/$sp/Ha412HOChr01.freq -g $allsnp/$sp -m 15 -n $sp -p $sp.pvalue.maf0.1.Rdata -v $sp.wavvar.maf0.1.Rdata -r $sp.rm_snp.Rdata"
   done ;;

*) echo "unknown step: $step (1-8)"; exit 1 ;;
esac
