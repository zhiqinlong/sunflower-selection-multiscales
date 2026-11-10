#!/bin/bash
# Variant filtering, imputation/phasing (Beagle) and ancestral-state polarisation (EST-SFS).
# Run on mrfox. Each step prints one command per chromosome (and taxon / outgroup pair); run them in parallel
# with ParaFly and wait for a step to finish before starting the next one:
#
#   bash filter_phase_polarize.sh <step> > com_<step>
#   ParaFly -c com_<step> -CPU 17
#
# Steps 1, 16 and 19 do their work directly (no ParaFly needed).
#
#   A. Filtering (all wild individuals together)          1-3
#   B. Imputation and phasing, per taxon                    4-8
#   C. Polarisation with two outgroups at a time            9-16
#   D. Phased haplotypes coded 0 = ancestral, 1 = derived  17-20 (input for rehh / iHS)
step=$1

# ---- inputs and tools (paths on mrfox) ---------------------------------------------------------------------
vcf_all=/SundanceScratch/Wild_HA412_VCF/w1506.allsite.snp2.vqsr90.vcf.gz      # joint VQSR callset (Todesco et al. 2020)
mappability=/TEMP_SCRATCH/zhiqin/rawdata/Ha412HOv2.0-20181130.fasta.k100e4.mappability.win1k
sample_lists=/TEMP_SCRATCH/zhiqin/usfs                                         # <sp>/<sp>.sample.lists
outgroup_vcf=/TEMP_SCRATCH/zhiqin/outgroup_chr_vcf                             # <chr>.outgroup.highmappability.noindel.annotate.vcf.gz
genetic_map=./map                                                              # Beagle genetic map per chromosome
beagle=/SundanceScratch/zhiqin/software/beagle.22Jul22.46e.jar
est_sfs=~/software/est-sfs-release-2.04/est-sfs
est_sfs_config=~/software/est-sfs-release-2.04/config-rate6.txt
python=/SundanceScratch/zhiqin/software/miniconda3/envs/python3.10.6/bin/python

taxa="ann arg pet fal"
chroms=$(echo Ha412HOChr{01..17})
focal=ann     # polarisation (part C) was run with H. annuus as the focal taxon; its calls are used for all four taxa

# file names used below
filtered() { echo filtered_vcf/$1.mappability.nomulti.tranche90.vcf.gz; }                         # <chr>
unphased() { echo beagle_input/$1/$1.$2.vqsr90.pass.highmappability.biallelic.fbeagle.noindel.mis.2; }  # <sp> <chr>
phased()   { echo beagle_output/$1/$1.$2.vqsr90.pass.highmappability.biallelic.fbeagle.mis.2.beagled; } # <sp> <chr>
mkdir -p filtered_vcf polarization hapmap new_hapmap
for sp in $taxa; do mkdir -p beagle_input/$sp beagle_output/$sp hapmap/$sp; done

# outgroups: H. giganteus (664647_GIG), H. decapetalus (DEC_1895), H. divaricatus (DIV_1956),
# H. grosseserratus (GRO_2043), one name per line in outgroup.sample.list. EST-SFS takes two outgroups, so all
# 6 pairs are run: polarization/pair_<i>_<j> holds the two sample names.
mapfile -t og < outgroup.sample.list
for ((i = 0; i < ${#og[@]}; i++)); do for ((j = i + 1; j < ${#og[@]}; j++)); do
  printf "%s\n%s\n" "${og[i]}" "${og[j]}" > polarization/pair_${i}_${j}
done; done
pairs=$(ls polarization/pair_[0-9]_[0-9])

case $step in

# ============ A. Filtering =====================================================================================
# 1. High-mappability regions. Column 3 of the mappability file is, for each 1-kb window (column 2 = centre), the
#    average number of genome locations to which each 100-bp sequence in the window matches with up to 4
#    mismatches (1 = unique). Windows with an average below 500 are kept (78.5% of the genome).
1) awk '$3 < 500 {print $1"\t"$2-500"\t"$2+500"\t"$3}' $mappability > filtered_vcf/high.mappability.bed ;;

# 2. Keep variants inside the high-mappability windows.
2) for chr in $chroms; do
     grep -w $chr filtered_vcf/high.mappability.bed > filtered_vcf/$chr.high.mappability.bed
     echo "bcftools view -R filtered_vcf/$chr.high.mappability.bed -Oz -o filtered_vcf/$chr.mappability.vcf.gz $vcf_all"
   done ;;

# 3. Remove multi-allelic sites and keep sites that pass VQSR (90% tranche).
3) for chr in $chroms; do
     echo "bcftools filter -e 'N_ALT>1' -Ou filtered_vcf/$chr.mappability.vcf.gz | bcftools view -f PASS -Oz -o $(filtered $chr)"
   done ;;

# ============ B. Imputation and phasing, per taxon =============================================================
# 4. Split into the four taxa.
4) for sp in $taxa; do for chr in $chroms; do
     echo "bcftools view -S $sample_lists/$sp/$sp.sample.lists -Oz -o $(unphased $sp $chr | sed 's/.noindel.mis.2$//').vcf.gz $(filtered $chr)"
   done; done ;;

# 5. SNPs only, called in at least 80% of the individuals of the taxon.
5) for sp in $taxa; do for chr in $chroms; do
     echo "vcftools --gzvcf $(unphased $sp $chr | sed 's/.noindel.mis.2$//').vcf.gz --max-missing 0.8 --remove-indels --recode --recode-INFO-all --out $(unphased $sp $chr)"
   done; done ;;

# 6. Impute and phase with Beagle v22Jul22.46e (Ne = 1,000,000, with the genetic map). The sed writes missing
#    genotypes as "./." so that Beagle reads them.
6) for sp in $taxa; do for chr in $chroms; do
     echo "sed -i 's/\t\.:/\t.\/.:/g' $(unphased $sp $chr).recode.vcf && java -Xmx10g -jar $beagle gt=$(unphased $sp $chr).recode.vcf ne=1000000 map=$genetic_map/$chr nthreads=4 out=$(phased $sp $chr)"
   done; done ;;

# 7. Phased VCF: set SNP IDs to CHROM:POS and index. Used for wavelet, iHS and haplotype analyses.
7) for sp in $taxa; do for chr in $chroms; do
     echo "tabix $(phased $sp $chr).vcf.gz && bcftools annotate --set-id +'%CHROM\:%POS' -Oz -o $(phased $sp $chr).annotated.vcf.gz $(phased $sp $chr).vcf.gz && tabix $(phased $sp $chr).annotated.vcf.gz"
   done; done ;;

# 8. Unphased (not imputed) VCF: same IDs. Used for site-frequency analyses (SweepFinder2, DMC, polarisation).
8) for sp in $taxa; do for chr in $chroms; do
     echo "bgzip $(unphased $sp $chr).recode.vcf && bcftools annotate --set-id +'%CHROM\:%POS' -Oz -o $(unphased $sp $chr).annotate.vcf.gz $(unphased $sp $chr).recode.vcf.gz"
   done; done ;;

# ============ C. Polarisation (EST-SFS, focal taxon = ann, each pair of outgroups) =============================
# 9. Outgroup sites with no missing genotype in either outgroup of the pair.
9) for p in $pairs; do for chr in $chroms; do
     echo "bcftools view -S $p -Ou $outgroup_vcf/$chr.outgroup.highmappability.noindel.annotate.vcf.gz | bcftools view -g ^miss -Ou | bcftools sort -Oz -o $p.$chr.nomiss.sorted.vcf.gz && bcftools index -f $p.$chr.nomiss.sorted.vcf.gz"
   done; done ;;

# 10. SNPs present in both the focal taxon and the outgroups.
10) for p in $pairs; do for chr in $chroms; do
      echo "bcftools isec -c all -n=2 -w 1 $p.$chr.nomiss.sorted.vcf.gz $(unphased $focal $chr).annotate.vcf.gz | grep -v '#' | cut -f 3 > $p.$chr.intersection.snp.txt"
    done; done ;;

# 11. Subset the focal and outgroup VCFs to those shared SNPs.
11) for p in $pairs; do for chr in $chroms; do
      echo "vcftools --gzvcf $(unphased $focal $chr).annotate.vcf.gz --snps $p.$chr.intersection.snp.txt --recode --recode-INFO-all --out $p.$chr.intersection.sample"
      echo "vcftools --gzvcf $p.$chr.nomiss.sorted.vcf.gz --snps $p.$chr.intersection.snp.txt --recode --recode-INFO-all --out $p.$chr.intersection.out"
    done; done ;;

# 12. Genotype tables (one column per individual).
12) for p in $pairs; do for chr in $chroms; do
      echo -e "CHROM\tPOS\tREF\tALT\t$(paste -sd '\t' $sample_lists/$focal/$focal.sample.lists)" > $p.$chr.samples.count
      echo -e "CHROM\tPOS\tREF\tALT\t$(paste -sd '\t' $p)" > $p.$chr.out.count
      echo "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' $p.$chr.intersection.sample.recode.vcf >> $p.$chr.samples.count"
      echo "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' $p.$chr.intersection.out.recode.vcf >> $p.$chr.out.count"
    done; done ;;

# 13. EST-SFS input: focal A,C,G,T counts, plus one allele per outgroup (a heterozygous outgroup contributes one
#     randomly chosen allele). Also records which allele is the major allele in the focal taxon.
13) for p in $pairs; do for chr in $chroms; do
      echo "$python ac_to_estsfs_format.py -com $p.$chr"
    done; done ;;

# 14. Run EST-SFS v2.04 (rate-6 model), one run per outgroup pair over all chromosomes.
14) for p in $pairs; do
      for chr in $chroms; do cat $p.${chr}_data_file.txt; done > ${p}_data_file.txt
      echo $RANDOM > $p.seed
      echo "$est_sfs $est_sfs_config ${p}_data_file.txt $p.seed $p.est_sfs.output.txt $p.est_sfs.output.p_anc.txt"
    done ;;

# 15. Call the ancestral allele: the major allele is ancestral if P(major is ancestral) > 0.5, otherwise the minor
#     allele. If both alleles are equally common, the first in A<C<G<T order is treated as the major allele.
15) for p in $pairs; do
      echo -e "CHROM\tPOS\tREF\tALT\tgeno_sum_s\tnonmissing\tmajor" > ${p}_major.txt
      for chr in $chroms; do grep -v nonmissing $p.${chr}_major.txt; done >> ${p}_major.txt
      echo "$python extract_derived.py -com $p"
    done ;;

# 16. Combine the 6 outgroup pairs: all SNPs polarised in at least one pair (13,682,118 SNPs), and the SNPs at
#     which the ALT allele is ancestral. Column 9 of *.derived_out.txt is the ancestral allele (REF or ALT).
16) cat polarization/*.derived_out.txt | grep -v CHROM | awk '$9 != "" {print $1"\t"$2}' | sort -k1,1 -k2,2n | uniq | awk '{print $1":"$2}' > all_polarized_snps.txt
    cat polarization/*.derived_out.txt | grep -v CHROM | awk '$9 == "ALT" {print $1"\t"$2}' | sort -k1,1 -k2,2n | uniq | awk '{print $1":"$2}' > alt_ancestral.txt
    wc -l all_polarized_snps.txt alt_ancestral.txt ;;

# ============ D. Haplotypes coded as ancestral / derived (input for rehh) =====================================
# 17. Keep the polarised SNPs in each taxon's phased VCF.
17) for sp in $taxa; do for chr in $chroms; do
      echo "vcftools --gzvcf $(phased $sp $chr).annotated.vcf.gz --snps all_polarized_snps.txt --recode --recode-INFO-all --out hapmap/$chr.$sp.polarized && bgzip hapmap/$chr.$sp.polarized.recode.vcf && tabix hapmap/$chr.$sp.polarized.recode.vcf.gz"
    done; done ;;

# 18. Convert to hap / legend / sample files.
18) for sp in $taxa; do for chr in $chroms; do
      echo "bcftools convert --haplegendsample hapmap/$sp/$chr.phase hapmap/$chr.$sp.polarized.recode.vcf.gz"
    done; done ;;

# 19. Unzip and write the map file (chromosome, position, position in Mb, position).
19) for sp in $taxa; do for chr in $chroms; do
      gzip -df hapmap/$sp/$chr.phase.legend.gz hapmap/$sp/$chr.phase.hap.gz
      grep -v id hapmap/$sp/$chr.phase.legend | awk '{split($1, a, ":"); print a[1]" "$2" "$2/1000000" "$2}' > hapmap/$sp/$chr.phase.map
    done; done ;;

# 20. Recode haplotypes so that 0 = ancestral and 1 = derived (flip SNPs in alt_ancestral.txt). Output: new_hapmap/<sp>/
20) for sp in $taxa; do for chr in $chroms; do
      echo "$python hap_transform_derived.py -sp $sp -chr $chr"
    done; done ;;

*) echo "unknown step: $step (1-20)"; exit 1 ;;
esac
