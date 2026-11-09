#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#
# H. annuus RNA-seq pipeline: raw fastq -> expression (TPM) -> co-expression network.
# Run on rorqual (Alliance Canada) from the project folder, which contains raw_data/ (180 fastq pairs from HeliaGene)
# and ref/ (RSEM/STAR index of HA412-HOv2, prefix ref/Ha412HOv2.0). Submit one step at a time; wait for each to finish.
#
#   bash   run_transcriptome.sh lists                                                       # step 0 (seconds, login node)
#   sbatch --array=1-180 --cpus-per-task=4  --mem-per-cpu=3G  run_transcriptome.sh fastp     # step 1
#   sbatch --array=1-45  --cpus-per-task=1  --mem-per-cpu=2G  run_transcriptome.sh merge     # step 2
#   sbatch --array=1-135 --cpus-per-task=4  --mem-per-cpu=20G run_transcriptome.sh rsem      # step 3
#   sbatch --cpus-per-task=40 --mem-per-cpu=4G run_transcriptome.sh wgcna 5                  # step 4 (top 5% edges)
#   sbatch --cpus-per-task=40 --mem-per-cpu=4G run_transcriptome.sh wgcna 1                  # step 4 (top 1% edges)
#
# Steps 5-6 run locally in R (paths inside the scripts):
#   Rscript tissue_specificity.R      # tau per gene; network gene IDs matched to positions
#   Rscript pleiotropy_level_plot.R   # Fig. 4A-C, Figs. S15-S17
set -euo pipefail
mkdir -p logs filter_fq merged_fq rsem_out qcout expression_values
step=$1
i=${SLURM_ARRAY_TASK_ID:-1}

case $step in

# 0. Sample lists. The 180 fastq pairs are 135 libraries: 45 libraries were sequenced twice, and the second
#    (top-up) run has "-b" after the library name. Library name = the text between "SM:" and "_LB:".
lists)
  ls raw_data | grep 'R1_001.fastq.gz$' > R1_fq
  sed 's/.*SM://; s/_LB:.*//' R1_fq | grep -v -- '-b$' | sort > all_libs.txt   # 135 libraries
  sed 's/.*SM://; s/_LB:.*//' R1_fq | grep -- '-b$' | sed 's/-b$//' | sort > merge_libs.txt   # 45 with two runs
  wc -l R1_fq all_libs.txt merge_libs.txt   # expect 180, 135, 45
  ;;

# 1. Trim and filter each sequencing run with fastp, then FastQC.
#    -f/-t 15: cut 15 bases from the front of R1/R2; -5/-3 20 with --cut_mean_quality 20: sliding-window trimming at both
#    ends; -q 20: base quality threshold. Output: filter_fq/<library>R1.filtered.fastq.gz (and R2).
#    (One library, ANN1430_leaf_HS_1, was re-run with the same settings because its R1 output was incomplete.)
fastp)
  module load fastp/0.24.0
  fq1=$(sed -n "${i}p" R1_fq)
  fq2=${fq1%R1_001.fastq.gz}R2_001.fastq.gz
  name=${fq1%%_LB:*}; name=${name##*SM:}
  fastp -i raw_data/$fq1 -I raw_data/$fq2 \
        -o filter_fq/${name}R1.filtered.fastq.gz -O filter_fq/${name}R2.filtered.fastq.gz \
        -w 4 -f 15 -t 15 -5 20 -3 20 --cut_mean_quality 20 -q 20 \
        -j logs/fastp_${name}.json -h logs/fastp_${name}.html
  module load StdEnv/2023 fastqc/0.12.1
  fastqc -o qcout/ -t 4 filter_fq/${name}R1.filtered.fastq.gz filter_fq/${name}R2.filtered.fastq.gz
  ;;

# 2. One fastq pair per library. For the 45 libraries with two runs, join run 1 and run 2 (gzip files can simply be
#    concatenated) and check that read numbers add up. The other 90 libraries are linked unchanged (done by task 1).
merge)
  if [ "$i" = 1 ]; then
    for l in $(comm -23 all_libs.txt merge_libs.txt); do
      for r in R1 R2; do ln -sf ../filter_fq/${l}${r}.filtered.fastq.gz merged_fq/; done
    done
  fi
  l=$(sed -n "${i}p" merge_libs.txt)
  for r in R1 R2; do
    cat filter_fq/${l}${r}.filtered.fastq.gz filter_fq/${l}-b${r}.filtered.fastq.gz > merged_fq/${l}${r}.filtered.fastq.gz
  done
  n() { zcat "$1" | awk 'END{print NR/4}'; }
  m1=$(n merged_fq/${l}R1.filtered.fastq.gz); m2=$(n merged_fq/${l}R2.filtered.fastq.gz)
  a=$(n filter_fq/${l}R1.filtered.fastq.gz);  b=$(n filter_fq/${l}-bR1.filtered.fastq.gz)
  echo -e "$l\trun1=$a\trun2=$b\tmerged_R1=$m1\tmerged_R2=$m2"
  [ "$m1" = "$m2" ] && [ "$m1" = $((a + b)) ] || { echo "read counts do not add up for $l"; exit 1; }
  ;;

# 3. Align with STAR and quantify with RSEM against HA412-HOv2. Output: rsem_out/<library>.genes.results
#    (expected counts and TPM per gene).
rsem)
  module load StdEnv/2023 star/2.7.11b rsem/1.3.3
  l=$(sed -n "${i}p" all_libs.txt)
  rsem-calculate-expression --paired-end --star --star-gzipped-read-file -p 4 \
        merged_fq/${l}R1.filtered.fastq.gz merged_fq/${l}R2.filtered.fastq.gz ref/Ha412HOv2.0 rsem_out/$l
  ;;

# 4. Co-expression network (WGCNA.R): DESeq2 VST of the 135 libraries, unsigned WGCNA adjacency (soft threshold from
#    pickSoftThreshold; beta = 8), keep the top X% of edges, then igraph degree, strength, betweenness and harmonic
#    centrality per gene. Argument = X (5 or 1). Output: adjacency_interactino_top<X>.txt
wgcna)
  module load StdEnv/2023 r/4.5.0
  export R_LIBS=~/.local/R/$EBVERSIONR/
  cp rsem_out/*.genes.results expression_values/
  Rscript WGCNA.R "$2"
  ;;

*) echo "unknown step: $step (use lists, fastp, merge, rsem or wgcna)"; exit 1 ;;
esac
