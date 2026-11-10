#!/bin/bash
# Gene-duplication events per orthogroup with OrthoFinder v3 (3.1.3), run on salal in /data/zlong/duplications.
#   bash orthofinder.sh 1   # rename sequences
#   bash orthofinder.sh 2   # keep one protein (longest isoform) per gene
#   bash orthofinder.sh 3   # OrthoFinder (~15 min on 80 threads)
#
# Proteomes (predicted proteins of three reference genomes):
#   Hann.fasta  H. annuus HA412-HOv2 (44,544 proteins, one per gene; headers ">mRNA:<gene ID>")
#   Harg.fasta  H. argophyllus ARG2202, haplotype 1 annotation (= Harg2202_H1.protein.fasta; 52,916 proteins)
#   Hpet.fasta  H. petiolaris subsp. neglectus, haplotype 1 annotation (= Hneg2702_H1.protein.fasta; 50,642 proteins).
#               No reference genome exists for subsp. petiolaris or fallax, so neglectus represents H. petiolaris.
#
# Output used downstream (run Results_Mar10_1; earlier runs Results_Mar06 / Results_Mar10 were not used):
#   Orthogroups/Orthogroups.txt                                  genes in each orthogroup
#   Comparative_Genomics_Statistics/Duplications_per_Orthogroup.tsv   number of duplication events per orthogroup
# These are read by duplications_per_gene.R.
step=$1
cd /data/zlong/duplications

if [[ $step == 1 ]]; then
  # Species prefix on the H. argophyllus / H. petiolaris IDs (both annotations use g1.t1, g2.t1, ...);
  # drop "mRNA:" from the HA412-HOv2 IDs so they match the gene IDs.
  sed 's/>/>Hpet./g' Hpet.fasta > Hpet.pep.fasta
  sed 's/>/>Harg./g' Harg.fasta > Harg.pep.fasta
  sed 's/>mRNA:/>/g' Hann.fasta > Hann.pep.fasta

elif [[ $step == 2 ]]; then
  # Longest isoform per gene (OrthoFinder's tools/primary_transcript.py; transcript IDs are <gene>.t<n>).
  # Writes primary_transcripts/<name>.pep.fasta: 44,544 / 46,913 / 44,786 proteins for Hann / Harg / Hpet.
  for sp in Hpet Harg Hann; do
    python primary_transcript.py $sp.pep.fasta last_dot_before_first_space
  done

elif [[ $step == 3 ]]; then
  # Orthogroups, gene trees, species tree, and duplication events from gene-tree / species-tree reconciliation.
  orthofinder -f primary_transcripts -t 80
fi
