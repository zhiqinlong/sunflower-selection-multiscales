# Gene duplication (OrthoFinder)

Number of gene-duplication events in the gene family (orthogroup) of each *H. annuus* HA412-HOv2 gene, used as one
proxy for pleiotropic potential (Fig. 4B, Figs. S15–S17).

| Step | File | Where | What it does |
|---|---|---|---|
| 1 | `orthofinder.sh` | salal (`/data/zlong/duplications`) | Renames proteins, keeps the longest isoform per gene, runs OrthoFinder v3.1.3 |
| 2 | `duplications_per_gene.R` | local | Assigns each HA412-HOv2 gene the duplication count of its orthogroup → `gene_duplications.csv` |

## Proteomes

| Label | Genome | Proteins (longest isoform) |
|---|---|---|
| `Hann` | *H. annuus* HA412-HOv2 | 44,544 |
| `Harg` | *H. argophyllus* ARG2202, haplotype 1 | 46,913 |
| `Hpet` | *H. petiolaris* subsp. *neglectus*, haplotype 1 (Hneg2702) | 44,786 |

No reference genome exists for subsp. *petiolaris* or *fallax*, so *neglectus* represents *H. petiolaris*.

## OrthoFinder result (run `Results_Mar10_1`)

- 136,243 proteins; 91.9% assigned to 37,595 orthogroups.
- 20,289 duplication events (all with ≥ 50% support).
- `Duplications_per_Orthogroup.tsv` gives the count per orthogroup ("Duplications (50% support)").

## Genes with a duplication count

OrthoFinder finds duplication events by building a gene tree for each orthogroup and comparing it with the species
tree. A gene tree cannot be built from fewer than 4 sequences, so OrthoFinder reports duplications only for
orthogroups with at least 4 genes (8,461 orthogroups). With three species, most smaller orthogroups contain one gene
per species.

Genes in orthogroups without a duplication count were set to NA. The duplication comparisons therefore use the
**14,460** HA412-HOv2 genes in orthogroups with at least 4 genes (median 2, mean 3.3, max 136 duplication events), and
both the adaptive and the background gene sets are drawn from these genes.

The protein IDs (`Ha412HOChr1g…`) differ from the gene IDs used elsewhere (`Ha412HOChr01g…`), so
`duplications_per_gene.R` matches them by genomic position. The same table is built inside
`../transcriptome/pleiotropy_level_plot.R` (`ann_dup`); the two are identical.
