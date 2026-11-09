# *H. annuus* RNA-seq: tissue specificity (τ) and co-expression network

135 RNA-seq libraries (180 sequencing runs) from 5 tissues: anther, leaf, nectary, ovary, pistil.

| File | Where | What it does |
|---|---|---|
| `run_transcriptome.sh` | rorqual (SLURM) | Steps 0–4: sample lists → fastp → merge the two runs of the same library → STAR + RSEM → WGCNA network. The submit command for each step and what each step does are explained at the top of the script and above each step. |
| `WGCNA.R` | rorqual (called by step 4) | DESeq2 VST, WGCNA adjacency, top 5% (or 1%) of edges, igraph centralities |
| `tissue_specificity.R` | local | τ per gene from mean TPM per tissue |
| `pleiotropy_level_plot.R` | local | Compares τ, centrality and duplication across gene sets (Fig. 4A–C, Figs. S15–S17) |

Software: fastp 0.24.0, STAR 2.7.11b, RSEM 1.3.3, R 4.5.0, DESeq2 1.50.2, WGCNA 1.74, igraph 2.3.3.
Paths are hard-coded; change them before running.
