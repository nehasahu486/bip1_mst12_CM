## Description

This repository contains the analyses, scripts and files used for _M. oryzae_ Bip1 project.


## Table of Contents

- [Copy numbers analysis](analysis/001_orthologs.md) - Comparative genomics to check orthologs of Bip1 (MGG_08118) and Mst12 (MGG_12958) across 42 fungal species
- [RNA Seq Analysis](analysis/002_rna_seq.md) - Differentially expressed genes (DEGs) in bip1 and mst12 null mutants when compared to wilde-type (Guy11)
- [Co-expression of bip1 and mst12](analysis/003_co_expression.md) - Comparing DEGs in bip1 and mst12 null mutants when compared to wilde-type (Guy11)
- [GO enrichment](analysis/004_go_enrichment.md) - Significantly enriched GO terms in various comparisons
- [Network analysis](analysis/005_gene_networks.md) - Networks based on gene function annotations showing DEGs over time course in bip1 and mst12 vs Guy11 (WT)
- [WGS confirmation](analysis/006_wgs_knockout.md) - Confirmation of bip1 gene knockout
- [Phosphosite conservation](analysis/007_phossite_consv.md) - Conservation of phosphorylation sites in Bip1 in other fungi


### Data Structure

- `analysis/` - Detailed Markdown files
- `data/processed/` - Processed and intermediate files from analysis on the HPC or other sources
- `data/results/` - Analysis results, plots, and summary tables
- `scripts/` - Scripts and/or workflows used
- `bin/` - Miscellaneous stuff

## Pipeline Execution Environment

The scripts developed and executed on the High-Performance Computing (HPC) cluster at the Sainsbury Laboratory (TSL) and Norwich BioScience Institute (NBI)

Please note that these bash scripts (.sh file extensions) within the pipeline specify usage allocations tailored to the SLURM batch-queue environment at TSL and NBI. These may require adjustments for different computing platforms.

For questions about the analytical methods, code implementation, or adapting these workflows for your data, please contact: <neha.sahu.tsl@gmail.com>

---

*Computational framework developed for Bip1 project in Magnaporthe oryzae.*


