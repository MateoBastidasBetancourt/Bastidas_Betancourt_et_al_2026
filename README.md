# Transcriptomic and epigenomic analysis of cerebral organoids and tissue (bulkRNA-seq, scRNA-seq, ChIP-seq and CUT&RUN-seq) from Bastidas-Betancourt et al. (2026)

This repository contains the R and bash scripts used to process and analyze bulkRNA-seq, scRNA-seq, ChIP-seq and CUT&RUN-seq datasets generated from human and rhesus macaque fetal cortical tissue and cerebral organoids in Bastidas-Betancourt et al. (2026).
## Overview

The repository includes workflows for comparative transcriptomic and epigenomic analysis between human and rhesus macaque aiming to identify novel candidate genes and networks in ape brain evolution.

## Repository structure
```text
├── preprocessing/
│   ├── CnR-seq/
│   │   ├── 01_Slurm_qc_trim_mapping.sh
│   │   ├── 02_Slurm_peak_calling.sh
│   │   ├── 03_Slurm_HOMER_repeatmasker.sh
│   │   └── 04_Slurm_HOMER_TE_analysis.sh
│   │
│   └── Bulk-RNAseq/
│       ├── 01_Slurm_index_liftoff.sh
│       └── 02_Slurm_fastqc_trim.sh
│
├── ChIP-seq reanalysis/
│   ├── 01_H3K27Ac_H3K4Me2_human_rhesus.R
│   └── 02_H3K27Ac_bed_genome_tracks.R
├── sc-RNAseq reanalysis/
│   ├── 01_Pseudobulk_Plots_human_organoids_Tynianskia.R
├── bulk-RNAseq/
│   ├── 01_Human_npc_neu.R
│   ├── 02_Human_rhesus_npc_neu.R
│   └── 03_Rhesus_electroporations.R
├── CnR-seq/
│   ├── 01_Peak_annotation.R
│   ├── 02_TE_annotation.R
│   └── 03_Plotting.R
├── Metadata/
│   ├── pheno_human_npcs.tab
│   ├── pheno_human_rhesus_npcs.tab
│   └── pheno_rhesus_electroporations.tab
│
├── renv.lock
└── README.md
```
## System Requirements
## Software

R 4.5.1

renv (for package version control)

edgeR

Seurat

ggplot2

dplyr

Package versions are recorded in renv.lock.

Operating System: Linux (HPC cluster environment and local testing)

## Hardware
Downstream R analyses can run on a standard workstation (≥16 GB RAM recommended).

Preprocessing scripts require an HPC environment with a SLURM job scheduler.

## Installation Guide

Clone the repository and restore the R environment:
```text
install.packages("renv")
renv::activate()
renv::restore() ##make sure to run after setting working directory where renv.lock file is
```
Estimated installation time: 15 minutes on a standard desktop computer with internet access.


## Instructions for Use
Use raw reads and/or processed data stored in Zenodo
- Bulk RNA-seq analysis

Use raw data as indicated in Table S1.

Run preprocessing scripts in ./preprocessing/Bulk-RNAseq/*

Run R scripts in ./Bulk-RNAseq/ in numerical order

- ChIP-seq reanalysis

Use processed data as indicated in Table S1.

Run scripts contained in ./ChIP-seq reanalysis/*

- scRNA-seq reanalysis

Use processed data from Tynianskaia et al., 2026

Run plotting scripts contained in ./sc-RNAseq reanalysis/01_Plots_human_organoids_Tynianskia.R

- CUNT&RUN-seq analysis

Run preprocessing scripts in ./preprocessing/CnR-seq/*

Run R scripts in ./CnR-seq/ in numerical order

