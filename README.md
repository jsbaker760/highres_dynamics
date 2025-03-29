Highly-resolved within-species dynamics in the human facial skin microbiome
=======================

This repository contains all code necessary to reproduce the analyses described in:

> __Jacob S. Baker, Evan Qu, Christopher P. Mancuso, A. Delphine Tripp, Arolyn Conwill, Tami D. Lieberman.__ Previously hidden intraspecies dynamics underlie the apparent stability of two important skin microbiome species
[BioRxiv](https://www.biorxiv.org/content/10.1101/2024.01.10.575018v2)

This project is part of my doctoral research in the [Lieberman Lab](http://lieberman.science).

# Overview
This repository is grouped into two major parts:

* [local_analysis](#all_figures)
  * local_analysis/all_figures.m calls all functions in local_analysis/scripts necessary to reproduce the figures in the manuscript, provided that pre-computed data structures are included.
  
* [snakemakes](#snakemake_raw_data)
  * These snakemakes were run in Slurm-based HPC environments to generate consolidated data structures from raw reads including genome assemblies and SNP tables, which were then analyzed locally (all_figures.m)
  
# data-availability
raw reads are available on the NCBI Sequence Read Archive under [Bioproject #PRJNA1052084](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1052084).

Pre-computed data structures small enough to be uploaded to GitHub are in local_analysis/data, including:
  * lineage coassemblies, annotations, and SNP tables (data/lineage_coassemblies)
  * lineage phylogenies (data/lineage_phylogenies)
  * Intraspecies abundances used in the manuscript (data/PHLAME_filtered, also in Table S3 in the paper)
  * bracken reports (data/metagenomics_bracken_reports)

The beta version of PHLAME used in this paper can be found [here](https://github.com/quevan/phlame_beta) (STAR Methods: Metagenomic lineage and phylotype abundances)
The current verion of PHLAME can be found [here](https://github.com/quevan/phlame) and is described in [this preprint](https://www.biorxiv.org/content/10.1101/2025.02.07.636498v1)

A full clone of data/ (36GB) is available on [DropBox](https://www.dropbox.com/scl/fo/2lartzjhtgubxvgamjcd7/AGw1K-xHXGwe99iwP_vWEZs?rlkey=dlor01adbbii9sgvctqu7w6uz&dl=0).

# all_figures

all_figures.m with [pre-computed data](#data-availability) will re-generate all figures in the manuscript:
    * Filters and analyzes metagenomics data (STAR Methods: Lineage diversity and sharing across subjects)
    * Filters and analyzes isolate data (STAR Methods: Filtering isolates)
    * clusters samples at the species level (STAR Methods: Species and genus-level microbiome composition)
    * generates distance matrices (STAR Methods: Calculating pairwise SNP distances for clustering lineage)
    * clusters isolate genomes into lineages with numerous figures not shown in the paper (STAR Methods: Clustering isolates into lineages)
    * identifies _de novo_ SNPs in each lineage and across homologous genes (STAR Methods: Within lineage SNPs)
    * generates lineage phylogenies (STAR Methods: Within lineage SNPs)
    * Measures multi-genotype transmissions, directionality, dMRCAs, dN/dS in lineage phylogenies, molecular clocks.
    * and more!!

# snakemake_raw_data
* "Snakemake 1"
  * 1a. Creates and annotates unfiltered genome assemblies from each isolate, which were used only to filter isolate genomes (STAR Methods:Isolating colonies and Filtering isolates)
  * 1b. Takes raw demultiplexed fastq reads from each isolate, performs basic filtering,
  * Aligns reads to Pacnes_C1 and SepidermidisATCC12228 reference genomes,
  * then, for both species, concatenates alignment data together into candidate mutation tables used in generating distance matrix (STAR Methods: Calculating pairwise SNP distances for clustering lineage)
    
* "Snakemake 2"
  * Takes raw demultiplexed fastq reads from each metagenomics samples and uses kraken2/bracken to estimate species and genus-level abundances for each sample (STAR Methods: Species and genus-level microbiome composition)
    
* "Snakemake 3"
  * Creates assemblies for each lineage (identified by clustering isolates, code in #all_figures)
  * Aligns reads from each isolate in each lineage to its respective reference genomes
  * Creates candidate mutation tables for each lineage
  * Annotates each co-assembly with Bakta (STAR Methods: Within lineage SNPs, Gene content within and across lineage assemblies)
