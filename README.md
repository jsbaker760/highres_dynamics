Highly-resolved within-species dynamics in the human facial skin microbiome
=======================
The manuscript corrosponding to this repo has not been peer reviewed, and this repository is likely to be ammended before publication
=======================

This repository contains all code necessary to reproduce the analyses described in:

> __Jacob S. Baker, Evan Qu, Christopher P. Mancuso, A. Delphine Tripp, Arolyn Conwill, Tami D. Lieberman.__ Highly-resolved within-species dynamics in the human facial skin microbiome

This project is part of my doctoral research in the [Lieberman Lab](http://lieberman.science).

# Overview
This repository is grouped into two major parts:

* [all_figures](#all_figures)
  * Running all_figures.m with [pre-computed data](#data-availability) will re-generate all figures in the manuscript:
    * Filters applied to metagenomics data
    * filters applied to isolates
    * distance matrix generation (/scripts/distance_matrix)
    * clustering isolate genomes into lineages (run_clustering.m)
    * identifying _de novo_ SNPs in each lineage and across homologues
    * generating phylogenies
    * clustering samples at the species level with UMAP
    * and more!
  * /PHLAME_filtered contains the within-species abundance information for metagenomics samples used in the manuscript
  * /PHLAME contains the (old) within-species abundance information used to filter out metagenomics samples
  * /metadata contains basic information on each sample
  
* [snakemake_raw_data](#snakemake_raw_data)
  * Running the Snakemakes in these folders with [raw data](#data-availability) will generate [pre-computed data](#data-availability) used in all_figures.
  
See below for information on [Data availability](#data-availability)


# all_figures
* all_figures.m calls all functions necessary for generated each figure in the manuscript,
  provided that pre-computed data structures are included.
* /scripts includes all functions called by all_figures.m


# snakemake_raw_data
* "Snakemake 1"
  * Takes raw demultiplexed fastq reads from each isolate, performs basic filtering
  * Aligns reads from each isolate to Pacnes_C1 and SepidermidisATCC12228 reference genomes
  * For both species, concatenates alignment data together into candidate mutation tables used in generating distance matrix
  * Creates and annotates genomes assemblies from each isolate to use in filtering
    
* "Snakemake 2"
  * Takes raw demultiplexed fastq reads from each metagenomics samples and performs basic filtering
  * Aligns reads from each metagenomics samples to Pacnes_C1 and SepidermidisATCC12228 reference genomes
  * For both species, concatenates alignment data together into candidate mutation tables used in generating distance matrix
  * Used kraken2/bracken to estimate species and genus-level abundances for each sample
    
* "Snakemake 3"
  * Creates assemblies for each lineage (identified by clustering isolates, code in #all_figures)
  * Aligns reads from each isolate in each lineage to its respective reference genomes
  * Creates candidate mutation tables for each lineage
  * Annotates each co-assembly with Bakta

# data-availability
Sequencing data for running Snakemakes is available on the NCBI Sequence Read Archive under Bioproject #SUB14054091.

Pre-computed MATLAB data structures necessary and sufficient to run all code in all_figures is available here.
