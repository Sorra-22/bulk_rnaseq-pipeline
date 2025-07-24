# RNA-seq Analysis Pipeline

A simple, educational RNA-seq analysis pipeline for bulk RNA sequencing data, designed for learning bioinformatics.

## Overview

This pipeline performs:
- Quality control (FastQC)
- Read trimming (Trimmomatic)
- Alignment (HISAT2)
- Quantification (featureCounts)
- Differential expression (DESeq2)

## Quick Start

### Prerequisites

- Ubuntu/Linux system
- Conda/Miniconda installed
- ~10GB free disk space
- 8GB+ RAM recommended

### Create conda environment:
conda env create -f environment.yml
conda activate rnaseq

### Install R packages:
Rscript scripts/install_r_packages.R

### Run Nextflow script:
-nextflow run

### Run the pipeline:
python scripts/rnaseq_pipeline.py

### Run differential expression:
python scripts/differential_expression.py
cd results/differential_expression
Rscript run_deseq2.R
