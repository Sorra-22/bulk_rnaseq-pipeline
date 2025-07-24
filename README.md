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

### Project Structure
rnaseq-pipeline/
├── README.md
├── environment.yml
├── requirements.txt
├── scripts/
│   ├── rnaseq_pipeline.py
│   ├── differential_expression.py
│   ├── explore_results.py
│   ├── pdf_report_generator.py
│   ├── download_test_data.sh
│   └── install_r_packages.R
├── data/
│   └── sample_info.csv
├── reference/
│   └── README.md
└── results/
    └── README.md

### Prerequisites

- Ubuntu/Linux system
- Conda/Miniconda installed
- ~10GB free disk space
- 8GB+ RAM recommended

### Create conda environment:
- conda env create -f environment.yml
- conda activate rnaseq

### Install R packages:
- Rscript scripts/install_r_packages.R

### Run the pipeline:
- python scripts/rnaseq_pipeline.py

### Run differential expression:
- python scripts/differential_expression.py
- cd results/differential_expression
- Rscript run_deseq2.R


