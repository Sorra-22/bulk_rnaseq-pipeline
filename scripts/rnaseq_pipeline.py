#!/usr/bin/env python3
"""
Simple RNA-seq Pipeline for Yeast Data
Author: Your Name
Date: 2024
"""

import os
import subprocess
import pandas as pd
import glob
from pathlib import Path
import logging
import json

# Set up logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(levelname)s - %(message)s')

class RNASeqPipeline:
    def __init__(self, project_dir):
        """Initialize pipeline with project directory structure"""
        self.project_dir = Path(project_dir)
        self.data_dir = self.project_dir / "data"
        self.results_dir = self.project_dir / "results"
        self.reference_dir = self.project_dir / "reference"
        
        # Create results subdirectories
        self.fastqc_dir = self.results_dir / "fastqc"
        self.trimmed_dir = self.results_dir / "trimmed"
        self.aligned_dir = self.results_dir / "aligned"
        self.counts_dir = self.results_dir / "counts"
        
        # Create directories if they don't exist
        for dir in [self.fastqc_dir, self.trimmed_dir, self.aligned_dir, self.counts_dir]:
            dir.mkdir(parents=True, exist_ok=True)
            
        # Reference files
        self.genome_fa = self.reference_dir / "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
        self.genome_gtf = self.reference_dir / "Saccharomyces_cerevisiae.R64-1-1.104.gtf"
        
    def run_command(self, cmd, description=""):
        """Run shell command and log output"""
        logging.info(f"Running: {description}")
        logging.info(f"Command: {cmd}")
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error(f"Error in {description}: {result.stderr}")
            raise RuntimeError(f"Command failed: {cmd}")
        
        return result
    
    def run_fastqc(self, sample_name, fastq_file):
        """Run FastQC quality control"""
        cmd = f"fastqc {fastq_file} -o {self.fastqc_dir}"
        self.run_command(cmd, f"FastQC for {sample_name}")
        
    def trim_reads(self, sample_name, fastq_file):
        """Trim adapters and low quality bases using Trimmomatic"""
        output_file = self.trimmed_dir / f"{sample_name}_trimmed.fastq.gz"
        
        # Basic trimming parameters for single-end reads
        cmd = f"""
        trimmomatic SE -threads 4 \
            {fastq_file} \
            {output_file} \
            ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:36
        """
        
        self.run_command(cmd, f"Trimming {sample_name}")
        return output_file
    
    def build_hisat2_index(self):
        """Build HISAT2 index for genome"""
        index_base = self.reference_dir / "yeast_index"
        
        # Check if index already exists
        if (index_base.parent / f"{index_base.name}.1.ht2").exists():
            logging.info("HISAT2 index already exists, skipping...")
            return str(index_base)
        
        cmd = f"hisat2-build {self.genome_fa} {index_base}"
        self.run_command(cmd, "Building HISAT2 index")
        return str(index_base)
    
    def check_samtools(self):
        """Check if samtools is working properly"""
        try:
            result = subprocess.run("samtools --version", shell=True, 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                logging.info(f"Samtools version: {result.stdout.split()[1]}")
                return True
            else:
                logging.error("Samtools not working properly")
                return False
        except Exception as e:
            logging.error(f"Error checking samtools: {e}")
            return False
    
    def align_reads(self, sample_name, trimmed_fastq):
        """Align reads using HISAT2"""
        # Check samtools first
        if not self.check_samtools():
            raise RuntimeError("Samtools is not working. Please reinstall it.")
            
        index_base = self.build_hisat2_index()
        sam_file = self.aligned_dir / f"{sample_name}.sam"
        bam_file = self.aligned_dir / f"{sample_name}.bam"
        sorted_bam = self.aligned_dir / f"{sample_name}_sorted.bam"
        
        # Run alignment
        cmd = f"""
        hisat2 -x {index_base} \
            -U {trimmed_fastq} \
            -S {sam_file} \
            --summary-file {self.aligned_dir}/{sample_name}_summary.txt
        """
        self.run_command(cmd, f"Aligning {sample_name}")
        
        # Convert SAM to BAM using pipe (more efficient)
        cmd = f"samtools view -bS {sam_file} -o {bam_file}"
        self.run_command(cmd, f"Converting SAM to BAM for {sample_name}")
        
        # Sort BAM file
        cmd = f"samtools sort {bam_file} -o {sorted_bam}"
        self.run_command(cmd, f"Sorting BAM for {sample_name}")
        
        # Index BAM file
        cmd = f"samtools index {sorted_bam}"
        self.run_command(cmd, f"Indexing BAM for {sample_name}")
        
        # Clean up intermediate files
        sam_file.unlink()
        bam_file.unlink()
        
        return sorted_bam
    
    def count_reads(self, sample_name, sorted_bam):
        """Count reads using featureCounts"""
        count_file = self.counts_dir / f"{sample_name}_counts.txt"
        
        cmd = f"""
        featureCounts -a {self.genome_gtf} \
            -o {count_file} \
            -t exon \
            -g gene_id \
            {sorted_bam}
        """
        
        self.run_command(cmd, f"Counting reads for {sample_name}")
        return count_file
    
    def merge_counts(self):
        """Merge all count files into a single matrix"""
        count_files = list(self.counts_dir.glob("*_counts.txt"))
        
        # Read first file to get gene names
        first_df = pd.read_csv(count_files[0], sep='\t', comment='#')
        count_matrix = first_df[['Geneid', 'Length']].copy()
        
        # Add counts from all samples
        for count_file in count_files:
            sample_name = count_file.stem.replace('_counts', '')
            df = pd.read_csv(count_file, sep='\t', comment='#')
            # The last column contains the counts
            count_matrix[sample_name] = df.iloc[:, -1]
        
        # Save merged counts
        output_file = self.counts_dir / "merged_counts.csv"
        count_matrix.to_csv(output_file, index=False)
        logging.info(f"Merged count matrix saved to {output_file}")
        
        return count_matrix
    
    def run_pipeline(self, sample_info_file):
        """Run complete pipeline for all samples"""
        # Read sample information
        samples_df = pd.read_csv(sample_info_file)
        
        # Process each sample
        for idx, row in samples_df.iterrows():
            sample_name = row['sample_id']
            fastq_file = self.data_dir / row['fastq_file']
            
            logging.info(f"\n{'='*50}")
            logging.info(f"Processing sample: {sample_name}")
            logging.info(f"{'='*50}")
            
            # 1. Quality control
            self.run_fastqc(sample_name, fastq_file)
            
            # 2. Trim reads
            trimmed_fastq = self.trim_reads(sample_name, fastq_file)
            
            # 3. Align reads
            sorted_bam = self.align_reads(sample_name, trimmed_fastq)
            
            # 4. Count reads
            count_file = self.count_reads(sample_name, sorted_bam)
            
        # 5. Merge all counts
        count_matrix = self.merge_counts()
        
        # 6. Run MultiQC
        cmd = f"multiqc {self.results_dir} -o {self.results_dir}/multiqc"
        self.run_command(cmd, "Running MultiQC")
        
        logging.info("\nPipeline completed successfully!")
        logging.info(f"Results are in: {self.results_dir}")
        
        return count_matrix

def main():
    """Main function to run the pipeline"""
    # Set up project directory
    project_dir = os.path.expanduser("~/rnaseq_project")
    
    # Initialize pipeline
    pipeline = RNASeqPipeline(project_dir)
    
    # Run pipeline
    sample_info = os.path.join(project_dir, "data", "sample_info.csv")
    count_matrix = pipeline.run_pipeline(sample_info)
    
    # Print summary statistics
    print("\nCount Matrix Summary:")
    print(f"Total genes: {len(count_matrix)}")
    print(f"Samples processed: {len(count_matrix.columns) - 2}")  # Minus Geneid and Length columns
    print("\nFirst few genes:")
    print(count_matrix.head())

if __name__ == "__main__":
    main()
