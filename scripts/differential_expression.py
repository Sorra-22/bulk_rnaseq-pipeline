#!/usr/bin/env python3
"""
Differential Expression Analysis using DESeq2 via rpy2
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

# Try to import rpy2, if not available, create R script instead
try:
    from rpy2.robjects import r, pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()
    USE_RPY2 = True
except ImportError:
    logging.warning("rpy2 not installed. Will create R script instead.")
    USE_RPY2 = False

class DifferentialExpression:
    def __init__(self, project_dir):
        self.project_dir = Path(project_dir)
        self.counts_dir = self.project_dir / "results" / "counts"
        self.de_dir = self.project_dir / "results" / "differential_expression"
        self.de_dir.mkdir(parents=True, exist_ok=True)
        
    def prepare_deseq2_input(self):
        """Prepare count matrix and sample info for DESeq2"""
        # Read count matrix
        count_matrix = pd.read_csv(self.counts_dir / "merged_counts.csv")
        
        # Read sample info
        sample_info = pd.read_csv(self.project_dir / "data" / "sample_info.csv")
        
        # Prepare count matrix (genes as rows, samples as columns)
        count_data = count_matrix.set_index('Geneid')
        count_data = count_data.drop('Length', axis=1)
        
        # Ensure column order matches sample info
        count_data = count_data[sample_info['sample_id'].tolist()]
        
        # Save prepared data
        count_data.to_csv(self.de_dir / "count_matrix_for_deseq2.csv")
        sample_info.to_csv(self.de_dir / "sample_info_for_deseq2.csv", index=False)
        
        return count_data, sample_info
    
    def create_r_script(self):
        """Create R script for DESeq2 analysis"""
        r_script = '''
# DESeq2 Analysis Script
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Set working directory
setwd("{de_dir}")

# Read data
count_data <- read.csv("count_matrix_for_deseq2.csv", row.names=1)
sample_info <- read.csv("sample_info_for_deseq2.csv")

# Ensure sample order matches
count_data <- count_data[, sample_info$sample_id]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_info,
    design = ~ condition
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "Mutant", "WT"))

# Add gene names and sort by p-value
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$pvalue),]

# Save results
write.csv(res_df, "deseq2_results.csv", row.names=FALSE)

# Create MA plot
pdf("MA_plot.pdf")
plotMA(res, ylim=c(-5,5))
dev.off()

# Create volcano plot
pdf("volcano_plot.pdf")
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")
ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue), color=significant)) +
    geom_point(alpha=0.6) +
    theme_minimal() +
    labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 P-value") +
    scale_color_manual(values=c("gray", "red"))
dev.off()

# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "normalized_counts.csv")

# Create heatmap of top 50 DE genes
top_genes <- head(res_df$gene[!is.na(res_df$padj) & res_df$padj < 0.05], 50)
if(length(top_genes) > 0) {{
    pdf("heatmap_top50_DE_genes.pdf", width=8, height=10)
    pheatmap(log2(normalized_counts[top_genes,] + 1),
             scale="row",
             clustering_distance_rows="correlation",
             clustering_distance_cols="correlation",
             show_rownames=TRUE,
             show_colnames=TRUE)
    dev.off()
}}

# Summary statistics
sink("deseq2_summary.txt")
cat("DESeq2 Analysis Summary\\n")
cat("======================\\n")
cat(sprintf("Total genes analyzed: %d\\n", nrow(dds)))
cat(sprintf("Significant genes (padj < 0.05): %d\\n", 
    sum(res_df$padj < 0.05, na.rm=TRUE)))
cat(sprintf("Upregulated genes: %d\\n", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm=TRUE)))
cat(sprintf("Downregulated genes: %d\\n", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm=TRUE)))
sink()

print("DESeq2 analysis completed!")
'''
        
        # Write R script
        r_script_formatted = r_script.format(de_dir=str(self.de_dir))
        r_script_path = self.de_dir / "run_deseq2.R"
        
        with open(r_script_path, 'w') as f:
            f.write(r_script_formatted)
            
        logging.info(f"R script created at: {r_script_path}")
        logging.info("Run it with: Rscript run_deseq2.R")
        
        return r_script_path
    
    def plot_results_python(self):
        """Create plots using Python (if DESeq2 results are available)"""
        results_file = self.de_dir / "deseq2_results.csv"
        
        if not results_file.exists():
            logging.warning("DESeq2 results not found. Run the R script first.")
            return
            
        # Read results
        results = pd.read_csv(results_file)
        
        # Create volcano plot
        plt.figure(figsize=(10, 8))
        
        # Define significance
        results['significant'] = (results['padj'] < 0.05) & (np.abs(results['log2FoldChange']) > 1)
        
        # Plot non-significant genes
        ns = results[~results['significant']]
        plt.scatter(ns['log2FoldChange'], -np.log10(ns['pvalue']), 
                   c='gray', alpha=0.5, s=20, label='Not significant')
        
        # Plot significant genes
        sig = results[results['significant']]
        plt.scatter(sig['log2FoldChange'], -np.log10(sig['pvalue']), 
                   c='red', alpha=0.8, s=30, label='Significant')
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10 P-value')
        plt.title('Volcano Plot')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Add threshold lines
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        
        plt.savefig(self.de_dir / 'volcano_plot_python.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create MA plot
        plt.figure(figsize=(10, 8))
        
        results['baseMean_log'] = np.log10(results['baseMean'] + 1)
        
        plt.scatter(results['baseMean_log'], results['log2FoldChange'], 
                   c='gray', alpha=0.5, s=20)
        plt.scatter(sig['baseMean_log'], sig['log2FoldChange'], 
                   c='red', alpha=0.8, s=30)
        
        plt.xlabel('Log10 Mean Expression')
        plt.ylabel('Log2 Fold Change')
        plt.title('MA Plot')
        plt.grid(True, alpha=0.3)
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        
        plt.savefig(self.de_dir / 'ma_plot_python.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Print summary
        print("\nDifferential Expression Summary:")
        print(f"Total genes: {len(results)}")
        print(f"Significant genes (padj < 0.05): {sum(results['padj'] < 0.05)}")
        print(f"Upregulated: {sum((results['padj'] < 0.05) & (results['log2FoldChange'] > 0))}")
        print(f"Downregulated: {sum((results['padj'] < 0.05) & (results['log2FoldChange'] < 0))}")
        
        # Top 10 DE genes
        top_genes = results.nsmallest(10, 'padj')[['gene', 'log2FoldChange', 'padj']]
        print("\nTop 10 differentially expressed genes:")
        print(top_genes.to_string(index=False))
        
    def run_analysis(self):
        """Run complete differential expression analysis"""
        # Prepare data
        count_data, sample_info = self.prepare_deseq2_input()
        
        # Create R script
        r_script = self.create_r_script()
        
        print(f"\nData prepared for DESeq2 analysis!")
        print(f"Count matrix shape: {count_data.shape}")
        print(f"\nTo run DESeq2 analysis, execute:")
        print(f"cd {self.de_dir}")
        print(f"Rscript run_deseq2.R")
        
        # If results exist, create Python plots
        self.plot_results_python()

def main():
    """Main function"""
    project_dir = Path.home() / "rnaseq_project"
    
    de_analysis = DifferentialExpression(project_dir)
    de_analysis.run_analysis()

if __name__ == "__main__":
    main()
