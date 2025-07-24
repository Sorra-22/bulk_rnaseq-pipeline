#!/usr/bin/env python3
"""
Explore RNA-seq Results - Interactive Analysis
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

class RNASeqExplorer:
    def __init__(self, project_dir):
        self.project_dir = Path(project_dir)
        self.results_dir = self.project_dir / "results"
        
    def check_alignment_stats(self):
        """Parse and visualize alignment statistics"""
        print("=== ALIGNMENT STATISTICS ===\n")
        
        align_dir = self.results_dir / "aligned"
        summary_files = list(align_dir.glob("*_summary.txt"))
        
        stats = {}
        for file in summary_files:
            sample = file.stem.replace("_summary", "")
            with open(file, 'r') as f:
                content = f.read()
                # Extract alignment rate (last line usually contains overall alignment rate)
                for line in content.split('\n'):
                    if 'overall alignment rate' in line:
                        rate = float(line.split('%')[0].split()[-1])
                        stats[sample] = rate
                        
        # Create bar plot
        if stats:
            plt.figure(figsize=(8, 6))
            samples = list(stats.keys())
            rates = list(stats.values())
            
            bars = plt.bar(samples, rates)
            plt.ylabel('Alignment Rate (%)')
            plt.title('Alignment Rates by Sample')
            plt.ylim(0, 100)
            
            # Color bars based on quality (>90% green, 80-90% yellow, <80% red)
            for bar, rate in zip(bars, rates):
                if rate > 90:
                    bar.set_color('green')
                elif rate > 80:
                    bar.set_color('orange')
                else:
                    bar.set_color('red')
                    
            # Add value labels on bars
            for bar, rate in zip(bars, rates):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f'{rate:.1f}%', ha='center', va='bottom')
                        
            plt.tight_layout()
            plt.savefig(self.results_dir / "alignment_rates.png", dpi=150)
            plt.show()
            
            # Print summary
            print(f"Average alignment rate: {np.mean(rates):.2f}%")
            print(f"Range: {min(rates):.2f}% - {max(rates):.2f}%")
            
    def explore_count_distribution(self):
        """Analyze count distribution"""
        print("\n=== COUNT DISTRIBUTION ===\n")
        
        # Read count matrix
        counts_file = self.results_dir / "counts" / "merged_counts.csv"
        if not counts_file.exists():
            print("Count file not found yet. Pipeline still running?")
            return
            
        counts = pd.read_csv(counts_file)
        count_cols = [col for col in counts.columns if col not in ['Geneid', 'Length']]
        
        # Calculate basic stats
        count_data = counts[count_cols]
        
        print("Library sizes (total counts per sample):")
        lib_sizes = count_data.sum()
        for sample, size in lib_sizes.items():
            print(f"  {sample}: {size:,} reads")
            
        # Plot count distributions
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Library sizes
        ax = axes[0, 0]
        lib_sizes.plot(kind='bar', ax=ax)
        ax.set_title('Library Sizes')
        ax.set_ylabel('Total Counts')
        ax.tick_params(axis='x', rotation=45)
        
        # 2. Count distribution (log scale)
        ax = axes[0, 1]
        log_counts = np.log10(count_data + 1)
        log_counts.melt().groupby('value').size().plot(ax=ax)
        ax.set_xlabel('Log10(Count + 1)')
        ax.set_ylabel('Frequency')
        ax.set_title('Count Distribution (Log Scale)')
        
        # 3. Zeros per sample
        ax = axes[1, 0]
        zeros = (count_data == 0).sum()
        zeros.plot(kind='bar', ax=ax)
        ax.set_title('Number of Zero-Count Genes')
        ax.set_ylabel('Number of Genes')
        ax.tick_params(axis='x', rotation=45)
        
        # 4. Sample correlation heatmap
        ax = axes[1, 1]
        corr = count_data.corr()
        sns.heatmap(corr, annot=True, fmt='.3f', cmap='coolwarm', ax=ax)
        ax.set_title('Sample Correlation')
        
        plt.tight_layout()
        plt.savefig(self.results_dir / "count_exploration.png", dpi=150)
        plt.show()
        
        # Identify potential issues
        print("\n=== QUALITY CHECKS ===")
        
        # Check for outlier samples
        mean_corr = corr.mean()
        for sample in mean_corr.index:
            if mean_corr[sample] < 0.9:
                print(f"⚠️  {sample} has low correlation with other samples ({mean_corr[sample]:.3f})")
                
        # Check for library size imbalances
        max_size = lib_sizes.max()
        min_size = lib_sizes.min()
        if max_size / min_size > 2:
            print(f"⚠️  Large library size imbalance: {max_size/min_size:.1f}x difference")
            
    def preview_de_results(self):
        """Preview differential expression results if available"""
        print("\n=== DIFFERENTIAL EXPRESSION PREVIEW ===\n")
        
        de_file = self.results_dir / "differential_expression" / "deseq2_results.csv"
        if not de_file.exists():
            print("DE results not available yet. Run differential expression analysis first.")
            return
            
        # Read results
        de_results = pd.read_csv(de_file)
        
        # Summary statistics
        sig_genes = de_results[de_results['padj'] < 0.05]
        print(f"Total genes tested: {len(de_results)}")
        print(f"Significant genes (padj < 0.05): {len(sig_genes)}")
        print(f"  - Upregulated: {len(sig_genes[sig_genes['log2FoldChange'] > 0])}")
        print(f"  - Downregulated: {len(sig_genes[sig_genes['log2FoldChange'] < 0])}")
        
        # Top genes
        print("\nTop 10 most significant genes:")
        top_genes = de_results.nsmallest(10, 'padj')[['gene', 'log2FoldChange', 'pvalue', 'padj']]
        print(top_genes.to_string(index=False))
        
        # Quick volcano plot
        plt.figure(figsize=(8, 6))
        
        # Define colors
        de_results['color'] = 'gray'
        de_results.loc[(de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 1), 'color'] = 'red'
        de_results.loc[(de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < -1), 'color'] = 'blue'
        
        # Plot
        for color, data in de_results.groupby('color'):
            plt.scatter(data['log2FoldChange'], -np.log10(data['pvalue']), 
                       c=color, alpha=0.5, s=20)
                       
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10 P-value')
        plt.title('Volcano Plot')
        plt.axvline(x=1, color='black', linestyle='--', alpha=0.3)
        plt.axvline(x=-1, color='black', linestyle='--', alpha=0.3)
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / "quick_volcano.png", dpi=150)
        plt.show()
        
    def generate_report(self):
        """Generate a simple HTML report"""
        print("\n=== GENERATING SUMMARY REPORT ===\n")
        
        html_content = """
        <html>
        <head>
            <title>RNA-seq Analysis Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                h1 { color: #2c3e50; }
                h2 { color: #34495e; border-bottom: 2px solid #3498db; }
                .metric { background: #ecf0f1; padding: 10px; margin: 10px 0; }
                img { max-width: 600px; margin: 20px 0; }
            </style>
        </head>
        <body>
            <h1>RNA-seq Analysis Report</h1>
            <p>Generated: {date}</p>
            
            <h2>Analysis Overview</h2>
            <div class="metric">
                <strong>Project Directory:</strong> {project_dir}
            </div>
            
            <h2>Quality Control</h2>
            <p>See MultiQC report for detailed metrics: 
               <a href="multiqc/multiqc_report.html">Open MultiQC Report</a></p>
            
            <h2>Alignment Statistics</h2>
            <img src="alignment_rates.png" alt="Alignment Rates">
            
            <h2>Count Distribution</h2>
            <img src="count_exploration.png" alt="Count Analysis">
            
            <h2>Differential Expression</h2>
            <img src="quick_volcano.png" alt="Volcano Plot">
        </body>
        </html>
        """
        
        from datetime import datetime
        
        report_content = html_content.format(
            date=datetime.now().strftime("%Y-%m-%d %H:%M"),
            project_dir=str(self.project_dir)
        )
        
        report_path = self.results_dir / "analysis_report.html"
        with open(report_path, 'w') as f:
            f.write(report_content)
            
        print(f"Report saved to: {report_path}")
        print(f"Open with: firefox {report_path}")

def main():
    """Run exploration"""
    import sys
    
    project_dir = Path.home() / "rnaseq_project"
    explorer = RNASeqExplorer(project_dir)
    
    print("RNA-seq Results Explorer")
    print("=" * 50)
    
    # Check what's available
    if (project_dir / "results" / "aligned").exists():
        explorer.check_alignment_stats()
        
    if (project_dir / "results" / "counts" / "merged_counts.csv").exists():
        explorer.explore_count_distribution()
        
    if (project_dir / "results" / "differential_expression" / "deseq2_results.csv").exists():
        explorer.preview_de_results()
        
    # Generate report if results exist
    if (project_dir / "results").exists():
        explorer.generate_report()

if __name__ == "__main__":
    main()
