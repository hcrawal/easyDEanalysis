# easyDEanalysis

`easyDEanalysis` is a streamlined, user-friendly R package designed to automate RNA-Seq differential expression workflows. It serves as a unified wrapper around the two industry-standard methods: **DESeq2** and **edgeR**.

Instead of writing hundreds of lines of repetitive code, researchers can execute a complete differential expression analysis, generate tabular data summaries, and export diagnostic plots using a single function call.

## Key Features
- **Dual Pipeline Integration:** Seamlessly toggle between `DESeq2` and `edgeR` workflows using a single argument.
- **Automated Text Summaries:** Automatically generates structured text files summarizing the analysis parameters and mapping details.
- **Unified Visualizations:** Renders quality control histograms and diagnostic plots directly into a clean, ready-to-publish PDF report.

## Installation

You can install the development version of `easyDEanalysis` directly from GitHub using the following R commands:

```R
# If you don't have remotes installed:
install.packages("remotes")

# Install easyDEanalysis
remotes::install_github("your_github_username/easyDEanalysis")
```

## Quick Start Guide

To run a complete analysis, simply load your raw count matrix, specify the number of samples (in pair), and choose your preferred method:

```R
library(easyDEanalysis)

# Run the pipeline using DESeq2 or edgeR
DE.analysis(M = "DESeq2", IN="count_matrix.txt", S="3:3", OUT="DE_analysis", LF="1", FD="0.05")
#S="3:3" means that in the count matrix file, there are 3 samples from condition 1 (eg. control) and 3 samples from condition 2 (treatment or stressed/patient samples); LF is the log2Fold change threshold, FD is the threshold for adjusted p-value or FDR
DE.analysis(M = "edgeR", IN="count_matrix.txt", S="3:3", OUT="DE_analysis", LF="1", FD="0.05")

DE.analysis(M = "edgeR", IN="count_matrix.txt", S="3:3")
DE.analysis(M = "DESeq2", IN="count_matrix.txt", S="3:3")

```

### Expected Outputs
The package will automatically output three files directly to your working directory:
1. 'DE_analysis_DESEeq2.txt' or ,DE_analysis_edgeR.txt' : A complete text file with calculated differential expression of all input IDs/genes/sequence from input count matrix.
2. 'DE_analysis_DESEeq2_LFC_X_and_padj_y.txt' or 'DE_analysis_edgeR_LFC_X_and_FDR_y.txt': A result file with differential expression values of input IDs/genes/sequence satisfying the provided threshold.
3. 'DE_analysis_DESEeq2_plots.txt' and 'DE_analysis_edgeR_plots.txt': A unified PDF containing the diagnostic histograms and differential expression plots.
