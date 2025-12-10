################################################################################
# Ecklonia radiata RNA-seq Differential Expression Analysis
# 
# Author: Olivia Wynn
# Date: 2024-2025
# 
# Purpose: 
# Complete pipeline for RNA-seq analysis of Ecklonia radiata responses to 
# marine heatwaves. Includes data import, quality control, differential 
# expression analysis, and GO enrichment.
#
# Workflow:
# 1. Data import (Salmon quantification + metadata)
# 2. Data preparation and quality control
# 3. Differential expression analysis (DESeq2)
# 4. PCA and exploratory analysis
# 5. GO enrichment analysis
# 6. Visualization
#
################################################################################

# ==============================================================================
# SECTION 1: LOAD REQUIRED PACKAGES
# ==============================================================================

# General data handling
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(reshape2)
library(readxl)
library(writexl)
library(tximport)

# Bioconductor packages
library(BiocManager)
library(DESeq2)
library(edgeR)
library(apeglm)
library(ashr)

# Quality control & visualization
library(PCAtools)
library(ggridges)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)

# GO enrichment
library(goseq)
library(fgsea)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)

# Other utilities
library(vegan)
library(rstatix)
library(patchwork)

# ==============================================================================
# SECTION 2: DATA IMPORT
# ==============================================================================

## 2.1 Import sample metadata
# Load sample information including treatment groups and replicates
sampleinfo.df <- read_tsv("samples.txt")

# Data preparation: Set up proper column structure
headers <- colnames(sampleinfo.df)
sampleinfo.df <- rbind(headers, sampleinfo.df)
colnames(sampleinfo.df) <- paste0("V", seq_len(ncol(sampleinfo.df)))
colnames(sampleinfo.df)[1:2] <- c("group", "sample")
sampleinfo.df <- sampleinfo.df[, 1:(ncol(sampleinfo.df) - 2)]

print(sampleinfo.df)

## 2.2 Import transcript-to-gene mapping
# Load TransDecoder gene-transcript mapping file
tx2gene <- read.table("transdecoder_gene_trans_map.txt", 
                      header = FALSE, 
                      col.names = c("transcript_id", "gene_id"))

# Clean gene IDs by removing isoform suffixes (_i*)
tx2gene <- tx2gene %>%
  mutate(transcript_id = gene_id,
         gene_id = gsub("_i[0-9]+", "", gene_id))

head(tx2gene)

## 2.3 Import Salmon quantification files
# Define file paths for Salmon output
files <- file.path("salmon", sampleinfo.df$sample, "quant.sf")
names(files) <- sampleinfo.df$sample

# Import gene-level counts using tximport
heatwave_gene_counts.df <- tximport(files, 
                                    type = "salmon", 
                                    tx2gene = tx2gene, 
                                    countsFromAbundance = "no")

# Convert to data frame for downstream analysis
heatwave_gene_counts.df <- as.data.frame(heatwave_gene_counts.df)

# ==============================================================================
# SECTION 3: EXPLORATORY DATA ANALYSIS - Count Distribution
# ==============================================================================

## 3.1 Visualize count distribution for one sample
# Check if data follows expected RNA-seq distribution
ggplot(heatwave_gene_counts.df) +
  geom_histogram(aes(x = `abundance.C1.2`), stat = "bin", bins = 200) +
  xlim(-5, 500) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  theme_bw()

## 3.2 Mean-variance relationship
# Calculate mean and variance across all samples
mean_counts <- apply(heatwave_gene_counts.df[, 2:31], 1, mean)
variance_counts <- apply(heatwave_gene_counts.df[, 2:31], 1, var)
check.df <- data.frame(mean_counts, variance_counts)

# Plot mean vs variance to determine appropriate statistical model
ggplot(check.df) +
  geom_point(aes(x = mean_counts, y = variance_counts)) + 
  geom_line(aes(x = mean_counts, y = mean_counts, color = "red")) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()

# Interpretation: Variance > mean indicates negative binomial distribution
# This confirms DESeq2 is appropriate (assumes negative binomial)

# ==============================================================================
# SECTION 4: DATA PREPARATION FOR DESeq2
# ==============================================================================

## 4.1 Prepare metadata
sampleinfo.df <- as.data.frame(sampleinfo.df)
rownames(sampleinfo.df) <- sampleinfo.df$sample

# Add experimental design factors
sampleinfo.df$timepoint <- factor(
  c(rep('during', 5), rep('post', 5),
    rep('during', 5), rep('post', 5),
    rep('during', 5), rep('post', 5)),
  levels = c("during", "post")
)

sampleinfo.df$heatwave <- factor(
  c(rep('control', 10),
    rep('single', 10),
    rep('double', 10)),
  levels = c("control", "single", "double")
)

## 4.2 Extract count matrix
# Separate raw counts from Salmon output
counts_matrix <- heatwave_gene_counts.df[, grep("counts", colnames(heatwave_gene_counts.df))]
colnames(counts_matrix) <- gsub("counts\\.", "", colnames(counts_matrix))
colnames(counts_matrix) <- gsub("\\.", "-", colnames(counts_matrix))

# Remove non-numeric column and round to integers
counts_matrix <- counts_matrix %>%
  dplyr::select(-`countsFromAbundance`)
counts_matrix <- round(as.matrix(counts_matrix))
counts_matrix <- as.data.frame(counts_matrix)

# Verify sample order matches metadata
all(colnames(counts_matrix) == rownames(sampleinfo.df))

## 4.3 Extract abundance matrix (for reference)
abundance_matrix <- heatwave_gene_counts.df[, grep("abundance", colnames(heatwave_gene_counts.df))]
colnames(abundance_matrix) <- gsub("abundance\\.", "", colnames(abundance_matrix))
colnames(abundance_matrix) <- gsub("\\.", "-", colnames(abundance_matrix))
abundance_matrix <- as.matrix(abundance_matrix)

## 4.4 Extract length matrix (for reference)
length_matrix <- heatwave_gene_counts.df[, grep("length", colnames(heatwave_gene_counts.df))]
colnames(length_matrix) <- gsub("length\\.", "", colnames(length_matrix))
colnames(length_matrix) <- gsub("\\.", "-", colnames(length_matrix))

# ==============================================================================
# SECTION 5: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS
# ==============================================================================

## 5.1 Create DESeqDataSet object
# Design formula includes main effects and interaction
heatwave.dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sampleinfo.df,
  design = ~ timepoint + heatwave + timepoint:heatwave
)

## 5.2 Normalization
# Apply median-of-ratios normalization
heatwave.dds <- estimateSizeFactors(heatwave.dds)

## 5.3 Filtering low-count genes
# Keep genes with at least 5 counts in at least 2 samples
# This removes noise and improves multiple testing correction
keep <- rowSums(counts(heatwave.dds) >= 5) >= 2
heatwave.dds_filtered <- heatwave.dds[keep, ]

cat("Genes after filtering:", nrow(heatwave.dds_filtered), "\n")

## 5.4 Variance stabilizing transformation
# VST for visualization and clustering (not for DE testing)
heatwave.vsd <- vst(heatwave.dds_filtered, blind = FALSE)
heatwave.vsd.mat <- assay(heatwave.vsd)

# ==============================================================================
# SECTION 6: QUALITY CONTROL - PCA ANALYSIS
# ==============================================================================

## 6.1 Perform PCA using PCAtools
heatwave.PCA_tools <- PCAtools::pca(
  heatwave.vsd.mat, 
  metadata = sampleinfo.df, 
  removeVar = 0.1
)

# Extract explained variance
explained_variance <- heatwave.PCA_tools$variance
print(explained_variance)

## 6.2 Determine optimal number of PCs
# Use Horn's parallel analysis and elbow method
heatwave_horn <- parallelPCA(heatwave.vsd.mat)
heatwave_elbow <- findElbowPoint(heatwave.PCA_tools$variance)

# Create scree plot
heatwave.PCA_tools.screeplot <- screeplot(
  heatwave.PCA_tools,
  components = getComponents(heatwave.PCA_tools, 1:20),
  vline = c(heatwave_horn$n, heatwave_elbow)
) +
  geom_label(aes(x = heatwave_horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 12)) +
  geom_label(aes(x = heatwave_elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 12)) +
  theme_bw()

print(heatwave.PCA_tools.screeplot)

## 6.3 Create PCA biplot
# Prepare data for plotting
pca_data_for_plotting <- as.data.frame(heatwave.PCA_tools$rotated)
pca_data_for_plotting$heatwave <- heatwave.PCA_tools$metadata$heatwave
pca_data_for_plotting$timepoint <- heatwave.PCA_tools$metadata$timepoint

# Define color palette (blue = cool, red = heat)
heatwave_colors <- c(
  "control" = "#1F78B4",
  "single"  = "#FF7F00",
  "double"  = "#E31A1C"
)

# Generate biplot
heatwave.biplot <- ggplot(pca_data_for_plotting, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(group = heatwave, fill = heatwave),
               geom = "polygon", alpha = 0.3, level = 0.95, type = 't') +
  geom_point(aes(color = heatwave, shape = timepoint), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.25) +
  scale_fill_manual(values = heatwave_colors, name = "Heatwave Treatment") +
  scale_colour_manual(values = heatwave_colors, name = "Heatwave Treatment") +
  scale_shape_manual(values = c(16, 17), name = "Timepoint",
                     labels = c("During", "Post")) +
  labs(x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
       y = paste0("PC2 (", round(explained_variance[2], 1), "%)")) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

print(heatwave.biplot)

# ==============================================================================
# SECTION 7: DIFFERENTIAL EXPRESSION TESTING
# ==============================================================================

## 7.1 Run DESeq2 analysis
# This performs dispersion estimation and statistical testing
heatwave.dds_filtered <- DESeq(heatwave.dds_filtered)

## 7.2 Extract results for specific contrasts
# Example: Single heatwave vs Control during heatwave exposure
results_single_vs_control_during <- results(
  heatwave.dds_filtered,
  contrast = c("heatwave", "single", "control"),
  alpha = 0.05
)

# Shrink log2 fold changes for visualization
results_single_vs_control_during_shrunk <- lfcShrink(
  heatwave.dds_filtered,
  contrast = c("heatwave", "single", "control"),
  res = results_single_vs_control_during,
  type = "ashr"
)

## 7.3 Summary of results
summary(results_single_vs_control_during_shrunk)

# Extract significant genes (padj < 0.05)
sig_genes <- results_single_vs_control_during_shrunk[
  !is.na(results_single_vs_control_during_shrunk$padj) & 
    results_single_vs_control_during_shrunk$padj < 0.05, 
]

cat("Number of significant genes:", nrow(sig_genes), "\n")

# ==============================================================================
# SECTION 8: VISUALIZATION - VOLCANO PLOT
# ==============================================================================

volcano_plot <- EnhancedVolcano(
  results_single_vs_control_during_shrunk,
  lab = rownames(results_single_vs_control_during_shrunk),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Single Heatwave vs Control',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.5,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0
)

print(volcano_plot)


################################################################################
# END OF PIPELINE
################################################################################