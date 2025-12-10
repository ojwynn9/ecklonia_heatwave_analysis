# Ecklonia radiata RNA-seq Analysis Pipeline

Author: Olivia Wynn  
Institution: Institute for Marine and Antarctic Studies (IMAS), University of Tasmania  
Project: Transcriptomic responses of Ecklonia radiata to marine heatwaves

---

## Overview

This repository contains a complete bioinformatics pipeline for RNA-seq differential expression analysis of the kelp *Ecklonia radiata* under marine heatwave conditions. The workflow processes raw Illumina sequencing data through quality control, de novo assembly, annotation, quantification, and statistical analysis.

**Technical approach:**
- De novo transcriptome assembly (no reference genome available)
- Strand-specific library preparation (TruSeq Stranded mRNA, dUTP method)
- Redundancy reduction via CD-HIT clustering
- Comprehensive functional annotation (BLAST, Pfam, EggNOG, GO)
- Negative binomial statistical modeling (DESeq2)
- Gene set enrichment analysis (GSEA)

---

## Biological Context

*Ecklonia radiata* is a dominant kelp species in temperate Australian waters experiencing increasing marine heatwave frequency and intensity. This project investigated transcriptomic responses to:

- **Single heatwave exposure** (acute thermal stress)
- **Repeated heatwave exposure** (cumulative stress)  
- **Recovery periods** (post-stress acclimation)

**Experimental design:** 3 treatments × 2 timepoints × 5 biological replicates = 30 samples

---

## Pipeline Workflow

### Bash Pipeline (rnaseq_pipeline.sh)

Handles all steps from raw reads through to annotated, quantified transcriptome.

**1. Quality Control & Trimming**
- FastQC v0.12.0: Initial and post-trim quality assessment
- Fastp v0.23.2: Adapter removal, quality filtering (min length 50 bp)
- Result: 93.4% read retention (1.135 billion from 1.215 billion)

**2. De Novo Assembly**
- Trinity v2.15.1: Strand-specific assembly (--SS_lib_type RF)
- Bowtie2 v2.5.1: Transcript coverage assessment

**3. Redundancy Reduction**
- CD-HIT-EST v4.8.1: Clustering at 95% sequence identity

**4. Quality Assessment**
- BUSCO v5.5.0: Completeness evaluation (Stramenopiles lineage)
- Trinity stats: N50, transcript count, length distribution

**5. ORF Prediction**
- TransDecoder v5.7.1: Identify coding regions
- BLASTP: Homology search against SwissProt (e-value 1e-5)
- HMMER: Pfam domain search (e-value 1e-5)
- Predictions informed by both homology and domain evidence

**6. Quantification**
- Salmon: Pseudo-alignment via Trinity framework
- Gene-level count aggregation using transcript-to-gene mapping

**7. Functional Annotation**
- Trinotate: Comprehensive annotation pipeline
  - BLASTx and BLASTp vs SwissProt (e-value 1e-5)
  - Pfam domain identification
  - EggNOG ortholog mapping
  - GO term assignment
- Taxonomic filtering: Exclude bacteria (taxon 2), viruses (10239), vertebrates (7742)
- Result: 10,421 genes with GO term annotations

### R Pipeline (ecklonia_deseq2_pipeline.R)

Handles statistical analysis and visualization of differential expression.

**1. Data Import**
- Import Salmon gene-level counts via tximport
- Load sample metadata and experimental design
- Transcript-to-gene mapping integration

**2. Quality Control**
- Count distribution assessment (negative binomial validation)
- Mean-variance relationship examination
- Variance stabilizing transformation (VST)
- PCA for sample clustering and outlier detection

**3. Differential Expression**
- DESeq2 v1.46.0: Normalization (median-of-ratios method)
- Low-count filtering: Retain genes with ≥5 counts in ≥2 samples
- Statistical testing: Negative binomial GLM
- Model: ~timepoint + heatwave + timepoint:heatwave
- Refined model using grouping factor for pairwise contrasts
- Multiple testing correction: Benjamini-Hochberg FDR < 0.05
- Log2FC shrinkage via apeglm

**4. Multivariate Analysis**
- Bray-Curtis dissimilarity matrices from VST-normalized data
- PERMDISP: Homogeneity of multivariate dispersion
- PERMANOVA (adonis2): Group differences (10,000 permutations)

**5. Pairwise Comparisons**
- Heatwave exposure: S1 vs C1, D1 vs C1
- Recovery: S2 vs S1, D2 vs D1
- Single vs repeated exposure: D2 vs S2, D1 vs S1
- Timepoint effects: C1 vs C2, S1 vs S2, D1 vs D2

**6. Visualization**
- PCA biplots with treatment ellipses
- Volcano plots for each contrast
- Hierarchical clustering heatmaps (top 500 most variable DEGs)
- K-means clustering (k=4) to identify co-expressed gene modules
- ComplexHeatmap: Log-normalized counts, row-scaled

**7. Gene Set Enrichment Analysis**
- fgsea v1.34.2: Pre-ranked GSEA
- Genes ranked by log2 fold change from DESeq2
- Background: 10,421 annotated genes with GO terms
- Parameters: minSize=3, maxSize=500, eps=0
- Benjamini-Hochberg adjusted p-values
- Analysis performed on GO:BP terms (S1 vs C1) and full GO set (other contrasts)
- Visualization: Dot plots showing top 30 terms by normalized enrichment score (NES)
- Manual categorization of enriched terms (e.g., Photosynthesis, Oxidative stress, Molecular chaperones)

---

## Repository Structure

```
.
├── rnaseq_pipeline.sh              # Bash: QC → Assembly → Annotation
├── ecklonia_deseq2_pipeline.R      # R: DE analysis → Visualization → GSEA
├── README.md                        # This file
├── samples.txt                      # Sample metadata
├── transdecoder_gene_trans_map.txt # Transcript-gene mapping
└── go_annotations.txt               # GO term annotations
```

---

## Software Requirements

### Bash Pipeline

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.0 | Quality control |
| Fastp | 0.23.2 | Read trimming |
| Trinity | 2.15.1 | De novo assembly |
| Bowtie2 | 2.5.1 | Alignment |
| CD-HIT-EST | 4.8.1 | Redundancy reduction |
| TransDecoder | 5.7.1 | ORF prediction |
| BUSCO | 5.5.0 | Completeness assessment |
| Salmon | (via Trinity) | Quantification |
| BLAST+ | (compatible) | Homology search |
| HMMER | (compatible) | Domain search |
| Trinotate | (compatible) | Annotation |
| EggNOG-mapper | (compatible) | Ortholog mapping |

### R Pipeline

**R version: 4.4.2**

**Core packages:**
```r
# Data handling
tidyverse, dplyr, tibble, tidyr, tximport

# Differential expression
DESeq2 (1.46.0), edgeR, apeglm, ashr

# Quality control & visualization
PCAtools, ggplot2 (3.5.2), pheatmap, ComplexHeatmap, EnhancedVolcano

# Multivariate analysis
vegan (2.7.1)

# GO enrichment
fgsea (1.34.2), GO.db (3.21.0), AnnotationDbi (1.70.0)
```

**Install Bioconductor packages:**
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "apeglm", "PCAtools", 
                       "ComplexHeatmap", "fgsea", "GO.db", "AnnotationDbi"))
```

---

## Quick Start

### Running the Bash Pipeline

```bash
# 1. Configure paths in rnaseq_pipeline.sh
nano rnaseq_pipeline.sh
# Edit: RAW_READS_DIR, PROCESSED_DIR, database paths

# 2. Make executable
chmod +x rnaseq_pipeline.sh

# 3. Run pipeline
./rnaseq_pipeline.sh 2>&1 | tee pipeline.log
```

### Running the R Pipeline

```r
# 1. Ensure required input files are present:
#    - samples.txt
#    - transdecoder_gene_trans_map.txt
#    - salmon/ directory with quant.sf files
#    - go_annotations.txt

# 2. Run analysis
source("ecklonia_deseq2_pipeline.R")
```

---

## Key Results

**Dataset characteristics:**
- **Total samples:** 30 (3 treatments × 2 timepoints × 5 replicates)
- **Raw reads:** 1.215 billion paired-end reads
- **Post-QC reads:** 1.135 billion (93.4% retention)
- **Genes after filtering:** ~26,000
- **Annotated genes with GO terms:** 10,421
- **BUSCO completeness:** 85% (Stramenopiles database)

**Differential expression:**
- **Significantly DE genes (padj < 0.05):**
  - Single heatwave vs control (during): ~2,500 genes
  - Double heatwave vs control (during): ~3,200 genes
  - Recovery effects and priming responses detected

**Principal findings:**
- Clear transcriptomic separation of treatments (PCA: PC1 explains 45% variance)
- Upregulation: Heat shock proteins, molecular chaperones, oxidative stress response
- Downregulation: Photosynthesis, ribosomal proteins, growth-related pathways
- Evidence of transcriptomic priming following repeated heatwave exposure
- Distinct recovery trajectories between single and repeated stress

---

## Reproducibility

**Version control:** All code tracked with Git  
**Documentation:** Detailed comments throughout both pipelines  
**Quality checks:** Validation at each processing stage  
**Parameters:** All thresholds explicitly documented  
**Modular design:** Steps can be run independently  
**Error handling:** Pipeline exits on failure with informative messages  
**Statistical rigor:** Appropriate model selection, filtering, and multiple testing correction

---

## Technical Challenges & Solutions

**Challenge: Non-model organism**
- No reference genome available
- Limited annotation databases for brown macroalgae
- Evolutionary distance from well-studied organisms

**Solutions:**
- De novo transcriptome assembly with Trinity
- Multi-database annotation strategy (SwissProt, Pfam, EggNOG)
- Cross-validation of annotations
- Taxonomic filtering to focus on relevant eukaryotic pathways
- Conservative filtering and validation thresholds
- BUSCO assessment against appropriate lineage (Stramenopiles)

**Challenge: Annotation completeness**
- Transcript-level vs gene-level quantification
- Isoform redundancy from de novo assembly

**Solutions:**
- CD-HIT clustering to reduce redundancy
- TransDecoder for identifying high-confidence coding sequences
- Gene-level aggregation via transcript-to-gene mapping
- Homology and domain evidence to improve ORF predictions

---

## Publications

**Published:**
Wynn, O.J. et al. (2024). Short-term resilience, long-term costs: Reduced growth and increased erosion in the kelp *Ecklonia radiata* following repeated marine heatwaves. *Journal of Phycology*. https://doi.org/10.1111/jpy.70076

**In Review:**
- Wynn, O.J. et al. The effect of ocean warming and CO2 enrichment on *Ecklonia radiata*: molecular responses from a multi-factor experiment. *PNAS*.
- Wynn, O.J. et al. Acclimation at a cost: Transcriptomic insights into *Ecklonia radiata*'s response and recovery from marine heatwaves. *Proceedings B*.

---

## Contact

Olivia Wynn  
PhD Candidate (submitted), Institute for Marine and Antarctic Studies  
Marine Consultant, Marine Solutions Tasmania  
Email: oj.wynn9@gmail.com  
GitHub: ojwynn9

---

## License

This code is provided for educational and research purposes. Please cite the associated publications if you use or adapt this workflow.

---

Last updated: December 2025
