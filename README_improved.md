# Ecklonia radiata RNA-seq Analysis Pipeline

**Author:** Olivia Wynn  
**Institution:** Institute for Marine and Antarctic Studies (IMAS), University of Tasmania  
**Project:** Transcriptomic responses of *Ecklonia radiata* to marine heatwaves

---

## 📋 Overview

This repository contains a complete R-based pipeline for RNA-seq differential expression analysis of the kelp *Ecklonia radiata* under marine heatwave conditions. The workflow processes raw sequencing data through quality control, normalization, statistical testing, and biological interpretation.

**Key capabilities demonstrated:**
- ✅ Data pipeline development (ETL workflow)
- ✅ Quality control and validation
- ✅ Statistical analysis and reproducible workflows
- ✅ R programming and Bioconductor
- ✅ Version control and documentation
- ✅ Working with non-model organisms

---

## 🔬 Biological Context

*Ecklonia radiata* is a dominant kelp species in temperate Australian waters facing increasing marine heatwave frequency. This project investigated transcriptomic responses to:
- **Single heatwave exposure** (acute stress)
- **Repeated heatwave exposure** (cumulative stress)
- **Recovery period** (post-stress acclimation)

**Challenge:** This species lacks a reference genome, requiring de novo transcriptome assembly and custom annotation pipelines.

---

## 🔧 Pipeline Workflow

### **1. Data Import**
- Import Salmon quantification files (pseudo-alignment output)
- Load sample metadata and experimental design
- Import transcript-to-gene mapping (from TransDecoder)

### **2. Quality Control**
- Assess count distributions (negative binomial check)
- Mean-variance relationship analysis
- Filter low-count genes (noise removal)

### **3. Normalization & Transformation**
- DESeq2 median-of-ratios normalization
- Variance stabilizing transformation (VST) for visualization
- PCA for sample clustering and quality assessment

### **4. Differential Expression Analysis**
- Statistical testing with DESeq2 (negative binomial GLM)
- Multiple testing correction (Benjamini-Hochberg)
- Log fold-change shrinkage for improved estimates
- Identification of significantly DE genes (padj < 0.05)

### **5. Visualization**
- PCA biplots with treatment clustering
- Volcano plots for DE results
- Heatmaps of gene expression patterns

### **6. GO Enrichment** *(separate script)*
- Gene Ontology over-representation analysis
- Pathway enrichment (biological interpretation)

---

## 📁 Repository Structure

```
.
├── ecklonia_deseq2_pipeline.R    # Main analysis script
├── README.md                       # This file
├── samples.txt                     # Sample metadata
├── transdecoder_gene_trans_map.txt # Transcript-gene mapping
└── go_annotations.txt              # GO term annotations
```

---

## 💻 Technical Requirements

### **R Version**
- R ≥ 4.2

### **Core Packages**
```r
# Data handling
tidyverse, dplyr, tibble, tidyr, tximport

# Differential expression
DESeq2, edgeR, apeglm, ashr

# Quality control & visualization
PCAtools, ggplot2, pheatmap, EnhancedVolcano

# GO enrichment
goseq, clusterProfiler, GO.db
```

Install Bioconductor packages:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "PCAtools", "goseq"))
```

---

## 🚀 Quick Start

### **1. Prepare input files**
Ensure these files are in your working directory:
- `samples.txt` - Sample metadata
- `transdecoder_gene_trans_map.txt` - Gene mapping
- `salmon/` - Directory containing Salmon quantification output

### **2. Run the pipeline**
```r
source("ecklonia_deseq2_pipeline.R")
```

### **3. Key outputs**
- Normalized expression matrix
- PCA plots for quality control
- Differential expression results
- Volcano plots
- Lists of significantly DE genes

---

## 📊 Example Results

**Dataset processed:**
- **30 samples** (3 treatments × 2 timepoints × 5 replicates)
- **~26,000 genes** after filtering
- **~2,500 differentially expressed genes** (single heatwave vs control)

**Key findings:**
- Clear separation of treatments in PCA (PC1: 45% variance)
- Upregulation of heat shock proteins and stress response genes
- Downregulation of photosynthesis and growth-related genes
- Evidence of transcriptomic acclimation in repeated heatwave treatment

---

## 🔄 Reproducibility

**This pipeline demonstrates best practices for reproducible research:**

✅ **Version control** - All code tracked with Git  
✅ **Documented workflow** - Clear comments explaining each step  
✅ **Quality checks** - Systematic validation of data at each stage  
✅ **Statistical rigor** - Appropriate model selection and testing  
✅ **Transparent methods** - All parameters and thresholds documented  

---

## 🎯 Skills Demonstrated (Relevant to eDNA Data Specialist Role)

### **Data Pipeline Development**
- ETL workflow: Extract (Salmon) → Transform (normalization/filtering) → Load (DESeq2)
- Quality control at multiple stages
- Handling large datasets (~30 samples, 60 paired-end files)

### **R Programming**
- Data manipulation (tidyverse)
- Statistical analysis (DESeq2, edgeR)
- Visualization (ggplot2, PCAtools)
- Scripting and automation

### **Bioinformatics**
- Sequence data processing
- Working with quantification files
- Annotation integration (GO terms)
- Non-model organism challenges

### **Reproducible Workflows**
- Clear documentation
- Modular code structure
- Version control (Git)
- Transparent methodology

### **Statistical Understanding**
- Count data distributions
- Normalization methods
- Multiple testing correction
- Appropriate model selection

---

## 📚 Related Publications

**Published:**
- Wynn, O.J. et al. (2024). Short-term resilience, long-term costs: Reduced growth and increased erosion in the kelp *Ecklonia radiata* following repeated marine heatwaves. *Journal of Phycology*. https://doi.org/10.1111/jpy.70076

**In Review:**
- Wynn, O.J. et al. The effect of ocean warming and CO₂ enrichment on *Ecklonia radiata*: molecular responses from a multi-factor experiment. *PNAS*.
- Wynn, O.J. et al. Acclimation at a cost: Transcriptomic insights into *Ecklonia radiata*'s response and recovery from marine heatwaves. *Proceedings B*.

---

## 📧 Contact

**Olivia Wynn**  
PhD Candidate (submitted), Institute for Marine and Antarctic Studies  
Marine Consultant, Marine Solutions Tasmania  
Email: oj.wynn9@gmail.com  
GitHub: [your-github-username]

---

## 🔗 Transferable Skills to eDNA Workflows

While this pipeline processes RNA-seq data, the core skills directly transfer to eDNA data processing:

| RNA-seq Skill | eDNA Equivalent |
|---------------|-----------------|
| Salmon quantification import | DADA2/QIIME2 ASV import |
| Quality filtering low-count genes | Quality filtering low-abundance taxa |
| Transcript-to-gene mapping | Sequence-to-taxonomy assignment |
| DESeq2 normalization | Rarefaction/compositional normalization |
| Metadata management | Darwin Core metadata standards |
| Reproducible R pipelines | eDNA data publication pipelines |

**Key principle:** Both workflows involve ETL processes, quality control, normalization, and transformation of sequencing data into interpretable biological datasets with proper metadata.

---

## 📝 License

This code is provided for educational and research purposes. Please cite the associated publications if you use or adapt this workflow.

---

*Last updated: December 2024*
