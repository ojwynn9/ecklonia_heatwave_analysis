#Packages used for analysis
#General packages
library(ggplot2)
library(stringr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(reshape2)
library(vegan)
library(readxl)
library(writexl)
library(BiocManager)
library(MASS)
library(tximport)
library(devtools)
library(trinotateR)
#DGE analysis packages
library(DESeq2)
library(edgeR)
library(glmpca)
library(apeglm)
library(ashr)
library(VennDiagram)
#ORA packages
library(goseq)
#Plotting packages
library(EnhancedVolcano)
library(gridExtra)
library(PCAtools)
library(DEGreport)
library(RColorBrewer)
library(ggpubr)
library(ggridges)
library(gtable)
library(plotly)
library(VennDiagram)
library(pheatmap)
#Other packages
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)
library(magick)
library(rstatix)
library(hexbin)
library(methods)
library(rmarkdown)
library(org.Mm.eg.db)
library(apeglm)
library(GSEABase)
library(fgsea)
library(DEXSeq)
library(msigdbr)
library(vsn)
library(clusterProfiler)
library(NbClust)
library(pheatmap)
library(enrichplot)
library(ggnewscale)
library(patchwork)
library(knitr)
library(testthat)
library(yaml)
library(biomaRt)
library(GO.db)
library(AnnotationDbi)
library(data.table)
#Import gene count matrices genereated from Salmon and associated metadata
#Import metadata
sampleinfo.df <- read_tsv("samples.txt")
sampleinfo.df
headers <- colnames(sampleinfo.df)
sampleinfo.df <- rbind(headers, sampleinfo.df)

colnames(sampleinfo.df) <- paste0("V", seq_len(ncol(sampleinfo.df))) # Reset column names to generic

colnames(sampleinfo.df)[1:2] <- c("group", "sample") # Rename columns

sampleinfo.df <- sampleinfo.df[, 1:(ncol(sampleinfo.df) - 2)] # Remove the last two columns

print(sampleinfo.df) # View the modified data frame

#Import the gene count matrix resulting from mapping of sequences via pseudo-alignment using Salmon
files <- file.path("C:/Users/oowynn/OneDrive - University of Tasmania/Documents/heatwave_transcriptomics/ecklonia_heatwave_transcriptomics/Ecklonia_heatwave_transcriptomics/Heat_transcriptomics_Ch4/salmon", sampleinfo.df$sample, "quant.sf")
names(files) <- sampleinfo.df$sample

# Assuming the file is tab-delimited, with no header
tx2gene <- read.table("transdecoder_gene_trans_map.txt", header = FALSE, col.names = c("transcript_id", "gene_id"))
tx2gene <- tx2gene %>%
  mutate(transcript_id = gene_id,  # Copy gene_id to transcript_id
         gene_id = gsub("_i[0-9]+", "", gene_id))  # Remove _i* from gene_id
head(tx2gene)

# Now use tx2gene in tximport
heatwave_gene_counts.df <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="no")
head(heatwave_gene_counts.df)
str(heatwave_gene_counts.df)
#convert it to a df
heatwave_gene_counts.df <- as.data.frame(heatwave_gene_counts.df)
#Check data types
str(heatwave_gene_counts.df)


# Decide on appropriate statistical model for downstream analysis 
#Visualise RNA-seq count distribution for one sample
ggplot(heatwave_gene_counts.df) +
  geom_histogram(aes(x = `abundance.C1.2`), stat = "bin", bins = 200) +
  xlim(-5, 500) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#Calculate mean and variance of count data
mean_counts <- apply(heatwave_gene_counts.df[, 2:31], 1, mean)
variance_counts <- apply(heatwave_gene_counts.df[, 2:31], 1, var)
check.df <- data.frame(mean_counts, variance_counts)
head(check.df)

#Plot mean versus variance
ggplot(check.df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()
#Overall, the variance appears to be greater than mean of counts
#This is typical with what you would expect for RNASeq data that follows a negative binomial
#This rules out use of Poisson distribution where mean=variance
#Rather, a negative binomial model is more appropriate to this analysis (in this case DESeq2 will be used)


# Data preparation
## Preparation of gene count matrices for downstream analysis
#Modify the heatwave_metadata.df table, adding  sample names as row names (needed for some of the DESeq2 functions)
sampleinfo.df <- as.data.frame(sampleinfo.df)
rownames(sampleinfo.df) <- sampleinfo.df$SampleName

# Set row names to the "sample" column
rownames(sampleinfo.df) <- sampleinfo.df$sample
# Add co2 and temp column info to sampleinfo.df
sampleinfo.df$timepoint <- c(
  'during', 'during', 'during', 'during', 'during',
  'post', 'post', 'post', 'post', 'post',
  'during', 'during', 'during', 'during', 'during',
  'post', 'post', 'post', 'post', 'post',
  'during', 'during', 'during', 'during', 'during',
  'post', 'post', 'post', 'post', 'post'
)
sampleinfo.df$heatwave <- c(
  'control', 'control', 'control', 'control', 'control',
  'control', 'control', 'control', 'control', 'control',
  'single', 'single', 'single', 'single', 'single',
  'single', 'single', 'single', 'single', 'single',
  'double', 'double', 'double', 'double', 'double',
  'double', 'double', 'double', 'double', 'double'
)
# Convert temp to a factor since it's representing treatment levels
sampleinfo.df$timepoint <- factor(sampleinfo.df$timepoint, levels = c("during", "post"))
sampleinfo.df$heatwave <- factor(sampleinfo.df$heatwave, levels = c("control", "single", "double"))
head(sampleinfo.df)


## Split heatwave_gene_counts.df into an abundance amtrix, count (raw) matrix and length matrix so column headers amtch rownames in sampleinfo.df
# Separate count matrix
counts_matrix <- heatwave_gene_counts.df[, grep("counts", colnames(heatwave_gene_counts.df))]
# Adjust column names of counts matrix to match metadata
colnames(counts_matrix) <- gsub("counts\\.", "", colnames(counts_matrix))
# Check if column names match with metadata sample names
all(colnames(counts_matrix) == rownames(metadata))  # This should return TRUE if matched
head(counts_matrix)
dim(counts_matrix)  # Should show the number of genes (rows) and the number of samples (columns)
dim(sampleinfo.df)  # Should show the number of samples (rows) and number of metadata variables (columns)
str(counts_matrix)
# remove non-numeric column
counts_matrix <- counts_matrix %>%
  dplyr::select(-`countsFromAbundance`)
counts_matrix <- round(as.matrix(counts_matrix))  # Round to nearest integer
counts_matrix <- as.data.frame(counts_matrix)
str(counts_matrix)
colnames(counts_matrix) <- gsub("\\.", "-", colnames(counts_matrix))  # Convert C1.2 → C1-2


## abundance
abundance_matrix <- heatwave_gene_counts.df[, grep("abundance", colnames(heatwave_gene_counts.df))]
colnames(abundance_matrix) <- gsub("abundance\\.", "", colnames(abundance_matrix))
all(colnames(abundance_matrix) == rownames(sampleinfo.df))  # This should return TRUE if matched
abundance_matrix <- as.matrix(abundance_matrix)
str(abundance_matrix)
colnames(abundance_matrix) <- gsub("\\.", "-", colnames(counts_matrix))  # Convert C1.2 → C1-2
## length
length_matrix <- heatwave_gene_counts.df[, grep("length", colnames(heatwave_gene_counts.df))]
colnames(length_matrix) <- gsub("length\\.", "", colnames(length_matrix))
all(colnames(length_matrix) == rownames(sampleinfo.df))  # This should return TRUE if matched
head(length_matrix)
colnames(length_matrix) <- gsub("\\.", "-", colnames(counts_matrix))  # Convert C1.2 → C1-2

# Convert temp and co2 to factors in sampleinfo.df
sampleinfo.df$timepoint <- as.factor(sampleinfo.df$timepoint)
sampleinfo.df$heatwave <- as.factor(sampleinfo.df$heatwave)

## Differential expression with DESeq2
# Now, create the DESeqDataSet
heatwave.dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                       colData = sampleinfo.df,
                                       design = ~ timepoint + heatwave + timepoint:heatwave)
#Perform the median of ratios method of normalization (from DESeq2) then filter genes
heatwave.dds <- estimateSizeFactors(heatwave.dds)
#heatwave.dds_filt <- rowSums(counts(heatwave.dds, normalized = TRUE) >= 5) >= 2
#table(heatwave.dds_filt)
#heatwave.dds_filt2 <- heatwave.dds[heatwave.dds_filt,]
#Number of genes (27,280) for heatwave analysis after filtering
#nrow(heatwave.dds_filt2)

## removing fitlering as low genes
heatwave.dds_filt2 <- heatwave.dds
nrow(heatwave.dds_filt2) # 35407

# Check levels
levels(sampleinfo.df$timepoint)  # Should be 'during' and 'post'
levels(sampleinfo.df$heatwave)   # Should be 'control', 'single', 'double'

# VST transformation
heatwave.vsd <- vst(heatwave.dds_filt2, blind = FALSE)
heatwave.vsd.mat <- assay(heatwave.vsd)

# Perform PCA using PCAtools
heatwave.PCA_tools <- PCAtools::pca(heatwave.vsd.mat, metadata = sampleinfo.df, removeVar = 0.1)

explained_variance <- heatwave.PCA_tools$variance
print(explained_variance)
# After running PCAtools::pca(...)
var_pc <- heatwave.PCA_tools$variance  # % variance per PC

# Get the percentages for PC1 and PC2 (works whether the vector is named or not)
pc1 <- if (!is.null(names(var_pc))) unname(var_pc["PC1"]) else var_pc[1]
pc2 <- if (!is.null(names(var_pc))) unname(var_pc["PC2"]) else var_pc[2]

cat(sprintf("PC1 = %.2f%%, PC2 = %.2f%%\n", pc1, pc2))


# Screeplot to identify the optimal number of PCs
heatwave_horn <- parallelPCA(heatwave.vsd.mat)
heatwave_elbow <- findElbowPoint(heatwave.PCA_tools$variance)

heatwave.PCA_tools.screeplot <- screeplot(heatwave.PCA_tools,
                                          components = getComponents(heatwave.PCA_tools, 1:20),
                                          vline = c(heatwave_horn$n, heatwave_elbow)) +
  geom_label(aes(x = heatwave_horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 12)) +
  geom_label(aes(x = heatwave_elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 12))
heatwave.PCA_tools.screeplot

# Biplot with ellipses by heatwave treatment
heatwave.heatwave.PCA_tools.biplot <- biplot(heatwave.PCA_tools,
                                             x = "PC1",
                                             y = "PC2",
                                             lab = NA,
                                             colby = 'heatwave',
                                             shape = 'timepoint',  # Use timepoint as shape
                                             pointSize = 3,
                                             hline = 0, hlineWidth = 0.25,
                                             vline = 0, vlineWidth = 0.25,
                                             xlim = NULL,  # Adjust if needed
                                             ylim = NULL) + # Adjust if needed
  stat_ellipse(aes(group = sampleinfo.df$heatwave, fill = sampleinfo.df$heatwave),
               geom = "polygon", alpha = 0.25, level = 0.95, type = 't') +
  scale_fill_manual(
    values = c("#56B4E9", "#009E73", "#F0E442"),
    name = "heatwave Treatment"
  ) +
  scale_colour_manual(
    values = c("#56B4E9", "#009E73", "#F0E442"),
    name = "heatwave Treatment"
  ) +
  scale_shape_manual(
    values = c(16, 17),  # Different shapes for 'during' and 'post'
    name = "Timepoint",
    labels = c("During", "Post")
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.position = 'top',
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black", face = 'bold'),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = 'bold'),
        axis.text.y = element_text(size = 12, colour = "black")) +
  guides(
    fill = "none", 
    linetype = "none",
    colour = guide_legend(title.position = "top", nrow = 1, override.aes = list(size = 3), order = 1),
    shape = guide_legend(title.position = "top", override.aes = list(size = 3), order = 2)
  ) +
  labs(
    shape = 'Timepoint',
    x = "Principal Component 1 (22.77%)",
    y = "Principal Component 2 (13.96%)"
  )

heatwave.heatwave.PCA_tools.biplot
ggsave(
  filename = "heatwave_PCA_biplot.png", 
  plot = heatwave.heatwave.PCA_tools.biplot, 
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 300
)

# Ellipses for each heatwave + timepoint combination
heatwave.heatwave.PCA_tools.biplot <- biplot(heatwave.PCA_tools,
                                             x = "PC1",
                                             y = "PC2",
                                             lab = NA,
                                             colby = 'heatwave',
                                             shape = 'timepoint',  # Use timepoint as shape
                                             pointSize = 3,
                                             hline = 0, hlineWidth = 0.25,
                                             vline = 0, vlineWidth = 0.25,
                                             xlim = NULL,  # Adjust if needed
                                             ylim = NULL) + # Adjust if needed
  # Draw ellipses for each heatwave + timepoint combination
  stat_ellipse(aes(group = interaction(sampleinfo.df$heatwave, sampleinfo.df$timepoint),
                   fill = interaction(sampleinfo.df$heatwave, sampleinfo.df$timepoint)),
               geom = "polygon", alpha = 0.3, level = 0.95, type = 't') +
  scale_fill_manual(
    values = c(
      "#E69F00", "#D55E00",  # Control-During, Control-Post (orange, dark orange)
      "#009E73", "#0072B2",  # Single-During, Single-Post (green, blue)
      "#F0E442", "#CC79A7"   # Double-During, Double-Post (yellow, pink)
    ),
    name = "heatwave x Timepoint"
  ) +
  scale_colour_manual(
    values = c(
      "#E69F00", "#D55E00",  # Control-During, Control-Post
      "#009E73", "#0072B2",  # Single-During, Single-Post
      "#F0E442", "#CC79A7"   # Double-During, Double-Post
    ),
    name = "heatwave x Timepoint"
  ) +
  scale_shape_manual(
    values = c(16, 17),  # Different shapes for 'during' and 'post'
    name = "Timepoint",
    labels = c("During", "Post")
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.position = 'top',
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black", face = 'bold'),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = 'bold'),
        axis.text.y = element_text(size = 12, colour = "black")) +
  guides(
    fill = guide_legend(title.position = "top", override.aes = list(alpha = 0.5), order = 1),
    shape = guide_legend(title.position = "top", override.aes = list(size = 3), order = 2)
  ) +
  labs(
    shape = 'Timepoint',
    x = "Principal Component 1",
    y = "Principal Component 2"
  )

heatwave.heatwave.PCA_tools.biplot


# Eigencorplot using heatwave and timepoint as metavariables
heatwave.PCA_tools.eigenplot <- eigencorplot(heatwave.PCA_tools,
                                             metavars = c('heatwave', 'timepoint'))
heatwave.PCA_tools.eigenplot


## PERMANOVA
# on both factors heatwave and timepoint
set.seed(1234)
heatwave.adonis <- adonis2(t(heatwave.vsd.mat) ~ heatwave * timepoint,
                           permutations = 10000,
                           data = sampleinfo.df,
                           method = "bray")
print(heatwave.adonis)

# heatwave only
set.seed(1234)
heatwave.adonis_heatwave <- adonis2(t(heatwave.vsd.mat) ~ heatwave,
                                    permutations = 10000,
                                    data = sampleinfo.df,
                                    method = "bray")
print(heatwave.adonis_heatwave)

# Timepoint only
set.seed(1234)
heatwave.adonis_timepoint <- adonis2(t(heatwave.vsd.mat) ~ timepoint,
                                     permutations = 10000,
                                     data = sampleinfo.df,
                                     method = "bray")
print(heatwave.adonis_timepoint)

# NMDS 
# Calculate Bray-Curtis dissimilarity matrix
bray_dist <- vegdist(t(heatwave.vsd.mat), method = "bray")

# Run NMDS (Non-metric Multidimensional Scaling)
set.seed(1234)
nmds_result <- metaMDS(bray_dist, k = 2)  # k = 2 for 2D visualization

# Extract NMDS coordinates
nmds_coordinates <- as.data.frame(scores(nmds_result, display = "sites"))
library(viridis)
# Add metadata
nmds_df <- cbind(nmds_coordinates, 
                 heatwave = sampleinfo.df$heatwave, 
                 Timepoint = sampleinfo.df$timepoint,
                 TreatmentGroup = interaction(sampleinfo.df$heatwave, sampleinfo.df$timepoint))
# NMDS Plot for heatwave and timepoint
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = TreatmentGroup, fill = TreatmentGroup, shape = Timepoint)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", alpha = 0.3, level = 0.95, aes(group = TreatmentGroup)) +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Color-blind friendly colors for points
  scale_fill_viridis(discrete = TRUE, option = "D") +  # Color for ellipses
  scale_shape_manual(values = c(16, 17),  # Different shapes for timepoints
                     name = "Timepoint",
                     labels = c("During", "Post")) +
  theme_minimal() +
  ggtitle("NMDS of heatwave and Timepoint Treatments") +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black", face = 'bold'),
    axis.text.x = element_text(size = 12, colour = "black"), 
    axis.title.y = element_text(size = 14, colour = "black", face = 'bold'), 
    axis.text.y = element_text(size = 12, colour = "black")
  )

# Extract PCA-ready matrix from VST/rlog-transformed data
### Exploration of transformation approaches

#Normalized counts are transformed either using vst/rlog
#Setting blind to TRUE results in a transformation unbiased to sample condition information (i.e., unsupervised transformation)
heatwave.vsd2 <- vst(heatwave.dds_filt2, blind = T) 
heatwave.rld <- rlog(heatwave.dds_filt2, blind = T)
heatwave.dds_filt2 <- estimateSizeFactors(heatwave.dds_filt2)
heatwave_trfmd.df <- bind_rows(as_tibble(assay(heatwave.vsd2)[, 1:2]) %>% mutate(transformation = "vst"),
                               as_tibble(assay(heatwave.rld)[, 1:2]) %>% mutate(transformation = "rlog"))

#Scatter plots of transformed counts from two samples
colnames(heatwave_trfmd.df)[1:2] <- c("x", "y")
heatwave_lvls <- c("vst", "rlog")
heatwave_trfmd.df$transformation <- factor(heatwave_trfmd.df$transformation, levels=heatwave_lvls)
ggplot(heatwave_trfmd.df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) 

## DGE analysis and QC
# 1. Ensure that GROUP is a factor
sampleinfo.df$group <- factor(sampleinfo.df$group, levels = c("C1", "C2", "S1", "S2", "D1", "D2"))


# Create DESeqDataSet object with interaction term
heatwave.dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                       colData = sampleinfo.df,
                                       design = ~ group)

# Set reference level
heatwave.dds$group <- relevel(heatwave.dds$group, ref = "C1")
# 1. Filtering raw counts before DESeq2 analysis
#keep <- rowSums(counts(heatwave.dds) >= 5) >= 2 # removed filtering for now!!
#heatwave.dds <- heatwave.dds[keep, ]

# 2. Run DESeq2 pipeline
heatwave.dge <- DESeq(heatwave.dds)

# 3. Normalized counts for PCA/heatmaps (AFTER DESeq)
norm_counts <- counts(heatwave.dge, normalized = TRUE)

# 4. Background genes for ORA should be the genes you input into DESeq2
background_genes <- rownames(heatwave.dds)
length(background_genes)  # Check the number of background genes

# Dispersion plots of gene-wise estimates (QC of model fitting)
plotDispEsts(heatwave.dge)

# Data should scatter around the fitted curve, with dispersion decreasing as mean expression levels increase

### 
# Get initial results
heatwave.dge.results <- results(heatwave.dge)
head(heatwave.dge.results)

summary(heatwave.dge.results)
sum(heatwave.dge.results$padj < 0.05, na.rm=TRUE) # p-val 0.05 = 2184 when remove filtering increases to 2202
background_genes <- rownames(heatwave.dds)
length(background_genes)  # Check the number of background genes

# Set contrasts
#During
heatwave_contrast1 <- c("group", "S1", "C1")
heatwave_contrast2 <- c("group", "D1", "C1")
heatwave_contrast3 <- c("group", "D1", "S1")

#Post
heatwave_contrast4 <- c("group", "S2", "C2")
heatwave_contrast5 <- c("group", "D2", "C2")
heatwave_contrast6 <- c("group", "D2", "S2")

#heatwave
heatwave_contrast7 <- c("group", "C1", "C2")
heatwave_contrast8 <- c("group", "S1", "S2")
heatwave_contrast9 <- c("group", "D1", "D2")

# Run pairwise comparisons
res_contrast1 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast1, alpha = 0.05)
res_contrast2 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast2, alpha = 0.05)
res_contrast3 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast3, alpha = 0.05)
res_contrast4 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast4, alpha = 0.05)
res_contrast5 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast5, alpha = 0.05)
res_contrast6 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast6, alpha = 0.05)
res_contrast7 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast7, alpha = 0.05)
res_contrast8 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast8, alpha = 0.05)
res_contrast9 <- DESeq2::results(heatwave.dge, contrast = heatwave_contrast9, alpha = 0.05)

summary(res_contrast1)
summary(res_contrast2)
summary(res_contrast3)
summary(res_contrast4)
summary(res_contrast5)
summary(res_contrast6)
summary(res_contrast7)
summary(res_contrast8)
summary(res_contrast9)

# Apply LFC shrinkage to the results tables (reduces noise from low gene counts)
heatwave_res_Ash1 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast1, res = res_contrast1, type = "ashr")
heatwave_res_Ash2 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast2, res = res_contrast2, type = "ashr")
heatwave_res_Ash3 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast3, res = res_contrast3, type = "ashr")
heatwave_res_Ash4 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast4, res = res_contrast4, type = "ashr")
heatwave_res_Ash5 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast5, res = res_contrast5, type = "ashr")
heatwave_res_Ash6 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast6, res = res_contrast6, type = "ashr")
heatwave_res_Ash7 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast7, res = res_contrast7, type = "ashr")
heatwave_res_Ash8 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast8, res = res_contrast8, type = "ashr")
heatwave_res_Ash9 <- lfcShrink(heatwave.dge, contrast = heatwave_contrast9, res = res_contrast9, type = "ashr")

# Convert results into data frames and assign contrast labels
heatwave_res1 <- heatwave_res_Ash1 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group S1 vs C1')

heatwave_res2 <- heatwave_res_Ash2 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group D1 vs C1')

heatwave_res3 <- heatwave_res_Ash3 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group D1 vs S1')

heatwave_res4 <- heatwave_res_Ash4 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group S2 vs C2')

heatwave_res5 <- heatwave_res_Ash5 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group D2 vs C2')

heatwave_res6 <- heatwave_res_Ash6 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group D2 vs S2')

heatwave_res7 <- heatwave_res_Ash7 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group C1 vs C2')

heatwave_res8 <- heatwave_res_Ash8 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group S1 vs S2')

heatwave_res9 <- heatwave_res_Ash9 %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  mutate(contrast = 'group D1 vs D2')

# Inspect the meaning of each column in the shrinked results
mcols(heatwave_res_Ash1, use.names = TRUE)
mcols(heatwave_res_Ash2, use.names = TRUE) 
mcols(heatwave_res_Ash3, use.names = TRUE)
mcols(heatwave_res_Ash4, use.names = TRUE)  
mcols(heatwave_res_Ash5, use.names = TRUE)  
mcols(heatwave_res_Ash6, use.names = TRUE)  
mcols(heatwave_res_Ash7, use.names = TRUE)
mcols(heatwave_res_Ash8, use.names = TRUE)  
mcols(heatwave_res_Ash9, use.names = TRUE)  

# Set thresholds for filtering significant genes
padj.cutoff <- 0.05
lfc.cutoff <- 0.0  # Corresponds to a fold change of 1.5

# Filter for significant genes (adjusted p-value < 0.05)
# Optionally, you can also apply the log2 fold change cutoff by uncommenting the lfc.cutoff part
heatwave_signif_genes1 <- heatwave_res1 %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes2 <- heatwave_res2 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes3 <- heatwave_res3 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes4 <- heatwave_res4 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes5 <- heatwave_res5 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes6 <- heatwave_res6 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes7 <- heatwave_res7 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes8 <- heatwave_res8 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

heatwave_signif_genes9 <- heatwave_res9 %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

# List of data frames to iterate over
datasets <- list(
  heatwave_signif_genes1,
  heatwave_signif_genes2,
  heatwave_signif_genes3,
  heatwave_signif_genes4,
  heatwave_signif_genes5,
  heatwave_signif_genes6,
  heatwave_signif_genes7,
  heatwave_signif_genes8,
  heatwave_signif_genes9
)

# Count the number of significant genes for each contrast
nrow(heatwave_signif_genes1)
nrow(heatwave_signif_genes2)
nrow(heatwave_signif_genes3)
nrow(heatwave_signif_genes4)
nrow(heatwave_signif_genes5)
nrow(heatwave_signif_genes6)
nrow(heatwave_signif_genes7)
nrow(heatwave_signif_genes8)
nrow(heatwave_signif_genes9)

# Combine all differential expression results into one dataframe
heatwave.DEGs <- rbind(heatwave_res1, heatwave_res2) %>%
  rbind(heatwave_res3) %>%
  rbind(heatwave_res4) %>%
  rbind(heatwave_res5) %>%
  rbind(heatwave_res6) %>%
  rbind(heatwave_res7) %>%
  rbind(heatwave_res8) %>%
  rbind(heatwave_res9)

# Check structure of the combined results
head(heatwave.DEGs)
tail(heatwave.DEGs)
str(heatwave.DEGs)

# Total number of entries in the gene column
total_entries <- nrow(heatwave.DEGs)
print(total_entries)

# Number of unique genes in the gene column
unique_genes <- length(unique(heatwave.DEGs$gene))
print(unique_genes)


### Visualization of DE genes
## Basic volcano plots
# Update the 'diffexpressed' column for significance classification
heatwave_signif_genes1 <- heatwave_signif_genes1 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes2 <- heatwave_signif_genes2 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes3 <- heatwave_signif_genes3 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes4 <- heatwave_signif_genes4 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes5 <- heatwave_signif_genes5 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes6 <- heatwave_signif_genes6 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes7 <- heatwave_signif_genes7 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes8 <- heatwave_signif_genes8 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))
heatwave_signif_genes9 <- heatwave_signif_genes9 %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

## Create volcano plots for each contrast
heatwave_volcano_plots <- list()
contrasts <- c("S1 vs C1", "D1 vs C1", "D1 vs S1", "S2 vs C2", "D2 vs C2", "D2 vs S2", "C1 vs C2", "S1 vs S2", "D1 vs D2")
heatwave_signif_genes_list <- list(heatwave_signif_genes1, heatwave_signif_genes2, heatwave_signif_genes3,
                                   heatwave_signif_genes4, heatwave_signif_genes5, heatwave_signif_genes6,
                                   heatwave_signif_genes7, heatwave_signif_genes8, heatwave_signif_genes9)

for (i in 1:length(contrasts)) {
  heatwave_volcano_plots[[i]] <- EnhancedVolcano(heatwave_signif_genes_list[[i]],
                                                 lab = NA,
                                                 x = 'log2FoldChange',
                                                 y = 'padj',
                                                 xlim = c(-3, 3), ylim = c(0, 10),
                                                 title = '', subtitle = contrasts[i],
                                                 caption = paste0(nrow(heatwave_signif_genes_list[[i]]), " DEGs"),
                                                 legendPosition = "none", legendLabSize = 8, legendIconSize = 2,
                                                 pCutoff = 0.05, FCcutoff = 0,
                                                 cutoffLineWidth = 0.25, cutoffLineCol = "black",
                                                 cutoffLineType = "dashed",
                                                 gridlines.minor = F, gridlines.major = F,
                                                 labSize = 3.0, axisLabSize = 14, colAlpha = 1) + 
    theme(plot.margin=unit(c(-1,-0.25,-0.95,0), "cm")) +
    guides(colour = F)
}

## Arrange and display volcano plots
heatwave_DEG_volc_plots_ALL <- ggarrange(plotlist = heatwave_volcano_plots, ncol = 3, nrow = 3, align = 'hv')

# Add overarching axis labels
y_label <- text_grob(expression(bold(-log[10](italic(P)))), rot = 90, size = 14)
x_label <- text_grob(expression(bold(Log[2]~Fold~Change)), size = 14)

heatwave_DEG_volc_plots_ALL <- annotate_figure(heatwave_DEG_volc_plots_ALL, left = y_label, bottom = x_label)
heatwave_DEG_volc_plots_ALL

Fig.3 <- heatwave_DEG_volc_plots_ALL
# Display the adjusted figure with a specific aspect ratio
ggsave("heatwave_volcano_plots_all.png", heatwave_DEG_volc_plots_ALL, 
       width = 12, height = 8)  # Adjust width and height for thinner plots

Fig.3 <- heatwave_DEG_volc_plots_ALL
ggsave('Fig.3.Volcanoplots.final.png', Fig.3, width = 4, height = 3, units = "in")
ggsave(plot = Fig.3, filename = "Fig.3.pdf",
       path = "./",
       width = 12,
       height = 8,
       units = "in",
       dpi = 300,
       device = "pdf")
ggsave(plot = Fig.3, filename = "Fig.3.png",
       path = "./",
       width = 12,
       height = 8,
       units = "in",
       dpi = 600,  # High resolution for clear rendering
       device = "png",
       bg = "white")

# Heatmap
# Extract normalized counts
norm_counts <- counts(heatwave.dge, normalized = TRUE)

# Check dimensions of normalized counts
dim(norm_counts)  
head(norm_counts[, 1:5])  # View a few columns

# Extract significant genes
signif_genes <- rownames(subset(heatwave.dge.results, padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff))

# Filter normalized counts for significant genes
heatmap_counts <- norm_counts[signif_genes, ]

# Standardize across genes (rows)
heatmap_counts_scaled <- t(scale(t(heatmap_counts)))

# Ensure no NA values are introduced
heatmap_counts_scaled[is.na(heatmap_counts_scaled)] <- 0  

# Check dimensions
dim(heatmap_counts_scaled)  

# Load required packages
library(pheatmap)
library(RColorBrewer)

# Define color scale
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# Define **colorblind-friendly** and **distinct** group annotation colors
group_colors <- c(
  "C1" = "#EF4567",   # Orange
  "C2" = "#56B4E9",   # Sky Blue
  "S1" = "#009E73",   # Green
  "S2" = "#F0E442",   # Yellow
  "D1" = "#D55E00",   # Red-Orange
  "D2" = "#007000"    # Pink-Purple
)


# Sample annotation based on experimental groups
annotation_col <- data.frame(Group = sampleinfo.df$group)
rownames(annotation_col) <- colnames(heatmap_counts_scaled)

# Generate the heatmap
pheatmap(
  heatmap_counts_scaled, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_col = annotation_col, 
  show_rownames = FALSE,  
  show_colnames = TRUE,  
  color = heatmap_colors,
  fontsize = 10,
  main = "Heatmap of Significant DEGs (Heatwave Experiment)"
)

# Define **colorblind-friendly** and **distinct** group annotation colors
group_colors <- c(
  "C1" = "#E69F00",   # Orange
  "C2" = "#56B4E9",   # Sky Blue
  "S1" = "#009E73",   # Green
  "S2" = "#F0E442",   # Yellow
  "D1" = "#D55E00",   # Red-Orange
  "D2" = "#CC79A7"    # Pink-Purple
)

# Sample annotation based on experimental groups
annotation_col <- data.frame(Group = sampleinfo.df$group)
rownames(annotation_col) <- colnames(heatmap_counts_scaled)

# Assign colors to the annotation
ann_colors <- list(Group = group_colors)

# **Force groups to be together** by manually ordering columns
ordered_cols <- order(annotation_col$Group)  # Sort columns by group
heatmap_counts_ordered <- heatmap_counts_scaled[, ordered_cols]
annotation_col <- annotation_col[ordered_cols, , drop=FALSE]

# Generate the heatmap with **no clustering of columns**
pheatmap(
  heatmap_counts_ordered, 
  cluster_rows = TRUE, 
  cluster_cols = FALSE,  # Disable column clustering to keep groups together
  annotation_col = annotation_col, 
  annotation_colors = ann_colors, # Apply distinct colors to groups
  show_rownames = FALSE,  
  show_colnames = TRUE,  
  color = heatmap_colors,
  fontsize = 10,
  main = "Heatmap of Significant DEGs (Grouped by Condition)"
)

# Define file path and save the heatmap as a high-quality PNG
png(filename = "Heatmap_DEGs_Grouped.png", width = 2000, height = 2000, res = 300)

# Generate the heatmap with no clustering of columns
pheatmap(
  heatmap_counts_ordered, 
  cluster_rows = TRUE, 
  cluster_cols = FALSE,  # Disable column clustering to keep groups together
  annotation_col = annotation_col, 
  annotation_colors = ann_colors, # Apply distinct colors to groups
  show_rownames = FALSE,  
  show_colnames = TRUE,  
  color = heatmap_colors,
  fontsize = 10,
  main = "Heatmap of Significant DEGs (Grouped by Condition)"
)

# Close the PNG device to finalize the image file
dev.off()



# Annotation
# Import go annotations
library(readr)

# Import the GO annotations file (tab-separated)
go_annotations <- read_tsv("go_annotations_new.txt", col_names = c("gene_id", "go_terms"))

# View the first few rows
head(go_annotations)

# Read in Trinotate Report
trinotate_report <- read_trinotate("trinotate_annotation_report_new.xls")
head(trinotate_report)
str(trinotate_report)
# Check column names
colnames(trinotate_report)

# Get a summary of the data
summary_trinotate(trinotate_report)

## Splitting GO terms
# convert report to data table
trinotate_report_dt <- setDT(trinotate_report)

# defining split_GO() function
split_GO <- function(x, hit = "gene_ontology_blast") {
  # Filter out rows where the GO column is NA
  y <- x[!is.na(x[[hit]]), .(gene_ontology = x[[hit]], gene_id = x$gene_id, 
                             transcript_id = x$transcript_id, prot_id = x$prot_id)]
  
  # Split multiple annotations in backtick-delimited list
  z <- strsplit(y$gene_ontology, "`")
  n <- sapply(z, length)
  
  # Split each GO annotation into 3 components by caret (^) delimiter
  z <- strsplit(unlist(z), "\\^")
  
  # Create a data frame with separate columns for each component
  x1 <- data.frame(
    gene = rep(y$gene_id, n),
    transcript = rep(y$transcript_id, n),
    protein = rep(gsub(".*\\|", "", y$prot_id), n),
    go = sapply(z, "[", 1),
    ontology = sapply(z, "[", 2),
    name = sapply(z, "[", 3),
    stringsAsFactors = FALSE
  )
  
  # Display a message with the number of annotations processed
  message(nrow(x1), " ", hit, " annotations")
  
  # Return the final data.table
  data.table(x1)
}

###Pfam go term extraction
split_go_data_pfam <- split_GO(trinotate_report, hit = "gene_ontology_Pfam") %>%
  dplyr::select(-transcript,-protein) %>%
  distinct() %>% dplyr::select(gene, go) #77286 annotations
head(split_go_data_pfam)

# Extract and split EggNOG GO terms
split_GO_eggnog <- function(x, hit = "EggNM.GOs") {
  # Filter out rows where the EggNOG GO column is NA
  y <- x[!is.na(x[[hit]]), .(gene_ontology = x[[hit]], gene_id = x$gene_id)]
  
  # Split multiple GO terms by commas
  z <- strsplit(y$gene_ontology, ",")
  n <- sapply(z, length)
  
  # Create a data frame with each GO term as a separate row
  x1 <- data.frame(
    gene = rep(y$gene_id, n),
    go = unlist(z),
    stringsAsFactors = FALSE
  )
  
  # Display a message with the number of annotations processed
  message(nrow(x1), " ", hit, " annotations")
  
  # Return the final data.table
  data.table(x1)
}

# Apply function to extract EggNOG GO terms
split_go_data_eggnog <- split_GO_eggnog(trinotate_report, hit = "EggNM.GOs")

# Combine Pfam and EggNOG GO terms
split_go_data <- bind_rows(split_go_data_pfam, split_go_data_eggnog) %>%
  distinct()

## new
# split the go annotations object
library(data.table)
# Ensure go_annotations is a data.table
go_annotations <- as.data.table(go_annotations)
# Function to split GO terms from go_annotations
split_GO_annotations <- function(x) {
  # Filter out rows where GO terms are NA
  y <- x[!is.na(x$go_terms), .(gene_id, go_terms)]
  
  # Split multiple GO terms by commas
  z <- strsplit(y$go_terms, ",")
  n <- sapply(z, length)
  
  # Create a data frame with each GO term as a separate row
  x1 <- data.frame(
    gene = rep(y$gene_id, n),  # Repeat gene IDs based on number of GO terms
    go = unlist(z),            # Flatten the list of GO terms
    stringsAsFactors = FALSE
  )
  
  # Display a message with the number of annotations processed
  message(nrow(x1), " GO annotations processed.")
  
  # Return the final data.table
  data.table(x1)
}
# Apply the function to your go_annotations dataset
split_go_data <- split_GO_annotations(go_annotations)
# View the result
head(split_go_data)

head(heatwave.DEGs) # contains all DE results from each contrast


# Create GO.list from the combined data
heatwave_GO.list <- split(split_go_data$gene, split_go_data$go)
GO.list <- split(split_go_data$go, split_go_data$gene)
heatwave_GO.list <- split(split_go_data$go, split_go_data$gene)
GO.list <- GO.list[lengths(GO.list) > 0] # Remove genes without GO terms
# Filter for rows with non-NA GO terms
genes_with_go <- split_go_data %>% filter(!is.na(go))

# Count unique genes with GO terms
num_genes_with_go <- genes_with_go %>% distinct(gene) %>% nrow()

# Print the result
cat("Number of genes with associated GO terms:", num_genes_with_go, "\n")

str(GO.list)
head(GO.list)
head(split_go_data)

# Filter background genes to include only those with GO terms
background_genes_filtered <- background_genes[background_genes %in% names(GO.list)]

# Check the updated count of background genes
length(background_genes_filtered)


#GO.list ready to go to carry go enrichment section. 
genes_with_go

### Over-representation analysis
## length data
# check distribution of length data to observe potential outliers
# if large outliers, might be more accurate to use median
boxplot(length_matrix[1:50, ], main = "Gene Lengths Distribution", las = 2)
# Calculate the mean length for each gene across all samples
#gene_lengths <- rowMeans(length_matrix, na.rm = TRUE)
# Calculate the median length for each gene across all samples
gene_lengths <- apply(length_matrix, 1, mean, na.rm = TRUE)

# Set gene names to the calculated mean lengths
names(gene_lengths) <- rownames(length_matrix)

# Filter gene lengths for background genes
gene_length_background <- gene_lengths[names(gene_lengths) %in% background_genes_filtered]

# Check consistency
length(gene_length_background)  # Number of genes with lengths
all(background_genes_filtered %in% names(gene_length_background))  # Should return TRUE

# Final validation
is.numeric(gene_length_background)  # Should return TRUE
all(!is.na(names(gene_length_background)))  # Should return TRUE

# Inspect the final named numeric vector
str(gene_length_background)
head(gene_length_background)


### Prepare lists of background genes and DEGs - upregulated and downregualted DEGs analysed separately
length(background_genes)# 27891 genes
length(background_genes_filtered) #10491


# Convert the named numeric vector into a data frame
gene_length_background_df <- data.frame(
  gene = names(gene_length_background),
  length = as.numeric(gene_length_background) # Ensure numeric length values
)

# Inspect the resulting data frame
str(gene_length_background_df)
head(gene_length_background_df)

# Create a new object with just the gene names
background_genes.list <- names(gene_length_background)
head(background_genes.list)
str(background_genes.list)

# Get the unique entries in the 'contrast' column
unique_contrasts <- unique(heatwave.DEGs$contrast)
unique_contrasts

## new code 
# Filter DEGs to include only those with GO terms
heatwave.DEGs_filtered <- heatwave.DEGs %>%
  filter(gene %in% names(GO.list))

# Check how many DEGs are left after filtering
nrow(heatwave.DEGs_filtered)


# now changing heatwave.DEGs to heatwave.DEGs_filtered
#List of background genes resulting from the DGE analysis (this means genes with low counts are already filtered out)
heatwave_ALL.DGE_assayed <- heatwave.DEGs_filtered %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_ALL.DGE_assayed_list <- heatwave_ALL.DGE_assayed[['gene']]
str(heatwave_ALL.DGE_assayed_list)

#List of DEGs for each contrast (filter for significant genes i.e., padj < 0.05)
heatwave_DEGs_assayed.1 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group S1 vs C1') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.1 <- heatwave_DEGs_assayed.1[['gene']]

heatwave_DEGs_assayed.2 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group D1 vs C1') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.2 <- heatwave_DEGs_assayed.2[['gene']]

heatwave_DEGs_assayed.3 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group D1 vs S1') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.3 <- heatwave_DEGs_assayed.3[['gene']]

heatwave_DEGs_assayed.4 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group S2 vs C2') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.4 <- heatwave_DEGs_assayed.4[['gene']]

heatwave_DEGs_assayed.5 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group D2 vs C2') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.5 <- heatwave_DEGs_assayed.5[['gene']]

heatwave_DEGs_assayed.6 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group D2 vs S2') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.6 <- heatwave_DEGs_assayed.6[['gene']]

heatwave_DEGs_assayed.7 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group C1 vs C2') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.7 <- heatwave_DEGs_assayed.7[['gene']]

heatwave_DEGs_assayed.8 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group S1 vs S2') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.8 <- heatwave_DEGs_assayed.8[['gene']]

heatwave_DEGs_assayed.9 <- heatwave.DEGs_filtered %>%
  #filter(log2FoldChange > 0) %>%
  #filter(log2FoldChange < 0) %>%
  filter(contrast == 'group D1 vs D2') %>%
  filter(padj < 0.05) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  mutate_all(trimws)
heatwave_DEGs_assayed_list.9 <- heatwave_DEGs_assayed.9[['gene']]

str(heatwave_DEGs_assayed_list.1)
str(heatwave_DEGs_assayed_list.2)
str(heatwave_DEGs_assayed_list.3)
str(heatwave_DEGs_assayed_list.4)
str(heatwave_DEGs_assayed_list.5)
str(heatwave_DEGs_assayed_list.6)
str(heatwave_DEGs_assayed_list.7)
str(heatwave_DEGs_assayed_list.8)
str(heatwave_DEGs_assayed_list.9)
### Calculate probability weighting function using goseq
#Construct a new vector that adds a 0 next to every gene that is not in the DEG list and a 1 next to every gene that is in the DEG list
heatwave_gene.vector.1 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.1)
names(heatwave_gene.vector.1) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.2 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.2)
names(heatwave_gene.vector.2) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.3 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.3)
names(heatwave_gene.vector.3) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.4 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.4)
names(heatwave_gene.vector.4) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.5 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.5)
names(heatwave_gene.vector.5) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.6 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.6)
names(heatwave_gene.vector.6) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.7 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.7)
names(heatwave_gene.vector.7) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.8 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.8)
names(heatwave_gene.vector.8) = heatwave_ALL.DGE_assayed_list

heatwave_gene.vector.9 = as.integer(heatwave_ALL.DGE_assayed_list %in% heatwave_DEGs_assayed_list.9)
names(heatwave_gene.vector.9) = heatwave_ALL.DGE_assayed_list


gene_length_background_df <- gene_length_background_df %>%
  rename(gene = 1, length = 2) %>%
  distinct()
str(heatwave_ALL.DGE_assayed_list)

#Ensure that the list of transcript lengths matches heatwave_gene.vector
# Subset the gene_length_background_df using the character vector
heatwave_transcript.length.1 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
str(heatwave_transcript.length.1)
head(heatwave_transcript.length.1)
# Ensure the dimensions match expectations
length(unique(heatwave_transcript.length.1$gene))  # Check how many unique genes are retained
gene_length_heatwave_list.1 <- deframe(heatwave_transcript.length.1)
str(gene_length_heatwave_list.1)
names(gene_length_heatwave_list.1)

# Subset the gene_length_background_df using the character vector
heatwave_transcript.length.2 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
str(heatwave_transcript.length.2)
head(heatwave_transcript.length.2)
gene_length_heatwave_list.2 <- deframe(heatwave_transcript.length.2)

# Subset the gene_length_background_df using the character vector
heatwave_transcript.length.3 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.3$gene))  # Check how many unique genes are retained
gene_length_heatwave_list.3 <- deframe(heatwave_transcript.length.3)

# Subset the gene_length_background_df for contrast 4
heatwave_transcript.length.4 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.4$gene))  # Check how many unique genes are retained
gene_length_heatwave_list.4 <- deframe(heatwave_transcript.length.4)

heatwave_transcript.length.5 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.5$gene))
gene_length_heatwave_list.5 <- deframe(heatwave_transcript.length.5)

heatwave_transcript.length.6 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.6$gene))
gene_length_heatwave_list.6 <- deframe(heatwave_transcript.length.6)

heatwave_transcript.length.7 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.7$gene))
gene_length_heatwave_list.7 <- deframe(heatwave_transcript.length.7)

heatwave_transcript.length.8 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.8$gene))
gene_length_heatwave_list.8 <- deframe(heatwave_transcript.length.8)

heatwave_transcript.length.9 <- gene_length_background_df[gene_length_background_df$gene %in% heatwave_ALL.DGE_assayed_list, ]
length(unique(heatwave_transcript.length.9$gene))
gene_length_heatwave_list.9 <- deframe(heatwave_transcript.length.9)


#Ensure order of vectors are the same and prepare vectors for pwf() function 
heatwave_gene.vector.1 <- heatwave_gene.vector.1[order(names(heatwave_gene.vector.1))]
gene_length_heatwave_list.1 <- gene_length_heatwave_list.1[order(names(gene_length_heatwave_list.1))]
gene_length_heatwave_list.1<- unname(gene_length_heatwave_list.1)
gene_length_heatwave_list.1 <- as.integer(gene_length_heatwave_list.1)

heatwave_gene.vector.2 <- heatwave_gene.vector.2[order(names(heatwave_gene.vector.2))]
gene_length_heatwave_list.2 <- gene_length_heatwave_list.2[order(names(gene_length_heatwave_list.2))]
gene_length_heatwave_list.2 <- unname(gene_length_heatwave_list.2)
gene_length_heatwave_list.2 <- as.integer(gene_length_heatwave_list.2)

heatwave_gene.vector.3 <- heatwave_gene.vector.3[order(names(heatwave_gene.vector.3))]
gene_length_heatwave_list.3 <- gene_length_heatwave_list.3[order(names(gene_length_heatwave_list.3))]
gene_length_heatwave_list.3 <- unname(gene_length_heatwave_list.3)
gene_length_heatwave_list.3 <- as.integer(gene_length_heatwave_list.3)

heatwave_gene.vector.4 <- heatwave_gene.vector.4[order(names(heatwave_gene.vector.4))]
gene_length_heatwave_list.4 <- gene_length_heatwave_list.4[order(names(gene_length_heatwave_list.4))]
gene_length_heatwave_list.4 <- unname(gene_length_heatwave_list.4)
gene_length_heatwave_list.4 <- as.integer(gene_length_heatwave_list.4)

heatwave_gene.vector.5 <- heatwave_gene.vector.5[order(names(heatwave_gene.vector.5))]
gene_length_heatwave_list.5 <- gene_length_heatwave_list.5[order(names(gene_length_heatwave_list.5))]
gene_length_heatwave_list.5 <- unname(gene_length_heatwave_list.5)
gene_length_heatwave_list.5 <- as.integer(gene_length_heatwave_list.5)

heatwave_gene.vector.6 <- heatwave_gene.vector.6[order(names(heatwave_gene.vector.6))]
gene_length_heatwave_list.6 <- gene_length_heatwave_list.6[order(names(gene_length_heatwave_list.6))]
gene_length_heatwave_list.6 <- unname(gene_length_heatwave_list.6)
gene_length_heatwave_list.6 <- as.integer(gene_length_heatwave_list.6)

heatwave_gene.vector.7 <- heatwave_gene.vector.7[order(names(heatwave_gene.vector.7))]
gene_length_heatwave_list.7 <- gene_length_heatwave_list.7[order(names(gene_length_heatwave_list.7))]
gene_length_heatwave_list.7 <- unname(gene_length_heatwave_list.7)
gene_length_heatwave_list.7 <- as.integer(gene_length_heatwave_list.7)

heatwave_gene.vector.8 <- heatwave_gene.vector.8[order(names(heatwave_gene.vector.8))]
gene_length_heatwave_list.8 <- gene_length_heatwave_list.8[order(names(gene_length_heatwave_list.8))]
gene_length_heatwave_list.8 <- unname(gene_length_heatwave_list.8)
gene_length_heatwave_list.8 <- as.integer(gene_length_heatwave_list.8)

heatwave_gene.vector.9 <- heatwave_gene.vector.9[order(names(heatwave_gene.vector.9))]
gene_length_heatwave_list.9 <- gene_length_heatwave_list.9[order(names(gene_length_heatwave_list.9))]
gene_length_heatwave_list.9 <- unname(gene_length_heatwave_list.9)
gene_length_heatwave_list.9 <- as.integer(gene_length_heatwave_list.9)


#Weigh the gene vector by lengths of the transcripts
heatwave_pwf.1 = nullp(heatwave_gene.vector.1, bias.data = gene_length_heatwave_list.1)
heatwave_pwf.2 = nullp(heatwave_gene.vector.2, bias.data = gene_length_heatwave_list.2)
heatwave_pwf.3 = nullp(heatwave_gene.vector.3, bias.data = gene_length_heatwave_list.3)
heatwave_pwf.4 = nullp(heatwave_gene.vector.4, bias.data = gene_length_heatwave_list.4)
heatwave_pwf.5 = nullp(heatwave_gene.vector.5, bias.data = gene_length_heatwave_list.5)
heatwave_pwf.6 = nullp(heatwave_gene.vector.6, bias.data = gene_length_heatwave_list.6)
heatwave_pwf.7 = nullp(heatwave_gene.vector.7, bias.data = gene_length_heatwave_list.7)
heatwave_pwf.8 = nullp(heatwave_gene.vector.8, bias.data = gene_length_heatwave_list.8)
heatwave_pwf.9 = nullp(heatwave_gene.vector.9, bias.data = gene_length_heatwave_list.9)

### Get GO terms and carry out over-representation analysis
#Prepare list of GO terms for each gene
print(unique_contrasts)
length(GO.list)
nrow(GO.list)
#Carry out go GO enrichment and depletion (using goseq)
heatwave_goseq.df1 = goseq(heatwave_pwf.1, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group S1 vs C1')
heatwave_goseq.df2 = goseq(heatwave_pwf.2, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group D1 vs C1')
heatwave_goseq.df3 = goseq(heatwave_pwf.3, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group D1 vs S1')
heatwave_goseq.df4 = goseq(heatwave_pwf.4, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group S2 vs C2')
heatwave_goseq.df5 = goseq(heatwave_pwf.5, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group D2 vs C2')
heatwave_goseq.df6 = goseq(heatwave_pwf.6, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group D2 vs S2')
heatwave_goseq.df7 = goseq(heatwave_pwf.7, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group C1 vs C2')
heatwave_goseq.df8 = goseq(heatwave_pwf.8, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group S1 vs S2')
heatwave_goseq.df9 = goseq(heatwave_pwf.9, gene2cat = GO.list, method = "Wallenius", use_genes_without_cat = F) %>% mutate(contrast = 'group D1 vs D2')

# Convert GO.list to a data frame without incorrect row names
GO.df <- data.frame(
  gene = rep(names(GO.list), lengths(GO.list)), # Repeat gene names based on the number of GO terms
  GO_ID = unlist(GO.list),                     # Flatten the GO list into a single vector
  row.names = NULL                             # Reset row names
)

# Inspect the corrected data frame
str(GO.df)
head(GO.df)

### Pass gene IDs to enriched/depleted GO terms
#Pass IDs of DEGs to output from enrichment (gene column = gene IDs of numInCat)
#GO_DEG.df <- GO.df %>% filter(gene %in% heatwave_ALL.DGE_assayed$gene) # ythis was why everything from hre down was the same
# Filter GO terms for each contrast separately
heatwave_DEGs_contrast1 <- names(heatwave_gene.vector.1[heatwave_gene.vector.1 == 1])
heatwave_DEGs_contrast2 <- names(heatwave_gene.vector.2[heatwave_gene.vector.2 == 1])
heatwave_DEGs_contrast3 <- names(heatwave_gene.vector.3[heatwave_gene.vector.3 == 1])
heatwave_DEGs_contrast4 <- names(heatwave_gene.vector.4[heatwave_gene.vector.4 == 1])
heatwave_DEGs_contrast5 <- names(heatwave_gene.vector.5[heatwave_gene.vector.5 == 1])
heatwave_DEGs_contrast6 <- names(heatwave_gene.vector.6[heatwave_gene.vector.6 == 1])
heatwave_DEGs_contrast7 <- names(heatwave_gene.vector.7[heatwave_gene.vector.7 == 1])
heatwave_DEGs_contrast8 <- names(heatwave_gene.vector.8[heatwave_gene.vector.8 == 1])
heatwave_DEGs_contrast9 <- names(heatwave_gene.vector.9[heatwave_gene.vector.9 == 1])

str(heatwave_DEGs_contrast1) #666 #new 296
str(heatwave_DEGs_contrast2) #1123 #new 506
str(heatwave_DEGs_contrast3) #2 #new 1
str(heatwave_DEGs_contrast4) #132 #new 75
str(heatwave_DEGs_contrast5) #1668 #new 840
str(heatwave_DEGs_contrast6) #95 #new 48
str(heatwave_DEGs_contrast7) #1329 #new 735
str(heatwave_DEGs_contrast8) #2468 #new 1095
str(heatwave_DEGs_contrast9) #2896 #new 1431

GO_DEG.df1 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast1)
GO_DEG.df2 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast2)
GO_DEG.df3 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast3)
GO_DEG.df4 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast4)
GO_DEG.df5 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast5)
GO_DEG.df6 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast6)
GO_DEG.df7 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast7)
GO_DEG.df8 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast8)
GO_DEG.df9 <- GO.df %>% filter(gene %in% heatwave_DEGs_contrast9)

# Aggregate gene IDs for each contrast
GO.df1 <- aggregate(gene ~ GO_ID, GO_DEG.df1, FUN = function(x) paste(unique(x), collapse = "; ")) %>% rename('category' = 1)
GO.df2 <- aggregate(gene ~ GO_ID, GO_DEG.df2, FUN = function(x) paste(unique(x), collapse = "; ")) %>% rename('category' = 1)
GO.df3 <- aggregate(gene ~ GO_ID, GO_DEG.df3, FUN = function(x) paste(unique(x), collapse = "; ")) %>% rename('category' = 1)
GO.df4 <- aggregate(gene ~ GO_ID, GO_DEG.df4, FUN = function(x) paste(unique(x), collapse = "; ")) %>%
  rename('category' = 1)
GO.df5 <- aggregate(gene ~ GO_ID, GO_DEG.df5, FUN = function(x) paste(unique(x), collapse = "; ")) %>%
  rename('category' = 1)
GO.df6 <- aggregate(gene ~ GO_ID, GO_DEG.df6, FUN = function(x) paste(unique(x), collapse = "; ")) %>%
  rename('category' = 1)
GO.df7 <- aggregate(gene ~ GO_ID, GO_DEG.df7, FUN = function(x) paste(unique(x), collapse = "; ")) %>%
  rename('category' = 1)
GO.df8 <- aggregate(gene ~ GO_ID, GO_DEG.df8, FUN = function(x) paste(unique(x), collapse = "; ")) %>%
  rename('category' = 1)
GO.df9 <- aggregate(gene ~ GO_ID, GO_DEG.df9, FUN = function(x) paste(unique(x), collapse = "; ")) %>%
  rename('category' = 1)

head(GO.df1)
head(GO.df2)
head(GO.df3)
head(GO.df4)
head(GO.df5)
head(GO.df6)
head(GO.df7)
head(GO.df8)
head(GO.df9)

# Number of rows in each data frame
cat("Number of rows in GO.df1:", nrow(GO.df1), "\n") #4384 #4385
cat("Number of rows in GO.df2:", nrow(GO.df2), "\n") #5546 #5496
cat("Number of rows in GO.df3:", nrow(GO.df3), "\n") #277 #277
cat("Number of rows in GO.df4:", nrow(GO.df4), "\n") #1742 #1766
cat("Number of rows in GO.df5:", nrow(GO.df5), "\n") #6936 #6899
cat("Number of rows in GO.df6:", nrow(GO.df6), "\n") #1175 #1175
cat("Number of rows in GO.df7:", nrow(GO.df7), "\n") #6012 #5866
cat("Number of rows in GO.df8:", nrow(GO.df8), "\n") #7999 #7926
cat("Number of rows in GO.df9:", nrow(GO.df9), "\n") #8676 #8658

heatwave_goseq.df1 <- left_join(heatwave_goseq.df1, GO.df1, by = 'category')
heatwave_goseq.df2 <- left_join(heatwave_goseq.df2, GO.df2, by = 'category')
heatwave_goseq.df3 <- left_join(heatwave_goseq.df3, GO.df3, by = 'category')
heatwave_goseq.df4 <- left_join(heatwave_goseq.df4, GO.df4, by = 'category')
heatwave_goseq.df5 <- left_join(heatwave_goseq.df5, GO.df5, by = 'category')
heatwave_goseq.df6 <- left_join(heatwave_goseq.df6, GO.df6, by = 'category')
heatwave_goseq.df7 <- left_join(heatwave_goseq.df7, GO.df7, by = 'category')
heatwave_goseq.df8 <- left_join(heatwave_goseq.df8, GO.df8, by = 'category')
heatwave_goseq.df9 <- left_join(heatwave_goseq.df9, GO.df9, by = 'category')

# Define the keywords to remove
unwanted_terms <- c("interleukin", "neuron", "neural", "gastrulation", 
                    "synaptic", "adhesion", "response to starvation", "central nervous system", 
                    "synapse", "pancreas", "pollen", "floral", "bone", "exocrine", 
                    "placenta", "dendrite", "dendritic", "eye", "odorant", "heart", 
                    "blood", "pain", "glial", "endothelial", "behavior", "locomotor", 
                    "insulin", "tRNA", "T cell", "L-leucine", "myotube", "oligodendrocyte",
                    "macrolide", "SMAD", "Notch", "chromosome", "ecdysteroid", "innate",
                    "transesterification", "organ", "vitamin", "post-embryonic", "dosage",
                    "viral", "epithelial", "canonical Wnt", "reproductive", "immune", "integration",
                    "development", "estrogen", "axon", "virus", "symbiont-mediated",
                    "mucopolysaccharide", "spinal", "water", "spermatogenesis", "male", "female", "bile",
                    "osteoblast",  "wnt", "insect", "face", "muscle", "cardiac", "node", "stem cell")

# Function to filter out unwanted terms
remove_unwanted_terms <- function(data, unwanted_terms) {
  data %>%
    filter(!str_detect(term, paste(unwanted_terms, collapse = "|")))
}

# Filter each data frame to remove unwanted terms
heatwave_goseq.df1 <- remove_unwanted_terms(heatwave_goseq.df1, unwanted_terms)
heatwave_goseq.df2 <- remove_unwanted_terms(heatwave_goseq.df2, unwanted_terms)
heatwave_goseq.df3 <- remove_unwanted_terms(heatwave_goseq.df3, unwanted_terms)
heatwave_goseq.df4 <- remove_unwanted_terms(heatwave_goseq.df4, unwanted_terms)
heatwave_goseq.df5 <- remove_unwanted_terms(heatwave_goseq.df5, unwanted_terms)
heatwave_goseq.df6 <- remove_unwanted_terms(heatwave_goseq.df6, unwanted_terms)
heatwave_goseq.df7 <- remove_unwanted_terms(heatwave_goseq.df7, unwanted_terms)
heatwave_goseq.df8 <- remove_unwanted_terms(heatwave_goseq.df8, unwanted_terms)
heatwave_goseq.df9 <- remove_unwanted_terms(heatwave_goseq.df9, unwanted_terms)

### Keep statistically significant GO terms and check the GO IDs for missing annotatations 
heatwave_enriched_signif.1 <- heatwave_goseq.df1 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method="BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')

heatwave_enriched_signif.1 <- heatwave_goseq.df1 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )


heatwave_enriched_signif.2 <- heatwave_goseq.df2 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method="BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 3) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')
heatwave_enriched_signif.2 <- heatwave_goseq.df2 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 3) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )


heatwave_enriched_signif.3 <- heatwave_goseq.df3 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method="BH")) %>%
  filter(over_represented_pvalue < 0.1) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.1 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')
heatwave_enriched_signif.3 <- heatwave_goseq.df3 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )

heatwave_enriched_signif.4 <- heatwave_goseq.df4 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')
heatwave_enriched_signif.4 <- heatwave_goseq.df4 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )

heatwave_enriched_signif.5 <- heatwave_goseq.df5 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')
heatwave_enriched_signif.5 <- heatwave_goseq.df5 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )

heatwave_enriched_signif.6 <- heatwave_goseq.df6 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')

heatwave_enriched_signif.6 <- heatwave_goseq.df6 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )


heatwave_enriched_signif.7 <- heatwave_goseq.df7 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')

heatwave_enriched_signif.7 <- heatwave_goseq.df7 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )

heatwave_enriched_signif.8 <- heatwave_goseq.df8 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')

heatwave_enriched_signif.8 <- heatwave_goseq.df8 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )


heatwave_enriched_signif.9 <- heatwave_goseq.df9 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(over_represented_pvalue < 0.05) %>%
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(over_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = over_represented_pvalue, direction = 'enriched')

heatwave_enriched_signif.9 <- heatwave_goseq.df9 %>%
  mutate(padj_BH = p.adjust(over_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = over_represented_pvalue,  # Rename for easier reference
    direction = 'enriched'
  )


# Number of rows for each enriched significance data frame
nrow(heatwave_enriched_signif.1) #106 #112
head(heatwave_enriched_signif.1)
ora_significant1 <- heatwave_enriched_signif.1 %>%
  filter(signif == 'yes')
nrow(ora_significant1)

nrow(heatwave_enriched_signif.2) #306
#heatwave_enriched_signif.2 %>%
  filter(signif == 'yes') %>%
  head()  # Show the first few rows

nrow(heatwave_enriched_signif.3) #0
nrow(heatwave_enriched_signif.4) #17
nrow(heatwave_enriched_signif.5) #121
nrow(heatwave_enriched_signif.6) #5
nrow(heatwave_enriched_signif.7) #237
nrow(heatwave_enriched_signif.8) #270
nrow(heatwave_enriched_signif.9) #177

# Now repeat for depleted terms
heatwave_depleted_signif.1 <- heatwave_goseq.df1 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.1 <- heatwave_goseq.df1 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.2 <- heatwave_goseq.df2 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.2 <- heatwave_goseq.df2 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )



heatwave_depleted_signif.3 <- heatwave_goseq.df3 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.1) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.1 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.3 <- heatwave_goseq.df3 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.4 <- heatwave_goseq.df4 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.4 <- heatwave_goseq.df4 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.5 <- heatwave_goseq.df5 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.5 <- heatwave_goseq.df5 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.6 <- heatwave_goseq.df6 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.6 <- heatwave_goseq.df6 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.7 <- heatwave_goseq.df7 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method="BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.7 <- heatwave_goseq.df7 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.8 <- heatwave_goseq.df8 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.8 <- heatwave_goseq.df8 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

heatwave_depleted_signif.9 <- heatwave_goseq.df9 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(under_represented_pvalue < 0.05) %>% 
  filter(numDEInCat >= 5) %>%
  mutate(signif = case_when(under_represented_pvalue & padj_BH < 0.05 ~ 'yes',
                            TRUE ~ 'no')) %>%
  mutate(pval = under_represented_pvalue, direction = 'depleted')

heatwave_depleted_signif.9 <- heatwave_goseq.df9 %>%
  mutate(padj_BH = p.adjust(under_represented_pvalue, method = "BH")) %>%
  filter(numDEInCat >= 5) %>%  # Only keep GO terms with at least 5 DEGs
  mutate(
    signif = case_when(
      padj_BH < 0.05 ~ 'yes',                   # Adjusted p-value < 0.05 = significant
      #over_represented_pvalue < 0.05 ~ 'raw_signif',  # Raw p-value < 0.05 = exploratory significance
      TRUE ~ 'no'
    ),
    pval = under_represented_pvalue,  # Rename for easier reference
    direction = 'depleted'
  )

# Number of rows for each enriched significance data frame
nrow(heatwave_depleted_signif.1) #38
nrow(heatwave_depleted_signif.2) #69
nrow(heatwave_depleted_signif.3) #0
nrow(heatwave_depleted_signif.4) #8
nrow(heatwave_depleted_signif.5) #82
nrow(heatwave_depleted_signif.6) #9
nrow(heatwave_depleted_signif.7) #223
nrow(heatwave_depleted_signif.8) #121
nrow(heatwave_depleted_signif.9) #115

heatwave_GO_signif_annot <- rbind(
  heatwave_enriched_signif.1, heatwave_enriched_signif.2, heatwave_enriched_signif.3,
  heatwave_enriched_signif.4, heatwave_enriched_signif.5, heatwave_enriched_signif.6,
  heatwave_enriched_signif.7, heatwave_enriched_signif.8, heatwave_enriched_signif.9
) %>%
  rbind(
    heatwave_depleted_signif.1, heatwave_depleted_signif.2, heatwave_depleted_signif.3,
    heatwave_depleted_signif.4, heatwave_depleted_signif.5, heatwave_depleted_signif.6,
    heatwave_depleted_signif.7, heatwave_depleted_signif.8, heatwave_depleted_signif.9
  ) %>%
  dplyr::select(!c(over_represented_pvalue, under_represented_pvalue))

head(heatwave_GO_signif_annot)

# filter out significant terms
heatwave_GO_signif_annot_filtered <- heatwave_GO_signif_annot %>%
  filter(signif == "yes")

# View the filtered data
head(heatwave_GO_signif_annot_filtered)

#REVIGO
# Filter for the specific contrast and ontology (e.g., Biological Process - BP)
revigo_input_1 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group S1 vs C1"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_1)
write.table(revigo_input_1, file = "revigo_input_1.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_2 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D1 vs C1"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_2)
write.table(revigo_input_2, file = "revigo_input_2.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_3 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D1 vs S1"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_3)
write.table(revigo_input_3, file = "revigo_input_3.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_4 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group S2 vs C2"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_4)
write.table(revigo_input_4, file = "revigo_input_4.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_5 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D2 vs C2"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_5)
write.table(revigo_input_5, file = "revigo_input_5.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_6 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D2 vs S2"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_6)
write.table(revigo_input_6, file = "revigo_input_6.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_7 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group C1 vs C2"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_7) #133
write.table(revigo_input_7, file = "revigo_input_7.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_8 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D1 vs S2"#, ontology == "BP"
  ) %>% 
  dplyr::select(category, pval)
head(revigo_input_8)
write.table(revigo_input_8, file = "revigo_input_8.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

revigo_input_9 <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D1 vs D2"#, ontology == "BP"
         ) %>% 
  dplyr::select(category, pval)

# Check the first few rows
head(revigo_input_9) #11
write.table(revigo_input_9, file = "revigo_input_9.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

nrow(revigo_input_1)
nrow(revigo_input_2)
nrow(revigo_input_3)
nrow(revigo_input_4)
nrow(revigo_input_5)
nrow(revigo_input_6)
nrow(revigo_input_7)
nrow(revigo_input_8)
nrow(revigo_input_9)

# Check if any GO terms exist for that contrast without significance filtering
test_contrast <- heatwave_GO_signif_annot_filtered %>%
  filter(contrast == "group D1 vs D2")

# View the first few rows
head(test_contrast)


## ORA
# Function to extract top enriched terms for a given contrast
get_top_terms_no_genes_enriched <- function(data, contrast, regulation = "enriched", top_n = 15) {
  data %>%
    filter(contrast == !!contrast, direction == regulation) %>%  # Filter by contrast and direction
    arrange(padj_BH) %>%                                            # Sort by p-value (ascending)
    head(top_n) %>%                                              # Select the top N terms
    dplyr::select(term, ontology, pval, padj_BH, numDEInCat, numInCat)  # Select relevant columns
}

revigo_ora_1 <- heatwave_enriched_signif.1 %>% arrange(pval) %>%
  filter(pval < 0.05) %>%
  dplyr::select(category, pval)


# Apply function to each contrast
# Extract top enriched terms for each contrast
top_enriched_S1_vs_C1 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.1, "group S1 vs C1")
top_enriched_D1_vs_C1 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.2, "group D1 vs C1")
top_enriched_D1_vs_S1 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.3, "group D1 vs S1")
top_enriched_S2_vs_C2 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.4, "group S2 vs C2")
top_enriched_D2_vs_C2 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.5, "group D2 vs C2")
top_enriched_D2_vs_S2 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.6, "group D2 vs S2")
top_enriched_C1_vs_C2 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.7, "group C1 vs C2")
top_enriched_S1_vs_S2 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.8, "group S1 vs S2")
top_enriched_D1_vs_D2 <- get_top_terms_no_genes_enriched(heatwave_enriched_signif.9, "group D1 vs D2")

# Plotting enriched GO terms for each contrast separately
# Function to plot enriched GO terms for a specific contrast
plot_enriched_GO_terms <- function(data, contrast_label) {
  data %>%
    filter(ontology == "BP") %>%  # Filter for Biological Processes (adjust as needed)
    mutate(gene_ratio = numDEInCat / numInCat) %>%  # Calculate gene ratio
    ggplot(aes(x = gene_ratio, y = fct_reorder(term, gene_ratio), size = numDEInCat, color = -log10(pval))) +
    geom_point() +
    labs(
      title = paste("Enriched GO Terms -", contrast_label),
      x = "Gene Ratio",
      y = "GO Term",
      color = "-log10(p-value)",
      size = "DEGs in GO Term"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White plot background
      plot.background = element_rect(fill = "white", color = NA),   # White outer background
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
}


# Generate plots for all contrasts
plot_S1_vs_C1_enriched_bp <- plot_enriched_GO_terms(top_enriched_S1_vs_C1, "group S1 vs C1")
plot_D1_vs_C1_enriched_bp <- plot_enriched_GO_terms(top_enriched_D1_vs_C1, "group D1 vs C1")
plot_D1_vs_S1_enriched_bp <- plot_enriched_GO_terms(top_enriched_D1_vs_S1, "group D1 vs S1")
plot_S2_vs_C2_enriched_bp <- plot_enriched_GO_terms(top_enriched_S2_vs_C2, "group S2 vs C2")
plot_D2_vs_C2_enriched_bp <- plot_enriched_GO_terms(top_enriched_D2_vs_C2, "group D2 vs C2")
plot_D2_vs_S2_enriched_bp <- plot_enriched_GO_terms(top_enriched_D2_vs_S2, "group D2 vs S2")
plot_C1_vs_C2_enriched_bp <- plot_enriched_GO_terms(top_enriched_C1_vs_C2, "group C1 vs C2")
plot_S1_vs_S2_enriched_bp <- plot_enriched_GO_terms(top_enriched_S1_vs_S2, "group S1 vs S2")
plot_D1_vs_D2_enriched_bp <- plot_enriched_GO_terms(top_enriched_D1_vs_D2, "group D1 vs D2")

# Display all plots
print(plot_S1_vs_C1_enriched_bp)
print(plot_D1_vs_C1_enriched_bp)
print(plot_D1_vs_S1_enriched_bp)
print(plot_S2_vs_C2_enriched_bp)
print(plot_D2_vs_C2_enriched_bp)
print(plot_D2_vs_S2_enriched_bp)
print(plot_C1_vs_C2_enriched_bp)
print(plot_S1_vs_S2_enriched_bp)
print(plot_D1_vs_D2_enriched_bp)


# Save plots separately
ggsave("GO_terms_S1_vs_C1_enriched_bp.png", plot = plot_S1_vs_C1_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D1_vs_C1_enriched_bp.png", plot = plot_D1_vs_C1_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D1_vs_S1_enriched_bp.png", plot = plot_D1_vs_S1_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_S2_vs_C2_enriched_bp.png", plot = plot_S2_vs_C2_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D2_vs_C2_enriched_bp.png", plot = plot_D2_vs_C2_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D2_vs_S2_enriched_bp.png", plot = plot_D2_vs_S2_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_C1_vs_C2_enriched_bp.png", plot = plot_C1_vs_C2_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_S1_vs_S2_enriched_bp.png", plot = plot_S1_vs_S2_enriched_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D1_vs_D2_enriched_bp.png", plot = plot_D1_vs_D2_enriched_bp, width = 8, height = 6, dpi = 300)

# Depleted
# Function to extract top depleted terms for a given contrast
get_top_terms_no_genes_depleted <- function(data, contrast, regulation = "depleted", top_n = 15) {
  data %>%
    filter(contrast == !!contrast, direction == regulation) %>%  # Filter by contrast and direction
    arrange(pval) %>%                                            # Sort by p-value (ascending)
    head(top_n) %>%                                              # Select the top N terms
    dplyr::select(term, ontology, pval, padj_BH, numDEInCat, numInCat)  # Select relevant columns
}

# Apply function to each contrast
# Extract top depleted terms for each contrast
top_depleted_S1_vs_C1 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.1, "group S1 vs C1")
top_depleted_D1_vs_C1 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.2, "group D1 vs C1")
top_depleted_D1_vs_S1 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.3, "group D1 vs S1")
top_depleted_S2_vs_C2 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.4, "group S2 vs C2")
top_depleted_D2_vs_C2 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.5, "group D2 vs C2")
top_depleted_D2_vs_S2 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.6, "group D2 vs S2")
top_depleted_C1_vs_C2 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.7, "group C1 vs C2")
top_depleted_S1_vs_S2 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.8, "group S1 vs S2")
top_depleted_D1_vs_D2 <- get_top_terms_no_genes_depleted(heatwave_depleted_signif.9, "group D1 vs D2")

# Plotting depleted GO terms for each contrast separately
# Function to plot depleted GO terms for a specific contrast
plot_depleted_GO_terms <- function(data, contrast_label) {
  data %>%
    filter(ontology == "BP") %>%  # Filter for Biological Processes (adjust as needed)
    mutate(gene_ratio = numDEInCat / numInCat) %>%  # Calculate gene ratio
    ggplot(aes(x = gene_ratio, y = fct_reorder(term, gene_ratio), size = numDEInCat, color = -log10(pval))) +
    geom_point() +
    labs(
      title = paste("depleted GO Terms -", contrast_label),
      x = "Gene Ratio",
      y = "GO Term",
      color = "-log10(p-value)",
      size = "DEGs in GO Term"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White plot background
      plot.background = element_rect(fill = "white", color = NA),   # White outer background
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
}


# Generate plots for all contrasts
plot_S1_vs_C1_depleted_bp <- plot_depleted_GO_terms(top_depleted_S1_vs_C1, "group S1 vs C1")
plot_D1_vs_C1_depleted_bp <- plot_depleted_GO_terms(top_depleted_D1_vs_C1, "group D1 vs C1")
plot_D1_vs_S1_depleted_bp <- plot_depleted_GO_terms(top_depleted_D1_vs_S1, "group D1 vs S1")
plot_S2_vs_C2_depleted_bp <- plot_depleted_GO_terms(top_depleted_S2_vs_C2, "group S2 vs C2")
plot_D2_vs_C2_depleted_bp <- plot_depleted_GO_terms(top_depleted_D2_vs_C2, "group D2 vs C2")
plot_D2_vs_S2_depleted_bp <- plot_depleted_GO_terms(top_depleted_D2_vs_S2, "group D2 vs S2")
plot_C1_vs_C2_depleted_bp <- plot_depleted_GO_terms(top_depleted_C1_vs_C2, "group C1 vs C2")
plot_S1_vs_S2_depleted_bp <- plot_depleted_GO_terms(top_depleted_S1_vs_S2, "group S1 vs S2")
plot_D1_vs_D2_depleted_bp <- plot_depleted_GO_terms(top_depleted_D1_vs_D2, "group D1 vs D2")

# Display all plots
print(plot_S1_vs_C1_depleted_bp)
print(plot_D1_vs_C1_depleted_bp)
print(plot_D1_vs_S1_depleted_bp)
print(plot_S2_vs_C2_depleted_bp)
print(plot_D2_vs_C2_depleted_bp)
print(plot_D2_vs_S2_depleted_bp)
print(plot_C1_vs_C2_depleted_bp)
print(plot_S1_vs_S2_depleted_bp)
print(plot_D1_vs_D2_depleted_bp)


# Save plots separately
ggsave("GO_terms_S1_vs_C1_depleted_bp.png", plot = plot_S1_vs_C1_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D1_vs_C1_depleted_bp.png", plot = plot_D1_vs_C1_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D1_vs_S1_depleted_bp.png", plot = plot_D1_vs_S1_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_S2_vs_C2_depleted_bp.png", plot = plot_S2_vs_C2_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D2_vs_C2_depleted_bp.png", plot = plot_D2_vs_C2_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D2_vs_S2_depleted_bp.png", plot = plot_D2_vs_S2_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_C1_vs_C2_depleted_bp.png", plot = plot_C1_vs_C2_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_S1_vs_S2_depleted_bp.png", plot = plot_S1_vs_S2_depleted_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_terms_D1_vs_D2_depleted_bp.png", plot = plot_D1_vs_D2_depleted_bp, width = 8, height = 6, dpi = 300)

# Combine plots for each contrast
# Function to combine enriched and depleted plots for each contrast

## new function
combine_GO_plots <- function(enriched_plot, depleted_plot, contrast_label) {
  enriched_plot + depleted_plot +
    plot_layout(ncol = 2, widths = c(1, 1)) +  # Adjusted widths to give more space
    plot_annotation(title = paste("Enriched and Depleted GO Terms:", contrast_label)) &
    theme(
      axis.text.y = element_text(size = 10),  # Adjust font size for better readability
      plot.margin = margin(10, 10, 10, 10)    # Add more space around plots
    )
}


# Combine plots for all contrasts
combined_plot_S1_vs_C1 <- combine_GO_plots(plot_S1_vs_C1_enriched_bp, plot_S1_vs_C1_depleted_bp, "group S1 vs C1")
combined_plot_D1_vs_C1 <- combine_GO_plots(plot_D1_vs_C1_enriched_bp, plot_D1_vs_C1_depleted_bp, "group D1 vs C1")
combined_plot_D1_vs_S1 <- combine_GO_plots(plot_D1_vs_S1_enriched_bp, plot_D1_vs_S1_depleted_bp, "group D1 vs S1")
combined_plot_S2_vs_C2 <- combine_GO_plots(plot_S2_vs_C2_enriched_bp, plot_S2_vs_C2_depleted_bp, "group S2 vs C2")
combined_plot_D2_vs_C2 <- combine_GO_plots(plot_D2_vs_C2_enriched_bp, plot_D2_vs_C2_depleted_bp, "group D2 vs C2")
combined_plot_D2_vs_S2 <- combine_GO_plots(plot_D2_vs_S2_enriched_bp, plot_D2_vs_S2_depleted_bp, "group D2 vs S2")
combined_plot_C1_vs_C2 <- combine_GO_plots(plot_C1_vs_C2_enriched_bp, plot_C1_vs_C2_depleted_bp, "group C1 vs C2")
combined_plot_S1_vs_S2 <- combine_GO_plots(plot_S1_vs_S2_enriched_bp, plot_S1_vs_S2_depleted_bp, "group S1 vs S2")
combined_plot_D1_vs_D2 <- combine_GO_plots(plot_D1_vs_D2_enriched_bp, plot_D1_vs_D2_depleted_bp, "group D1 vs D2")

# Display combined plots
print(combined_plot_S1_vs_C1)
print(combined_plot_D1_vs_C1)
print(combined_plot_D1_vs_S1)
print(combined_plot_S2_vs_C2)
print(combined_plot_D2_vs_C2)
print(combined_plot_D2_vs_S2)
print(combined_plot_C1_vs_C2)
print(combined_plot_S1_vs_S2)
print(combined_plot_D1_vs_D2)

## This now works for combined plots rather than the function above 4/2/2025
# Combine plots for S1 vs C1
combined_plotBP_S1_vs_C1 <- plot_S1_vs_C1_enriched_bp + plot_S1_vs_C1_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): S1 vs C1") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_S1_vs_C1)
png("combined_plotBP_S1_vs_C1.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_S1_vs_C1)
dev.off()
# Export the combined plot for 26C vs 16C
ggsave("combined_plotBP_S1_vs_C1.png", 
       plot = combined_plotBP_S1_vs_C1, 
       width = 10, height = 5, units = "in", dpi = 300)

# Combine plots for D1 vs C1
combined_plotBP_D1_vs_C1 <- plot_D1_vs_C1_enriched_bp + plot_D1_vs_C1_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): D1 vs C1") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_D1_vs_C1)
png("combined_plotBP_D1_vs_C1.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_D1_vs_C1)
dev.off()
ggsave("combined_plotBP_D1_vs_C1.png", 
       plot = combined_plotBP_D1_vs_C1, 
       width = 10, height = 5, units = "in", dpi = 300)

# Combine plots for D1 vs S1
combined_plotBP_D1_vs_S1 <- plot_D1_vs_S1_enriched_bp + plot_D1_vs_S1_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): D1 vs S1") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_D1_vs_S1)
png("combined_plotBP_D1_vs_S1.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_D1_vs_S1)
dev.off()
ggsave("combined_plotBP_D1_vs_S1.png", 
       plot = combined_plotBP_D1_vs_S1, 
       width = 10, height = 5, units = "in", dpi = 300)


# Combine plots for S2 vs C2
combined_plotBP_S2_vs_C2 <- plot_S2_vs_C2_enriched_bp + plot_S2_vs_C2_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): S2 vs C2") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_S2_vs_C2)
png("combined_plotBP_S2_vs_C2.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_S2_vs_C2)
dev.off()
ggsave("combined_plotBP_S2_vs_C2.png", 
       plot = combined_plotBP_S2_vs_C2, 
       width = 10, height = 5, units = "in", dpi = 300)


# Combine plots for D2 vs C2
combined_plotBP_D2_vs_C2 <- plot_D2_vs_C2_enriched_bp + plot_D2_vs_C2_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): D2 vs C2") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_D2_vs_C2)
png("combined_plotBP_D2_vs_C2.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_D2_vs_C2)
dev.off()
ggsave("combined_plotBP_D2_vs_C2.png", 
       plot = combined_plotBP_D2_vs_C2, 
       width = 10, height = 5, units = "in", dpi = 300)

# Combine plots for D2 vs S2
combined_plotBP_D2_vs_S2 <- plot_D2_vs_S2_enriched_bp + plot_D2_vs_S2_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): D2 vs S2") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_D2_vs_S2)
png("combined_plotBP_D2_vs_S2.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_D2_vs_S2)
dev.off()
ggsave("combined_plotBP_D2_vs_S2.png", 
       plot = combined_plotBP_D2_vs_S2, 
       width = 10, height = 5, units = "in", dpi = 300)


# Combine plots for C1 vs C2
combined_plotBP_C1_vs_C2 <- plot_C1_vs_C2_enriched_bp + plot_C1_vs_C2_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): C1 vs C2") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_C1_vs_C2)
png("combined_plotBP_C1_vs_C2.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_C1_vs_C2)
dev.off()
ggsave("combined_plotBP_C1_vs_C2.png", 
       plot = combined_plotBP_C1_vs_C2, 
       width = 10, height = 5, units = "in", dpi = 300)

# Combine plots for S1 vs S2
combined_plotBP_S1_vs_S2 <- plot_S1_vs_S2_enriched_bp + plot_S1_vs_S2_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): S1 vs S2") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_S1_vs_S2)
png("combined_plotBP_S1_vs_S2.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_S1_vs_S2)
dev.off()
ggsave("combined_plotBP_S1_vs_S2.png", 
       plot = combined_plotBP_S1_vs_S2, 
       width = 10, height = 5, units = "in", dpi = 300)

# Combine plots for D1 vs D2
combined_plotBP_D1_vs_D2 <- plot_D1_vs_D2_enriched_bp + plot_D1_vs_D2_depleted_bp +
  plot_layout(ncol = 2, widths = c(1, 1))  +
  plot_annotation(title = "Enriched and Depleted GO Terms (Biological Process): D1 vs D2") &
  theme(axis.text.y = element_text(size = 10))
print(combined_plotBP_D1_vs_D2)
png("combined_plotBP_D1_vs_D2.png", width = 15, height = 10, units = "in", res = 600)
print(combined_plotBP_D1_vs_D2)
dev.off()
ggsave("combined_plotBP_D1_vs_D2.png", 
       plot = combined_plotBP_D1_vs_D2, 
       width = 10, height = 5, units = "in", dpi = 300)

# Save combined plots
ggsave("combined_GO_terms_S1_vs_C1.png", plot = combined_plot_S1_vs_C1, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_D1_vs_C1.png", plot = combined_plot_D1_vs_C1, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_D1_vs_S1.png", plot = combined_plot_D1_vs_S1, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_S2_vs_C2.png", plot = combined_plot_S2_vs_C2, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_D2_vs_C2.png", plot = combined_plot_D2_vs_C2, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_D2_vs_S2.png", plot = combined_plot_D2_vs_S2, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_C1_vs_C2.png", plot = combined_plot_C1_vs_C2, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_S1_vs_S2.png", plot = combined_plot_S1_vs_S2, width = 12, height = 8, dpi = 400)
ggsave("combined_GO_terms_D1_vs_D2.png", plot = combined_plot_D1_vs_D2, width = 12, height = 8, dpi = 400)

## GSEA
## GSEA for Heatwave Project
# Create ranked lists for each contrast
heatwave.DEGs.gsea_groupS1vsC1 <- heatwave.DEGs %>%
  filter(contrast == 'group S1 vs C1') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupD1vsC1 <- heatwave.DEGs %>%
  filter(contrast == 'group D1 vs C1') %>%
  distinct(gene, .keep_all = TRUE) %>% 
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupD1vsS1 <- heatwave.DEGs %>%
  filter(contrast == 'group D1 vs S1') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupS2vsC2 <- heatwave.DEGs %>%
  filter(contrast == 'group S2 vs C2') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupD2vsC2 <- heatwave.DEGs %>%
  filter(contrast == 'group D2 vs C2') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupD2vsS2 <- heatwave.DEGs %>%
  filter(contrast == 'group D2 vs S2') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupC1vsC2 <- heatwave.DEGs %>%
  filter(contrast == 'group C1 vs C2') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupS1vsS2 <- heatwave.DEGs %>%
  filter(contrast == 'group S1 vs S2') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

heatwave.DEGs.gsea_groupD1vsD2 <- heatwave.DEGs %>%
  filter(contrast == 'group D1 vs D2') %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, log2FoldChange) %>%
  arrange(-log2FoldChange) %>%
  pull(log2FoldChange, name = gene)

# Introduce minor variations to avoid duplicate rankings
tied_values_S1vsC1 <- which(duplicated(heatwave.DEGs.gsea_groupS1vsC1) | duplicated(heatwave.DEGs.gsea_groupS1vsC1, fromLast = TRUE))
tied_values_D1vsC1 <- which(duplicated(heatwave.DEGs.gsea_groupD1vsC1) | duplicated(heatwave.DEGs.gsea_groupD1vsC1, fromLast = TRUE))
tied_values_D1vsS1 <- which(duplicated(heatwave.DEGs.gsea_groupD1vsS1) | duplicated(heatwave.DEGs.gsea_groupD1vsS1, fromLast = TRUE))
tied_values_S2vsC2 <- which(duplicated(heatwave.DEGs.gsea_groupS2vsC2) | duplicated(heatwave.DEGs.gsea_groupS2vsC2, fromLast = TRUE))
tied_values_D2vsC2 <- which(duplicated(heatwave.DEGs.gsea_groupD2vsC2) | duplicated(heatwave.DEGs.gsea_groupD2vsC2, fromLast = TRUE))
tied_values_D2vsS2 <- which(duplicated(heatwave.DEGs.gsea_groupD2vsS2) | duplicated(heatwave.DEGs.gsea_groupD2vsS2, fromLast = TRUE))
tied_values_C1vsC2 <- which(duplicated(heatwave.DEGs.gsea_groupC1vsC2) | duplicated(heatwave.DEGs.gsea_groupC1vsC2, fromLast = TRUE))
tied_values_S1vsS2 <- which(duplicated(heatwave.DEGs.gsea_groupS1vsS2) | duplicated(heatwave.DEGs.gsea_groupS1vsS2, fromLast = TRUE))
tied_values_D1vsD2 <- which(duplicated(heatwave.DEGs.gsea_groupD1vsD2) | duplicated(heatwave.DEGs.gsea_groupD1vsD2, fromLast = TRUE))

heatwave.DEGs.gsea_groupS1vsC1[tied_values_S1vsC1] <- heatwave.DEGs.gsea_groupS1vsC1[tied_values_S1vsC1] + runif(length(tied_values_S1vsC1), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupD1vsC1[tied_values_D1vsC1] <- heatwave.DEGs.gsea_groupD1vsC1[tied_values_D1vsC1] + runif(length(tied_values_D1vsC1), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupD1vsS1[tied_values_D1vsS1] <- heatwave.DEGs.gsea_groupD1vsS1[tied_values_D1vsS1] + runif(length(tied_values_D1vsS1), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupS2vsC2[tied_values_S2vsC2] <- heatwave.DEGs.gsea_groupS2vsC2[tied_values_S2vsC2] + runif(length(tied_values_S2vsC2), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupD2vsC2[tied_values_D2vsC2] <- heatwave.DEGs.gsea_groupD2vsC2[tied_values_D2vsC2] + runif(length(tied_values_D2vsC2), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupD2vsS2[tied_values_D2vsS2] <- heatwave.DEGs.gsea_groupD2vsS2[tied_values_D2vsS2] + runif(length(tied_values_D2vsS2), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupC1vsC2[tied_values_C1vsC2] <- heatwave.DEGs.gsea_groupC1vsC2[tied_values_C1vsC2] + runif(length(tied_values_C1vsC2), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupS1vsS2[tied_values_S1vsS2] <- heatwave.DEGs.gsea_groupS1vsS2[tied_values_S1vsS2] + runif(length(tied_values_S1vsS2), min = 0, max = 0.001)
heatwave.DEGs.gsea_groupD1vsD2[tied_values_D1vsD2] <- heatwave.DEGs.gsea_groupD1vsD2[tied_values_D1vsD2] + runif(length(tied_values_D1vsD2), min = 0, max = 0.001)

# Prepare GO list using split_go_data
heatwave_GO.list <- split(split_go_data$gene, split_go_data$go)
head(heatwave_GO.list)
library(stats)

library(GO.db)
library(AnnotationDbi)

# Example: Your GO terms from heatwave_GO.list
go_ids <- names(heatwave_GO.list)  # Extract GO terms from your list

# Query GO.db for only these GO terms
go_ontology_filtered <- AnnotationDbi::select(GO.db, 
                                              keys = go_ids, 
                                              columns = "ONTOLOGY", 
                                              keytype = "GOID")

# View results
head(go_ontology_filtered)

# filter for BP
bp_go_ids <- go_ontology_filtered %>%
  filter(ONTOLOGY == "BP") %>%
  pull(GOID)  # Extract only BP GO IDs

# Filter heatwave_GO.list to keep only BP terms
GO_BP.list <- heatwave_GO.list[names(heatwave_GO.list) %in% bp_go_ids]

# Check structure
head(GO_BP.list)
str(GO_BP.list)

# Count the number of GO terms in the original full list
num_go_terms_all <- length(heatwave_GO.list)

# Count the number of GO terms after filtering for Biological Process (BP)
num_go_terms_bp <- length(GO_BP.list)

# Print the results
cat("Total GO terms in heatwave_GO.list:", num_go_terms_all, "\n")
cat("Total GO terms in GO_BP.list (Biological Process only):", num_go_terms_bp, "\n")

# Run GSEA analysis for each contrast
# 10/02 filtered for BP using GO_BP.list now instead of heatwave_GO.list
fgseaRes_S1vsC1 <- fgsea::fgsea(
  pathways = GO_BP.list, 
  stats = heatwave.DEGs.gsea_groupS1vsC1,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_D1vsC1 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupD1vsC1,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_D1vsS1 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupD1vsS1,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_S2vsC2 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupS2vsC2,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_D2vsC2 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupD2vsC2,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_D2vsS2 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupD2vsS2,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_C1vsC2 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupC1vsC2,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_S1vsS2 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupS1vsS2,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

fgseaRes_D1vsD2 <- fgsea::fgsea(
  pathways = heatwave_GO.list, 
  stats = heatwave.DEGs.gsea_groupD1vsD2,
  eps = 0.0, 
  minSize = 3, 
  maxSize = 500,
  nproc = 1
) %>%
  as_tibble() %>%
  arrange(padj)

# Example summaries
summary(fgseaRes_S1vsC1$pval)
summary(fgseaRes_S1vsC1$padj)
head(fgseaRes_S1vsC1)
head(heatwave.DEGs.gsea_groupS1vsC1)
summary(heatwave.DEGs.gsea_groupS1vsC1)

summary(fgseaRes_D1vsC1$pval)
summary(fgseaRes_D1vsC1$padj)
head(fgseaRes_D1vsC1)
head(heatwave.DEGs.gsea_groupD1vsC1)
summary(heatwave.DEGs.gsea_groupD1vsC1)

heatwave_GO_signif_down <- heatwave_GO_signif_annot
heatwave_GO_signif_down2 <- heatwave_GO_signif_down %>%
  mutate(direction2 = 'down') %>%
  filter(!(term %>%str_detect("brain|heart|blood|chemokine|lipid")))

heatwave_GO_signif_up <- heatwave_GO_signif_annot
heatwave_GO_signif_up2 <- heatwave_GO_signif_up %>%
  mutate(direction2 = 'up') %>%
  filter(!(term %>%str_detect('angiogenesis|blood|costamere|node|chemokine|lipid')))

### Identify significant DEGs associated with enriched GO terms (focus on significant GO terms with more than 5 DEGs)
heatwave_GO_signif_DEGs <- rbind(heatwave_GO_signif_up2, heatwave_GO_signif_down2) %>%
  filter(direction == 'enriched') %>%
  dplyr::select(category, term, ontology, contrast, gene) %>%
  separate_rows(gene, sep = '; ')

# Create a lookup table for GO terms and their names
go_terms_lookup <- heatwave_GO_signif_DEGs %>%
  distinct(category, term) %>%
  rename(GO_ID = category, Term_Name = term)

# Add GO Term Names to GSEA Results
fgseaRes_S1vsC1 <- fgseaRes_S1vsC1 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_D1vsC1 <- fgseaRes_D1vsC1 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_D1vsS1 <- fgseaRes_D1vsS1 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_S2vsC2 <- fgseaRes_S2vsC2 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_D2vsC2 <- fgseaRes_D2vsC2 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_D2vsS2 <- fgseaRes_D2vsS2 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_C1vsC2 <- fgseaRes_C1vsC2 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_S1vsS2 <- fgseaRes_S1vsS2 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

fgseaRes_D1vsD2 <- fgseaRes_D1vsD2 %>%
  left_join(go_terms_lookup, by = c("pathway" = "GO_ID"))

# Removing obsolete GO terms
fgseaRes_S1vsC1 <- fgseaRes_S1vsC1 %>% filter(!is.na(Term_Name))
fgseaRes_D1vsC1 <- fgseaRes_D1vsC1 %>% filter(!is.na(Term_Name))
fgseaRes_D1vsS1 <- fgseaRes_D1vsS1 %>% filter(!is.na(Term_Name))
fgseaRes_S2vsC2 <- fgseaRes_S2vsC2 %>% filter(!is.na(Term_Name))
fgseaRes_D2vsC2 <- fgseaRes_D2vsC2 %>% filter(!is.na(Term_Name))
fgseaRes_D2vsS2 <- fgseaRes_D2vsS2 %>% filter(!is.na(Term_Name))
fgseaRes_C1vsC2 <- fgseaRes_C1vsC2 %>% filter(!is.na(Term_Name))
fgseaRes_S1vsS2 <- fgseaRes_S1vsS2 %>% filter(!is.na(Term_Name))
fgseaRes_D1vsD2 <- fgseaRes_D1vsD2 %>% filter(!is.na(Term_Name))

# Print the heads of the updated objects
head(fgseaRes_S1vsC1)
head(fgseaRes_D1vsC1)
head(fgseaRes_D1vsS1)
head(fgseaRes_S2vsC2)
head(fgseaRes_D2vsC2)
head(fgseaRes_D2vsS2)
head(fgseaRes_C1vsC2)
head(fgseaRes_S1vsS2)
head(fgseaRes_D1vsD2)

# Removing obsolete GO terms
fgseaRes_S1vsC1 <- fgseaRes_S1vsC1 %>% dplyr::select(-Term_Name.x, -Term_Name.y) %>% filter(!is.na(Term_Name))
fgseaRes_D1vsC1 <- fgseaRes_D1vsC1 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_D1vsS1 <- fgseaRes_D1vsS1 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_S2vsC2 <- fgseaRes_S2vsC2 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_D2vsC2 <- fgseaRes_D2vsC2 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_D2vsS2 <- fgseaRes_D2vsS2 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_C1vsC2 <- fgseaRes_C1vsC2 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_S1vsS2 <- fgseaRes_S1vsS2 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))
fgseaRes_D1vsD2 <- fgseaRes_D1vsD2 %>% dplyr::select(-Term_Name.x, -Term_Name.y)# %>% filter(!is.na(Term_Name))


## new bit of code - REVIGO
revigo_gsea_input_S1vsS2 <- fgseaRes_S1vsS2 %>% 
  dplyr::select(pathway, padj) %>% 
  arrange(padj)  # Sort by significance

revigo_gsea_input_D1vsD2 <- fgseaRes_D1vsD2 %>% 
  dplyr::select(pathway, padj) %>% 
  arrange(padj)  # Sort by significance

# Save for REVIGO
write.table(revigo_gsea_input_S1vsS2, "revigo_gsea_S1vsS2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(revigo_gsea_input_D1vsD2, "revigo_gsea_D1vsD2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

revigo_results_S1vsS2 <- read.table("revigo_results_S1vsS2.tsv", sep = "\t", header = TRUE)
revigo_results_D1vsD2 <- read.table("revigo_results_D1vsD2.tsv", sep = "\t", header = TRUE)

# merge revigo results and gsea results
filtered_gsea_S1vsS2 <- fgseaRes_S1vsS2 %>% 
  filter(pathway %in% revigo_results_S1vsS2$TermID)

filtered_gsea_D1vsD2 <- fgseaRes_D1vsD2 %>% 
  filter(pathway %in% revigo_results_D1vsD2$TermID)

## plot
library(ggplot2)

ggplot(filtered_gsea_S1vsS2, aes(x = NES, y = reorder(Term_Name, NES), size = size, color = -log10(padj))) +
  geom_point() +
  labs(
    title = "GSEA Results (S1 vs S2) - Non-Redundant GO Terms",
    x = "Normalized Enrichment Score (NES)",
    y = "GO Term",
    color = "-log10(padj)",
    size = "Gene Set Size"
  ) +
  theme_minimal()
#### END OF NEW CODE


library(GO.db)

# Combine and print results
gsea_results <- list(
  Contrast_S1vsC1 = fgseaRes_S1vsC1,
  Contrast_D1vsC1 = fgseaRes_D1vsC1,
  Contrast_D1vsS1 = fgseaRes_D1vsS1,
  Contrast_S2vsC2 = fgseaRes_S2vsC2,
  Contrast_D2vsC2 = fgseaRes_D2vsC2,
  Contrast_D2vsS2 = fgseaRes_D2vsS2,
  Contrast_C1vsC2 = fgseaRes_C1vsC2,
  Contrast_S1vsS2 = fgseaRes_S1vsS2,
  Contrast_D1vsD2 = fgseaRes_D1vsD2
)


# Prepare GO list using split_go_data
heatwave_GO.list <- split(split_go_data$gene, split_go_data$go)
head(heatwave_GO.list)
library(stats)

# Define the keywords to remove
unwanted_terms <- c("interleukin", "neuron", "neural", "cytokine", "gastrulation", "synaptic", "response to starvation", "central nervous system",
                    "synapse", "pancreas", "pollen", "floral", "bone", "exocrine", "placenta", "dendrite", "dendritic", 
                    "eye", "odorant", "heart", "neurotransmitter", "blood", "pain", "glial", "endothelial", "behavior", "locomotor", "insulin",
                    "T cell", "L-leucine", "myotube", "oligodendrocyte", "muscle", "cardiac")

unwanted_terms <- c("interleukin", "neuron", "neural", "gastrulation", 
  "synaptic", "adhesion", "response to starvation", "central nervous system", 
  "synapse", "pancreas", "pollen", "floral", "bone", "exocrine", 
  "placenta", "dendrite", "dendritic", "eye", "odorant", "heart", 
  "blood", "pain", "glial", "endothelial", "behavior", "locomotor", 
  "insulin", "tRNA", "T cell", "L-leucine", "myotube", "oligodendrocyte",
  "macrolide", "SMAD", "Notch", "chromosome", "ecdysteroid", "innate",
  "transesterification", "organ", "vitamin", "post-embryonic", "dosage",
  "viral", "epithelial", "canonical Wnt", "reproductive", "immune", "integration",
  "development", "estrogen", "axon", "virus", "symbiont-mediated",
  "mucopolysaccharide", "spinal", "water", "spermatogenesis", "male", "female", "bile",
  "osteoblast", "wnt", "body", "mechanical", "body", "fluid", "sound",
  "neurotransmitter", "motile", "nucleophagy", "microtubule", "intraciliary",
  "insect", "face", "muscle", "cardiac", "node", "stem cell", "alcohol", "azole",
  "sulfatase", "ciliary", "N-acetylgalactosamine-4-sulfatase", "brush", "azurophil",
  "V-type", "snoRNA", "neurogenesis", "telomere")


# Function to filter out unwanted terms
remove_unwanted_terms <- function(data, unwanted_terms) {
  data %>%
    dplyr::filter(!str_detect(Term_Name, paste(unwanted_terms, collapse = "|")))
}

# Apply to S1 vs C1 contrast
fgseaRes_S1vsC1 <- remove_unwanted_terms(fgseaRes_S1vsC1, unwanted_terms)

# Group terms for S1 vs C1
fgseaRes_S1vsC1_grouped <- fgseaRes_S1vsC1 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)lipid|carbohydrate|metabolic|polysaccharide|metabolism|metabolism|acyltransferase") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)methylation|phosphorylation") ~ 'Post-translational modification',
    str_detect(Term_Name, "(?i)localization|transmembrane|ubiquitination") ~ 'Protein processing and transport',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|translation|recombination|amino|transcription|nucleotide") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|acidification|proton|ions") ~ 'pH and Ion Homeostasis',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_S1vsC1_grouped <- fgseaRes_S1vsC1_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_S1vsC1_grouped$group))

# Filter significant terms for S1 vs C1
fgseaRes_S1vsC1_signif <- fgseaRes_S1vsC1_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_S1vsC1 <- fgseaRes_S1vsC1_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_S1vsC1 <- top_terms_S1vsC1 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_S1vsC1)

# Inspect the grouped terms in top_terms_S1vsC1
head(top_terms_S1vsC1)

# Filter terms in the group "Other processes"
other_processes_terms_S1vsC1 <- top_terms_S1vsC1 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_S1vsC1)

# If you want to see the count of terms in "Other processes"
n_other_processes_S1vsC1 <- nrow(other_processes_terms_S1vsC1)
cat("Number of terms in 'Other processes':", n_other_processes_S1vsC1, "\n")

# Reorder groups for better organization
top_terms_S1vsC1$group <- factor(top_terms_S1vsC1$group, levels = c(
  'Oxidative stress markers',
  'Post-translational modification',
  'Protein processing and transport',
  'Stress and repair mechanisms',
  'Molecular chaperones',
  'Growth processes',
  'Nitrogen transport',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for S1 vs C1
gsea_enriched_S1vsC1 <- ggplot(top_terms_S1vsC1, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Oxidative stress markers' = "#E69F00",
    'Post-translational modification' = "#56B4E9",
    'Protein processing and transport' = "#009E73",
    'Growth processes' = "#0072B2",
    'Stress and repair mechanisms' = "#D55E00",
    'Molecular chaperones' = "#CC79A7",
    'Nitrogen transport' = "#FC8D62",
    'DNA and RNA processing' = "#8DA0CB",
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for S1 vs C1"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_S1vsC1 <- gsea_enriched_S1vsC1
Fig_S1vsC1

ggsave(plot = Fig_S1vsC1, filename = "Fig_S1vsC1.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_S1vsC1, filename = "Fig_S1vsC1.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")

#### D1 vs C1
# Apply to D1 vs C1 contrast
fgseaRes_D1vsC1 <- remove_unwanted_terms(fgseaRes_D1vsC1, unwanted_terms)

# Group terms for D1 vs C1
fgseaRes_D1vsC1_grouped <- fgseaRes_D1vsC1 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide reductase|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)lipid|carbohydrate|metabolic|thiamine|polysaccharide|metabolism|amino acid metabolism|acyltransferase activity|carboxyl- or carbamoyltransferase activity|transaminase activity") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)apoptosis|apoptotic|death") ~ 'Apoptosis processes',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|gene|translation|recombination|amino|transcription|nucleotide") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)lipase") ~ 'Lipid metabolism',
    str_detect(Term_Name, "(?i)endoplasmic reticulum|ERAD") ~ 'Protein processing',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|acidification|proton|ions|chloride") ~ 'pH and Ion Homeostasis',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_D1vsC1_grouped <- fgseaRes_D1vsC1_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_D1vsC1_grouped$group))

# Filter significant terms for D1 vs C1
fgseaRes_D1vsC1_signif <- fgseaRes_D1vsC1_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_D1vsC1 <- fgseaRes_D1vsC1_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_D1vsC1 <- top_terms_D1vsC1 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_D1vsC1)

# Inspect the grouped terms in top_terms_D1vsC1
head(top_terms_D1vsC1)

# Filter terms in the group "Other processes"
other_processes_terms_D1vsC1 <- top_terms_D1vsC1 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_D1vsC1)

# If you want to see the count of terms in "Other processes"
n_other_processes_D1vsC1 <- nrow(other_processes_terms_D1vsC1)
cat("Number of terms in 'Other processes':", n_other_processes_D1vsC1, "\n")

# Reorder groups for better organization
top_terms_D1vsC1$group <- factor(top_terms_D1vsC1$group, levels = c(
  'Oxidative stress markers',
  'Lipid metabolism',
  'Protein processing',
  'Photosynthesis',
  'Stress and repair mechanisms',
  'Molecular chaperones',
  'Growth processes',
  'Nitrogen transport',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for D1 vs C1
gsea_enriched_D1vsC1 <- ggplot(top_terms_D1vsC1, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Lipid metabolism' = "#E69F00",
    'DNA and RNA processing' = "#56B4E9",
    'Growth processes' = "#009E73",
    'Photosynthesis' = "#0072B2",
    'Stress and repair mechanisms' = "#D55E00",
    'Molecular chaperones' = "#CC79A7",
    'Nitrogen transport' = "#FC8D62",
    'Protein processing' = "#8DA0CB",
    'pH and Ion Homeostasis' = "#CCCC43"
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for D1 vs C1"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_D1vsC1 <- gsea_enriched_D1vsC1
Fig_D1vsC1

ggsave(plot = Fig_D1vsC1, filename = "Fig_D1vsC1.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_D1vsC1, filename = "Fig_D1vsC1.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")

#S2 vs C2
# Apply to S2 vs C2 contrast
fgseaRes_S2vsC2 <- remove_unwanted_terms(fgseaRes_S2vsC2, unwanted_terms)

# Group terms for S2 vs C2
fgseaRes_S2vsC2_grouped <- fgseaRes_S2vsC2 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide reductase|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)lipid|carbohydrate|metabolic|thiamine|polysaccharide|metabolism|amino acid metabolism|acyltransferase activity|carboxyl- or carbamoyltransferase activity|transaminase activity") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)apoptosis|apoptotic|death") ~ 'Apoptosis processes',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|gene|translation|recombination|amino|transcription|nucleotide") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|electron|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)lipase") ~ 'Lipid metabolism',
    str_detect(Term_Name, "(?i)endoplasmic reticulum|ERAD|endopeptidase|peptidase") ~ 'Protein processing',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|O-methyltransferase|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|pH|ion|acidification|proton|chloride") ~ 'pH and Ion Homeostasis',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_S2vsC2_grouped <- fgseaRes_S2vsC2_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_S2vsC2_grouped$group))

# Filter significant terms for S2 vs C2
fgseaRes_S2vsC2_signif <- fgseaRes_S2vsC2_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_S2vsC2 <- fgseaRes_S2vsC2_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_S2vsC2 <- top_terms_S2vsC2 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_S2vsC2)

# Inspect the grouped terms in top_terms_S2vsC2
head(top_terms_S2vsC2)

# Filter terms in the group "Other processes"
other_processes_terms_S2vsC2 <- top_terms_S2vsC2 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_S2vsC2)

# If you want to see the count of terms in "Other processes"
n_other_processes_S2vsC2 <- nrow(other_processes_terms_S2vsC2)
cat("Number of terms in 'Other processes':", n_other_processes_S2vsC2, "\n")

# Reorder groups for better organization
top_terms_S2vsC2$group <- factor(top_terms_S2vsC2$group, levels = c(
  'Oxidative stress markers',
  'Lipid metabolism',
  'Protein processing',
  'Photosynthesis',
  'Stress and repair mechanisms',
  'Molecular chaperones',
  'Growth processes',
  'Nitrogen transport',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for S2 vs C2
gsea_enriched_S2vsC2 <- ggplot(top_terms_S2vsC2, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Growth processes' = "#E69F00",
    'pH and Ion Homeostasis' = "#56B4E9",
    'DNA and RNA processing' = "#009E73",
    'Photosynthesis' = "#0072B2",
    'Stress and repair mechanisms' = "#D55E00",
    'Molecular chaperones' = "#CC79A7",
    'Nitrogen transport' = "#FC8D62",
    'Protein processing' = "#8DA0CB",
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for S2 vs C2"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_S2vsC2 <- gsea_enriched_S2vsC2
Fig_S2vsC2

ggsave(plot = Fig_S2vsC2, filename = "Fig_S2vsC2.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_S2vsC2, filename = "Fig_S2vsC2.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")


### D2 vs C2
# Apply to D2 vs C2 contrast
fgseaRes_D2vsC2 <- remove_unwanted_terms(fgseaRes_D2vsC2, unwanted_terms)

# Group terms for D2 vs C2
fgseaRes_D2vsC2_grouped <- fgseaRes_D2vsC2 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide reductase|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)lipid|carbohydrate|thiamine|polysaccharide|metabolism|amino acid metabolism|acyltransferase activity|carboxyl- or carbamoyltransferase activity|transaminase activity") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)apoptosis|apoptotic|death") ~ 'Apoptosis processes',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|gene|translation|recombination|amino|transcription|nucleotide") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)lipase") ~ 'Lipid metabolism',
    str_detect(Term_Name, "(?i)endoplasmic reticulum|ERAD|endopeptidase|peptidase") ~ 'Protein processing',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|O-methyltransferase|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|pH|ion|acidification|proton|chloride") ~ 'pH and Ion Homeostasis',
    str_detect(Term_Name, "(?i)iron-sulfur|cluster|electron") ~ 'Metal Cofactors & Redox Reactions',
    str_detect(Term_Name, "(?i)carotenoid|tetraterpenoid|pigment|molecule") ~ 'Pigment and Secondary Metabolism',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_D2vsC2_grouped <- fgseaRes_D2vsC2_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_D2vsC2_grouped$group))

# Filter significant terms for D2 vs C2
fgseaRes_D2vsC2_signif <- fgseaRes_D2vsC2_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_D2vsC2 <- fgseaRes_D2vsC2_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_D2vsC2 <- top_terms_D2vsC2 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_D2vsC2)

# Inspect the grouped terms in top_terms_D2vsC2
head(top_terms_D2vsC2)

# Filter terms in the group "Other processes"
other_processes_terms_D2vsC2 <- top_terms_D2vsC2 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_D2vsC2)

# If you want to see the count of terms in "Other processes"
n_other_processes_D2vsC2 <- nrow(other_processes_terms_D2vsC2)
cat("Number of terms in 'Other processes':", n_other_processes_D2vsC2, "\n")

# Reorder groups for better organization
top_terms_D2vsC2$group <- factor(top_terms_D2vsC2$group, levels = c(
  'Oxidative stress markers',
  'Lipid metabolism',
  'Protein processing',
  'Photosynthesis',
  'Pigment and Secondary Metabolism',
  'Metal Cofactors & Redox Reactions',
  'Growth processes',
  'Nitrogen transport',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for D2 vs C2
gsea_enriched_D2vsC2 <- ggplot(top_terms_D2vsC2, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Oxidative stress markers' = "#E69F00",
    'Pigment and Secondary Metabolism' = "#56B4E9",
    'Metal Cofactors & Redox Reactions' = "#009E73",
    'Photosynthesis' = "#0072B2",
    'Growth processes' = "#D55E00",
    'DNA and RNA processing' = "#CC79A7",
    'Nitrogen transport' = "#FC8D62",
    'Protein processing' = "#8DA0CB",
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for D2 vs C2"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_D2vsC2 <- gsea_enriched_D2vsC2
Fig_D2vsC2

ggsave(plot = Fig_D2vsC2, filename = "Fig_D2vsC2.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_D2vsC2, filename = "Fig_D2vsC2.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")

### D2 vs S2
# Apply to D2 vs S2 contrast
fgseaRes_D2vsS2 <- remove_unwanted_terms(fgseaRes_D2vsS2, unwanted_terms)

# Group terms for D2 vs S2
fgseaRes_D2vsS2_grouped <- fgseaRes_D2vsS2 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide reductase|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)lipid|carbohydrate|catabolic|metabolic|thiamine|polysaccharide|metabolism|amino acid metabolism|acyltransferase activity|carboxyl- or carbamoyltransferase activity|transaminase activity") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)apoptosis|apoptotic|death") ~ 'Apoptosis processes',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|gene|translation|recombination|amino|transcription|nucleotide") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)lipase") ~ 'Lipid metabolism',
    str_detect(Term_Name, "(?i)endoplasmic reticulum|ERAD|endopeptidase|peptidase|cysteine") ~ 'Protein processing',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|O-methyltransferase|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|pH|ion|proton|chloride") ~ 'pH and Ion Homeostasis',
    str_detect(Term_Name, "(?i)iron-sulfur|cluster|electron|oxidoreductase") ~ 'Metal Cofactors & Redox Reactions',
    str_detect(Term_Name, "(?i)carotenoid|plastoglobule|tetraterpenoid|pigment") ~ 'Pigment and Secondary Metabolism',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_D2vsS2_grouped <- fgseaRes_D2vsS2_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_D2vsS2_grouped$group))

# Filter significant terms for D2 vs S2
fgseaRes_D2vsS2_signif <- fgseaRes_D2vsS2_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_D2vsS2 <- fgseaRes_D2vsS2_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_D2vsS2 <- top_terms_D2vsS2 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_D2vsS2)

# Inspect the grouped terms in top_terms_D2vsS2
head(top_terms_D2vsS2)

# Filter terms in the group "Other processes"
other_processes_terms_D2vsS2 <- top_terms_D2vsS2 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_D2vsS2)

# If you want to see the count of terms in "Other processes"
n_other_processes_D2vsS2 <- nrow(other_processes_terms_D2vsS2)
cat("Number of terms in 'Other processes':", n_other_processes_D2vsS2, "\n")

# Reorder groups for better organization
top_terms_D2vsS2$group <- factor(top_terms_D2vsS2$group, levels = c(
  'Oxidative stress markers',
  'Lipid metabolism',
  'Protein processing',
  'Photosynthesis',
  'Pigment and Secondary Metabolism',
  'Metal Cofactors & Redox Reactions',
  'Growth processes',
  'Metabolism',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for D2 vs S2
gsea_enriched_D2vsS2 <- ggplot(top_terms_D2vsS2, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Oxidative stress markers' = "#E69F00",
    'Metabolism' = "#56B4E9",
    'Pigment and Secondary Metabolism' = "#009E73",
    'Photosynthesis' = "#0072B2",
    'Metal Cofactors & Redox Reactions' = "#D55E00",
    'Growth Processes' = "#CC79A7",
    'DNA and RNA processing' = "#FC8D62",
    'Protein processing' = "#8DA0CB",
    'pH and Ion Homeostasis' = "#CCCC56",
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for D2 vs S2"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_D2vsS2 <- gsea_enriched_D2vsS2
Fig_D2vsS2

ggsave(plot = Fig_D2vsS2, filename = "Fig_D2vsS2.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_D2vsS2, filename = "Fig_D2vsS2.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")

#S1 vs S2
# Apply to S1 vs S2 contrast
fgseaRes_S1vsS2 <- remove_unwanted_terms(fgseaRes_S1vsS2, unwanted_terms)

# Group terms for S1 vs S2
fgseaRes_S1vsS2_grouped <- fgseaRes_S1vsS2 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide reductase|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)carbohydrate|catabolic|metabolic|thiamine|polysaccharide|metabolism|amino acid metabolism|transaminase activity") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)apoptosis|apoptotic|death") ~ 'Apoptosis processes',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)CoA") ~ 'Lipid metabolism',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|gene|translation|recombination|amino|transcription|nucleotide") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)dehydrogenase|ligase|transferase|hydrolase") ~ 'Enzymes',
    str_detect(Term_Name, "(?i)endoplasmic reticulum|ERAD|endopeptidase|peptidase|cysteine|protein") ~ 'Protein processing',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|O-methyltransferase|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|pH|ion|proton|chloride") ~ 'pH and Ion Homeostasis',
    str_detect(Term_Name, "(?i)iron-sulfur|cluster|electron|oxidoreductase") ~ 'Metal Cofactors & Redox Reactions',
    str_detect(Term_Name, "(?i)carotenoid|plastoglobule|tetraterpenoid|pigment") ~ 'Pigment and Secondary Metabolism',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_S1vsS2_grouped <- fgseaRes_S1vsS2_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_S1vsS2_grouped$group))

# Filter significant terms for S1 vs S2
fgseaRes_S1vsS2_signif <- fgseaRes_S1vsS2_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_S1vsS2 <- fgseaRes_S1vsS2_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_S1vsS2 <- top_terms_S1vsS2 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_S1vsS2)

# Inspect the grouped terms in top_terms_S1vsS2
head(top_terms_S1vsS2)

# Filter terms in the group "Other processes"
other_processes_terms_S1vsS2 <- top_terms_S1vsS2 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_S1vsS2)

# If you want to see the count of terms in "Other processes"
n_other_processes_S1vsS2 <- nrow(other_processes_terms_S1vsS2)
cat("Number of terms in 'Other processes':", n_other_processes_S1vsS2, "\n")

# Reorder groups for better organization
top_terms_S1vsS2$group <- factor(top_terms_S1vsS2$group, levels = c(
  'Molecular chaperones',
  'Enzymes',
  'Protein processing',
  'Photosynthesis',
  'Lipid metabolism',
  'Metal Cofactors & Redox Reactions',
  'Growth processes',
  'Metabolism',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for S1 vs S2
gsea_enriched_S1vsS2 <- ggplot(top_terms_S1vsS2, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Enzymes' = "#E69F00",
    'Metabolism' = "#56B4E9",
    'Lipid metabolism' = "#009E73",
    'Photosynthesis' = "#0072B2",
    'Metal Cofactors & Redox Reactions' = "#D55E00",
    'Molecular chaperones' = "#CC79A7",
    'Growth processes' = "#FC8D62",
    'Protein processing' = "#8DA0CB",
    'DNA and RNA processing' = "#CCCC45",
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for S1 vs S2"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_S1vsS2 <- gsea_enriched_S1vsS2
Fig_S1vsS2

ggsave(plot = Fig_S1vsS2, filename = "Fig_S1vsS2.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_S1vsS2, filename = "Fig_S1vsS2.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")

## D1 vs D2
# Apply to D1 vs D2 contrast
fgseaRes_D1vsD2 <- remove_unwanted_terms(fgseaRes_D1vsD2, unwanted_terms)

# Group terms for D1 vs D2
fgseaRes_D1vsD2_grouped <- fgseaRes_D1vsD2 %>%
  mutate(group = case_when(
    str_detect(Term_Name, "(?i)ascorbate|catalase|glutathione|peroxidase|peroxiredoxin|detoxification|superoxide|s-oxide reductase|methionine sulfoxide reductase|antioxidant") ~ 'Oxidative stress markers',
    str_detect(Term_Name, "(?i)carbohydrate|catabolic|metabolic|thiamine|polysaccharide|metabolism|amino acid metabolism|transaminase activity") ~ 'Metabolism',
    str_detect(Term_Name, "(?i)apoptosis|apoptotic|death") ~ 'Apoptosis processes',
    str_detect(Term_Name, "(?i)TOR|cell cycle") ~ 'Growth processes',
    str_detect(Term_Name, "(?i)CoA") ~ 'Lipid metabolism',
    str_detect(Term_Name, "(?i)DNA|RNA|mRNA|gene|translation|recombination|amino|transcription|nucleotide|spliceosome") ~ 'DNA and RNA processing',
    str_detect(Term_Name, "(?i)photosynthesis|photosynthetic|photosystem|chlorophyll|chloroplast|light|thylakoid|plastid|porphyrin|tetrapyrrole") ~ 'Photosynthesis',
    str_detect(Term_Name, "(?i)chaperone|dnaj|heat shock|hsp70|hsp90") ~ 'Molecular chaperones',
    str_detect(Term_Name, "(?i)ammonium|nitrate|nitrite|nitrogen") ~ 'Nitrogen transport',
    str_detect(Term_Name, "(?i)fructose|gapdh|glyceraldehyde|glycolytic|glucose|phosphoglycerate") ~ 'Glycolysis',
    str_detect(Term_Name, "(?i)dehydrogenase|ligase|transferase|hydrolase") ~ 'Enzymes',
    str_detect(Term_Name, "(?i)endoplasmic reticulum|ERAD|endopeptidase|peptidase|cysteine|protein|preribosome") ~ 'Protein processing',
    str_detect(Term_Name, "(?i)proteolysis|damage|repair|O-methyltransferase|deubiquitination|stress|repair|peroxisome|unfolded|refolding|refolded|folding|autophagy|autophagosome|microautophagy|abscisic acid") ~ 'Stress and repair mechanisms',
    str_detect(Term_Name, "(?i)proton-transporting|pH|ion|proton|chloride") ~ 'pH and Ion Homeostasis',
    str_detect(Term_Name, "(?i)iron-sulfur|cluster|electron|oxidoreductase") ~ 'Metal Cofactors & Redox Reactions',
    str_detect(Term_Name, "(?i)carotenoid|plastoglobule|tetraterpenoid|pigment") ~ 'Pigment and Secondary Metabolism',
    TRUE ~ 'Other processes'
  ))

# Ensure no "NA" group appears
fgseaRes_D1vsD2_grouped <- fgseaRes_D1vsD2_grouped %>%
  mutate(group = replace_na(group, "Other processes"))

# Check for NA values
sum(is.na(fgseaRes_D1vsD2_grouped$group))

# Filter significant terms for D1 vs D2
fgseaRes_D1vsD2_signif <- fgseaRes_D1vsD2_grouped %>% filter(padj < 0.05)

# Filter to include top 30 NES scores
top_terms_D1vsD2 <- fgseaRes_D1vsD2_grouped %>%
  filter(!is.na(group)) %>%
  filter(group != "NA") %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

# Count the occurrences of each group in the filtered top terms
group_summary_D1vsD2 <- top_terms_D1vsD2 %>%
  dplyr::count(group, sort = TRUE)

# View the grouped summary to understand group representation
print(group_summary_D1vsD2)

# Inspect the grouped terms in top_terms_D1vsD2
head(top_terms_D1vsD2)

# Filter terms in the group "Other processes"
other_processes_terms_D1vsD2 <- top_terms_D1vsD2 %>%
  filter(group == "Other processes") %>%
  dplyr::select(Term_Name, NES, padj)

# View the terms in "Other processes"
print(other_processes_terms_D1vsD2)

# If you want to see the count of terms in "Other processes"
n_other_processes_D1vsD2 <- nrow(other_processes_terms_D1vsD2)
cat("Number of terms in 'Other processes':", n_other_processes_D1vsD2, "\n")

# Reorder groups for better organization
top_terms_D1vsD2$group <- factor(top_terms_D1vsD2$group, levels = c(
  'Molecular chaperones',
  'Enzymes',
  'Protein processing',
  'Photosynthesis',
  'Lipid metabolism',
  'Metal Cofactors & Redox Reactions',
  'Growth processes',
  'Metabolism',
  'pH and Ion Homeostasis',
  'DNA and RNA processing',
  'Other processes'
))

# Dot plot for D1 vs D2
gsea_enriched_D1vsD2 <- ggplot(top_terms_D1vsD2, aes(x = NES, y = reorder(Term_Name, NES), size = -log10(pval), color = group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    'Enzymes' = "#E69F00",
    'Metabolism' = "#56B4E9",
    'Growth processes' = "#009E73",
    'Photosynthesis' = "#0072B2",
    'DNA and RNA processing' = "#D55E00",
    'Molecular chaperones' = "#CC79A7",
    'pH and Ion Homeostasis' = "#FC8D62",
    'Protein processing' = "#8DA0CB",
    'Other processes' = "gray"
  )) +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x = "Normalized Enrichment Score (NES)", 
    y = "GO Term", 
    size = "-log10(p-value)",
    title = "Enriched GO Terms for D1 vs D2"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), 
        axis.text.y = element_text(size = 14),  # **Increase Y-axis text size**
        axis.text.x = element_text(size = 16),  # **Increase X-axis text size**
        legend.position = "top",  
        legend.text = element_text(size = 14),  # **Increase legend text size**
        legend.title = element_text(size = 16, face = "bold"),  # **Increase legend title size and bold**
        strip.text.y = element_text(size = 16, face = "bold", angle = 0),  # **Increase facet labels & bold**
        panel.grid.major.x = element_blank(),
        panel.spacing.y = unit(0.5, "lines")
  ) +
  guides(color = guide_legend(ncol = 2), size = guide_legend(override.aes = list(alpha = 1)))

# Save the plot
Fig_D1vsD2 <- gsea_enriched_D1vsD2
Fig_D1vsD2

ggsave(plot = Fig_D1vsD2, filename = "Fig_D1vsD2.dotplot.enriched.pdf",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "pdf")

ggsave(plot = Fig_D1vsD2, filename = "Fig_D1vsD2.dotplot.enriched.png",
       path = "./",
       width = 15,
       height = 17.5,
       units = "in",
       dpi = 400,
       device = "png")

# View the DEGs in D1 vs S1 contrast
ecklonia_DEGs_D1_vs_S1 <- heatwave.DEGs %>% filter(contrast == "group D1 vs S1") %>% filter(padj < 0.05)


# Print the selected DEGs
print(ecklonia_DEGs_D1_vs_S1)

# Filter Trinotate report for the two genes
trinotate_results <- trinotate_report %>%
  filter(gene_id %in% c("TRINITY_DN13575_c0_g1", "TRINITY_DN19337_c0_g1")) %>%
  dplyr::select(gene_id, EggNM.PFAMs)  # Adjust column names if needed

# Print results
print(trinotate_results)


trinotate_results1 <- trinotate_report %>%
  filter(gene_id %in% c("TRINITY_DN19337_c0_g1"))
# Print results
print(trinotate_results1)



# looking at top DEGs in some contrasts
# Filter DEGs for the desired contrast
unique_contrasts
heatwave_DEGs_D1_vs_D2 <- heatwave.DEGs_filtered %>%
  filter(contrast == "group D1 vs D2")

trinotate_results <- trinotate_report %>%
  dplyr::select(gene_id, EggNM.PFAMs, EggNM.Description)  # Adjust column names if needed

# Ensure column names match
colnames(trinotate_results)

# Merge the filtered DEGs with Trinotate annotations
heatwave_DEGs_D1_vs_D2_annotated <- heatwave_DEGs_D1_vs_D2 %>%
  left_join(trinotate_results, by = c("gene" = "gene_id"))

# View the first few rows
head(heatwave_DEGs_D1_vs_D2_annotated)

# Count how many genes have annotations
sum(!is.na(heatwave_DEGs_D1_vs_D2_annotated$EggNM.PFAMs))


heatwave_DEGs_D1_vs_D2_significant <- heatwave_DEGs_D1_vs_D2_annotated %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > 0.58)

nrow(heatwave_DEGs_D1_vs_D2_significant)

heatwave_DEGs_D1_vs_D2_significant_ranked <- heatwave_DEGs_D1_vs_D2_significant %>%
  arrange(padj)
head(heatwave_DEGs_D1_vs_D2_significant_ranked)
heatwave_DEGs_D1_vs_D2_significant_ranked <- heatwave_DEGs_D1_vs_D2_significant_ranked %>%
  filter(!is.na(EggNM.PFAMs))
head(heatwave_DEGs_D1_vs_D2_significant_ranked)
view(heatwave_DEGs_D1_vs_D2_significant_ranked)

# Split into upregulated and downregulated based on log2FoldChange
heatwave_DEGs_D1_vs_D2_upregulated <- heatwave_DEGs_D1_vs_D2_significant_ranked %>%
  filter(log2FoldChange > 0.58)

heatwave_DEGs_D1_vs_D2_downregulated <- heatwave_DEGs_D1_vs_D2_significant_ranked %>%
  filter(log2FoldChange < -0.58)

# Merge the filtered DEGs with Trinotate annotations
heatwave_DEGs_D1_vs_D2_upregulated <- heatwave_DEGs_D1_vs_D2_upregulated %>%
  left_join(trinotate_results, by = c("gene" = "gene_id"))

heatwave_DEGs_D1_vs_D2_downregulated <- heatwave_DEGs_D1_vs_D2_downregulated %>%
  left_join(trinotate_results, by = c("gene" = "gene_id"))

# View the first few rows of each
head(heatwave_DEGs_D1_vs_D2_upregulated)
head(heatwave_DEGs_D1_vs_D2_downregulated)
# Remove duplicate genes, keeping the first occurrence
heatwave_DEGs_D1_vs_D2_upregulated <- heatwave_DEGs_D1_vs_D2_upregulated %>%
  distinct(gene, .keep_all = TRUE)

# View first few rows to confirm
head(heatwave_DEGs_D1_vs_D2_upregulated)

# Remove duplicate genes, keeping the first occurrence
heatwave_DEGs_D1_vs_D2_downregulated <- heatwave_DEGs_D1_vs_D2_downregulated %>%
  distinct(gene, .keep_all = TRUE)

# View first few rows to confirm
head(heatwave_DEGs_D1_vs_D2_downregulated)

view(heatwave_DEGs_D1_vs_D2_downregulated)


# Optionally, save them as separate files
write.csv(heatwave_DEGs_D1_vs_D2_upregulated, "Upregulated_DEGs_D1_vs_D2.csv", row.names = FALSE)
write.csv(heatwave_DEGs_D1_vs_D2_downregulated, "Downregulated_DEGs_D1_vs_D2.csv", row.names = FALSE)

view(heatwave_DEGs_D1_vs_D2_upregulated)

