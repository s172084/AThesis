
library(GEOquery)
library(limma)
library(Biobase)
library(tidyverse)
library(viridis)
library(pheatmap)
library(ggpubr)
library(edgeR)
library(umap)
library(dplyr)
library(ggplot2)
library(factoextra)
library(estimate)
library(ggthemes)
library(broom)
library(cowplot)
library(ggrepel)
library(fs)
library(hrbrthemes)

# Soft file is first downloaded from : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3189
# Open a compressed file with the GEO Accession GSE3189 
# gse3189 <- getGEO(filename='Desktop/GSE3189_family.soft.gz')
# 
# #Look at the metadata 
# Meta(gse3189)
# 
# # Look at the overall design 
# Meta(gse3189)$overall_design
# 
# # Look at the sample IDs
# Meta(gse3189)$sample_id
# 
# # Names of all the GSM objects contained in the GSE
# names(GSMList(gse3189))
# 
# # Get the name of the GPL object. 
# names(GPLList(gse3189))

# -------------------------- Easier  ---------------------------------
# Getting GSE Series Matrix files as an ExpressionSet (based on Bioconductor instructions)

# Extract expression matrix
new_gse3189 <- getGEO('gse3189',GSEMatrix=TRUE)
show(new_gse3189)


# -------------------- Viewing Sample Information ----------------------
# number of samples 
ncol(new_gse3189[[1]]) # There are 70 samples. 

# Detailed sample information from the phenotype data.
pData(new_gse3189[[1]])
mel_pheno <-pData(new_gse3189[[1]])

# View specific columns of sample data
pData(new_gse3189[[1]])[1:40, c("title", "source_name_ch1", "characteristics_ch1")]

# -------------------- Viewing Gene Information ----------------------
# Number of genes/probes
nrow(new_gse3189[[1]])

# Gene/probe names
featureNames(new_gse3189[[1]])

# Feature/gene annotation data
fData(new_gse3189[[1]])

# View first few rows of feature data * this has the GO information inside of it. 
head(fData(new_gse3189[[1]]))
# -------------------- Viewing Expression Data Information ----------------------

# Get the expression matrix
expr_matrix <- exprs(new_gse3189[[1]])
expr_matrix

# Dimensions of expression matrix (genes x samples)
dim(expr_matrix)

# View first few rows and columns
expr_matrix[1:5, 1:5]

# Summary statistics
summary(expr_matrix)

# Basic info
cat("Number of genes:", nrow(new_gse3189[[1]]), "\n")
cat("Number of samples:", ncol(new_gse3189[[1]]), "\n")

# -------------------- Viewing Expression Data Information ----------------------

# Create boxplot of expression values
#boxplot(expr_matrix,
#        col = blues9, main = "Gene Expression Distribution")

# For density plot of a specific gene
df <- data.frame(gene = expr_matrix[1,])  # First gene across all samples

ggplot(df, mapping = aes(x = gene)) +
  geom_density() +
  labs(title= "Expression Distribution of First Gene")

# -------------------- Viewing Expression Data Information ----------------------

# Extract expression matrix and phenotype data
expr_matrix <- exprs(new_gse3189[[1]])

mel_pheno <- pData(new_gse3189[[1]])

# -------------------- Sample Selection ----------------------
# First, let's examine the phenotype data to understand sample types
print("Examining phenotype data:")
print(table(mel_pheno$characteristics_ch1))
print(table(mel_pheno$description))

# Function to identify sample types from both columns
identify_sample_type <- function(char_ch1, description) {
  # Convert to lowercase for easier matching
  char_lower <- tolower(char_ch1)
  desc_lower <- tolower(description)
  
  # Check for melanoma
  if(grepl("melanoma", char_lower) | grepl("melanoma", desc_lower)) {
    return("Melanoma")
  }
  # Check for nevus
  else if(grepl("nevus|nevi", char_lower) | grepl("nevus|nevi", desc_lower)) {
    return("Nevus")
  }
  # Check for normal
  else if(grepl("normal", char_lower) | grepl("normal", desc_lower)) {
    return("Normal")
  }
  else {
    return("Unknown")
  }
}

# Apply function to identify sample types
mel_pheno$sample_type <- mapply(identify_sample_type, 
                                mel_pheno$characteristics_ch1, 
                                mel_pheno$description)


# -------------------- Extract Specific Numbers of Samples ----------------------
set.seed(123)  # For reproducible sampling

# Extract samples
normal_samples <- mel_pheno[mel_pheno$sample_type == "Normal", ]
nevus_samples <- mel_pheno[mel_pheno$sample_type == "Nevus", ]
melanoma_samples <- mel_pheno[mel_pheno$sample_type == "Melanoma", ]

normal_samples

# Sample the required numbers (if available)
extract_samples <- function(df, n_samples, sample_name) {
  if(nrow(df) >= n_samples) {
    selected <- df[sample(nrow(df), n_samples), ]
    cat(paste("Selected", n_samples, sample_name, "samples\n"))
    return(selected)
  } else {
    cat(paste("Warning: Only", nrow(df), sample_name, "samples available, need", n_samples, "\n"))
    return(df)
  }
}

# Extract required samples
selected_normal <- extract_samples(normal_samples, 7, "Normal")
selected_nevus <- extract_samples(nevus_samples, 14, "Nevus")
selected_melanoma <- extract_samples(melanoma_samples, 21, "Melanoma")

# Combine selected samples for analysis (Melanoma vs Nevus only)
analysis_samples <- rbind(selected_melanoma, selected_nevus)
analysis_samples$sample_type <- factor(analysis_samples$sample_type, levels = c("Nevus", "Melanoma"))

print("Selected samples for analysis:")
print(table(analysis_samples$sample_type))

# -------------------- Prepare Expression Data ----------------------
# Extract expression data for selected samples
sample_ids <- rownames(analysis_samples)
analysis_expr <- expr_matrix[, sample_ids]
analysis_expr

# Check dimensions
cat("Expression matrix dimensions:", dim(analysis_expr), "\n")
cat("Number of samples in analysis:", ncol(analysis_expr), "\n")


# -------------------- Data Preprocessing ----------------------
# Remove genes with low expression (optional but recommended)
# Calculate mean expression for each gene
gene_means <- rowMeans(analysis_expr)

# Filter genes with very low expression (you can adjust this threshold)
expressed_genes <- gene_means > quantile(gene_means, 0.25)  # Keep top 75% expressed genes
analysis_expr_filtered <- analysis_expr[expressed_genes, ]

cat("Genes after filtering:", nrow(analysis_expr_filtered), "\n")

# -------------------- Limma-Voom Analysis ----------------------
# Create design matrix
design <- model.matrix(~ sample_type, data = analysis_samples)
design

colnames(design) <- c("Intercept", "Melanoma_vs_Nevus")

print("Design matrix:")
print(design)

# Convert to DGEList object for voom
dge <- edgeR::DGEList(counts = analysis_expr_filtered)
dge

# Calculate normalization factors
dge <- calcNormFactors(dge)
dge

# Apply limma voom transformation
v <- voom(dge, design, plot = TRUE)
v

# Fit linear model
fit <- lmFit(v, design)
fit

# Apply empirical Bayes
fit <- eBayes(fit)
fit

# Mean-variance trend
plotSA(fit, main = "Mean variance trend, GSE3189")

# -------------------- Results ----------------------
# Get differential expression results
results <- topTable(fit, coef = "Melanoma_vs_Nevus", number = Inf, sort.by = "P")
head(results,  n = 20)


# Summary of results
cat("Differential Expression Results Summary:\n")
cat("Total genes analyzed:", nrow(results), "\n")

cat("Significantly upregulated genes (adj.P.Val < 0.05, logFC > 1):", 
    sum(results$adj.P.Val < 0.05 & results$logFC > 1), "\n")

cat("Significantly downregulated genes (adj.P.Val < 0.05, logFC < -1):", 
    sum(results$adj.P.Val < 0.05 & results$logFC < -1), "\n")

# View top results
print("Top 20 differentially expressed genes:")
print(results[1:20, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])

#--------- A column for gene names. 
results


# Add gene symbols from the ExpressionSet
# Using the $ operator with backticks
results$gene_symbol <- fData(new_gse3189[[1]])[rownames(results), ]$`Gene Symbol`


# View results with gene symbols
head(results, 20)

# Summary with gene symbols
cat("Top 20 differentially expressed genes:\n")
print(results[1:20, c("gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])

# ----------------------Looking at the extremes ----------------------------

# Add gene symbols first
results$gene_symbol <- fData(new_gse3189[[1]])[rownames(results), ]$`Gene Symbol`

# Get top 20 UPREGULATED genes (highest positive logFC)
top_upregulated <- results %>%
  filter(logFC > 0) %>%  # Only positive fold changes
  arrange(desc(logFC)) %>%  # Sort by highest logFC first
  head(20)

top_upregulated

# Get top 20 DOWNREGULATED genes (most negative logFC)
top_downregulated <- results %>%
  filter(logFC < 0) %>%  # Only negative fold changes
  arrange(logFC) %>%  # Sort by most negative logFC first
  head(20)

# Display results
cat("TOP 20 UPREGULATED GENES (highest fold change):\n")
print(top_upregulated[, c("gene_symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val")])

cat("\nTOP 20 DOWNREGULATED GENES (lowest fold change):\n")
print(top_downregulated[, c("gene_symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val")])


# Optional: Get significantly upregulated/downregulated genes
cat("\nTOP 20 SIGNIFICANTLY UPREGULATED GENES (adj.P.Val < 0.05):\n")
sig_upregulated <- results %>%
  filter(adj.P.Val < 0.05 & logFC > 0) %>%
  arrange(desc(logFC)) %>%
  head(20)
print(sig_upregulated[, c("gene_symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val")])


cat("\nTOP 20 SIGNIFICANTLY DOWNREGULATED GENES (adj.P.Val < 0.05):\n")
sig_downregulated <- results %>%
  filter(adj.P.Val < 0.05 & logFC < 0) %>%
  arrange(logFC) %>%
  head(20)
print(sig_downregulated[, c("gene_symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val")])

# -------------------- Visualization ----------------------

# Volcano plot
volcano_data <- data.frame(
  logFC = results$logFC,
  negLog10P = -log10(results$P.Value),
  significant = results$adj.P.Val < 0.05 & abs(results$logFC) > 1
)

volcano_plot <- ggplot(volcano_data, 
                       aes(x = logFC, y = negLog10P, fill = significant)) + # padj. less than 0.05. 
  geom_point(alpha = 0.6, pch = 21, size = 1) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "#69b3a2")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkblue", alpha = 0.75) +
  labs(
    title = "Volcano Plot: Melanoma vs Nevus",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Significant\n(adj.P < 0.05 & |logFC| > 1)"
  ) +
  theme_classic(base_family = 'Times')
  
print(volcano_plot)

# ------------------------------------------------------------



#----
# MA plot
ma_plot <- ggplot(results,
                  mapping = aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.5, pch = 21, fill = ifelse(results$adj.P.Val < 0.05, "#69b3a2", "gray")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkblue",  alpha = 0.75) +
  geom_smooth(method = 'gam', color = "black") + 
  labs(
    title = "MA Plot: Melanoma vs Nevus",
    x = "Average Expression",
    y = "Log2 Fold Change"
  ) +
  theme_classic(base_family = 'Times')
  
print(ma_plot)



# -------------------- Save Results ----------------------
# Save results to file
write.csv(results, "melanoma_vs_nevus_DE_results.csv", row.names = TRUE)

# Save significant genes only
significant_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
write.csv(significant_genes, "significant_melanoma_vs_nevus_genes.csv", row.names = TRUE)

cat("Analysis complete! Results saved to CSV files.\n")
cat("Significant genes (adj.P < 0.05, |logFC| > 1):", nrow(significant_genes), "\n")

# -------------------- Additional Quality Control ----------------------
# Sample correlation heatmap

# Calculate sample correlations
sample_cor <- cor(analysis_expr_filtered)
sample_cor


# Create annotation for heatmap
annotation_col <- data.frame(
  Sample_Type = analysis_samples$sample_type,
  row.names = colnames(sample_cor)
)

# Plot correlation heatmap
 
pheatmap(sample_cor,
      annotation_col = annotation_col,
      main = "Sample Correlation Heatmap",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_col = 8,             
      fontsize_row = 8,
      family = "Times",
      color = viridis::cividis(28),
      angle_col = 45)

# --- PCA plot ------- 
pca_result <- prcomp(t(analysis_expr_filtered), scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Sample_Type = analysis_samples$sample_type
)

pca_plot <- ggplot(pca_data, 
                   mapping = aes(x = PC1, y = PC2,stroke = 1.2, fill = Sample_Type, shape=Sample_Type)) +
  geom_point(size = 6, pch = 21, alpha = 0.75) +
  labs(
    title = "PCA Plot \n Melanoma vs Nevus Samples",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
  ) +
  scale_fill_manual(
    values = c("Nevus" = "#FF7420", "Melanoma" = "#69b3a2")  # Custom bold colors
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),  # Bold centered title
    axis.title = element_text(size = 16, face = "bold"),  # Bold axis titles
    axis.text = element_text(size = 14),  # Larger axis tick labels
    legend.title = element_text(size = 14, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 12),  # Larger legend text
    axis.line = element_line(size = 1),  # Thicker axis lines
    axis.ticks = element_line(size = 1)  # Thicker axis ticks
  )

print(pca_plot)

# Scree Plot
factoextra::fviz_eig(pca_result, addlabels = TRUE,
                     main = 'Explained variance(%)for each Principal Component', 
                     col = "#69b3a2")

# My own scree plot
pca_result %>%
  broom::tidy(matrix = "eigenvalues") %>%
  # Make a bar plot of the principal components.
  ggplot(mapping = aes(x = PC, y = percent)) +
  # make a bar plot
  geom_col(colour = "darkblue", fill = "#9933ff", alpha = 0.3) +
  # Add a line and then some points.
  geom_line()+
  geom_point(shape=21, color="black", fill="darkslategray4", size=5) +
  # Adjust the x axis and the y axis.
  scale_x_continuous(breaks = 1:12, limits = c(0,10)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  # Add a grid
  theme_minimal_hgrid(12)+
  theme_classic(base_family = "Helvetica",
                base_size = 12) +
  # Add some labels.
  labs(
    y = "Explained variance(%) by each Principal Component",
    x = "The Principal Component",
    title = "Principal Component Analysis Plot",
    subtitle = " "
    ,
    caption = "Figure")

summary(pca_result)










# UMAP
ump <-umap(t(analysis_expr_filtered), n_neighbors = 15, random_state = 123)

plot(ump$layout, 
     main = "UMAP plot, nbrs=15", xlab = "", ylab = "", pch = 20, cex = 1.5)
text(ump$layout, 
     labels = rownames(ump$layout), cex = 0.45, pos = 3)


# Create UMAP data frame with sample type information
umap_data <- data.frame(
  UMAP1 = ump$layout[,1],
  UMAP2 = ump$layout[,2],
  Sample_Type = analysis_samples$sample_type,
  Sample_ID = rownames(analysis_samples)
)

umap_plot <- ggplot(umap_data, 
                    mapping = aes(x = UMAP1, y = UMAP2,
                                  color = Sample_Type, shape = Sample_Type)) +
  geom_point(size = 6, alpha = 0.75) +
  #scale_color_brewer(type = "qual") + 
  scale_color_manual(values = c("Nevus" = "#cc5c00", "Melanoma" = "#4c8575"))+
  ggrepel::geom_label_repel(aes(label = Sample_ID), 
                  colour = "black",
                  size = 3, 
                  max.overlaps = Inf,
                  box.padding = 0.6,
                  point.padding = 0.3) +
  labs(
    title = "UMAP Plot: Melanoma vs Nevus Samples",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_classic()

print(umap_plot)
# ===============================================================================
# ===============================================================================

# ClusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# A DEG result dataframe called results
# Filter significant DEGs (adjusted p < 0.05, |logFC| > 1)
sig_genes <- results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)

# 2. Get unique gene symbols from filtered results
gene_symbols <- unique(sig_genes$gene_symbol)
gene_symbols

# 3. Map gene symbols to Entrez IDs using org.Hs.eg.db
gene_entrez <-clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez

# 4. Keep only genes that mapped successfully - Entrez IDs for significant DEGs
sig_genes_entrez <- gene_entrez$ENTREZID
sig_genes_entrez

sig_genes_mapped <- sig_genes %>%
  inner_join(gene_entrez, by = c("gene_symbol" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

gene_list <- sig_genes_mapped$logFC
gene_list 

names(gene_list) <- sig_genes_mapped$ENTREZID

# Sort decreasingly (important for some enrichment analyses and plotting)
gene_list <- sort(gene_list, decreasing=TRUE)
gene_list

# GO enrichment (Biological Process)
ego <- clusterProfiler::enrichGO(gene = sig_genes_entrez,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Visualize top enriched terms
clusterProfiler::dotplot(ego, showCategory=15) + ggtitle("Top GO Biological Processes")

# KEGG Pathway enrichment
ekegg <-  clusterProfiler::enrichKEGG(gene = sig_genes_entrez,
                    organism = 'hsa',
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)

clusterProfiler::dotplot(ekegg, showCategory=15) + ggtitle("Top KEGG Pathways")


clusterProfiler::cnetplot(ego, categorySize="geneNum", foldChange=gene_list) + 
  ggtitle("Gene-Pathway Network")

#=================================================================================
#==================================================================================

# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Create volcano plot data with gene symbols
volcano_data <- data.frame(
  logFC = results$logFC,
  negLog10P = -log10(results$P.Value),
  adj_P_Val = results$adj.P.Val,
  gene_symbol = results$gene_symbol,
  probe_id = rownames(results),
  significant = results$adj.P.Val < 0.05 & abs(results$logFC) > 1
)

# Get probe IDs for top genes to label
top_up_probes <- rownames(top_upregulated)
top_down_probes <- rownames(top_downregulated)

# Create labelling column - only label top up/down regulated genes
volcano_data$label <- ifelse(
  volcano_data$probe_id %in% c(top_up_probes, top_down_probes), 
  volcano_data$gene_symbol, 
  ""
)

# Remove empty gene symbols from labels
volcano_data$label <- ifelse(
  is.na(volcano_data$label) | volcano_data$label == "", 
  "", 
  volcano_data$label
)

# Create the volcano plot with labels
volcano_plot_labeled <- ggplot(volcano_data, 
                               aes(x = logFC, y = negLog10P, 
                                   fill = significant)) +

  geom_point(alpha = 0.6, pch = 21, size = 1) +
  # Add gene labels for top up/down regulated genes only
  ggrepel::geom_label_repel(
    aes(label = label),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    fill = NA,        # Remove background fill
    color = "black",  # Make labels black for readability
    fontface = "bold"
  ) +
  
  # Highlight the labeled points
  geom_point(
    data = volcano_data[volcano_data$label != "", ],
    aes(x = logFC, y = negLog10P),
    color = "darkblue",
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "#69b3a2")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkblue", alpha = 0.75) +
  
  labs(
    title = "Volcano Plot: Melanoma vs Nevus",
    subtitle = "Top 20 up/down-regulated genes labeled",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Significant\n(adj.P < 0.05 & |logFC| > 1)"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Display the plot
print(volcano_plot_labeled)
# this is just for the labels. 



library(GSVA)
library(ggplot2)
library(factoextra)
library(pheatmap)
library(dplyr)
library(viridis)
# =============================================================================
# ssGSEA Analysis for Top 20 Up/Down Regulated Genes (Melanoma vs Nevus)
# =============================================================================

# Extract your top and bottom 20 genes from the differential expression results
sig_upregulated
sig_downregulated
# Get top 20 UPREGULATED genes (highest positive logFC)
top_upregulated <- sig_upregulated

# Get top 20 DOWNREGULATED genes (most negative logFC)
top_downregulated <- sig_downregulated

# Print the genes we're using
cat("TOP 20 UPREGULATED GENES for ssGSEA:\n")
print(top_upregulated[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val")])

cat("\nTOP 20 DOWNREGULATED GENES for ssGSEA:\n")
print(top_downregulated[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val")])

# =============================================================================
# Prepare Gene Sets for GSVA
# =============================================================================

# Create gene sets using probe IDs (which match your expression matrix rownames)
gene_sets <- list(
  "Melanoma_Upregulated_Top20" = rownames(top_upregulated),
  "Melanoma_Downregulated_Top20" = rownames(top_downregulated)
)

# Optional: Create combined gene set
gene_sets[["All_Top_DE_Genes"]] <- c(rownames(top_upregulated), rownames(top_downregulated))

# Display gene sets
cat("\nGene sets created:\n")
for(i in 1:length(gene_sets)) {
  cat(names(gene_sets)[i], ":", length(gene_sets[[i]]), "genes\n")
}

# =============================================================================
# Prepare Expression Matrix for ssGSEA
# =============================================================================

# Use your filtered expression matrix from the differential expression analysis
expr_matrix_ssgsea <- analysis_expr_filtered

# Ensure we have all the genes in our gene sets
missing_genes <- setdiff(unlist(gene_sets), rownames(expr_matrix_ssgsea))
if(length(missing_genes) > 0) {
  cat("Warning: Missing genes from expression matrix:", missing_genes, "\n")
  # Remove missing genes from gene sets
  gene_sets <- lapply(gene_sets, function(x) x[x %in% rownames(expr_matrix_ssgsea)])
}

cat("Expression matrix dimensions:", dim(expr_matrix_ssgsea), "\n")
cat("Samples:", colnames(expr_matrix_ssgsea), "\n")

# =============================================================================
# Run ssGSEA Analysis
# =============================================================================

# Perform ssGSEA
# Note: Since you're using microarray data (from GEO), use kcdf="Gaussian"
# For count data, you would use kcdf="Poisson"
ssGSEA_enrichments <- GSVA::gsva(expr = expr_matrix_ssgsea,
                           gset.idx.list = gene_sets,
                           method = "ssgsea", 
                           kcdf = "Gaussian",        # Appropriate for microarray data
                           min.sz = 5,               # Minimum gene set size
                           max.sz = 500,             # Maximum gene set size
                           parallel.sz = 1,         # Number of cores
                           verbose = TRUE)

# Check results
cat("ssGSEA results dimensions:", dim(ssGSEA_enrichments), "\n")
print("ssGSEA enrichment scores:")
print(ssGSEA_enrichments)

# =============================================================================
# Visualization and Analysis
# =============================================================================

# 1. Create a data frame for plotting
enrichment_df <- as.data.frame(t(ssGSEA_enrichments))
enrichment_df$Sample_ID <- rownames(enrichment_df)
enrichment_df$Sample_Type <- analysis_samples[rownames(enrichment_df), "sample_type"]

# View the enrichment scores
print("Enrichment scores by sample:")
print(enrichment_df)

# 2. Box plots comparing enrichment between sample types
library(reshape2)
library(ggridges)
enrichment_long <- melt(enrichment_df, 
                        id.vars = c("Sample_ID", "Sample_Type"),
                        variable.name = "Gene_Set", 
                        value.name = "Enrichment_Score")


# Create box plots
boxplot_enrichment <- ggplot(enrichment_long, 
                             aes(x = Sample_Type, y = Enrichment_Score, fill = Sample_Type)) +
  geom_boxplot(colour = "darkblue", alpha = 0.7, outlier.alpha = 0) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
  facet_wrap(~Gene_Set, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Nevus" = "#cc5c00", "Melanoma" = "#4c8575"))+
  #scale_fill_brewer(palette="Accent")+
  labs(title = "Boxplots of ssGSEA Enrichment Scores: Melanoma vs Nevus",
       subtitle = "Top 20 Up/Down Regulated Gene Sets",
       x = "Sample Type",
       y = "Enrichment Score",
       fill = "Sample Type",
       caption = "GEO Dataset") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"))

print(boxplot_enrichment)

# 3. Heatmap of enrichment scores
annotation_col <- data.frame(
  Sample_Type = analysis_samples[colnames(ssGSEA_enrichments), "sample_type"],
  row.names = colnames(ssGSEA_enrichments)
)

pheatmap(ssGSEA_enrichments,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = viridis::viridis(50),
         main = "ssGSEA Enrichment Scores Heatmap",
         fontsize_row = 10,
         fontsize_col = 8,
         angle_col = 90)

# 4. PCA on ssGSEA enrichment scores
if(nrow(ssGSEA_enrichments) > 1) {  # Only do PCA if we have multiple gene sets
  pca_ssGSEA <- prcomp(t(ssGSEA_enrichments), scale. = TRUE, center = TRUE)
  
  # Calculate variance explained
  var_explained <- round(100 * pca_ssGSEA$sdev^2 / sum(pca_ssGSEA$sdev^2), 2)
  
  # Create PCA data frame
  pc_ssGSEA <- data.frame(
    Sample_ID = rownames(pca_ssGSEA$x),
    PC1 = pca_ssGSEA$x[, 1], 
    PC2 = pca_ssGSEA$x[, 2],
    Sample_Type = analysis_samples[rownames(pca_ssGSEA$x), "sample_type"]
  )
  
  # PCA plot
  pca_plot <- ggplot(pc_ssGSEA, 
                     mapping = aes(x = PC1, y = PC2, 
                                   fill = Sample_Type, shape = Sample_Type, stroke = 1.2)) +
    geom_point(colour = "darkblue",pch = 21, size = 6, alpha = 0.8) +
    scale_fill_manual(values = c("Nevus" = "#cc5c00", "Melanoma" = "#4c8575"))+
    labs(title = "PCA - ssGSEA Enrichment Scores",
         subtitle = "Based on Top 20 Up/Down Regulated Gene Sets",
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         color = "Sample Type",
         shape = "Sample Type") +
    theme_classic() +
    theme(legend.position = "bottom")
  
  print(pca_plot)
}

# =============================================================================
# Statistical Testing
# =============================================================================

# Perform t-tests for each gene set between sample types
stats_results <- data.frame()

for(gene_set in rownames(ssGSEA_enrichments)) {
  nevus_scores <- ssGSEA_enrichments[gene_set, analysis_samples$sample_type == "Nevus"]
  melanoma_scores <- ssGSEA_enrichments[gene_set, analysis_samples$sample_type == "Melanoma"]
  
  # Perform t-test
  t_test <- t.test(melanoma_scores, nevus_scores)
  
  # Store results
  stats_results <- rbind(stats_results, data.frame(
    Gene_Set = gene_set,
    Melanoma_Mean = mean(melanoma_scores),
    Nevus_Mean = mean(nevus_scores),
    Mean_Difference = mean(melanoma_scores) - mean(nevus_scores),
    T_Statistic = t_test$statistic,
    P_Value = t_test$p.value,
    CI_Lower = t_test$conf.int[1],
    CI_Upper = t_test$conf.int[2]
  ))
}

# Add adjusted p-values
stats_results$Adj_P_Value <- p.adjust(stats_results$P_Value, method = "BH")

# Print statistical results
cat("\nStatistical comparison of enrichment scores between Melanoma and Nevus:\n")
print(stats_results)


cat("\n=============================================================================\n")
cat("SUMMARY OF ssGSEA ANALYSIS\n")
cat("=============================================================================\n")

for(i in 1:nrow(stats_results)) {
  gene_set <- stats_results$Gene_Set[i]
  p_val <- stats_results$P_Value[i]
  adj_p_val <- stats_results$Adj_P_Value[i]
  mean_diff <- stats_results$Mean_Difference[i]
  
  cat("\n", gene_set, ":\n")
  cat("  Mean difference (Melanoma - Nevus):", round(mean_diff, 4), "\n")
  cat("  P-value:", format(p_val, scientific = TRUE), "\n")
  cat("  Adjusted P-value:", format(adj_p_val, scientific = TRUE), "\n")
  
  if(adj_p_val < 0.05) {
    direction <- ifelse(mean_diff > 0, "higher", "lower")
    cat("  ** SIGNIFICANT: Melanoma samples have", direction, "enrichment **\n")
  } else {
    cat("  Not significant at adj.P < 0.05\n")
  }
}

# =============================================================================
# Save Results
# =============================================================================

# Save enrichment scores
write.csv(enrichment_df, "ssGSEA_enrichment_scores_top20genes.csv", row.names = FALSE)

# Save statistical results
write.csv(stats_results, "ssGSEA_statistical_results_top20genes.csv", row.names = FALSE)

# Save gene sets used
gene_sets_df <- data.frame(
  Gene_Set = rep(names(gene_sets), sapply(gene_sets, length)),
  Probe_ID = unlist(gene_sets),
  Gene_Symbol = results[unlist(gene_sets), "gene_symbol"]
)
write.csv(gene_sets_df, "gene_sets_used_ssGSEA.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved to CSV files.\n")
cat("Files created:\n")
cat("- ssGSEA_enrichment_scores_top20genes.csv\n")
cat("- ssGSEA_statistical_results_top20genes.csv\n")
cat("- gene_sets_used_ssGSEA.csv\n")

# =================================================

# =============================================================================
# Section 4: Tumor Purity Estimation for Melanoma vs Nevus Samples
# =============================================================================
# To Install it : 
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)

# ESTIMATE Analysis for Tumor Purity Calculation
# Melanoma vs Nevus samples from GEO dataset
# Simplified ESTIMATE Analysis Workflow
library(estimate)
library(ggplot2)
library(dplyr)
library(reshape2)


cat("=== ESTIMATE Analysis: Tumor Purity ===\n")

# Step 1: Prepare Expression Data

# Extract gene symbols from your dataset
gene_symbols <- fData(new_gse3189[[1]])[rownames(analysis_expr_filtered), ]$`Gene Symbol`


# Create expression matrix and assign gene symbols as rownames
expr_for_estimate <- analysis_expr_filtered
rownames(expr_for_estimate) <- gene_symbols

# Remove rows with missing or invalid gene symbols
expr_for_estimate <- expr_for_estimate[!is.na(rownames(expr_for_estimate)) & 
                                         rownames(expr_for_estimate) != "" & 
                                         rownames(expr_for_estimate) != "---", ]

# Handle duplicate gene symbols by keeping the row with the highest mean expression
if (any(duplicated(rownames(expr_for_estimate)))) {
  expr_for_estimate <- expr_for_estimate %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    group_by(gene_symbol) %>%
    summarize(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames("gene_symbol")
}

cat("Expression matrix finalized with", nrow(expr_for_estimate), "unique genes.\n")

# Step 2–5: Define and Run the ESTIMATE Pipeline
cat("Running ESTIMATE pipeline...\n")

# Function to handle the ESTIMATE analysis
ESTIMATE_fun <- function(expr) {
  # Write the expression data to a GCT-compatible file
  write.table(as.data.frame(expr), file = "rma.data.gct", quote = F, sep = '\t')
  
  # Filter common genes and calculate ESTIMATE scores
  filterCommonGenes(input.f = 'rma.data.gct', output.f = "RMA_10412.gct", id = "GeneSymbol")
  estimateScore("RMA_10412.gct", "estimate_score.gct", platform = "affymetrix")
  
  # Read the calculated ESTIMATE scores
  estimate <- read.table("estimate_score.gct", sep = '\t', row.names = 1, header = T, skip = 2)
  estimate <- estimate[, -1]  # Remove the description column
  
  # Clean up temporary files
  file.remove("rma.data.gct")
  file.remove("RMA_10412.gct")
  file.remove("estimate_score.gct")
  
  # Process scores and return a formatted data frame
  estimate_scores <- as.data.frame(t(estimate))
  colnames(estimate_scores) <- c('StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity')
  return(estimate_scores)
}

# Run the ESTIMATE function
scores <- ESTIMATE_fun(expr_for_estimate)
scores

# Add sample type information to the scores
scores$Sample_Type <- analysis_samples[rownames(scores), "sample_type"]
scores$Sample_Type 

# Step 6: Statistical Comparisons
cat("Performing statistical comparisons...\n")

# Perform group-based statistical testing
comparison_results <- data.frame(
  Score = c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"),
  Melanoma_Mean = NA,
  Nevus_Mean = NA,
  P_Value = NA
)

for (i in 1:nrow(comparison_results)) {
  score <- comparison_results$Score[i]
  
  # Extract scores for each group
  melanoma_vals <- scores[scores$Sample_Type == "Melanoma", score]
  nevus_vals <- scores[scores$Sample_Type == "Nevus", score]
  
  # Calculate group means and p-value using a t-test
  comparison_results$Melanoma_Mean[i] <- round(mean(melanoma_vals), 3)
  comparison_results$Nevus_Mean[i] <- round(mean(nevus_vals), 3)
  comparison_results$P_Value[i] <- round(t.test(melanoma_vals, nevus_vals)$p.value, 4)
}

cat("Statistical comparisons complete. Results:\n")
print(comparison_results)

# Step 7: Visualization


# Melt data for ggplot
plot_data <- melt(scores[, c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity", "Sample_Type")],
                  id.vars = "Sample_Type")

# Visualize ESTIMATE scores by sample type
p1 <- ggplot(plot_data, 
             mapping = aes(x = Sample_Type, y = value, fill = Sample_Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  theme_classic(base_family = "Times") +
  labs(title = "ESTIMATE Scores: Melanoma vs Nevus",
       x = "Sample Type", y = "Score") +
  scale_fill_manual(values = c("Nevus" = "#cc5c00", "Melanoma" = "#4c8575"))+
  theme(strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

print(p1)

# Step 8: Save Results and Clean Up
cat("Analysis complete. Saving results...\n")

# Save scores and comparison results for later use
write.table(scores, file = "estimate_scores.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(comparison_results, file = "comparison_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Results saved to 'estimate_scores.txt' and 'comparison_results.txt'.\n")

# There are 2 samples that have a low tumour purity in melanoma. 
# Subset melanoma samples
melanoma_scores <- scores[scores$Sample_Type == "Melanoma", ]
# Sort by TumorPurity ascending and get the top 2
lowest_purity <- melanoma_scores[order(melanoma_scores$TumorPurity), ][1:2, ]

lowest_purity

# There are nevus samples that also have low tumour purity (Subset nevus samples)
nevus_scores <- scores[scores$Sample_Type == "Nevus", ]

# Sort by TumorPurity ascending
low_purity_nevus <- nevus_scores[order(nevus_scores$TumorPurity), ]

# View the lowest ones
head(low_purity_nevus, 3) 
