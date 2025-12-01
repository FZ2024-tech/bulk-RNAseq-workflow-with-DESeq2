
# Bulk RNA-Seq Analysis: Alpha Mutant vs WildType (Demo)


# This script demonstrates bulk RNA-Seq analysis workflow using a subset of publicly available samples.
# Only 4 samples were selected for demonstration purposes.

# Data source / citation:
# Lee HK, Knabl L, Knabl Sr L, Wieser M, Mur A, Zabernigg A, Schumacher J, Kapferer S, Kaiser N, Furth PA, Hennighausen L. 
# Immune transcriptome analysis of COVID-19 patients infected with SARS-CoV-2 variants carrying the E484K escape mutation identifies a distinct gene module. 
# Scientific Reports. 2022;12:2784. https://doi.org/10.1038/s41598-022-06985-3


# 1. Install and Load Required Packages

# Bioconductor packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "sva", "edgeR"))

# CRAN packages
install.packages(c("tidyverse", "pheatmap", "ggrepel", "R.utils"))

# Load libraries
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(data.table)
library(pheatmap)
library(ggrepel)
library(R.utils)


# 2. Set Working Directory and List Files


setwd("E:/Bulk RNA_Seq-2025/Raw_Data")
getwd()
list.files()



# Example: read each sample individually
library(data.table)

# List all your file names
files <- c( "GSM5728881_Alpha_E484K_S39_1st.txt.gz",
  "GSM5728884_Alpha_E484K_S40_1st.txt.gz",
  "GSM5728887_Alpha_E484K_S41_1st.txt.gz",
  "GSM5728890_Alpha_E484K_S45_1st.txt.gz",
  "GSM5728892_Alpha_E484K_S51_1st.txt.gz",
  "GSM5728894_Alpha_E484K_S52_1st.txt.gz",
  "GSM5728911_Alpha_S08_1st.txt.gz",
  "GSM5728913_Alpha_S10_1st.txt.gz",
  "GSM5728915_Alpha_S11_1st.txt.gz",
  "GSM5728918_Alpha_S19_1st.txt.gz",
  "GSM5728920_Alpha_S20_1st.txt.gz",
  "GSM5728922_Alpha_S21_1st.txt.gz")

# Read all samples into a list
samples_list <- lapply(files, fread)


#  for all samples dynamically:
count_files <- list.files(pattern = "*_1st.txt.gz")
count_list <- list()
gene_ids <- fread(count_files[1])[[1]]  # Assume first column contains gene IDs

for (file in count_files) {
  sample_data <- fread(file)
  sample_name <- gsub("_1st.txt.gz", "", file)
  count_list[[sample_name]] <- sample_data[[2]]  # Second column: counts
  }

# Create count matrix
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_ids


# 4. Map Gene IDs to Symbols

# Test common ID types
gene_samples <- head(rownames(count_matrix), 20)
id_types <- c("ENTREZID", "ENSEMBL", "REFSEQ", "SYMBOL")
success_rates <- sapply(id_types, function(keytype) {
  mapped <- tryCatch(mapIds(org.Hs.eg.db, keys = gene_samples, column = "SYMBOL", keytype = keytype, multiVals = "first"),
                     error = function(e) rep(NA, length(gene_samples)))
  mean(!is.na(mapped))
})

best_type <- id_types[which.max(success_rates)]
gene_symbols <- if(best_type == "SYMBOL") rownames(count_matrix) else 
  mapIds(org.Hs.eg.db, keys = rownames(count_matrix), column = "SYMBOL", keytype = best_type, multiVals = "first")

# Replace IDs with gene symbols, handle duplicates
symbols_available <- !is.na(gene_symbols)
final_gene_names <- ifelse(symbols_available, gene_symbols, rownames(count_matrix))
duplicate_genes <- duplicated(final_gene_names) | duplicated(final_gene_names, fromLast = TRUE)
final_gene_names[duplicate_genes] <- paste0(final_gene_names[duplicate_genes], "_", rownames(count_matrix)[duplicate_genes])
rownames(count_matrix) <- final_gene_names


# Prepare Sample Metadata


samples <- data.frame(
  sample_id = c("GSM5728881", "GSM5728884", "GSM5728887", "GSM5728890"),
  group = c("Mutant", "Mutant", "Mutant", "WildType"),
  stringsAsFactors = FALSE)

rownames(samples) <- samples$sample_id
samples$group <- factor(samples$group, levels = c("WildType", "Mutant"))


# Subset Count Data


common_samples <- intersect(rownames(samples), colnames(count_matrix))
count_sub <- count_matrix[, common_samples, drop = FALSE]
samples_sub <- samples[common_samples, , drop = FALSE]
stopifnot(all(rownames(samples_sub) == colnames(count_sub)))

#  DESeq2 Analysis

# Create DESeq2 object
dds_sub <- DESeqDataSetFromMatrix(countData = count_sub, colData = samples_sub, design = ~ group)

#filter low counts
dds_sub <- dds_sub[rowSums(counts(dds_sub)) > 10, ]  # filter low counts

# Run DESeq2
dds_sub <- DESeq(dds_sub)
res_sub <- results(dds_sub, contrast = c("group", "Mutant", "WildType"))
print(paste("After filtering:", nrow(dds_sub), "genes remaining"))
summary(res_sub)

resultsNames(dds_sub)

# LFC shrinkage
resLFC_sub <- lfcShrink(dds_sub, contrast = c("group", "WildType", "Mutant"), type = "normal")# Create results dataframe


# Create results dataframe
res_df <- as.data.frame(resLFC_sub) %>% 
  tibble::rownames_to_column("gene") %>% 
  arrange(padj)

print("Top 10 significant genes:")
print(res_df %>% filter(padj < 0.05) %>% head(10))


# Save Results


dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, "results/tables/DESeq2_Full.csv", row.names = FALSE)

sig_DEGs <- res_df %>% filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig_DEGs, "results/tables/DESeq2_Significant_DEGs.csv", row.names = FALSE)


#  Visualization

vsd_sub <- vst(dds_sub, blind = FALSE)

# PCA
print("Creating visualizations...")

# Variance-stabilizing transformation
vsd_sub <- vst(dds_sub, blind = FALSE)

# PCA Plot
pcaData <- plotPCA(vsd_sub, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 10) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  ggtitle("PCA: Alpha E484K vs Alpha WT") +
  theme_minimal()

ggsave("results/figures/PCA_Alpha_Mutant_vs_WT.png", pca_plot, width = 8, height = 6, dpi = 300)



# Volcano plot 
volcano_data <- res_df %>% 
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1,
         label = ifelse(significant & abs(log2FoldChange) > 1.5, gene, NA))

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15, box.padding = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot: Mutant vs WildType",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "none")


# HeatMap
topN <- 20
top_genes <- res_df %>% 
  filter(!is.na(padj)) %>% 
  arrange(padj) %>% 
  pull(gene) %>% 
  head(topN)

print("Top genes for heatmap:")
print(top_genes)

# Create heatmap matrix
mat_top <- assay(vsd_sub)[top_genes, , drop = FALSE]
mat_top_z <- t(scale(t(mat_top)))

# Create annotation
annotation_col <- data.frame(Group = samples_sub$group)
rownames(annotation_col) <- colnames(mat_top_z)

# Create heatmap
pheatmap(mat_top_z,
         annotation_col = annotation_col,
         main = paste("Top", topN, "DEGs: Alpha Mutant vs WT"),
         filename = "results/figures/Heatmap_top20_DEGs.png",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         width = 8,
         height = 6)
print("Heatmap saved"


# Create annotation
annotation_col <- data.frame(Group = samples_sub$group)
rownames(annotation_col) <- colnames(mat_top_z)


# Heatmap Top 20 DEGs
pheatmap(mat_top_z,
         annotation_col = annotation_col,
         main = paste("Top", topN, "DEGs: Alpha Mutant vs WT"),
         filename = "results/figures/Heatmap_top20_DEGs.png",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         width = 8,
         height = 6)
print("Heatmap saved")
