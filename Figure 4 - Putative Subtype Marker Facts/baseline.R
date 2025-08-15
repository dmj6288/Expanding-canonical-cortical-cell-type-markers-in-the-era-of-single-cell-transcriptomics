rm(list=ls())

suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dendextend)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))

set.seed(42) 

cell_names <- c("Astro", "Excite", "Inhibit", "Micro", "Oligo", "OPC")

rename_list <- list("Astro", "Astro", "Excite", "Inhibit", "Micro", "Micro", "Oligo", "OPC")
names(rename_list) <- c("Astro", "Astrocyte", "Excite", "Inhibit", "Microglia", "Micro", "Oligo", "OPC")

region_names <- tstrsplit(list.files("source/agg_matrix_RN011/"), "_1.rds", fixed=TRUE)[[1]]

all_canonical_gene_list <- readRDS("source/canonical_marker_lists.rds")
all_canonical_genes <- Reduce(union, all_canonical_gene_list)

counts_to_tpm <- function(counts_matrix, gene_lengths_kb) {
  # Step 1: Compute RPK (Reads Per Kilobase)
  rpk <- counts_matrix / gene_lengths_kb
  
  # Step 2: Compute the column sums (scaling factor per sample)
  scaling_factor <- colSums(rpk)
  
  # Step 3: Normalize by the scaling factor and multiply by 10^6
  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  
  return(tpm)
}


cohen_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

# # Compute Cohen's d for each row
# cohen_d_values <- sapply(common_rows, function(r) {
#   cohen_d(A[r, ], B[r, ])
# })


DESeq2_normalization <- function(countData, colData){
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ 1)
  dds <- DESeq(dds)
  
  normalized_counts <- counts(dds, normalized = TRUE)
  
  return(normalized_counts)
}

cell_specificity_index <- function(normalized_data_matrix){
  
  # Find the maximal values for each column
  #normalized_data_matrix <- normalized_data
  max_values <- apply(normalized_data_matrix, 1, max)
  
  # Divide each column by its maximal value
  zero_to_one_normalized_matrix <- sweep(normalized_data_matrix, 1, max_values, FUN="/")
  
  #zero_to_one_normalized_matrix <- normalized_data_matrix/apply(normalized_data_matrix, 2, max)
  specificity_values            <- (length(colnames(zero_to_one_normalized_matrix)) - rowSums(zero_to_one_normalized_matrix))/
    (length(colnames(zero_to_one_normalized_matrix)) - 1)
  
  combined_spec_exp <- cbind(zero_to_one_normalized_matrix, specificity_values)
  
  return(combined_spec_exp)
  
}

