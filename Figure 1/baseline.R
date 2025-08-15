rm(list=ls())

suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(openxlsx)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dendextend)))

cell_names <- c("Astro", "Excite", "Inhibit", "Micro", "Oligo", "OPC")

rename_list <- list("Astro", "Astro", "Excite", "Inhibit", "Micro", "Micro", "Oligo", "OPC")
names(rename_list) <- c("Astro", "Astrocyte", "Excite", "Inhibit", "Microglia", "Micro", "Oligo", "OPC")

all_canonical_gene_list <- readRDS("interm/canonical_marker_lists.rds")

all_canonical_gene_list[["Excite"]] <- c("CAMK2A", "SLC17A7", "SATB2")
all_canonical_gene_list[["Inhibit"]] <- c("GAD1", "GAD2")

source_file_path       <- "source/coldata_RN012/"
source_file_path_preQC <- "source/pre_QC_metadata_RN012/"

counts_to_tpm <- function(counts_matrix, gene_lengths_kb) {
  # Step 1: Compute RPK (Reads Per Kilobase)
  rpk <- counts_matrix / gene_lengths_kb
  
  # Step 2: Compute the column sums (scaling factor per sample)
  scaling_factor <- colSums(rpk)
  
  # Step 3: Normalize by the scaling factor and multiply by 10^6
  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  
  return(tpm)
}

get_second_to_last <- function(x) {
  parts <- unlist(strsplit(x, "_"))
  if (length(parts) >= 2) {
    return(parts[length(parts) - 1])
  } else {
    return(NA)  # or handle differently if fewer than 2 parts
  }
}


DESeq2_normalization <- function(countData, colData){
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ 1)
  dds <- DESeq(dds)
  
  normalized_counts <- counts(dds, normalized = TRUE)
  
  return(list(normalized_counts, dds))
}

df_row <- c()

for (cell in names(all_canonical_gene_list)){
  
  a <- cbind(all_canonical_gene_list[[cell]], rep(paste0(cell, " marker"), length(all_canonical_gene_list[[cell]])))
  df_row <- rbind(df_row, a)
  
}
rownames(df_row) <- df_row[, 1]
df_row <- as.data.frame(df_row)
df_row <- df_row[-c(1)]
colnames(df_row) <- "color"

inverse_gene_cell_hash <- list()

for (key in names(all_canonical_gene_list)){
  for (gene in all_canonical_gene_list[[key]]){
    inverse_gene_cell_hash[[gene]] <- key
  }
}


colors <- c("red", "green", "blue", "orange", "violet", "brown")#brewer.pal(6, "Set1")

# Define your colors as a named vector
color_vector <- c(
  "Astro marker"   = colors[1],
  "Excite marker"  = colors[2],
  "Inhibit marker" = colors[3],
  "Micro marker"   = colors[4],
  "Oligo marker"   = colors[5],
  "OPC marker"     = colors[6]
)

# Define your colors as a named vector
cell_vector <- c(
  "Astro"   = colors[1],
  "Excite"  = colors[2],
  "Inhibit" = colors[3],
  "Micro"   = colors[4],
  "Oligo"   = colors[5],
  "OPC"     = colors[6]
)


annotation_colors <- list()

annotation_colors$color       <- color_vector
annotation_colors$`Cell Type` <- cell_vector
