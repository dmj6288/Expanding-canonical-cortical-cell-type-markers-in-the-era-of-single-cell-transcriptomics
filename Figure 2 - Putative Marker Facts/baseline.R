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
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(ggrepel)))

set.seed(42) 

cell_names <- c("Astro", "Excite", "Inhibit", "Micro", "Oligo", "OPC")

rename_list <- list("Astro", "Astro", "Excite", "Inhibit", "Micro", "Micro", "Oligo", "OPC")
names(rename_list) <- c("Astro", "Astrocyte", "Excite", "Inhibit", "Microglia", "Micro", "Oligo", "OPC")

regions_to_avoid <- readRDS("source/regions_to_avoid.rds")

region_names <- tstrsplit(list.files("source/agg_matrices_RN012/"), "_1.rds", fixed=TRUE)[[1]]

all_canonical_gene_list    <- readRDS("source/canonical_marker_lists.rds")
all_canonical_gene_list[["Excite"]] <- c("CAMK2A", "SLC17A7", "SATB2")
all_canonical_gene_list[["Inhibit"]] <- c("GAD1", "GAD2")

saveRDS(all_canonical_gene_list, "interm/pan_canonical_list.rds")

all_canonical_genes <- Reduce(union, all_canonical_gene_list)
studies_to_keep_threshold <- readRDS("source/studies_to_keep_threshold.rds")

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

cell_specific_generate_gene_sets <- function(
    fidelity_thresh,           # now: named list or vector per cell type
    background_thresh,         # now: named list or vector per cell type
    padj_thresh = 0.05,
    log2FC_thresh = 1,
    cohen_d_thresh = cohen_d_thresh,
    most_frequent = 5,
    filtered_DESeq_matrix,
    df_merged,
    cell_names,
    regions_to_avoid
) {
  
  # Filter significant genes *per cell type* based on their thresholds
  df_merged_sig <- df_merged %>%
    rowwise() %>%
    filter(
      fidelity > fidelity_thresh[[Cell]],
      background < background_thresh[[Cell]],
      padj < padj_thresh,
      log2FoldChange > log2FC_thresh
    ) %>%
    ungroup()
  # fidelity_thresh [1] 10 20
  # Remove regions to avoid
  df_merged_sig <- df_merged_sig[df_merged_sig$region %in% tstrsplit(regions_to_avoid, "_")[[1]] == FALSE, ]
  
  # Top genes per Cell
  top_genes <- df_merged_sig %>%
    group_by(Cell, gene) %>%
    summarise(freq = n(), .groups = "drop_last") %>%
    arrange(Cell, desc(freq)) %>%
    slice_max(order_by = freq, n = most_frequent, with_ties = FALSE)
  
  # Build subset count matrix
  subset_count_matrix <- filtered_DESeq_matrix[top_genes$gene, ]
  global_gene_set_results <- list()
  putative  <- data.frame()
  canonical <- data.frame()
  
  # Loop through cells
  for (cell in cell_names) {
    log2_filtered_TPM_matrix <- log2(subset_count_matrix + 1)
    current_count_matrix <- log2_filtered_TPM_matrix
    genes <- rownames(current_count_matrix)

    current_count_matrix_cell <- current_count_matrix[, tstrsplit(colnames(current_count_matrix), "_", fixed=TRUE)[[1]] == cell]
    current_count_matrix_BG   <- current_count_matrix[, tstrsplit(colnames(current_count_matrix), "_", fixed=TRUE)[[1]] != cell]

    # Mann-Whitney U test
    p_values_mann <- sapply(genes, function(r) {
      test_result <- wilcox.test(current_count_matrix_cell[r, ], current_count_matrix_BG[r, ])
      test_result$p.value
    })

    adjusted_p_mann <- p.adjust(p_values_mann, method = "BH")

    # Cohen's d
    cohen_d_values <- sapply(genes, function(r) {
      cohen_d(current_count_matrix_cell[r, ], current_count_matrix_BG[r, ])
    })

    summary_table <- data.frame(
      row       = genes,
      mean_A    = sapply(genes, function(r) mean(current_count_matrix_cell[r, ])),
      mean_B    = sapply(genes, function(r) mean(current_count_matrix_BG[r, ])),
      p_value   = p_values_mann,
      adjusted_p = adjusted_p_mann,
      cohen_d   = cohen_d_values
    )

    summary_table <- summary_table[summary_table$row %in% top_genes[top_genes$Cell == cell, ]$gene, ]
    summary_table$cell <- cell

    # Apply Cohen's d filter
    summary_table <- summary_table[summary_table$cohen_d > cohen_d_thresh[[cell]], ]
    global_gene_set_results[[cell]] <- summary_table$row
    
    # Collect putative & canonical
    putative  <- rbind(putative,
                       pct_metrics[pct_metrics$gene %in% global_gene_set_results[[cell]] & pct_metrics$cluster == cell, ])
    canonical <- rbind(canonical,
                       pct_metrics[pct_metrics$gene %in% all_canonical_gene_list[[cell]] & pct_metrics$cluster == cell, ])
  }
  
  return(list(global_gene_set_results, putative, canonical))
}


run_GLM_analysis <- function(global_gene_set_results, putative, canonical, cell_counts) {
  
  results_list_fid <- list()
  results_list_bg  <- list()
  
  for (cluster in names(global_gene_set_results)) {
    
    # Filter for putative and canonical cell types
    putative_cell_type  <- putative  %>% filter(cluster == !!cluster) %>% mutate(group = "putative")
    canonical_cell_type <- canonical %>% filter(cluster == !!cluster) %>% mutate(group = "canonical")
    
    # Combine data
    test_data_combined <- bind_rows(putative_cell_type, canonical_cell_type) %>%
      mutate(region_key = region)
    
    # Cell counts lookup
    cell_counts_lookup <- cell_counts %>%
      dplyr::select(region_name = 1, cell_count = {{cluster}}) %>%
      filter(!is.na(cell_count)) %>%
      mutate(region_name = as.character(region_name))
    
    test_data_combined <- test_data_combined %>%
      left_join(cell_counts_lookup, by = c("region_key" = "region_name")) %>%
      dplyr::rename(cell_total = cluster)
    
    # Fidelity (pct1)
    test_data_combined_fid <- test_data_combined %>%
      mutate(pct1 = as.numeric(pct1),
             expressing = as.numeric(round(pct1 / 100 * as.numeric(cell_count))),
             non_expressing = as.numeric(cell_count) - expressing) %>%
      filter(!is.na(expressing), !is.na(non_expressing), !is.na(group))
    
    if (nrow(test_data_combined_fid) > 0 && length(unique(test_data_combined_fid$group)) > 1) {
      glm_fit <- glm(cbind(expressing, non_expressing) ~ group, family = binomial, data = test_data_combined_fid)
      coef_table <- summary(glm_fit)$coefficients
      
      if ("groupputative" %in% rownames(coef_table)) {
        res <- data.frame(
          cluster     = cluster,
          estimate    = coef_table["groupputative", "Estimate"],
          std.error   = coef_table["groupputative", "Std. Error"],
          z.value     = coef_table["groupputative", "z value"],
          p.value     = coef_table["groupputative", "Pr(>|z|)"],
          odds_ratio  = exp(coef_table["groupputative", "Estimate"])
        )
        results_list_fid[[cluster]] <- res
      }
    }
    
    # Background (pct2)
    test_data_combined_bg <- test_data_combined %>%
      mutate(pct2 = as.numeric(pct2),
             expressing = as.numeric(round(pct2 / 100 * as.numeric(cell_count))),
             non_expressing = as.numeric(cell_count) - expressing) %>%
      filter(!is.na(expressing), !is.na(non_expressing), !is.na(group))
    
    if (nrow(test_data_combined_bg) > 0 && length(unique(test_data_combined_bg$group)) > 1) {
      glm_fit <- glm(cbind(expressing, non_expressing) ~ group, family = binomial, data = test_data_combined_bg)
      coef_table <- summary(glm_fit)$coefficients
      
      if ("groupputative" %in% rownames(coef_table)) {
        res <- data.frame(
          cluster     = cluster,
          estimate    = coef_table["groupputative", "Estimate"],
          std.error   = coef_table["groupputative", "Std. Error"],
          z.value     = coef_table["groupputative", "z value"],
          p.value     = coef_table["groupputative", "Pr(>|z|)"],
          odds_ratio  = exp(coef_table["groupputative", "Estimate"])
        )
        results_list_bg[[cluster]] <- res
      }
    }
  }
  
  results_df_fid <- bind_rows(results_list_fid)
  results_df_bg  <- bind_rows(results_list_bg)
  
  GLM_results <- rbind(results_df_fid, results_df_bg)
  return(GLM_results)
}

generate_gene_sets <- function(
    fidelity_thresh = 80,
    background_thresh = 10,
    padj_thresh = 0.05,
    log2FC_thresh = 1,
    cohen_d_thresh = 3,
    most_frequent = 5,
    filtered_DESeq_matrix,
    df_merged,
    cell_names,
    regions_to_avoid
) {
  min_occurrences <- 5   # threshold for minimum occurrences
  
  # Filter significant genes
  df_merged_sig <- df_merged %>%
    filter(
      fidelity > fidelity_thresh,
      background < background_thresh,
      padj < padj_thresh,
      log2FoldChange > log2FC_thresh
    )
  
  # Remove regions to avoid
  df_merged_sig <- df_merged_sig[df_merged_sig$region %in% tstrsplit(regions_to_avoid, "_")[[1]] == FALSE, ]
  
  # Top genes per Cell
  top_genes <- df_merged_sig %>%
    group_by(Cell, gene) %>%
    summarise(freq = n(), .groups = "drop_last") %>%
    arrange(Cell, desc(freq)) %>%
    slice_max(order_by = freq, n = most_frequent, with_ties = FALSE)
  
  # Build subset count matrix
  subset_count_matrix <- filtered_DESeq_matrix[top_genes$gene, ]
  global_gene_set_results <- list()
  putative  <- data.frame()
  canonical <- data.frame()
  
  # Loop through cells
  for (cell in cell_names) {
    log2_filtered_TPM_matrix <- log2(subset_count_matrix + 1)
    current_count_matrix <- log2_filtered_TPM_matrix
    
    genes <- rownames(current_count_matrix)
    
    current_count_matrix_cell <- current_count_matrix[, tstrsplit(colnames(current_count_matrix), "_", fixed=TRUE)[[1]] == cell]
    current_count_matrix_BG   <- current_count_matrix[, tstrsplit(colnames(current_count_matrix), "_", fixed=TRUE)[[1]] != cell]
    
    # Mann-Whitney U test
    p_values_mann <- sapply(genes, function(r) {
      test_result <- wilcox.test(current_count_matrix_cell[r, ], current_count_matrix_BG[r, ])
      test_result$p.value
    })
    
    adjusted_p_mann <- p.adjust(p_values_mann, method = "BH")
    
    # Cohen's d
    cohen_d_values <- sapply(genes, function(r) {
      cohen_d(current_count_matrix_cell[r, ], current_count_matrix_BG[r, ])
    })
    
    # Summary table
    summary_table <- data.frame(
      row       = genes,
      mean_A    = sapply(genes, function(r) mean(current_count_matrix_cell[r, ])),
      mean_B    = sapply(genes, function(r) mean(current_count_matrix_BG[r, ])),
      p_value   = p_values_mann,
      adjusted_p = adjusted_p_mann,
      cohen_d   = cohen_d_values
    )
    
    # Keep only top_genes for the cell
    summary_table <- summary_table[summary_table$row %in% top_genes[top_genes$Cell == cell, ]$gene, ]
    summary_table$cell <- cell
    
    # Filter by Cohen's d threshold
    summary_table <- summary_table[summary_table$cohen_d > cohen_d_thresh, ]
    global_gene_set_results[[cell]] <- summary_table$row
    
    # Collect putative & canonical
    putative  <- rbind(putative,
                       pct_metrics[pct_metrics$gene %in% global_gene_set_results[[cell]] & pct_metrics$cluster == cell, ])
    
    canonical <- rbind(canonical,
                       pct_metrics[pct_metrics$gene %in% all_canonical_gene_list[[cell]] & pct_metrics$cluster == cell, ])
  }
  
  return(list(
    global_gene_set_results = global_gene_set_results,
    putative = putative,
    canonical = canonical
  ))
}

# Example usage:
# results <- generate_gene_sets(80, 10, 0.05, 1, 3, filtered_DESeq_matrix, df_merged)
# global_gene_set_results <- results$global_gene_set_results
# putative <- results$putative
# canonical <- results$canonical


# Example usage:
# GLM_results <- run_GLM_analysis(global_gene_set_results, putative, canonical, cell_counts)
