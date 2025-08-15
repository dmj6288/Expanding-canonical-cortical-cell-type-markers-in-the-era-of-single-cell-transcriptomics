suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(DESeq2)))

# INPUT
pseudo_bulk_ready_data_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/DEX_data/RN012/"
DESeq_coldata_path          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/DEX_cols/RN012/"

# OUTPUT
output_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/pbres/RN012/"

# Slurm array task index
file_idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
dex_to_analyze <- list.files(pseudo_bulk_ready_data_path)

# Safety check
if (file_idx > length(dex_to_analyze) | file_idx < 1) {
    stop("SLURM_ARRAY_TASK_ID is out of range.")
}

file_name <- dex_to_analyze[file_idx]
cat("Processing file:", file_name, "\n")

DE_matrix  <- readRDS(paste0(pseudo_bulk_ready_data_path, file_name))
DE_coldata <- readRDS(paste0(DESeq_coldata_path, file_name))

contrast <- "compare_B_vs_A"
result_table <- list()

for (cell in names(DE_matrix)) {
    cat("Running DE for cell type:", cell, "\n")

    count_matrix <- as.matrix(DE_matrix[[cell]])
    counts_filtered <- count_matrix[rowSums(count_matrix >= 10) >= (ncol(count_matrix) / 2), ]

    if (nrow(counts_filtered) >= 2) {  # Add a safeguard to avoid empty matrices
        dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                      colData = DE_coldata[[cell]],
                                      design = ~ compare)
        dds <- DESeq(dds)
        res_shrink <- lfcShrink(dds, coef="compare_B_vs_A", type="ashr")
        result_table[[cell]] <- res_shrink
    } else {
        cat("Skipping cell type:", cell, " (not enough genes)\n")
    }
}

saved_file <- tstrsplit(file_name, ".rds", fixed = TRUE)[[1]]
saveRDS(result_table, paste0(output_path, saved_file, ".rds"))

cat("✅ Finished:", file_name, "\n")
