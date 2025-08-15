suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(S4Vectors)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(comprehenr)))

support_functions_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

# INPUT
sce_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/sces/RN012/"

# OUTPUT
pseudo_bulk_ready_data <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/DEX_data/RN012/"
DESeq_coldata          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/DEX_cols/RN012/"
coldata_path           <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/coldata/RN012/"
aggmatrix_path         <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/Aggregate_matrices/RN012/"

PC_genes <- readRDS("rds_files/human_PC_genes.rds")

source(paste0(support_functions_path, "aggregate_function_Matrix.Utils.R"))

# Get array index
file_idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sce_folders_to_analyze <- list.files(sce_path)

if (file_idx > length(sce_folders_to_analyze) | file_idx < 1) {
    stop("SLURM_ARRAY_TASK_ID is out of range.")
}

sce_folder <- sce_folders_to_analyze[file_idx]
cat("Processing file:", sce_folder, "\n")

sce <- readRDS(paste0(sce_path, sce_folder))
saveRDS(colData(sce), paste0(coldata_path, sce_folder))

groups <- colData(sce)[, c("cluster_id", "sample_id")]
groups$cluster_id <- as.factor(groups$cluster_id)
groups$sample_id <- as.factor(groups$sample_id)

levels(groups$cluster_id) <- sort(levels(groups$cluster_id))
cluster_names <- levels(groups$cluster_id)
sample_names <- levels(groups$sample_id)

aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")
aggr_counts <- t(aggr_counts)

pseudo_groups <- as.data.frame(sub("^[^_]+_", "", colnames(aggr_counts)))

#pseudo_groups <- as.data.frame(sub(".*_#?", "", colnames(aggr_counts)))

colnames(pseudo_groups) <- c("Sample")
rownames(pseudo_groups) <- colnames(aggr_counts)

saveRDS(aggr_counts, paste0(aggmatrix_path, sce_folder))

counts_ls <- list()
metadata_ls <- list()

for (cell in cluster_names) {

    #current_cluster <- aggr_counts[rownames(aggr_counts) %in% PC_genes,
    #                               tstrsplit(colnames(aggr_counts), "_#", fixed = TRUE)[[1]] %in% cell]

    current_cluster <- aggr_counts[rownames(aggr_counts) %in% PC_genes,
                                   tstrsplit(colnames(aggr_counts), "_", fixed = TRUE)[[1]] %in% cell]

    if (dim(as.matrix(current_cluster))[2] > 1) {

        #background_clusters <- aggr_counts[rownames(aggr_counts) %in% PC_genes,
        #                                   !(tstrsplit(colnames(aggr_counts), "_#", fixed = TRUE)[[1]] %in% cell)]

	background_clusters <- aggr_counts[rownames(aggr_counts) %in% PC_genes,
                                           !(tstrsplit(colnames(aggr_counts), "_", fixed = TRUE)[[1]] %in% cell)]

        background <- pseudo_groups[!(rownames(pseudo_groups) %in% colnames(current_cluster)), , drop = FALSE]

        colnames(background) <- "Sample"
        rownames(background) <- colnames(background_clusters)

        aggr_counts_background <- t(aggregate.Matrix(t(background_clusters), groupings = background, fun = "sum"))
        DE_matrix <- cbind(current_cluster, aggr_counts_background)

        metadata <- data.frame(Sample = c(colnames(current_cluster), colnames(aggr_counts_background)),
                               compare = c(rep("B", ncol(current_cluster)), rep("A", ncol(aggr_counts_background))))

        counts_ls[[cell]] <- DE_matrix
        metadata_ls[[cell]] <- metadata
    }
}

file_name <- tstrsplit(sce_folder, ".rds", fixed = TRUE)[[1]]

saveRDS(counts_ls,   paste0(pseudo_bulk_ready_data, file_name, ".rds"))
saveRDS(metadata_ls, paste0(DESeq_coldata, file_name, ".rds"))

cat("✅ Finished:", sce_folder, "\n")
