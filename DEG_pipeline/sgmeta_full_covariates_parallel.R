suppressMessages(suppressWarnings(library(stringr)))
print("STRINGR IMPORTED")
suppressMessages(suppressWarnings(library(Seurat)))
#suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(hash)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

support_functions <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

metadata_path                <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/rds_files/RN012/"
seurat_objects_path          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/QCfiltered_data/RN012/"
sce_path                     <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/sces/RN012/"

PC_genes                     <- readRDS("rds_files/human_PC_genes.rds")
seurat_objects               <- list.files(metadata_path)

# Get array index
file_idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
seurat_objects_to_analyze <- seurat_objects

if (file_idx > length(seurat_objects_to_analyze) | file_idx < 1) {
    stop("SLURM_ARRAY_TASK_ID is out of range.")
}
 
study_metadata <- seurat_objects_to_analyze[file_idx]
cat("Processing file:", study_metadata, "\n")

filtered_seurat_object <- readRDS(paste0(seurat_objects_path, study_metadata))

filtered_seurat_object <- UpdateSeuratObject(filtered_seurat_object)

print(str(filtered_seurat_object))

print(paste(metadata_path, "/", study_metadata, sep =""))


df                           <- readRDS(paste(metadata_path, "/", study_metadata, sep =""))

filtered_seurat_object_with_identified_cells <- subset(filtered_seurat_object, cells = rownames(df))

df                           <- df[colnames(filtered_seurat_object_with_identified_cells), ]

print(identical(colnames(filtered_seurat_object_with_identified_cells), rownames(df)))

filtered_seurat_object_with_identified_cells@meta.data$orig.ident <- df$sample_id
filtered_seurat_object_with_identified_cells@meta.data$human_age  <- df$age
filtered_seurat_object_with_identified_cells@meta.data$sex        <- df$sex

filtered_seurat_object_with_identified_cells <- SetIdent(filtered_seurat_object_with_identified_cells, value = df$trad)

#all.pct <- FindAllMarkers(
#  object       = filtered_seurat_object_with_identified_cells,
#  assay        = "RNA",
#  slot         = "data",          # or "counts" if you really want raw
#  only.pos     = FALSE,           # so that it returns both pos+neg
#  min.pct      = 0,               # include every gene
#  logfc.threshold = 0             # include every gene
#)
# assume seurat_obj is already loaded

# 1) Pull out the raw counts matrix
counts_mat <- GetAssayData(filtered_seurat_object_with_identified_cells, assay = "RNA", slot = "counts")

# 2) Grab your cluster assignments
clusters <- Idents(filtered_seurat_object_with_identified_cells)

# 3) For each cluster, compute pct1/pct2 from counts>0
pct_list <- lapply(levels(clusters), function(cl) {
  # cells in vs. out of this cluster
  cells_in  <- WhichCells(filtered_seurat_object_with_identified_cells, idents = cl)
  cells_out <- setdiff(colnames(counts_mat), cells_in)
  
  # logical matrix: TRUE if count > 0
  expr_pos <- counts_mat > 0
  
  # fraction of cells with count>0 in vs. out (×100 for percent)
  pct1 <- Matrix::rowSums(expr_pos[, cells_in, drop = FALSE])  / length(cells_in)  * 100
  pct2 <- Matrix::rowSums(expr_pos[, cells_out, drop = FALSE]) / length(cells_out) * 100
  
  data.frame(
    gene    = rownames(counts_mat),
    cluster = cl,
    pct1    = pct1,
    pct2    = pct2,
    row.names = NULL
  )
})

# 4) Combine into one data.frame
pct_df <- do.call(rbind, pct_list)

# quick look
head(pct_df)

saveRDS(pct_df, paste0("fidelity_metrics/RN012/", study_metadata))

# print(str(filtered_seurat_object_with_identified_cells))

# FILTERING THE DATA

counts                       <- filtered_seurat_object_with_identified_cells@assays$RNA@layers$counts

colnames(counts)             <- colnames(filtered_seurat_object_with_identified_cells)
rownames(counts)             <- rownames(filtered_seurat_object_with_identified_cells)

PC_counts                    <- counts[rownames(counts) %in% PC_genes, ]
metadata                     <- filtered_seurat_object_with_identified_cells@meta.data


# print("------FILTERED SEURAT STRUCTURE------")
# print(str(filtered_seurat_object_with_identified_cells))


metadata$cluster_id          <- df[colnames(PC_counts),"trad"]
metadata$sample_id           <- df[colnames(PC_counts),"sample_id"]

sce                          <- SingleCellExperiment(assays = list(counts = PC_counts), colData = metadata)

file_name <- tstrsplit(study_metadata, ".rds", fixed=TRUE)[[1]]
saveRDS(sce, paste(sce_path, file_name, "_1.rds", sep =""))
