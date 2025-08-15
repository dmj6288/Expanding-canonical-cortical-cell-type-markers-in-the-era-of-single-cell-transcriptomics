suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(S4Vectors)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(comprehenr)))

support_functions_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

source(paste0(support_functions_path, "counts_to_tpm.R"))
# INPUT

sce_path               <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/sces/RN012/"
considered_genes       <- readRDS("rds_files/common_genes_agg_RN012.rds")
#gene_lengths           <- readRDS("rds_files/gene_lengths_RN009.rds")

sce_folders_to_analyze   <- list.files(sce_path)
full_sc_metadata         <- c()
counts_list              <- list()

for (sce_folder in sce_folders_to_analyze){

	sce <- readRDS(paste(sce_path, sce_folder, sep =""))

	print("----READ COMPLETE----")
	
	groups <- colData(sce)[, c("cluster_id", "sample_id", "human_age", "sex")]
	groups$region <- tstrsplit(sce_folder, "_1")[[1]]

	full_sc_metadata <- rbind(full_sc_metadata, groups)
	current_counts   <- counts(sce)	

	counts_list[[sce_folder]] <- current_counts[considered_genes, ]
	
}

all_counts <- do.call(cbind, counts_list)

print(dim(all_counts))

#subset_gene_length <- gene_lengths[gene_lengths$hgnc_symbol %in% rownames(all_counts), ]
#all_counts <- all_counts[subset_gene_length$hgnc_symbol, ]

#full_TPM_matrix <- counts_to_tpm(all_counts, subset_gene_length$Length)

saveRDS(all_counts      , "single_cell_global/all_counts_RN012.rds")
#saveRDS(full_TPM_matrix,  "single_cell_global/all_counts_TPM.rds")
saveRDS(full_sc_metadata, "single_cell_global/all_metadata_RN012.rds")
