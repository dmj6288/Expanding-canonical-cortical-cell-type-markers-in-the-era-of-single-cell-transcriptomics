suppressWarnings(suppressMessages(library(data.table)))

# INPUT

res_table_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/pbres/RN012/"

# OUTPUT

gene_list_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/MarkerPaper/gene_lists/RN012/"


PC_genes       <- readRDS("rds_files/human_PC_genes.rds")

all_files                <- list.files(res_table_path)
all_file_list            <- list()
all_file_list_considered <- list()
all_result_tables        <- list()

for (file_name in all_files){

	current_file                      <- readRDS(paste0(res_table_path, file_name))
	current_file_result_table         <- list()
	current_file_gene_list            <- list()
	current_file_considered_gene_list <- list()

	for (cell in names(current_file)){

		current_gene_set                          <- na.omit(current_file[[cell]])
		current_file_result_table[[cell]]         <- current_gene_set
		current_file_considered_gene_list[[cell]] <- rownames(current_gene_set)
		current_file_gene_list[[cell]]            <- rownames(current_gene_set[!is.na(current_gene_set$padj) & current_gene_set$padj < 0.05 & current_gene_set$log2FoldChange > 0, ])
		print("ROWNAMES")
		current_file_considered_gene_list[[cell]] <- current_file_considered_gene_list[[cell]][current_file_considered_gene_list[[cell]] %in% PC_genes]
                current_file_gene_list[[cell]]            <- current_file_gene_list[[cell]][current_file_gene_list[[cell]] %in% PC_genes]

	}

	region <- tstrsplit(file_name, ".rds", fixed=TRUE)[[1]]
	all_file_list[[region]]            <- current_file_gene_list
	all_file_list_considered[[region]] <- current_file_considered_gene_list
	all_result_tables[[region]]        <- current_file_result_table

}

#saveRDS(all_file_list,            paste0(gene_list_path, "DE_genesFC/DE_genesFC_level0.rds"))
saveRDS(all_result_tables,        paste0(gene_list_path, "DE_results_RN012.rds"))

#saveRDS(all_file_list_considered, paste0(gene_list_path, "Considered_FC/", region, ".rds"))
