
library(dplyr)


can_class <- readRDS("results/RN012/classification_pan_canonical_RN012.rds")
put_class <- readRDS("results/RN012/classification_putative_frequent_RN012.rds")

can_cluster <- readRDS("results/RN012/clustering_pan_canonical_RN012.rds")
put_cluster <- readRDS("results/RN012/clustering_putative_frequent_RN012.rds")

# Put all tables into a list
tables <- list(put_class, can_class, put_cluster, can_cluster)

# Merge by study
merged_df <- Reduce(function(x, y) merge(x, y, by = "study", all = TRUE), tables)

colnames(merged_df) <- c("study",  "putative_accuracy", "putative_F1", "canonical_accuracy", "canonical_F1", "putative_ari", "canonical_ari")

saveRDS(merged_df, "results/RN012/total_pan_results_RN012.rds")
