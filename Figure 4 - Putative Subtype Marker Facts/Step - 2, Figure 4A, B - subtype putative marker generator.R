rm(list=ls())

suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(ggplot2)))


full_count_matrix          <- readRDS("interm/sample_wise40_count_matrix_RAW_RN011.rds")
DE_tables_with_specificity <- readRDS("interm/final_DE_table_with_DESeq2_specificity_and_fidelity_RN011.rds")

macro_names_table  <- read.csv("source/extended_subtype_names.csv")

macro_names        <- macro_names_table$macro
names(macro_names) <- macro_names_table$x

matrix_columns             <- tstrsplit(colnames(full_count_matrix), "_", fixed=TRUE)[[1]]
n_clusters                 <- length(unique(matrix_columns))

all_gene_lists             <- c()

for (region in names(DE_tables_with_specificity)){
  
  for (cell_name in names(DE_tables_with_specificity[[region]])){
    
    df <- DE_tables_with_specificity[[region]][[cell_name]]
    df <- df[df$padj < 0.05, ]
    #df <- head(df[order(df$fidelity, decreasing = TRUE), ], 1000)
    df <- df[df$background < 15 & df$fidelity > 80, ]
    #df <- df[df$specificity_values > 0.65, ]
    df <- df[df$log2FoldChange > 2, ]
    
    df <- df[order(df$fidelity, decreasing = TRUE), ]
    df <- head(df, 15)
    
    if (dim(df)[1] > 0){
      
      df$DEGs <- rownames(df)
      df$Region <- region
      
      all_gene_lists <- rbind(all_gene_lists, df)
      
    }
    
  }
  
}



# > identical(subtype_DEG_table_top15$cluster, subtype_DEG_table_top15_without_specificity$cluster)
# [1] TRUE
# > identical(subtype_DEG_table_top15$fidelity, subtype_DEG_table_top15_without_specificity$fidelity)
# [1] TRUE
# > identical(subtype_DEG_table_top15$background, subtype_DEG_table_top15_without_specificity$background)
# [1] TRUE
# > identical(subtype_DEG_table_top15$specificity_values, subtype_DEG_table_top15_without_specificity$specificity_values)
# [1] TRUE
# > identical(subtype_DEG_table_top15$log2FoldChange, subtype_DEG_table_top15_without_specificity$log2FoldChange)
# [1] TRUE
# > identical(subtype_DEG_table_top15$padj, subtype_DEG_table_top15_without_specificity$padj)
# [1] TRUE
# > identical(subtype_DEG_table_top15$DEGs, subtype_DEG_table_top15_without_specificity$DEGs)
# [1] TRUE

all_gene_lists$macro <- macro_names[as.character(all_gene_lists$cluster)]

all_gene_lists$Region <- tstrsplit(all_gene_lists$Region, "_")[[1]]
#all_gene_lists <- all_gene_lists[all_gene_lists$macro %in% c("Inhibit", "Excite"), ]

write.csv(all_gene_lists, "interm/subtype_DEG_table_top15_without_specificity.csv")


library(dplyr)

summary_df <- all_gene_lists %>%
  group_by(Region, macro) %>%
  summarize(
    fidelity_mean = mean(fidelity, na.rm = TRUE),
    fidelity_sd = sd(fidelity, na.rm = TRUE),
    background_mean = mean(background, na.rm = TRUE),
    background_sd = sd(background, na.rm = TRUE),
  #  log2FC_mean = mean(log2FoldChange, na.rm = TRUE),
  #  log2FC_sd = sd(log2FoldChange, na.rm = TRUE),
    .groups = 'drop'
  )

library(tidyr)


plot_df <- summary_df %>%
  pivot_longer(cols = ends_with("_mean"), 
               names_to = "metric", 
               values_to = "mean") %>%
  mutate(
    metric = gsub("_mean", "", metric),
    sd = case_when(
      metric == "fidelity" ~ fidelity_sd,
      metric == "background" ~ background_sd,
#      metric == "log2FC" ~ log2FC_sd
    )
  )


fidelity_df <- plot_df[plot_df$metric == "fidelity", ]
p <- ggplot(data = fidelity_df, aes(x = macro, y = mean, fill = macro)) + geom_boxplot(linewidth = 0.5)  + theme_bw() +
  ylab("") +
  xlab("Metric") +
  theme(
    axis.text = element_text(size = 7),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("red", "green", "blue", "orange", "violet", "brown"))

ggsave(filename = "figures/Region_and_subtype_averaged_fidelity.png", plot = p, width = 960, height = 720, units = "px", device = "png")

bg_df <- plot_df[plot_df$metric == "background", ]
p <- ggplot(data = bg_df, aes(x = macro, y = mean, fill = macro)) + geom_boxplot(linewidth = 0.5)  + theme_bw() +
  ylab("") +
  xlab("Metric") +
  theme(
    axis.text = element_text(size = 7),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("red", "green", "blue", "orange", "violet", "brown"))

ggsave(filename = "figures/Region_and_subtype_averaged_background.png", plot = p, width = 960, height = 720, units = "px", device = "png")


# library(ggplot2)
# 
# plot_df_list <- list()
# all_fidelity_metrics <- c()
# 
# for (cell in unique(plot_df$macro)){
#   
#   current_plot_df <- plot_df[plot_df$macro %in% cell, ]
# 
#   p <- ggplot(current_plot_df, aes(x = Region, y = mean, fill = metric)) +
#     geom_bar(stat = "identity", position = "dodge", width = 0.5) +
#     geom_errorbar(aes(ymin = mean - 0.5*sd, ymax = mean + 0.5*sd), 
#                   position = position_dodge(0.45), 
#                   width = 0.3) +
#     scale_fill_manual(values = c("darkblue", "orange"))        +
#     labs(y = "Mean Value") + theme_bw() +  guides(fill="none") + 
#     ylim(c(0, 100)) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15), 
#           axis.text.y = element_text(size = 15), 
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank())
#   p  
#   
#   plot_df_list[[cell]] <- current_plot_df
#   
#   ggsave(filename = paste0("figures/", cell, "_subtype_fidelity_plot_top15.png"), plot = p, height = 1400, width = 3000, device = "png", units = "px")
#   
# }
# 
# 





