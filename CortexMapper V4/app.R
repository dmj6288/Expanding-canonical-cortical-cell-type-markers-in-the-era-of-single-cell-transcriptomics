
options(shiny.fullstacktrace = TRUE)
options(shiny.maxRequestSize = 160*1024^2)
options(shiny.error = traceback)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(dendextend)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(DT)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(shiny)))
suppressMessages(suppressWarnings(library(dplyr)))

cell_names       <- c("Astro", "Excite", "Inhibit", "Micro", "Oligo", "OPC")
regions_to_avoid <- c("ACCK", "ACCV", "M1")

# Define the UI (User Interface)
ui <- fluidPage(
  titlePanel("Yi lab marker generator"),
  
  
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("pvalue", "Maximum Adjusted p-value", min = 0, max = 0.25, value = 1, step = 0.01),
      sliderInput("FC", "Minimum Fold Change", min = 0, max = 4, value = 0, step = 0.1),
      sliderInput("spec", "Minimum Specificity value", min = 0, max = 1, value = 0, step = 0.01),
      textInput("Astro_markers", "Your choice of canonical Astrocyte markers", value = c("ALDH1L1, AQP4, GFAP, GJA1, SLC1A2, SLC4A4"), placeholder = "AQP4"),
      textInput("Excite_markers", "Your choice of canonical Excitatory Neuron markers", value = c("CAMK2A, SLC17A7, SATB2"), placeholder = "CAMK2A"),
      textInput("Inhibit_markers", "Your choice of canonical Inhibitory Neuron markers", value = c("GAD1, GAD2"), placeholder = "GAD1"),
      textInput("Micro_markers", "Your choice of canonical Microglia markers", value = c("AIF1, APBB1IP, CSF1R, CX3CR1, DOCK8, HLA-DRA, P2RY12, PTPRC, TYROBP"), placeholder = "APBB1IP"),
      textInput("Oligo_markers", "Your choice of canonical Oligodendrocyte markers", value = c("MBP, MOBP, MOG, OPALIN, PLP1, ST18"), placeholder = "MBP"),
      textInput("OPC_markers", "Your choice of canonical OPC markers", value = c("CSPG4, PCDH15, PDGFRA, VCAN"), placeholder = "PDGFRA"),
      numericInput("Frequency", "The number of times genes were observed", value = 5, min = 0),
      DTOutput("my_datatable"),
      actionButton("go",label = "Plot Data"),
      uiOutput("radio_matrix")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap",
                 uiOutput("plotContainer")),
         tabPanel("Marker lists",
                  div(
                    class = "table-container",
                    DTOutput("gene_table")
                  )),
        tabPanel("PCA plots",
                 div(
                   class = "table-container",
                   # Adding a text input for the user to enter markers:
                   textInput("markersText", "Enter markers (comma separated):", 
                             value = "GAD1, GAD2"),
                   plotOutput("PCAplot", width = "300px", height = "300px")
                 )),
        tabPanel("Regional gene rank",
                 div(
                   class = "table-container",
                   textInput("rank_genes", "Gene", value = "PDGFRA", placeholder = "PDGFRA"),
                   textInput("rank_cell_type", "Cell Type", value = "OPC", placeholder = "OPC"),
                   DTOutput("ranking_table")
                 ))
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  pct_metrics <- reactive({
    readRDS("source/all_pct_metrics_RN012.rds")
  })
  
  cell_counts <- reactive({
    read_excel("source/Supplementary Table 2B - postQC.xlsx")
  })
  
  filtered_DESeq_matrix <- reactive({
    readRDS("source/DESeq2_normalized_counts.rds")
  })
  
  df_merged <- reactive({
    readRDS("source/merged_gene_info.rds")
  })
  
  all_canonical_gene_list <- reactive({
    readRDS("source/pan_canonical_marker_list.rds")
  })
  
  vst_mat <- reactive({
    readRDS("source/vst_mat.rds")
  })
  
  PCA_seurat <- reactive({
    readRDS("source/final_seurat_object_RN012.rds")
  })
  
  v <- reactiveValues(data = {
    df <- data.frame(
      Astrocyte       = c(80, 15, 2.25),
      Inhibitory      = c(80, 15, 0),
      Excitatory      = c(90,  5, 2),
      Microglia       = c(95,  5, 4.5),
      Oligodendrocyte = c(95, 10, 1),
      OPC             = c(90,  5, 2)
    )
    
    # fidelity_thresh            <- list(80  , 80, 90  , 95  , 95, 90)
    # background_thresh          <- list(15  , 15, 5  , 5  , 10, 5)
    # cohen_d_thresh             <- list(2.25 , 0, 2 , 4.5, 1 , 2)
    
    row.names(df) <- c("Fidelity", "Background", "Cohen distance")
    df
  })
  
  
  #output the datatable based on the dataframe (and make it editable)
  output$my_datatable <- renderDT({
    DT::datatable(v$data, editable = TRUE)
  })
  
  
  
# Define a reactive expression that returns multiple input values as a list
user_inputs <- reactive({

  specificity_vector = input$spec
  log2FC_thresh      = input$FC
  padj_thresh        = input$pvalue

  frequency          = input$Frequency
  
  Astro_markers   = as.character(tstrsplit(input$Astro_markers, ", "))
  Excite_markers  = as.character(tstrsplit(input$Excite_markers, ", "))
  Inhibit_markers = as.character(tstrsplit(input$Inhibit_markers, ", "))
  Micro_markers   = as.character(tstrsplit(input$Micro_markers, ", "))
  Oligo_markers   = as.character(tstrsplit(input$Oligo_markers, ", "))
  OPC_markers     = as.character(tstrsplit(input$OPC_markers, ", "))
  
  rank_genes     = input$rank_genes
  rank_cell_type = input$rank_cell_type
  
  all_canonical_gene_list = list(Astro_markers, Excite_markers, Inhibit_markers, Micro_markers, Oligo_markers, OPC_markers)
  names(all_canonical_gene_list) = cell_names
  
  #print(as.character(tstrsplit(input$Astro_markers, ", ")))
  

  fidelity_thresh = v$data["Fidelity", ]
  background_thresh = v$data["Background", ]
  cohen_d_thresh = v$data["Cohen distance", ]
  
  list(
    specificity_vector = specificity_vector,
    log2FC_thresh      = log2FC_thresh,
    rank_genes         = rank_genes,
    rank_cell_type     = rank_cell_type,
    padj_thresh        = padj_thresh,
    frequency          = frequency,
    fidelity_thresh    = fidelity_thresh,
    background_thresh  = background_thresh,
    cohen_d_thresh     = cohen_d_thresh,
    all_canonical_gene_list = all_canonical_gene_list
  )

})

observeEvent(input$my_datatable_cell_edit, {
  info <- input$my_datatable_cell_edit
  v$data[info$row, info$col] <- as.numeric(info$value)
})

  
cohen_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / pooled_sd
}

server <- function(input, output, session) {
  output$checkbox_matrix <- renderUI({
    rows <- lapply(1:4, function(i) {
      fluidRow(
        lapply(1:5, function(j) {
          column(
            width = 2,  # to fit 5 checkboxes per row (12/5 ≈ 2)
            checkboxInput(
              inputId = paste0("check_", i, "_", j),
              label = paste("C", i, j),
              value = FALSE
            )
          )
        })
      )
    })
    tagList(rows)
  })
}

cell_specific_generate_gene_sets <- reactive({
  
    df_merged               <- df_merged()
    pct_metrics             <- pct_metrics()
    count_matrix            <- filtered_DESeq_matrix()
    
    #all_canonical_gene_list <- all_canonical_gene_list()
    
    ui_vals <- user_inputs()
    
    all_canonical_gene_list <- ui_vals$all_canonical_gene_list
    print(all_canonical_gene_list)
    fidelity_thresh   <- as.numeric(ui_vals$fidelity_thresh)
    background_thresh <- as.numeric(ui_vals$background_thresh)
    padj_thresh       <- as.numeric(ui_vals$padj_thresh)
    log2FC_thresh     <- ui_vals$log2FC_thresh
    cohen_d_thresh    <- ui_vals$cohen_d_thresh
    frequency         <- as.numeric(ui_vals$frequency)  # Frequency is a string from textInput
    
    
    names(fidelity_thresh)     <- cell_names
    names(background_thresh)   <- cell_names
    names(cohen_d_thresh)      <- cell_names
    
    # print(background_thresh[["Astro"]])
    # print(fidelity_thresh[["Astro"]])
    
    # Filter significant genes *per cell type* based on their thresholds
    df_merged_sig <- df_merged %>% rowwise() %>% filter(fidelity       > fidelity_thresh[Cell],
                                                        background     < background_thresh[Cell],
                                                        padj           < padj_thresh,
                                                        log2FoldChange > log2FC_thresh) %>% ungroup()

    df_merged_sig <- df_merged_sig[df_merged_sig$region %in% regions_to_avoid == FALSE, ]
    
    # print(nrow(df_merged_sig))
    # 
    # # Top genes per Cell
    if (nrow(df_merged_sig) == 0) {
      top_genes <- data.frame(Cell = character(), gene = character(), freq = numeric())
    } else {
      top_genes <- df_merged_sig %>% group_by(Cell, gene) %>% summarise(freq = n(), .groups = "drop_last") %>% arrange(Cell, desc(freq)) %>% slice_max(order_by = freq, n = frequency, with_ties = FALSE)
    }
    
    # print(top_genes)

    # # Build subset count 
    subset_count_matrix     <- count_matrix[top_genes$gene, ]
    global_gene_set_results <- list()
    putative                <- data.frame()
    canonical               <- data.frame()
    
    # print(dim(subset_count_matrix))
    
    # Loop through cells
    for (cell in cell_names) {
      
      log2_filtered_TPM_matrix  <- log2(subset_count_matrix + 1)
      current_count_matrix      <- log2_filtered_TPM_matrix
      genes                     <- rownames(current_count_matrix)

      current_count_matrix_cell <- current_count_matrix[, tstrsplit(colnames(current_count_matrix), "_", fixed=TRUE)[[1]] == cell]
      current_count_matrix_BG   <- current_count_matrix[, tstrsplit(colnames(current_count_matrix), "_", fixed=TRUE)[[1]] != cell]
      
      # print(dim(current_count_matrix_cell))
      # print(dim(current_count_matrix_BG))

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

      summary_table      <- summary_table[summary_table$row %in% top_genes[top_genes$Cell == cell, ]$gene, ]
      summary_table$cell <- cell
      
      # Apply Cohen's d filter
      summary_table <- summary_table[summary_table$cohen_d > cohen_d_thresh[[cell]], ]
      global_gene_set_results[[cell]] <- summary_table$row

      #print(head(pct_metrics))
      
      # Collect putative & canonical
      putative  <- rbind(putative,  pct_metrics[pct_metrics$gene %in% global_gene_set_results[[cell]] & pct_metrics$cluster == cell, ])
      canonical <- rbind(canonical, pct_metrics[pct_metrics$gene %in% all_canonical_gene_list[[cell]] & pct_metrics$cluster == cell, ])
    }

    #print(str(gene_list_by_cell))
    
    list(global_gene_set_results = global_gene_set_results, putative = putative, canonical = canonical)

})


output$heatmapPlot <- renderPlot({

  gene_name_lists_by_cell <- cell_specific_generate_gene_sets()
  
  # mat <- matrix(rnorm(36), nrow = 6)
  # rownames(mat) <- paste0("Gene", 1:6)
  # colnames(mat) <- paste0("Sample", 1:6)
  # 
  # # Plot a simple heatmap using base R
  # heatmap(mat, 
  #         main = "Debug Heatmap", 
  #         Colv = NA, 
  #         Rowv = NA, 
  #         scale = "row", 
  #         col = colorRampPalette(c("blue", "white", "red"))(50))
  
  # gene_name_lists_by_cell <- cell_specific_generate_gene_sets()

  # By index:
  global_gene_set_results <- gene_name_lists_by_cell$global_gene_set_results
  all_markers             <- unique(unlist(global_gene_set_results))
  
  vst_mat                <- vst_mat()
  common_gene_set        <- unique(intersect(rownames(vst_mat), all_markers))
  
  plot_matrix            <- vst_mat[common_gene_set, ]
  print(dim(plot_matrix))
  df_row                 <- data.frame(common_gene_set)
  rownames(df_row)       <- common_gene_set
  colnames(df_row)       <- "Marker"
  
  print(dim(df_row))
  
  inverse_gene_cell_hash <- list()

  for (key in names(global_gene_set_results)){
    for (gene in global_gene_set_results[[key]]){
      inverse_gene_cell_hash[[gene]] <- key
    }
  }

  print(length(names(inverse_gene_cell_hash)))
  
  for(marker in df_row$Marker){
    df_row[marker, "color"] <- paste0(inverse_gene_cell_hash[[marker]], " marker")
  }
  df_row <- df_row[-c(1)]

  df_col                                  <- data.frame(tstrsplit(colnames(plot_matrix), "_")[[1]])
  rownames(df_col)                        <- colnames(plot_matrix)
  colnames(df_col)                        <- "Cell Type"

  # # Load RColorBrewer for color palettes
  #
  colors <- c("red", "green", "blue", "orange", "violet", "brown")#brewer.pal(6, "Set1")

  # Define your colors as a named vector
  color_vector <- c(
    "Astro marker"   = colors[1],
    "Excite marker"  = colors[2],
    "Inhibit marker" = colors[3],
    "Micro marker"   = colors[4],
    "Oligo marker"   = colors[5],
    "OPC marker"     = colors[6])

  # Define your colors as a named vector
  cell_vector <- c(
    "Astro"   = colors[1],
    "Excite"  = colors[2],
    "Inhibit" = colors[3],
    "Micro"   = colors[4],
    "Oligo"   = colors[5],
    "OPC"     = colors[6])

  annotation_colors <- list()

  annotation_colors$color       <- as.factor(color_vector)
  annotation_colors$`Cell Type` <- as.factor(cell_vector)
  grid.newpage()  # Start a new graphical page

  pheatmap_obj <- pheatmap(
    plot_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = df_col,
    annotation_row = df_row,
    annotation_colors = annotation_colors,
    fontsize_col = 7,
    fontsize_row = 10,
    fontsize = 25,
    cellwidth = 7,
    cellheight = 10,
    clustering_method = "ward.D2",
    legend = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    angle_col = 90)

  grid.draw(pheatmap_obj)

})


output$plotContainer <- renderUI({
  
  gene_name_lists_by_cell <- cell_specific_generate_gene_sets()
  global_gene_set_results <- gene_name_lists_by_cell$global_gene_set_results
  all_markers             <- unique(unlist(global_gene_set_results))
  
  plotOutput("heatmapPlot", height = paste0(as.character(300 + 10*length(all_markers)), "px"), width = "2800px")
  #plotOutput("heatmapPlot", height = "500px", width = "500px")
  
})

# server.R (or inside your server function)
output$PCAplot <- renderPlot({
    # 1. grab your Seurat object
    so <- PCA_seurat()

    # 2. split & trim your comma‐separated input
    genes <- trimws( unlist(strsplit(input$markersText, ",")) )

    # coords <- Embeddings(so, "pca")[, 1:2]
    # xlim_shared <- range(coords[,1])
    # ylim_shared <- range(coords[,2])

    # 3. make one FeaturePlot per gene
    gene_plots <- lapply(genes, function(g) {
      FeaturePlot(object   = so, features = g, reduction = "pca", pt.size = 1.5, order = TRUE) +
        theme(plot.margin = unit(c(0,0,0,0), "cm"))
    })

    # 4. add your PCA clustering plot
    colors <- c("green", "blue", "red", "violet", "orange", "brown")
    p_cluster <- DimPlot(so, reduction = "pca", group.by  = "cluster_id", pt.size   = 1.5, cols = colors)

    all_plots <- c(gene_plots, list(p_cluster))

    # 5. arrange in a grid
    plot_grid(plotlist = all_plots, ncol = 4)
  },
  # DYNAMIC WIDTH: recalc whenever markersText changes
  width = 1350,
  height = function(){
    n <- length(trimws(unlist(strsplit(input$markersText, ","))))
    # e.g. allocate 600px per column (adjust to taste)
    (round((n)/4) + 1)*250
  }
)

# server.R (or inside your server function)
output$barplot <- renderPlot({
  
  GLM <- run_GLM_analysis()
  
  GLM_results <- GLM$GLM_results
  
  ggplot(GLM_results, aes(cluster, odds_ratio, fill = metric)) + geom_bar(stat="identity", position = "dodge", width = 0.65)       +
    labs(title="Multiple Bar plots")                    +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "tomato", linewidth = 1.5)       +
    scale_fill_manual(values = c("darkblue", "orange")) +
    xlab("Cell Type") + ylab("Odds Ratio")              +
    theme_bw()                                          +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))},
# DYNAMIC WIDTH: recalc whenever markersText changes
    width = 1200, height = 600
)


table_data <- reactive({

  gene_name_lists_by_cell <- cell_specific_generate_gene_sets()
  global_gene_set_results <- gene_name_lists_by_cell$global_gene_set_results
  all_markers             <- unique(unlist(global_gene_set_results))

  df <- data.frame(
    Cell_Type = names(global_gene_set_results),
    Genes = sapply(global_gene_set_results, paste, collapse = ", "),
    Count = sapply(global_gene_set_results, length)# Convert to comma-separated strings
  )

  return(df)  # Return the dataframe
})

canonical_gene_table_data <- reactive({
  
  ui_vals <- user_inputs()
  
  canonical_gene_set <- ui_vals$all_canonical_gene_list
  
  print(canonical_gene_set)
  #global_gene_set_results <- gene_name_lists_by_cell$global_gene_set_results
  #all_markers             <- unique(unlist(global_gene_set_results))
  
  df <- data.frame(
    Cell_Type = names(canonical_gene_set),
    Genes = sapply(canonical_gene_set, paste, collapse = ", "),
    Count = sapply(canonical_gene_set, length)# Convert to comma-separated strings
  )
  
  return(df)  # Return the dataframe
})


rank_table <- reactive({
  
  ui_vals        <- user_inputs()
  
  canonical_gene_set <- ui_vals$all_canonical_gene_list
  print(canonical_gene_set)
  
  
  df_merged      <- df_merged()
  
  #print(ui_vals)
  
  rank_genes     <- ui_vals$rank_genes
  rank_cell_type <- ui_vals$rank_cell_type
  
  print(rank_genes)
  print(rank_cell_type)
  
  genes <- as.character(tstrsplit(rank_genes, ", "))
  print(genes)
  
  display_table <- df_merged[df_merged$gene %in% genes & df_merged$cluster == rank_cell_type, ]
  display_table <- display_table[order(display_table$gene, decreasing = FALSE), ]
  
  return(display_table)  # Return the dataframe
})


output$gene_table <- renderDT({

  datatable(table_data(),
            options = list(autoWidth = TRUE, scrollX = TRUE, pageLength = 100),
            escape = FALSE) %>%
    formatStyle(
      columns = "Genes", # Replace with your column name
      whiteSpace = "normal",
      wordWrap = "break-word"
    )
})

output$Canonical_marker_table <- renderDT({
  
  datatable(canonical_gene_table_data(),
            options = list(autoWidth = TRUE, scrollX = TRUE),
            escape = FALSE) %>%
    formatStyle(
      columns = "Genes", # Replace with your column name
      whiteSpace = "normal",
      wordWrap = "break-word"
    )
})

output$ranking_table <- renderDT({
  
  print(rank_table())
  datatable(rank_table(),
            options = list(autoWidth = TRUE, scrollX = TRUE, pageLength = 100),
            escape = FALSE)
})


}

shinyApp(ui = ui, server = server)
