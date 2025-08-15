suppressMessages(suppressWarnings(library(FNN)))      # for fast kNN
suppressMessages(suppressWarnings(library(igraph)))   # for community detection
suppressMessages(suppressWarnings(library(future)))
suppressMessages(suppressWarnings(library(future.apply)))
suppressMessages(suppressWarnings(library(data.table)))   # for rbindlist()
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(Matrix)))

plan(multisession, workers = 24)

#sc_mat  <- readRDS("../single_cell_global/all_counts.rds")

raw_counts <- readRDS("../single_cell_global/all_counts_RN012.rds")

# compute library-size normalized, log1p counts:
lib_sizes <- colSums(raw_counts)
lognorm   <- t( t(raw_counts) / lib_sizes * 1e4 )   # counts per 10k
lognorm   <- log2(lognorm + 1)

sc_mat    <- lognorm

sc_meta <- readRDS("../single_cell_global/all_metadata_RN012.rds")
markers <- Reduce(union, readRDS("../rds_files/validation_Gene_lists/putative_gene_set_region_frequent.rds"))

use_genes <- intersect(markers, rownames(sc_mat))
X         <- t(sc_mat[use_genes, ])

rm(sc_mat)

studies    <- unique(sc_meta$region)
n_clusters <- length(unique(sc_meta$cell_type))

loo_fold_louvain <- function(s) {

  print(s)
  # 1) split train/test
  train_idx <- sc_meta$region != s
  test_idx  <- !train_idx
  
  X_train   <- X[train_idx, , drop = FALSE]
  X_test    <- X[test_idx,  , drop = FALSE]

  # 2) build kNN graph on training cells
  knn_res   <- get.knn(X_train, k = 20)    # choose ~10–30 neighbors
  edges     <- do.call(rbind, lapply(seq_len(nrow(knn_res$nn.index)), function(i) {
    cbind(rep(i, ncol(knn_res$nn.index)), knn_res$nn.index[i,])
  }))
  g         <- graph_from_edgelist(edges, directed = FALSE)
  g         <- simplify(g)                       # collapse duplicate edges

  # 3) run Louvain community detection
  comm      <- cluster_louvain(g)
  cl_train  <- membership(comm)
  clusters  <- sort(unique(cl_train))
 
  # 4) compute centroids and assign held-out cells (same as before)
  centroids   <- t(sapply(clusters, function(cl) {colMeans( X_train[ cl_train == cl, , drop = FALSE ] )}))
  rownames(centroids) <- clusters
  assign_test <- apply(X_test, 1, function(x) which.min(colSums((centroids - x)^2)))

  # 5) ARI
  true_labels <- as.integer(factor(sc_meta$cluster_id[test_idx]))
  ari_val     <- adjustedRandIndex(true_labels, assign_test)

  data.table(study = s, ari = ari_val)
}

# then:

results_list <- future_lapply(studies, loo_fold_louvain, future.seed = TRUE)

sc_results   <- rbindlist(results_list)
saveRDS(sc_results, "results/RN012/clustering_putative_frequent_RN012.rds")


