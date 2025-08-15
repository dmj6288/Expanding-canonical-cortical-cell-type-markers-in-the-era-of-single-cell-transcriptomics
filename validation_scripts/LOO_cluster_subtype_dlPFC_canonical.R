suppressMessages(suppressWarnings(library(FNN)))      # for fast kNN
suppressMessages(suppressWarnings(library(igraph)))   # for community detection
suppressMessages(suppressWarnings(library(future)))
suppressMessages(suppressWarnings(library(future.apply)))
suppressMessages(suppressWarnings(library(data.table)))   # for rbindlist()
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(dbscan)))

#plan(multisession, workers = 16)

#sc_mat  <- readRDS("../single_cell_global/all_counts.rds")

sce_dlPFC  <- readRDS("../sces/RN011/dlPFC_1.rds")
#raw_counts <- readRDS("../single_cell_global/all_counts_RN012.rds")
raw_counts <- counts(sce_dlPFC)

# compute library-size normalized, log1p counts:
lib_sizes <- colSums(raw_counts)
lognorm   <- t( t(raw_counts) / lib_sizes * 1e4 )   # counts per 10k
lognorm   <- log2(lognorm + 1)

sc_mat    <- lognorm

#sc_meta <- readRDS("../single_cell_global/all_metadata_RN012.rds")
sc_meta <- colData(sce_dlPFC)
#markers <- Reduce(union, readRDS("../rds_files/validation_Gene_lists/fidelity_bf_informed_gene_list_RN012.rds"))
markers <- readRDS("../rds_files/validation_Gene_lists/canonical_dlPFC_subtypes_RN012.rds")

use_genes <- intersect(unlist(markers), rownames(sc_mat))
X         <- t(sc_mat[use_genes, ])

rm(sc_mat)

studies    <- unique(sc_meta$region)
n_clusters <- length(unique(sc_meta$cell_type))

X_boot   <- X[, use_genes]
# 1) split train/test
#train_idx <- sc_meta$region != s
train_idx <- sample(c(1:dim(X_boot)[1]), size = 0.7*dim(X_boot)[1], replace = FALSE)
test_idx  <- setdiff(c(1:dim(X_boot)[1]), train_idx)
#test_idx  <- !train_idx

X_train   <- X_boot[train_idx, , drop = FALSE]
X_test    <- X_boot[test_idx,  , drop = FALSE]

#knn_res   <- get.knn(X_train, k = 20)    # choose ~10–30 neighbors
#edges     <- do.call(rbind, lapply(seq_len(nrow(knn_res$nn.index)), function(i) {
#  cbind(rep(i, ncol(knn_res$nn.index)), knn_res$nn.index[i,])
#}))
#g         <- graph_from_edgelist(edges, directed = FALSE)
#g         <- simplify(g)                       # collapse duplicate edges

# 3) run Louvain community detection
#comm      <- cluster_louvain(g)
#cl_train  <- membership(comm)

pca_res <- prcomp(X_train, center=TRUE, scale.=TRUE)
X_pca10 <- pca_res$x[,1:20]

# 2) HDBSCAN
hdb <- hdbscan(X_pca10, minPts = 20)
cl_train <- hdb$cluster   # 0 = noise
clusters  <- sort(unique(cl_train))

# 4) compute centroids and assign held-out cells (same as before)
centroids   <- t(sapply(clusters, function(cl) {colMeans( X_train[ cl_train == cl, , drop = FALSE ] )}))
rownames(centroids) <- clusters
assign_test <- apply(X_test, 1, function(x) which.min(colSums((centroids - x)^2)))

# 5) ARI
true_labels <- as.integer(factor(sc_meta$cluster_id[test_idx]))
ari_val     <- adjustedRandIndex(true_labels, assign_test)

df <- data.table(ari = ari_val)

# then:

#results_list <- future_lapply(c(1:1000), loo_fold_louvain, future.seed = TRUE)

#sc_results   <- rbindlist(results_list)
saveRDS(df, "results/RN012/clustering_dlPFC_canonical_subtype_RN012.rds")


