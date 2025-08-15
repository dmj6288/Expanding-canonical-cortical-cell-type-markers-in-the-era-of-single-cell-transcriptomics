#─────────────────────────────────────────────────────────────────────────────
# 0) Libraries
#─────────────────────────────────────────────────────────────────────────────
library(data.table)    # fast tables
library(Matrix)        # for use_genes issue
library(glmnet)        # regularized (multi)logistic regression
library(mclust)        # for adjustedRandIndex() if you want ARI
library(future)
library(future.apply)
library(SingleCellExperiment)
#plan(multisession, workers = )


# (we’ll compute F1 by hand below)

#─────────────────────────────────────────────────────────────────────────────
# 1) Load data & marker list
#─────────────────────────────────────────────────────────────────────────────
sce_dlPFC  <- readRDS("../sces/RN011/dlPFC_1.rds")
sc_mat     <- counts(sce_dlPFC)
sc_meta    <- colData(sce_dlPFC)
markers    <- unique(readRDS("../rds_files/validation_Gene_lists/dlPFC_subtype_markers_top15_genes_RN012.rds"))

#─────────────────────────────────────────────────────────────────────────────
# 2) Normalize: CP10K → log1p → z-score
#─────────────────────────────────────────────────────────────────────────────
totals  <- colSums(sc_mat)                                 # library size per cell
cp10k   <- t(t(sc_mat) / totals * 1e4)                     # counts-per-10k
log_mat <- log1p(cp10k)                                    # log(CP10K + 1)

# subset to markers & transpose so rows=cells, cols=genes
use_genes <- intersect(markers, rownames(log_mat))
X_all     <- t(log_mat[use_genes, , drop=FALSE])           # now: (n_cells × n_markers)

# z-score each column (gene)
X_all <- scale(X_all)                                      # (center=TRUE, scale=TRUE)


train_ix <- sample(c(1:dim(X_all)[1]), size = 0.7*dim(X_all)[1], replace = FALSE)
test_ix  <- setdiff(c(1:dim(X_all)[1]), train_ix)

X_train  <- X_all[train_ix, , drop=FALSE]
X_test   <- X_all[test_ix,  , drop=FALSE]
y_train  <- sc_meta$cluster_id[train_ix]
y_test   <- sc_meta$cluster_id[test_ix]

alpha <- 0
seed_base <- 42

set.seed(seed_base)
        cvfit <- cv.glmnet(
          x      = X_train, y = y_train,
          family = "multinomial", alpha = alpha,
          type.multinomial = "ungrouped", parallel = FALSE
        )

preds <- as.character(predict(cvfit, X_test, s="lambda.min", type="class"))
acc <- mean(preds == y_test)
classes <- sort(unique(y_test))
f1_vec <- sapply(classes, function(cl) {
          tp <- sum(preds==cl & y_test==cl); fp <- sum(preds==cl & y_test!=cl)
          fn <- sum(preds!=cl & y_test==cl)
          prec <- if ((tp+fp)>0) tp/(tp+fp) else NA_real_
          rec  <- if ((tp+fn)>0) tp/(tp+fn) else NA_real_
          if (!is.na(prec) && !is.na(rec) && (prec+rec)>0) 2*prec*rec/(prec+rec) else NA_real_
        })
result <- data.table(accuracy = acc, macro_F1 = mean(f1_vec, na.rm=TRUE))

#─────────────────────────────────────────────────────────────────────────────
# 5) Collate & visualize
#─────────────────────────────────────────────────────────────────────────────

saveRDS(result, "LOO_classification_dlPFC_putative_RN012.rds")
