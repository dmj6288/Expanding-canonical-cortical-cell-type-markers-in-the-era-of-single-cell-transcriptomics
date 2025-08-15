#─────────────────────────────────────────────────────────────────────────────
# 0) Libraries
#─────────────────────────────────────────────────────────────────────────────
library(data.table)    # fast tables
library(Matrix)        # for use_genes issue
library(glmnet)        # regularized (multi)logistic regression
library(mclust)        # for adjustedRandIndex() if you want ARI
library(future)
library(future.apply)
plan(multisession, workers = 12)


# (we’ll compute F1 by hand below)

#─────────────────────────────────────────────────────────────────────────────
# 1) Load data & marker list
#─────────────────────────────────────────────────────────────────────────────
sc_mat   <- readRDS("../single_cell_global/all_counts_RN012.rds")       # rows=genes, cols=cells
sc_meta  <- readRDS("../single_cell_global/all_metadata_RN012.rds")     # must contain 'cell_type' & 'study'
markers  <- Reduce(union, readRDS("../rds_files/validation_Gene_lists/pan_canonical_marker_list.rds"))

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

#─────────────────────────────────────────────────────────────────────────────
# 3) Prepare for LOSO
#─────────────────────────────────────────────────────────────────────────────
studies     <- unique(sc_meta$region)
n_studies   <- length(studies)
results     <- vector("list", n_studies)

#─────────────────────────────────────────────────────────────────────────────
# 4) LOSO loop: train multinomial logistic & evaluate
#─────────────────────────────────────────────────────────────────────────────
for (i in seq_along(studies)) {
  s        <- studies[i]
  train_ix <- sc_meta$region != s
  test_ix  <- !train_ix

  X_train  <- X_all[train_ix, , drop=FALSE]
  X_test   <- X_all[test_ix,  , drop=FALSE]
  y_train  <- sc_meta$cluster_id[train_ix]
  y_test   <- sc_meta$cluster_id[test_ix]

  # — 4.1 Fit regularized multinomial logistic (alpha=0 for ridge; use alpha=1 for lasso)
  set.seed(42 + i)  # reproducibility
  cvfit <- cv.glmnet(
    x      = X_train,
    y      = y_train,
    family = "multinomial",
    alpha  = 0,
    type.multinomial = "ungrouped",
    parallel = FALSE         # turn TRUE if you have a parallel backend
  )

  # — 4.2 Predict on held-out fold
  preds <- predict(
    cvfit,
    newx = X_test,
    s    = "lambda.min",
    type = "class"
  )
  preds <- as.character(preds)  # a factor matrix → char vector

  # — 4.3 Compute accuracy
  acc <- mean(preds == y_test)

  # — 4.4 Compute per-class precision/recall/F1, then macro-F1
  classes <- sort(unique(y_test))
  f1_vec  <- sapply(classes, function(cl) {
    tp <- sum(preds == cl & y_test == cl)
    fp <- sum(preds == cl & y_test != cl)
    fn <- sum(preds != cl & y_test == cl)

    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) {
      2 * prec * rec / (prec + rec)
    } else {
      NA_real_
    }
  })
  macro_f1 <- mean(f1_vec, na.rm = TRUE)

  # — 4.5 Store results
  results[[i]] <- data.table(
    study    = s,
    accuracy = acc,
    macro_F1 = macro_f1
  )
}

#─────────────────────────────────────────────────────────────────────────────
# 5) Collate & visualize
#─────────────────────────────────────────────────────────────────────────────
df_results <- rbindlist(results)
print( df_results )
saveRDS(df_results, "results/RN012/classification_pan_canonical_RN012.rds")

# Simple boxplots across the 19 folds:
par(mfrow = c(1,2))
boxplot(df_results$accuracy, main = "LOSO Accuracy", ylab = "Accuracy")
boxplot(df_results$macro_F1, main = "LOSO Macro-F1", ylab = "Macro-F1")

# You can also run a paired test against another marker set, e.g.:
# wilcox.test(macro_F1 ~ marker_set, data = df_results, paired = TRUE)
