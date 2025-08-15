#─────────────────────────────────────────────────────────────────────────────
# 0) Libraries
#─────────────────────────────────────────────────────────────────────────────
library(data.table)    # fast tables
library(Matrix)        # for use_genes issue
library(glmnet)        # regularized (multi)logistic regression
library(mclust)        # for adjustedRandIndex() if you want ARI
library(future)
library(future.apply)
#plan(multisession, workers = )


# (we’ll compute F1 by hand below)

#─────────────────────────────────────────────────────────────────────────────
# 1) Load data & marker list
#─────────────────────────────────────────────────────────────────────────────
sc_mat   <- readRDS("../single_cell_global/all_counts_RN012.rds")       # rows=genes, cols=cells
sc_meta  <- readRDS("../single_cell_global/all_metadata_RN012.rds")     # must contain 'cell_type' & 'study'
markers  <- Reduce(union, readRDS("../rds_files/validation_Gene_lists/putative_gene_set_region_frequent.rds"))

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

print("Data Ready")

#─────────────────────────────────────────────────────────────────────────────
# 4) LOSO loop: train multinomial logistic & evaluate
#─────────────────────────────────────────────────────────────────────────────

#' Run LOSO multinomial‐logistic classification in parallel
#'
#' @param X_all   numeric matrix, rows=cells, cols=genes (z-scored)
#' @param sc_meta data.frame or data.table with columns 'region' and 'cluster_id'
#' @param alpha   elastic-net mixing parameter (0 = ridge, 1 = lasso)
#' @param seed_base integer, base to add study index for reproducibility
#' @param workers integer, number of parallel sessions
#' @return data.table with columns study, accuracy, macro_F1

run_loso_cv_future <- function(X_all, sc_meta, alpha = 0, seed_base = 42, workers = 36) {

  plan(multisession, workers = workers)
  
  studies <- unique(sc_meta$region)
  
  # progressor for one tick per study
  with_progress({
    p <- progressor(steps = length(studies))
    
    results_list <- future_lapply(
      seq_along(studies),
      function(i) {
        s <- studies[i]
        # tick the bar
        p(sprintf("Fold %d/%d (%s)", i, length(studies), s))
        
        # … your one_fold body inlined here …
        train_ix <- sc_meta$region != s
        test_ix  <- !train_ix
        X_train  <- X_all[train_ix, , drop=FALSE]
        X_test   <- X_all[test_ix,  , drop=FALSE]
        y_train  <- sc_meta$cluster_id[train_ix]
        y_test   <- sc_meta$cluster_id[test_ix]
        
        set.seed(seed_base + i)
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
        result <- data.table(study = s, accuracy = acc, macro_F1 = mean(f1_vec, na.rm=TRUE))
        rm(X_train, X_test, y_train, y_test, cvfit, preds, f1_vec, classes)
        gc()
        return(result)
      }
    )
  })
  
  rbindlist(results_list)
}

#─────────────────────────────────────────────────────────────────────────────
# 5) Collate & visualize
#─────────────────────────────────────────────────────────────────────────────

library(progressr)
handlers("txtprogressbar")
df_results <- run_loso_cv_future(X_all, sc_meta)

saveRDS(df_results, "results/RN012/classification_putative_frequent_RN012.rds")

print( df_results )

# Simple boxplots across the 19 folds:
par(mfrow = c(1,2))
boxplot(df_results$accuracy, main = "LOSO Accuracy", ylab = "Accuracy")
boxplot(df_results$macro_F1, main = "LOSO Macro-F1", ylab = "Macro-F1")

# You can also run a paired test against another marker set, e.g.:
# wilcox.test(macro_F1 ~ marker_set, data = df_results, paired = TRUE)
