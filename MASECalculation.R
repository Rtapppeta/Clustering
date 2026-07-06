# ==========================================================
# FULL SCRIPT: Calculate MASE for each gene in each scheme
# ==========================================================

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# Check what object(s) were loaded
print(ls())

# IMPORTANT:
# Replace "ImmuneCounts" below with the actual object name
# from ls() if needed
dat <- counts

# -------------------------------
# 3. Build expression-only matrix
# -------------------------------
# Assumes first column contains gene names
rownames(dat) <- dat[, 1]
dat_expr_only <- dat[, -1, drop = FALSE]

# Convert all expression columns to numeric
dat_expr_only <- as.data.frame(lapply(dat_expr_only, as.numeric))
rownames(dat_expr_only) <- rownames(dat)

# Convert to matrix
expr_mat <- as.matrix(dat_expr_only)

# Check result
cat("Expression matrix dimensions:\n")
print(dim(expr_mat))
print(head(expr_mat[, 1:min(5, ncol(expr_mat)), drop = FALSE]))

# -------------------------------
# 4. Check membership_df exists
# -------------------------------
if (!exists("membership_df")) {
  stop("membership_df was not found. Run your clustering script first so membership_df exists.")
}

# clustering scheme columns
scheme_cols <- c("KM_RAW", "KM_SMOOTH", "HC_RAW", "HC_SMOOTH", "POLY")

missing_cols <- scheme_cols[!scheme_cols %in% colnames(membership_df)]
if (length(missing_cols) > 0) {
  stop(paste("These scheme columns are missing from membership_df:",
             paste(missing_cols, collapse = ", ")))
}

# -------------------------------
# 5. Align genes between datasets
# -------------------------------
# If membership_df has rownames, align by gene name
if (!is.null(rownames(membership_df)) && !is.null(rownames(expr_mat))) {
  common_genes <- intersect(rownames(expr_mat), rownames(membership_df))
  
  if (length(common_genes) == 0) {
    stop("No matching gene names between expr_mat and membership_df rownames.")
  }
  
  expr_mat <- expr_mat[common_genes, , drop = FALSE]
  membership_df <- membership_df[common_genes, , drop = FALSE]
} else {
  # otherwise assume same row order
  if (nrow(expr_mat) != nrow(membership_df)) {
    stop("expr_mat and membership_df do not have the same number of rows, and rownames are missing.")
  }
}

cat("Aligned data dimensions:\n")
print(dim(expr_mat))
print(dim(membership_df))

# -------------------------------
# 6. Define MASE function
# -------------------------------
# Numerator:
#   mean absolute difference between gene profile and its cluster mean profile
# Denominator:
#   mean absolute difference between consecutive time points of that gene
compute_gene_mase <- function(gene_profile, cluster_profile) {
  numerator <- mean(abs(gene_profile - cluster_profile), na.rm = TRUE)
  denominator <- mean(abs(diff(gene_profile)), na.rm = TRUE)
  
  if (is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }
  
  numerator / denominator
}

# -------------------------------
# 7. Calculate MASE for each scheme
# -------------------------------
mase_list <- list()

for (scheme in scheme_cols) {
  cat("\nProcessing scheme:", scheme, "\n")
  
  clusters <- membership_df[[scheme]]
  clusters <- as.numeric(clusters)
  
  if (all(is.na(clusters))) {
    warning(paste("All cluster assignments are NA for scheme", scheme))
    next
  }
  
  clust_ids <- sort(unique(clusters[!is.na(clusters)]))
  
  # Build cluster mean profiles
  cluster_means_list <- list()
  
  for (cl in clust_ids) {
    idx <- which(clusters == cl)
    
    if (length(idx) == 1) {
      cluster_means_list[[as.character(cl)]] <- expr_mat[idx, ]
    } else {
      cluster_means_list[[as.character(cl)]] <- colMeans(expr_mat[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  
  # Compute gene-level MASE values
  mase_vals <- rep(NA_real_, nrow(expr_mat))
  
  for (i in seq_len(nrow(expr_mat))) {
    if (is.na(clusters[i])) next
    
    gene_profile <- expr_mat[i, ]
    cl <- as.character(clusters[i])
    cluster_profile <- cluster_means_list[[cl]]
    
    mase_vals[i] <- compute_gene_mase(gene_profile, cluster_profile)
  }
  
  mase_list[[scheme]] <- data.frame(
    Gene = rownames(expr_mat),
    Scheme = scheme,
    Cluster = clusters,
    MASE = mase_vals,
    stringsAsFactors = FALSE
  )
}

# Combine all schemes into one dataframe
mase_df <- do.call(rbind, mase_list)
rownames(mase_df) <- NULL

# -------------------------------
# 8. Inspect output
# -------------------------------
cat("\nFirst few rows of mase_df:\n")
print(head(mase_df))

cat("\nSummary of MASE values:\n")
print(summary(mase_df$MASE))

# -------------------------------
# 9. Save output
# -------------------------------
write.csv(mase_df, "gene_scheme_mase.csv", row.names = FALSE)

cat("\nSaved output to gene_scheme_mase.csv\n")