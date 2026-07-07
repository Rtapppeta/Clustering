############################################################
# FOLLOW-UP VISUALIZATIONS FOR ALL SCHEMES (k = 4)
#
# Makes:
#   1) gene_membership_heatmap.png
#   2) cluster_size_barplot.png
#   3) cluster_mean_trajectories_all_schemes.png
#
# Schemes:
#   - KM_RAW
#   - KM_SMOOTH
#   - HC_RAW
#   - HC_SMOOTH
#   - POLY
############################################################

set.seed(123)

# install once if needed:
# install.packages(c("ggplot2", "reshape2", "RColorBrewer"))

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# ----------------------------------------------------------
# 0) Load data
# ----------------------------------------------------------
setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ----------------------------------------------------------
# 1) Clean data
# ----------------------------------------------------------
df <- as.data.frame(counts, stringsAsFactors = FALSE)

first_col_numeric_test <- suppressWarnings(as.numeric(as.character(df[[1]])))
if (mean(is.na(first_col_numeric_test)) > 0.5) {
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
}

df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]

expr <- as.matrix(sapply(df, function(x) as.numeric(as.character(x))))
if (is.null(dim(expr))) expr <- matrix(expr, ncol = 1)

rownames(expr) <- rownames(df)
colnames(expr) <- colnames(df)
expr[is.na(expr)] <- 0

# ----------------------------------------------------------
# 2) Average A/B replicates by timepoint
# ----------------------------------------------------------
time_num <- as.numeric(sub("^([0-9]+).*", "\\1", colnames(expr)))
times <- sort(unique(time_num))

expr_avg <- sapply(times, function(t) rowMeans(expr[, time_num == t, drop = FALSE]))
expr_avg <- as.matrix(expr_avg)
rownames(expr_avg) <- rownames(expr)
colnames(expr_avg) <- paste0("T", times)

k <- 4

# ----------------------------------------------------------
# 3) Raw and smoothed matrices
# ----------------------------------------------------------
expr_log <- log1p(expr_avg)

X_raw <- t(scale(t(expr_log)))
X_raw[!is.finite(X_raw)] <- 0
rownames(X_raw) <- rownames(expr_avg)
colnames(X_raw) <- colnames(expr_avg)

smooth_gene <- function(y, x) {
  predict(smooth.spline(x = x, y = y, df = 5), x = x)$y
}

expr_smooth_log <- t(apply(expr_log, 1, smooth_gene, x = times))
expr_smooth_log <- as.matrix(expr_smooth_log)
rownames(expr_smooth_log) <- rownames(expr_avg)
colnames(expr_smooth_log) <- colnames(expr_avg)

X_smooth <- t(scale(t(expr_smooth_log)))
X_smooth[!is.finite(X_smooth)] <- 0
rownames(X_smooth) <- rownames(expr_avg)
colnames(X_smooth) <- colnames(expr_avg)

# ----------------------------------------------------------
# 4) Helpers
# ----------------------------------------------------------
get_means <- function(X, cl, labels) {
  out <- t(sapply(seq_along(labels), function(i) {
    idx <- which(cl == i)
    if (length(idx) == 0) {
      rep(NA_real_, ncol(X))
    } else {
      colMeans(X[idx, , drop = FALSE])
    }
  }))
  out <- as.matrix(out)
  rownames(out) <- labels
  colnames(out) <- colnames(X)
  out
}

safe_cor_hclust <- function(X, k = 4) {
  row_sd <- apply(X, 1, sd)
  good <- row_sd > 0
  
  X_good <- X[good, , drop = FALSE]
  dist_mat <- as.dist(1 - cor(t(X_good)))
  hc <- hclust(dist_mat, method = "average")
  cl_good <- cutree(hc, k = k)
  
  cl <- rep(NA_integer_, nrow(X))
  cl[good] <- cl_good
  cl[is.na(cl)] <- 1
  cl
}

all_perms <- function(v) {
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    rest <- v[-i]
    subp <- all_perms(rest)
    for (j in seq_along(subp)) {
      out[[length(out) + 1]] <- c(v[i], subp[[j]])
    }
  }
  out
}

match_clusters_and_means <- function(raw_means, smooth_means, smooth_clusters) {
  perms <- all_perms(seq_len(nrow(raw_means)))
  scores <- sapply(perms, function(p) {
    sum(sapply(seq_len(nrow(raw_means)), function(i) {
      sqrt(sum((raw_means[i, ] - smooth_means[p[i], ])^2, na.rm = TRUE))
    }))
  })
  best_p <- perms[[which.min(scores)]]
  
  means_matched <- smooth_means[best_p, , drop = FALSE]
  
  old_to_new <- integer(length(best_p))
  for (j in seq_along(best_p)) old_to_new[j] <- which(best_p == j)
  
  clusters_matched <- old_to_new[smooth_clusters]
  
  list(
    means = means_matched,
    clusters = clusters_matched
  )
}

classify_polynomial_gene <- function(y, x) {
  deg <- if (length(unique(x)) <= 3) 2 else 3
  
  if (deg == 2) {
    fit <- lm(y ~ x + I(x^2))
    co <- coef(fit)
    b0 <- unname(co[1])
    c1 <- if ("x" %in% names(co)) unname(co["x"]) else 0
    c2 <- if ("I(x^2)" %in% names(co)) unname(co["I(x^2)"]) else 0
    
    local_max_x <- NA
    local_min_x <- NA
    
    if (abs(c2) > 1e-12) {
      xcrit <- -c1 / (2 * c2)
      if (xcrit >= min(x) && xcrit <= max(x)) {
        if (c2 < 0) local_max_x <- xcrit
        if (c2 > 0) local_min_x <- xcrit
      }
    }
    
    cand_x <- c(min(x), max(x))
    if (!is.na(local_max_x)) cand_x <- c(cand_x, local_max_x)
    if (!is.na(local_min_x)) cand_x <- c(cand_x, local_min_x)
    cand_y <- b0 + c1 * cand_x + c2 * cand_x^2
    
  } else {
    fit <- lm(y ~ x + I(x^2) + I(x^3))
    co <- coef(fit)
    b0 <- unname(co[1])
    c1 <- if ("x" %in% names(co)) unname(co["x"]) else 0
    c2 <- if ("I(x^2)" %in% names(co)) unname(co["I(x^2)"]) else 0
    c3 <- if ("I(x^3)" %in% names(co)) unname(co["I(x^3)"]) else 0
    
    local_max_x <- NA
    local_min_x <- NA
    
    D <- c2^2 - 3 * c1 * c3
    crit_all <- numeric(0)
    
    if (abs(c3) > 1e-12 && D > 0) {
      crit_all <- c(
        (-c2 + sqrt(D)) / (3 * c3),
        (-c2 - sqrt(D)) / (3 * c3)
      )
      crit_all <- crit_all[crit_all >= min(x) & crit_all <= max(x)]
      
      if (length(crit_all) > 0) {
        second_deriv <- 2 * c2 + 6 * c3 * crit_all
        if (any(second_deriv < 0)) local_max_x <- crit_all[which(second_deriv < 0)[1]]
        if (any(second_deriv > 0)) local_min_x <- crit_all[which(second_deriv > 0)[1]]
      }
    }
    
    cand_x <- c(min(x), max(x), crit_all)
    cand_y <- b0 + c1 * cand_x + c2 * cand_x^2 + c3 * cand_x^3
  }
  
  global_max_x <- cand_x[which.max(cand_y)]
  
  group <- "complex_or_other"
  if (!is.na(local_max_x) && is.na(local_min_x)) {
    group <- "peak"
  } else if (is.na(local_max_x) && !is.na(local_min_x)) {
    group <- "valley"
  } else if (is.na(local_max_x) && is.na(local_min_x)) {
    if (global_max_x == max(x)) group <- "increasing_or_late_high"
    if (global_max_x == min(x)) group <- "decreasing_or_early_high"
  }
  
  group
}

# ----------------------------------------------------------
# 5) Build all cluster/group assignments
# ----------------------------------------------------------
# kmeans
km_raw <- kmeans(X_raw, centers = k, nstart = 50, iter.max = 1000)
km_smooth <- kmeans(X_smooth, centers = k, nstart = 50, iter.max = 1000)

means_km_raw <- get_means(X_raw, km_raw$cluster, paste0("KM_RAW_", 1:k))
means_km_smooth <- get_means(X_smooth, km_smooth$cluster, paste0("KM_SMOOTH_", 1:k))

km_match <- match_clusters_and_means(means_km_raw, means_km_smooth, km_smooth$cluster)
means_km_smooth <- km_match$means
rownames(means_km_smooth) <- paste0("KM_SMOOTH_", 1:k)
cl_km_smooth_matched <- km_match$clusters

# hclust
cl_hc_raw <- safe_cor_hclust(X_raw, k = k)
cl_hc_smooth <- safe_cor_hclust(X_smooth, k = k)

means_hc_raw <- get_means(X_raw, cl_hc_raw, paste0("HC_RAW_", 1:k))
means_hc_smooth <- get_means(X_smooth, cl_hc_smooth, paste0("HC_SMOOTH_", 1:k))

hc_match <- match_clusters_and_means(means_hc_raw, means_hc_smooth, cl_hc_smooth)
means_hc_smooth <- hc_match$means
rownames(means_hc_smooth) <- paste0("HC_SMOOTH_", 1:k)
cl_hc_smooth_matched <- hc_match$clusters

# polynomial groups -> numeric 1:4 for heatmap/barplot
poly_group <- sapply(seq_len(nrow(expr_avg)), function(i) {
  classify_polynomial_gene(expr_avg[i, ], times)
})

poly_map <- c(
  "peak" = 1,
  "valley" = 2,
  "increasing_or_late_high" = 3,
  "decreasing_or_early_high" = 4
)

cl_poly <- unname(poly_map[poly_group])

# genes not in the 4 main groups become NA
# remove them so all schemes compare the same genes
membership_df <- data.frame(
  gene = rownames(expr_avg),
  KM_RAW = km_raw$cluster,
  KM_SMOOTH = cl_km_smooth_matched,
  HC_RAW = cl_hc_raw,
  HC_SMOOTH = cl_hc_smooth_matched,
  POLY = cl_poly,
  stringsAsFactors = FALSE
)

membership_df <- membership_df[!is.na(membership_df$POLY), ]

# ----------------------------------------------------------
# 6) Heatmap of gene membership across schemes
# ----------------------------------------------------------
scheme_cols <- c("KM_RAW", "KM_SMOOTH", "HC_RAW", "HC_SMOOTH", "POLY")
membership_mat <- as.matrix(membership_df[, scheme_cols])

# order rows by similarity in membership patterns
row_order <- hclust(dist(membership_mat), method = "average")$order
membership_mat <- membership_mat[row_order, , drop = FALSE]

# optional column order by scheme similarity
col_order <- hclust(dist(t(membership_mat)), method = "average")$order
membership_mat <- membership_mat[, col_order, drop = FALSE]

heat_df <- melt(membership_mat)
colnames(heat_df) <- c("gene_index", "scheme", "cluster")
heat_df$gene_index <- factor(heat_df$gene_index, levels = unique(heat_df$gene_index))
heat_df$scheme <- factor(heat_df$scheme, levels = colnames(membership_mat))
heat_df$cluster <- factor(heat_df$cluster, levels = 1:4)

cluster_colors <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02",
  "3" = "#7570b3",
  "4" = "#e7298a"
)

p_heat <- ggplot(heat_df, aes(x = scheme, y = gene_index, fill = cluster)) +
  geom_tile() +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Gene membership across clustering schemes",
    x = "Scheme",
    y = "Genes (ordered by similarity)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("gene_membership_heatmap.png", p_heat, width = 8, height = 10, dpi = 300)

# ----------------------------------------------------------
# 7) Barplot of number of genes in each cluster in each scheme
# ----------------------------------------------------------
bar_df <- melt(membership_df[, scheme_cols])
colnames(bar_df) <- c("scheme", "cluster")
bar_df$scheme <- factor(bar_df$scheme, levels = scheme_cols)
bar_df$cluster <- factor(bar_df$cluster, levels = 1:4)

p_bar <- ggplot(bar_df, aes(x = scheme, fill = cluster)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Number of genes in each cluster across schemes",
    x = "Scheme",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("cluster_size_barplot.png", p_bar, width = 8, height = 5, dpi = 300)

# ----------------------------------------------------------
# 8) Average expression of each cluster as companion to PCA
# ----------------------------------------------------------
############################################################
# FOLLOW-UP VISUALIZATIONS FOR ALL SCHEMES (k = 4)
#
# Makes:
#   1) gene_membership_heatmap.png
#   2) cluster_size_barplot.png
#   3) cluster_mean_trajectories_all_schemes.png
#
# Schemes:
#   - KM_RAW
#   - KM_SMOOTH
#   - HC_RAW
#   - HC_SMOOTH
#   - POLY
############################################################

set.seed(123)

# install once if needed:
# install.packages(c("ggplot2", "reshape2", "RColorBrewer"))

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# ----------------------------------------------------------
# 0) Load data
# ----------------------------------------------------------
setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ----------------------------------------------------------
# 1) Clean data
# ----------------------------------------------------------
df <- as.data.frame(counts, stringsAsFactors = FALSE)

first_col_numeric_test <- suppressWarnings(as.numeric(as.character(df[[1]])))
if (mean(is.na(first_col_numeric_test)) > 0.5) {
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
}

df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]

expr <- as.matrix(sapply(df, function(x) as.numeric(as.character(x))))
if (is.null(dim(expr))) expr <- matrix(expr, ncol = 1)

rownames(expr) <- rownames(df)
colnames(expr) <- colnames(df)
expr[is.na(expr)] <- 0

# ----------------------------------------------------------
# 2) Average A/B replicates by timepoint
# ----------------------------------------------------------
time_num <- as.numeric(sub("^([0-9]+).*", "\\1", colnames(expr)))
times <- sort(unique(time_num))

expr_avg <- sapply(times, function(t) rowMeans(expr[, time_num == t, drop = FALSE]))
expr_avg <- as.matrix(expr_avg)
rownames(expr_avg) <- rownames(expr)
colnames(expr_avg) <- paste0("T", times)

k <- 4

# ----------------------------------------------------------
# 3) Raw and smoothed matrices
# ----------------------------------------------------------
expr_log <- log1p(expr_avg)

X_raw <- t(scale(t(expr_log)))
X_raw[!is.finite(X_raw)] <- 0
rownames(X_raw) <- rownames(expr_avg)
colnames(X_raw) <- colnames(expr_avg)

smooth_gene <- function(y, x) {
  predict(smooth.spline(x = x, y = y, df = 5), x = x)$y
}

expr_smooth_log <- t(apply(expr_log, 1, smooth_gene, x = times))
expr_smooth_log <- as.matrix(expr_smooth_log)
rownames(expr_smooth_log) <- rownames(expr_avg)
colnames(expr_smooth_log) <- colnames(expr_avg)

X_smooth <- t(scale(t(expr_smooth_log)))
X_smooth[!is.finite(X_smooth)] <- 0
rownames(X_smooth) <- rownames(expr_avg)
colnames(X_smooth) <- colnames(expr_avg)

# ----------------------------------------------------------
# 4) Helpers
# ----------------------------------------------------------
get_means <- function(X, cl, labels) {
  out <- t(sapply(seq_along(labels), function(i) {
    idx <- which(cl == i)
    if (length(idx) == 0) {
      rep(NA_real_, ncol(X))
    } else {
      colMeans(X[idx, , drop = FALSE])
    }
  }))
  out <- as.matrix(out)
  rownames(out) <- labels
  colnames(out) <- colnames(X)
  out
}

safe_cor_hclust <- function(X, k = 4) {
  row_sd <- apply(X, 1, sd)
  good <- row_sd > 0
  
  X_good <- X[good, , drop = FALSE]
  dist_mat <- as.dist(1 - cor(t(X_good)))
  hc <- hclust(dist_mat, method = "average")
  cl_good <- cutree(hc, k = k)
  
  cl <- rep(NA_integer_, nrow(X))
  cl[good] <- cl_good
  cl[is.na(cl)] <- 1
  cl
}

all_perms <- function(v) {
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    rest <- v[-i]
    subp <- all_perms(rest)
    for (j in seq_along(subp)) {
      out[[length(out) + 1]] <- c(v[i], subp[[j]])
    }
  }
  out
}

match_clusters_and_means <- function(raw_means, smooth_means, smooth_clusters) {
  perms <- all_perms(seq_len(nrow(raw_means)))
  scores <- sapply(perms, function(p) {
    sum(sapply(seq_len(nrow(raw_means)), function(i) {
      sqrt(sum((raw_means[i, ] - smooth_means[p[i], ])^2, na.rm = TRUE))
    }))
  })
  best_p <- perms[[which.min(scores)]]
  
  means_matched <- smooth_means[best_p, , drop = FALSE]
  
  old_to_new <- integer(length(best_p))
  for (j in seq_along(best_p)) old_to_new[j] <- which(best_p == j)
  
  clusters_matched <- old_to_new[smooth_clusters]
  
  list(
    means = means_matched,
    clusters = clusters_matched
  )
}

classify_polynomial_gene <- function(y, x) {
  deg <- if (length(unique(x)) <= 3) 2 else 3
  
  if (deg == 2) {
    fit <- lm(y ~ x + I(x^2))
    co <- coef(fit)
    b0 <- unname(co[1])
    c1 <- if ("x" %in% names(co)) unname(co["x"]) else 0
    c2 <- if ("I(x^2)" %in% names(co)) unname(co["I(x^2)"]) else 0
    
    local_max_x <- NA
    local_min_x <- NA
    
    if (abs(c2) > 1e-12) {
      xcrit <- -c1 / (2 * c2)
      if (xcrit >= min(x) && xcrit <= max(x)) {
        if (c2 < 0) local_max_x <- xcrit
        if (c2 > 0) local_min_x <- xcrit
      }
    }
    
    cand_x <- c(min(x), max(x))
    if (!is.na(local_max_x)) cand_x <- c(cand_x, local_max_x)
    if (!is.na(local_min_x)) cand_x <- c(cand_x, local_min_x)
    cand_y <- b0 + c1 * cand_x + c2 * cand_x^2
    
  } else {
    fit <- lm(y ~ x + I(x^2) + I(x^3))
    co <- coef(fit)
    b0 <- unname(co[1])
    c1 <- if ("x" %in% names(co)) unname(co["x"]) else 0
    c2 <- if ("I(x^2)" %in% names(co)) unname(co["I(x^2)"]) else 0
    c3 <- if ("I(x^3)" %in% names(co)) unname(co["I(x^3)"]) else 0
    
    local_max_x <- NA
    local_min_x <- NA
    
    D <- c2^2 - 3 * c1 * c3
    crit_all <- numeric(0)
    
    if (abs(c3) > 1e-12 && D > 0) {
      crit_all <- c(
        (-c2 + sqrt(D)) / (3 * c3),
        (-c2 - sqrt(D)) / (3 * c3)
      )
      crit_all <- crit_all[crit_all >= min(x) & crit_all <= max(x)]
      
      if (length(crit_all) > 0) {
        second_deriv <- 2 * c2 + 6 * c3 * crit_all
        if (any(second_deriv < 0)) local_max_x <- crit_all[which(second_deriv < 0)[1]]
        if (any(second_deriv > 0)) local_min_x <- crit_all[which(second_deriv > 0)[1]]
      }
    }
    
    cand_x <- c(min(x), max(x), crit_all)
    cand_y <- b0 + c1 * cand_x + c2 * cand_x^2 + c3 * cand_x^3
  }
  
  global_max_x <- cand_x[which.max(cand_y)]
  
  group <- "complex_or_other"
  if (!is.na(local_max_x) && is.na(local_min_x)) {
    group <- "peak"
  } else if (is.na(local_max_x) && !is.na(local_min_x)) {
    group <- "valley"
  } else if (is.na(local_max_x) && is.na(local_min_x)) {
    if (global_max_x == max(x)) group <- "increasing_or_late_high"
    if (global_max_x == min(x)) group <- "decreasing_or_early_high"
  }
  
  group
}

# ----------------------------------------------------------
# 5) Build all cluster/group assignments
# ----------------------------------------------------------
# kmeans
km_raw <- kmeans(X_raw, centers = k, nstart = 50, iter.max = 1000)
km_smooth <- kmeans(X_smooth, centers = k, nstart = 50, iter.max = 1000)

means_km_raw <- get_means(X_raw, km_raw$cluster, paste0("KM_RAW_", 1:k))
means_km_smooth <- get_means(X_smooth, km_smooth$cluster, paste0("KM_SMOOTH_", 1:k))

km_match <- match_clusters_and_means(means_km_raw, means_km_smooth, km_smooth$cluster)
means_km_smooth <- km_match$means
rownames(means_km_smooth) <- paste0("KM_SMOOTH_", 1:k)
cl_km_smooth_matched <- km_match$clusters

# hclust
cl_hc_raw <- safe_cor_hclust(X_raw, k = k)
cl_hc_smooth <- safe_cor_hclust(X_smooth, k = k)

means_hc_raw <- get_means(X_raw, cl_hc_raw, paste0("HC_RAW_", 1:k))
means_hc_smooth <- get_means(X_smooth, cl_hc_smooth, paste0("HC_SMOOTH_", 1:k))

hc_match <- match_clusters_and_means(means_hc_raw, means_hc_smooth, cl_hc_smooth)
means_hc_smooth <- hc_match$means
rownames(means_hc_smooth) <- paste0("HC_SMOOTH_", 1:k)
cl_hc_smooth_matched <- hc_match$clusters

# polynomial groups -> numeric 1:4 for heatmap/barplot
poly_group <- sapply(seq_len(nrow(expr_avg)), function(i) {
  classify_polynomial_gene(expr_avg[i, ], times)
})

poly_map <- c(
  "peak" = 1,
  "valley" = 2,
  "increasing_or_late_high" = 3,
  "decreasing_or_early_high" = 4
)

cl_poly <- unname(poly_map[poly_group])

# genes not in the 4 main groups become NA
# remove them so all schemes compare the same genes
membership_df <- data.frame(
  gene = rownames(expr_avg),
  KM_RAW = km_raw$cluster,
  KM_SMOOTH = cl_km_smooth_matched,
  HC_RAW = cl_hc_raw,
  HC_SMOOTH = cl_hc_smooth_matched,
  POLY = cl_poly,
  stringsAsFactors = FALSE
)

membership_df <- membership_df[!is.na(membership_df$POLY), ]

# ----------------------------------------------------------
## ----------------------------------------------------------
# Heatmap of gene membership across clustering schemes
# base R heatmap()
# ----------------------------------------------------------
scheme_cols <- c("KM_RAW", "KM_SMOOTH", "HC_RAW", "HC_SMOOTH", "POLY")

# build matrix of cluster memberships
membership_mat <- as.matrix(membership_df[, scheme_cols])

# make sure values are numeric
membership_mat <- apply(membership_mat, 2, as.numeric)

# optional: add gene names as row names if available
# rownames(membership_mat) <- membership_df$Gene

# colors for clusters 1, 2, 3, 4
cluster_colors <- c(
  "#1b9e77",  # cluster 1
  "#d95f02",  # cluster 2
  "#7570b3",  # cluster 3
  "#e7298a"   # cluster 4
)

# -----------------------------
# build clustering objects first
# -----------------------------
row_hc <- hclust(dist(membership_mat), method = "complete")
col_hc <- hclust(dist(t(membership_mat)), method = "complete")

# -----------------------------------------
# heatmap with more room for row dendrogram
# -----------------------------------------
png("gene_membership_heatmap_baseR_tall.png", width = 1200, height = 2200, res = 200)

heatmap(
  membership_mat,
  scale = "none",
  col = cluster_colors,
  breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),   # one color per cluster
  Rowv = as.dendrogram(row_hc),          # use explicit row tree
  Colv = as.dendrogram(col_hc),          # use explicit column tree
  labRow = NA,                           # hide gene labels
  labCol = colnames(membership_mat),
  margins = c(10, 8),
  lwid = c(3, 8),                        # more width for row dendrogram
  lhei = c(3, 8),
  main = "Gene membership across clustering schemes"
)

legend(
  "topright",
  legend = paste("Cluster", 1:4),
  fill = cluster_colors,
  border = NA,
  bty = "n",
  cex = 1
)

dev.off()

# -----------------------------------------
# row dendrogram alone: easier to inspect
# -----------------------------------------
pdf("gene_membership_row_dendrogram.pdf", width = 10, height = 18)
plot(
  row_hc,
  hang = -1,
  labels = FALSE,
  main = "Row dendrogram: gene membership patterns",
  xlab = "",
  sub = "",
  cex = 0.3
)
dev.off()

png("gene_membership_row_dendrogram.png", width = 1200, height = 2400, res = 200)
plot(
  row_hc,
  hang = -1,
  labels = FALSE,
  main = "Row dendrogram: gene membership patterns",
  xlab = "",
  sub = "",
  cex = 0.3
)
dev.off()

# ----------------------------------------------------------
# 7) Barplot of number of genes in each cluster in each scheme
# ----------------------------------------------------------
bar_df <- melt(membership_df[, scheme_cols])
colnames(bar_df) <- c("scheme", "cluster")
bar_df$scheme <- factor(bar_df$scheme, levels = scheme_cols)
bar_df$cluster <- factor(bar_df$cluster, levels = 1:4)

cluster_count_labels <- as.data.frame(table(bar_df$scheme, bar_df$cluster))
colnames(cluster_count_labels) <- c("scheme", "cluster", "gene_count")

p_bar <- ggplot(bar_df, aes(x = scheme, fill = cluster)) +
  geom_bar(position = "stack") +
  geom_text(
    data = cluster_count_labels,
    aes(x = scheme, y = gene_count, label = gene_count, group = cluster),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Number of genes in each cluster across schemes",
    x = "Scheme",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("cluster_size_barplot_labeled.png", p_bar, width = 8, height = 5, dpi = 300)

# ----------------------------------------------------------
# 8) Average expression of each cluster as companion to PCA
# ----------------------------------------------------------
keep_genes <- membership_df$gene

get_means_from_membership <- function(X, membership_vec, scheme_name) {
  out <- do.call(rbind, lapply(1:4, function(cl) {
    idx <- which(membership_vec == cl)
    if (length(idx) == 0) return(NULL)
    data.frame(
      scheme = scheme_name,
      cluster = paste0("Cluster ", cl),
      time = times,
      mean_expr = colMeans(X[idx, , drop = FALSE]),
      stringsAsFactors = FALSE
    )
  }))
  out
}

traj_df <- rbind(
  get_means_from_membership(X_raw[keep_genes, , drop = FALSE], membership_df$KM_RAW, "KM_RAW"),
  get_means_from_membership(X_smooth[keep_genes, , drop = FALSE], membership_df$KM_SMOOTH, "KM_SMOOTH"),
  get_means_from_membership(X_raw[keep_genes, , drop = FALSE], membership_df$HC_RAW, "HC_RAW"),
  get_means_from_membership(X_smooth[keep_genes, , drop = FALSE], membership_df$HC_SMOOTH, "HC_SMOOTH"),
  get_means_from_membership(X_raw[keep_genes, , drop = FALSE], membership_df$POLY, "POLY")
)

traj_df$scheme <- factor(traj_df$scheme, levels = scheme_cols)
traj_df$cluster <- factor(traj_df$cluster, levels = paste0("Cluster ", 1:4))

cluster_colors <- c(
  "Cluster 1" = "#1b9e77",
  "Cluster 2" = "#d95f02",
  "Cluster 3" = "#7570b3",
  "Cluster 4" = "#e7298a"
)

label_df <- traj_df[traj_df$time == max(traj_df$time), ]

p_traj <- ggplot(traj_df, aes(x = time, y = mean_expr, color = cluster, group = cluster)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ scheme, scales = "fixed") +
  scale_color_manual(values = cluster_colors) +
  geom_text(
    data = label_df,
    aes(label = cluster),
    hjust = -0.1,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "Average expression trajectories of each cluster across schemes",
    x = "Timepoint",
    y = "Mean scaled expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(min(times), max(times) + 1)) +
  guides(color = guide_legend(title = "Cluster"))

ggsave("cluster_mean_trajectories_all_schemes_labeled.png", p_traj, width = 12, height = 7, dpi = 300)
