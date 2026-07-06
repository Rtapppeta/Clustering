############################################################
# PCA comparison of all clustering/grouping schemes (k = 4)
#
# Schemes included:
#   1) kmeans raw
#   2) kmeans smooth (df = 5)
#   3) hclust raw (correlation distance + average linkage)
#   4) hclust smooth (correlation distance + average linkage)
#   5) polynomial-regression groups on raw data:
#        - peak
#        - valley
#        - increasing_or_late_high
#        - decreasing_or_early_high
#
# Main output:
#   PCA_all_schemes_k4_with_polynomial.png
#
# Optional output:
#   polynomial_group_assignments.csv
############################################################

set.seed(123)

# ----------------------------------------------------------
# 0) Load data
# ----------------------------------------------------------
setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ----------------------------------------------------------
# 1) Clean data safely
# ----------------------------------------------------------
df <- as.data.frame(counts, stringsAsFactors = FALSE)

# If first column looks like gene names, move it to rownames
first_col_numeric_test <- suppressWarnings(as.numeric(as.character(df[[1]])))
if (mean(is.na(first_col_numeric_test)) > 0.5) {
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
}

# Keep only timepoint columns like 1A, 21B
df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]

# Convert to numeric matrix
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
# 4) Helper functions
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

# all permutations for matching smooth to raw
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

match_means <- function(raw_means, smooth_means) {
  perms <- all_perms(seq_len(nrow(raw_means)))
  scores <- sapply(perms, function(p) {
    sum(sapply(seq_len(nrow(raw_means)), function(i) {
      sqrt(sum((raw_means[i, ] - smooth_means[p[i], ])^2, na.rm = TRUE))
    }))
  })
  smooth_means[perms[[which.min(scores)]], , drop = FALSE]
}

# ----------------------------------------------------------
# 5) kmeans raw / smooth (k = 4)
# ----------------------------------------------------------
km_raw <- kmeans(X_raw, centers = k, nstart = 50, iter.max = 1000)
km_smooth <- kmeans(X_smooth, centers = k, nstart = 50, iter.max = 1000)

means_km_raw <- get_means(X_raw, km_raw$cluster, paste0("KM_RAW_", 1:k))
means_km_smooth <- get_means(X_smooth, km_smooth$cluster, paste0("KM_SMOOTH_", 1:k))

# match smooth to raw for readability
means_km_smooth <- match_means(means_km_raw, means_km_smooth)
rownames(means_km_smooth) <- paste0("KM_SMOOTH_", 1:k)

# ----------------------------------------------------------
# 6) hclust raw / smooth (k = 4)
# ----------------------------------------------------------
cl_hc_raw <- safe_cor_hclust(X_raw, k = k)
cl_hc_smooth <- safe_cor_hclust(X_smooth, k = k)

means_hc_raw <- get_means(X_raw, cl_hc_raw, paste0("HC_RAW_", 1:k))
means_hc_smooth <- get_means(X_smooth, cl_hc_smooth, paste0("HC_SMOOTH_", 1:k))

# match smooth to raw for readability
means_hc_smooth <- match_means(means_hc_raw, means_hc_smooth)
rownames(means_hc_smooth) <- paste0("HC_SMOOTH_", 1:k)

# ----------------------------------------------------------
# 7) Polynomial-regression grouping on raw data
#    Uses raw polynomial coefficients, not poly()
# ----------------------------------------------------------
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

poly_group <- sapply(seq_len(nrow(expr_avg)), function(i) {
  classify_polynomial_gene(expr_avg[i, ], times)
})

poly_assignments <- data.frame(
  gene = rownames(expr_avg),
  poly_group = poly_group,
  stringsAsFactors = FALSE
)

write.csv(poly_assignments, "polynomial_group_assignments.csv", row.names = FALSE)

# Keep only the 4 core groups
poly_levels <- c(
  "peak",
  "valley",
  "increasing_or_late_high",
  "decreasing_or_early_high"
)

poly_present <- poly_levels[poly_levels %in% poly_group]

means_poly <- t(sapply(poly_present, function(g) {
  idx <- which(poly_group == g)
  colMeans(X_raw[idx, , drop = FALSE])
}))

means_poly <- as.matrix(means_poly)
rownames(means_poly) <- paste0("POLY_", poly_present)
colnames(means_poly) <- colnames(X_raw)

# ----------------------------------------------------------
# 8) Combine all schemes into one PCA
# ----------------------------------------------------------
png("PCA_all_schemes_k4_with_polynomial.png", width = 1100, height = 750)

plot(
  pc[,1], pc[,2],
  col = col_map[method],
  pch = 19,
  cex = 1.8,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA of cluster/group means: all schemes (k = 4 + polynomial groups)"
)

text(pc[,1], pc[,2], labels = rownames(all_means), pos = 3, cex = 0.75)

# ---- ADD LINES CONNECTING POINTS WITHIN EACH SCHEME ----
# ---- DRAW CLOSED "SQUARE-LIKE" SHAPES ----
unique_methods <- unique(method)

for (m in unique_methods) {
  idx <- which(method == m)
  
  pts <- pc[idx, , drop = FALSE]
  
  # Compute center
  center <- colMeans(pts)
  
  # Compute angles relative to center
  angles <- atan2(pts[,2] - center[2], pts[,1] - center[1])
  
  # Order points around center (circular order)
  ord <- order(angles)
  pts_ord <- pts[ord, , drop = FALSE]
  
  # Close the shape by repeating first point
  pts_closed <- rbind(pts_ord, pts_ord[1, ])
  
  lines(
    pts_closed[,1],
    pts_closed[,2],
    col = col_map[m],
    lwd = 2
  )
}

# ---- MOVE LEGEND TO LEFT ----
legend(
  "topleft",
  legend = names(col_map),
  col = unname(col_map),
  pch = 19,
  lwd = 2,
  bty = "n",
  title = "Scheme"
)

cat("=== KMeans Raw ===\n"); print(table(km_raw$cluster))
cat("=== KMeans Smooth ===\n"); print(table(km_smooth$cluster))
cat("=== HClust Raw ===\n"); print(table(cl_hc_raw))
cat("=== HClust Smooth ===\n"); print(table(cl_hc_smooth))
cat("=== Polynomial ===\n"); print(table(poly_group))

dev.off()