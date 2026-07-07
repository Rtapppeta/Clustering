############################################################
# SIMPLE CLUSTERING COMPARISON SCRIPT
# Makes only:
#   1) hclust_raw_vs_smooth.png
#   2) PCA_all_clustering_schemes_k3.png
############################################################

set.seed(123)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ----------------------------------------------------------
# 1) Clean data safely
# ----------------------------------------------------------
df <- as.data.frame(counts, stringsAsFactors = FALSE)

# If first column is gene names, use it as rownames
first_col_numeric_test <- suppressWarnings(as.numeric(as.character(df[[1]])))
if (mean(is.na(first_col_numeric_test)) > 0.5) {
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
}

# Keep only columns like 1A, 21B
df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]

# Convert safely to numeric matrix
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

k <- 3

# ----------------------------------------------------------
# 3) Build raw matrix
#    log first, then scale per gene
# ----------------------------------------------------------
expr_log <- log1p(expr_avg)

X_raw <- t(scale(t(expr_log)))
X_raw[!is.finite(X_raw)] <- 0
rownames(X_raw) <- rownames(expr_avg)
colnames(X_raw) <- colnames(expr_avg)

# ----------------------------------------------------------
# 4) Build smoothed matrix
#    IMPORTANT: smooth the log-transformed data, not raw counts
# ----------------------------------------------------------
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
# 5) kmeans
# ----------------------------------------------------------
km_raw <- kmeans(X_raw, centers = k, nstart = 50, iter.max = 1000)
km_smooth <- kmeans(X_smooth, centers = k, nstart = 50, iter.max = 1000)

means_km_raw <- t(sapply(1:k, function(cl) {
  colMeans(X_raw[km_raw$cluster == cl, , drop = FALSE])  #instead of looping through k as a number, loop through a vector of the names (peak,valley etc.)
}))
means_km_smooth <- t(sapply(1:k, function(cl) {
  colMeans(X_smooth[km_smooth$cluster == cl, , drop = FALSE])
}))

means_km_raw <- as.matrix(means_km_raw)
means_km_smooth <- as.matrix(means_km_smooth)

rownames(means_km_raw) <- paste0("KM_RAW_", 1:k)
rownames(means_km_smooth) <- paste0("KM_SMOOTH_", 1:k)
colnames(means_km_raw) <- colnames(X_raw)
colnames(means_km_smooth) <- colnames(X_smooth)

# ----------------------------------------------------------
# 6) hclust using correlation distance + average linkage
#    Need to handle zero-variance genes
# ----------------------------------------------------------
safe_cor_dist <- function(X) {
  row_sd <- apply(X, 1, sd)
  good <- row_sd > 0
  
  X_good <- X[good, , drop = FALSE]
  cor_mat <- cor(t(X_good))
  dist_mat <- as.dist(1 - cor_mat)
  
  list(dist = dist_mat, good = good)
}

raw_dist_info <- safe_cor_dist(X_raw)
smooth_dist_info <- safe_cor_dist(X_smooth)

hc_raw <- hclust(raw_dist_info$dist, method = "average")
hc_smooth <- hclust(smooth_dist_info$dist, method = "average")

cl_hc_raw_good <- cutree(hc_raw, k = k)
cl_hc_smooth_good <- cutree(hc_smooth, k = k)

# Assign zero-variance rows to NA first
cl_hc_raw <- rep(NA, nrow(X_raw))
cl_hc_smooth <- rep(NA, nrow(X_smooth))

cl_hc_raw[raw_dist_info$good] <- cl_hc_raw_good
cl_hc_smooth[smooth_dist_info$good] <- cl_hc_smooth_good

# For zero-variance rows, assign them to cluster 1 so code still runs
cl_hc_raw[is.na(cl_hc_raw)] <- 1
cl_hc_smooth[is.na(cl_hc_smooth)] <- 1

means_hc_raw <- t(sapply(1:k, function(cl) {
  colMeans(X_raw[cl_hc_raw == cl, , drop = FALSE])
}))
means_hc_smooth <- t(sapply(1:k, function(cl) {
  colMeans(X_smooth[cl_hc_smooth == cl, , drop = FALSE])
}))

means_hc_raw <- as.matrix(means_hc_raw)
means_hc_smooth <- as.matrix(means_hc_smooth)

rownames(means_hc_raw) <- paste0("HC_RAW_", 1:k)
rownames(means_hc_smooth) <- paste0("HC_SMOOTH_", 1:k)
colnames(means_hc_raw) <- colnames(X_raw)
colnames(means_hc_smooth) <- colnames(X_smooth)

# ----------------------------------------------------------
# 7) Match smooth clusters to raw clusters
# ----------------------------------------------------------
match_means <- function(raw_means, smooth_means) {
  perms <- list(
    c(1,2,3), c(1,3,2),
    c(2,1,3), c(2,3,1),
    c(3,1,2), c(3,2,1)
  )
  
  scores <- sapply(perms, function(p) {
    sum(sapply(1:k, function(i) {
      sqrt(sum((raw_means[i, ] - smooth_means[p[i], ])^2))
    }))
  })
  
  smooth_means[perms[[which.min(scores)]], , drop = FALSE]
}

means_km_smooth <- match_means(means_km_raw, means_km_smooth)
means_hc_smooth <- match_means(means_hc_raw, means_hc_smooth)

rownames(means_km_smooth) <- paste0("KM_SMOOTH_", 1:k)
rownames(means_hc_smooth) <- paste0("HC_SMOOTH_", 1:k)

# ----------------------------------------------------------
# 8) PNG: hclust raw vs smooth trajectories
# ----------------------------------------------------------
png("hclust_raw_vs_smooth.png", width = 1100, height = 450)
par(mfrow = c(1, k), mar = c(4, 4, 3, 1))

for (i in 1:k) {
  y_raw <- means_hc_raw[i, ]
  y_sm <- means_hc_smooth[i, ]
  yr <- range(c(y_raw, y_sm))
  
  plot(times, y_raw, type = "l", lwd = 2, lty = 1,
       xlab = "Timepoint", ylab = "Mean scaled expression",
       main = paste0("hclust Cluster ", i), ylim = yr)
  lines(times, y_sm, lwd = 2, lty = 2)
  legend("topleft", legend = c("RAW", "SMOOTH df=5"),
         lty = c(1, 2), lwd = 2, bty = "n", cex = 0.9)
}

dev.off()
par(mfrow = c(1, 1))

# ----------------------------------------------------------
# 9) PNG: one PCA with all clustering schemes
# ----------------------------------------------------------
all_means <- rbind(
  means_km_raw,
  means_km_smooth,
  means_hc_raw,
  means_hc_smooth  # add another row of means_poly (recluster first)
)

pca_all <- prcomp(all_means, center = TRUE, scale. = FALSE)
pc <- pca_all$x[, 1:2]

method <- c(rep("KM_RAW", 3), rep("KM_SMOOTH", 3), rep("HC_RAW", 3), rep("HC_SMOOTH", 3))
cluster_num <- rep(1:3, 4)

col_map <- c("1" = 1, "2" = 2, "3" = 3)
pch_map <- c("KM_RAW" = 1, "KM_SMOOTH" = 16, "HC_RAW" = 2, "HC_SMOOTH" = 17)

png("PCA_all_clustering_schemes_k3.png", width = 1000, height = 700)
plot(pc[,1], pc[,2],
     col = col_map[as.character(cluster_num)],
     pch = pch_map[method],
     cex = 2,
     xlab = "PC1", ylab = "PC2",
     main = "PCA of cluster means: kmeans vs hclust, raw vs smooth")

text(pc[,1], pc[,2], labels = rownames(all_means), pos = 3, cex = 0.75)

legend("topleft",
       legend = c("Cluster 1", "Cluster 2", "Cluster 3"),
       col = c(1, 2, 3), pch = 16, bty = "n", title = "Cluster")

legend("topright",
       legend = names(pch_map),
       pch = unname(pch_map), bty = "n", title = "Method")

dev.off()



