
set.seed(123)

# ---- 0) Load ----
setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ---- 1) Clean to numeric matrix (genes as rownames) ----
df <- as.data.frame(counts, stringsAsFactors = FALSE)
genes <- df[[1]]
df <- df[, -1, drop = FALSE]
rownames(df) <- genes
df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]

expr <- as.matrix(data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                             check.names = FALSE))
rownames(expr) <- rownames(df)
expr[is.na(expr)] <- 0

# ---- 2) Average A/B technical replicates by timepoint ----
cn <- colnames(expr)
time_num <- as.numeric(sub("^([0-9]+).*", "\\1", cn))
times <- sort(unique(time_num))

expr_avg <- sapply(times, function(t) rowMeans(expr[, time_num == t, drop = FALSE]))
expr_avg <- as.matrix(expr_avg)
rownames(expr_avg) <- rownames(expr)
colnames(expr_avg) <- paste0("T", times)

k <- 3

# (A) RAW clustering=
X_raw <- t(scale(t(log1p(expr_avg))))
X_raw[!is.finite(X_raw)] <- 0
rownames(X_raw) <- rownames(expr_avg)
colnames(X_raw) <- colnames(expr_avg)

km_raw <- kmeans(X_raw, centers = k, nstart = 50, iter.max = 1000)

clusters_raw <- data.frame(gene = rownames(X_raw), cluster_raw = as.integer(km_raw$cluster))
clusters_raw <- clusters_raw[order(clusters_raw$gene), ]
write.csv(clusters_raw, "clusters_raw_k3.csv", row.names = FALSE)

means_raw <- t(sapply(1:k, function(cl) colMeans(X_raw[km_raw$cluster == cl, , drop = FALSE])))
rownames(means_raw) <- paste0("Cluster", 1:k)
colnames(means_raw) <- colnames(X_raw)

pca_raw <- prcomp(means_raw, center = TRUE, scale. = FALSE)
pc_raw <- pca_raw$x[, 1:2, drop = FALSE]
png("PCA_cluster_means_raw.png", width = 900, height = 650)
plot(pc_raw[,1], pc_raw[,2], pch = 16, cex = 2, xlab = "PC1", ylab = "PC2",
     main = "PCA of cluster-average trajectories (RAW)")
text(pc_raw[,1], pc_raw[,2], labels = rownames(means_raw), pos = 3)
dev.off()


# (B) SPLINE(df=5) clustering
smooth_gene_df5 <- function(y, x) predict(smooth.spline(x = x, y = y, df = 5), x = x)$y

expr_smooth <- t(apply(expr_avg, 1, smooth_gene_df5, x = times))
rownames(expr_smooth) <- rownames(expr_avg)
colnames(expr_smooth) <- colnames(expr_avg)

X_smooth <- t(scale(t(log1p(expr_smooth))))
X_smooth[!is.finite(X_smooth)] <- 0
rownames(X_smooth) <- rownames(expr_smooth)
colnames(X_smooth) <- colnames(expr_smooth)

km_smooth <- kmeans(X_smooth, centers = k, nstart = 50, iter.max = 1000)

clusters_smooth <- data.frame(gene = rownames(X_smooth), cluster_smooth = as.integer(km_smooth$cluster))
clusters_smooth <- clusters_smooth[order(clusters_smooth$gene), ]
write.csv(clusters_smooth, "clusters_smooth_df5_k3.csv", row.names = FALSE)

means_smooth <- t(sapply(1:k, function(cl) colMeans(X_smooth[km_smooth$cluster == cl, , drop = FALSE])))
rownames(means_smooth) <- paste0("Cluster", 1:k)
colnames(means_smooth) <- colnames(X_smooth)

pca_smooth <- prcomp(means_smooth, center = TRUE, scale. = FALSE)
pc_smooth <- pca_smooth$x[, 1:2, drop = FALSE]
png("PCA_cluster_means_smooth_df5.png", width = 900, height = 650)
plot(pc_smooth[,1], pc_smooth[,2], pch = 16, cex = 2, xlab = "PC1", ylab = "PC2",
     main = "PCA of cluster-average trajectories (SPLINE df=5)")
text(pc_smooth[,1], pc_smooth[,2], labels = rownames(means_smooth), pos = 3)
dev.off()


# (C) MATCH SMOOTH clusters to RAW clusters (fix label switching)
D <- matrix(0, k, k)
for (i in 1:k) for (j in 1:k) {
  D[i, j] <- sqrt(sum((means_raw[i, ] - means_smooth[j, ])^2))
}

perms <- list(
  c(1,2,3), c(1,3,2),
  c(2,1,3), c(2,3,1),
  c(3,1,2), c(3,2,1)
)
scores <- sapply(perms, function(p) sum(D[cbind(1:k, p)]))
best_p <- perms[[which.min(scores)]]

# Reorder smooth means so row i matches raw Cluster i
means_smooth_m <- means_smooth[best_p, , drop = FALSE]
rownames(means_smooth_m) <- rownames(means_raw)

# Remap smooth gene cluster labels to matched numbering
# old smooth cluster j becomes new label = which(best_p == j)
map_old_to_new <- integer(k)
for (j in 1:k) map_old_to_new[j] <- which(best_p == j)

clusters_smooth$cluster_smooth_matched <- map_old_to_new[clusters_smooth$cluster_smooth]
write.csv(clusters_smooth, "clusters_smooth_df5_k3.csv", row.names = FALSE)

match_map <- data.frame(
  raw_cluster = paste0("Cluster", 1:k),
  smooth_cluster_original = paste0("Cluster", best_p),
  stringsAsFactors = FALSE
)
write.csv(match_map, "cluster_match_map_df5.csv", row.names = FALSE)


# (D) Combined PCA of 6 cluster means (RAW + MATCHED SMOOTH) + lines
means_combined <- rbind(means_raw, means_smooth_m)

method <- c(rep("RAW", k), rep("SMOOTH_df5", k))
cluster_lab <- rep(paste0("Cluster", 1:k), times = 2)
rownames(means_combined) <- paste(method, cluster_lab, sep = "_")

pca_comb <- prcomp(means_combined, center = TRUE, scale. = FALSE)
pc <- pca_comb$x[, 1:2, drop = FALSE]

col_map <- c(Cluster1 = 1, Cluster2 = 2, Cluster3 = 3)
pch_map <- c(RAW = 1, SMOOTH_df5 = 16)

pt_col <- col_map[cluster_lab]
pt_pch <- pch_map[method]

png("PCA_cluster_means_RAW_vs_SMOOTH_df5_matched.png", width = 950, height = 650)
plot(pc[,1], pc[,2], col = pt_col, pch = pt_pch, cex = 2,
     xlab = "PC1", ylab = "PC2",
     main = "Combined PCA of cluster means: RAW vs SPLINE(df=5) (matched)")
text(pc[,1], pc[,2], labels = rownames(means_combined), pos = 3, cex = 0.8)

for (i in 1:k) {
  i_raw <- which(method == "RAW" & cluster_lab == paste0("Cluster", i))
  i_sm  <- which(method == "SMOOTH_df5" & cluster_lab == paste0("Cluster", i))
  segments(pc[i_raw,1], pc[i_raw,2], pc[i_sm,1], pc[i_sm,2], lwd = 2)
}

legend("topleft", legend = names(col_map), col = unname(col_map), pch = 16,
       title = "Cluster", bty = "n")
legend("topright", legend = names(pch_map), pch = unname(pch_map),
       title = "Method", bty = "n")
dev.off()

# (E) Overlay cluster mean trajectories over time (RAW vs MATCHED SMOOTH)
png("Cluster_mean_trajectories_RAW_vs_SMOOTH_df5_matched.png", width = 1100, height = 450)
par(mfrow = c(1, k), mar = c(4, 4, 3, 1))

for (i in 1:k) {
  cl <- paste0("Cluster", i)
  y_raw <- means_raw[cl, ]
  y_sm  <- means_smooth_m[cl, ]
  rng <- range(c(y_raw, y_sm), finite = TRUE)
  
  plot(times, y_raw, type = "l", lwd = 2, lty = 1,
       xlab = "Timepoint", ylab = "Mean scaled expression",
       main = cl, ylim = rng)
  lines(times, y_sm, lwd = 2, lty = 2)
  legend("topleft", legend = c("RAW", "SPLINE df=5"),
         lty = c(1, 2), lwd = 2, bty = "n", cex = 0.9)
}

dev.off()
par(mfrow = c(1, 1))