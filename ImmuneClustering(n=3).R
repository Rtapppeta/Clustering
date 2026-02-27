#Immune Counts Clustering n=3

############################################################
# Gene time-course clustering with k-means (k = 3)
# Rows = genes, Columns = time points
# Output:
#   - gene_kmeans_clusters_k3.csv  (gene -> cluster assignment)
#   - cluster_mean_trajectories_k3.png (average trajectory per cluster)
############################################################

set.seed(123)

# ---- 0) Load data ----
# If you already loaded ImmuneCounts.RData earlier, you can comment these two lines out.
setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("ImmuneCounts.RData")

# ---- 1) Get expression matrix (genes x time points) ----
if (!exists("counts")) stop("Object `counts` not found. Did you load ImmuneCounts.RData?")
expr <- as.matrix(counts)

if (is.null(rownames(expr))) stop("`counts` has no rownames. Row names must be gene IDs/names.")
if (nrow(expr) < 2 || ncol(expr) < 2) stop("`counts` must be at least 2 genes x 2 time points.")
if (!is.numeric(expr)) storage.mode(expr) <- "numeric"

message("Counts dim: ", nrow(expr), " genes x ", ncol(expr), " time points")

# ---- 2) Preprocess: log-transform + per-gene scaling across time ----
# log1p stabilizes variance for count-like data
expr_log <- log1p(expr)

# Scale each gene across time points to cluster by SHAPE over time
expr_scaled <- t(scale(t(expr_log)))

# Replace NA/Inf (e.g., genes with zero variance across time) with 0
expr_scaled[!is.finite(expr_scaled)] <- 0

# ---- 3) K-means clustering across genes ----
k <- 3
km <- kmeans(expr_scaled, centers = k, nstart = 50, iter.max = 1000)

# ---- 4) Gene -> cluster assignment table ----
gene_clusters <- data.frame(
  gene = rownames(expr_scaled),
  cluster = as.integer(km$cluster),
  stringsAsFactors = FALSE
)

# Sort for readability
gene_clusters <- gene_clusters[order(gene_clusters$cluster, gene_clusters$gene), ]

# Save assignments
write.csv(gene_clusters, "gene_kmeans_clusters_k3.csv", row.names = FALSE)

# Print summary
cat("\nCluster sizes:\n")
print(table(gene_clusters$cluster))
cat("\nFirst few assignments:\n")
print(head(gene_clusters, 10))

# ---- 5) Plot mean trajectory per cluster (scaled expression) ----
# Compute mean scaled expression at each time point for each cluster
cluster_means <- matrix(NA_real_, nrow = k, ncol = ncol(expr_scaled))
for (cl in 1:k) {
  idx <- which(km$cluster == cl)
  if (length(idx) == 0) next
  cluster_means[cl, ] <- colMeans(expr_scaled[idx, , drop = FALSE])
}

# Use time point labels if they exist, else 1..T
tp <- colnames(expr_scaled)
if (is.null(tp)) tp <- seq_len(ncol(expr_scaled))
x <- seq_len(ncol(expr_scaled))

png("cluster_mean_trajectories_k3.png", width = 900, height = 600)
matplot(
  x, t(cluster_means),
  type = "l", lty = 1, lwd = 2,
  xaxt = "n",
  xlab = "Time point",
  ylab = "Mean scaled expression (per gene z-score)"
)
axis(1, at = x, labels = tp)
legend("topright", legend = paste("Cluster", 1:k), lty = 1, lwd = 2, bty = "n")
title("K-means clustering of gene expression trajectories (k = 3)")
dev.off()

cat("\nSaved:\n")
cat(" - gene_kmeans_clusters_k3.csv\n")
cat(" - cluster_mean_trajectories_k3.png\n")

# gene_clusters is the final in-memory object:
# every gene assigned to cluster 1..3
############################################################
