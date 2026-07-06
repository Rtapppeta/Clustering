set.seed(123)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

df <- as.data.frame(counts, stringsAsFactors = FALSE)

first_col_is_gene <- is.null(rownames(df)) || all(rownames(df) == as.character(seq_len(nrow(df))))


if (first_col_is_gene) {
  x1 <- df[[1]]
  suppressWarnings(x1_num <- as.numeric(as.character(x1)))
  if (mean(is.na(x1_num)) > 0.5) {
    rownames(df) <- as.character(df[[1]])
    df <- df[, -1, drop = FALSE]
  }
}

df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]
stopifnot(ncol(df) > 1, nrow(df) > 1)

# ---- Numeric matrix (preserve rownames) ----
expr <- as.matrix(data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                             check.names = FALSE))
rownames(expr) <- rownames(df)
expr[is.na(expr)] <- 0

# ---- Timepoints + average A/B replicates ----
cn <- colnames(expr)
time_num <- as.numeric(sub("^([0-9]+).*", "\\1", cn))
times <- sort(unique(time_num))

expr_avg <- sapply(times, function(t) rowMeans(expr[, time_num == t, drop = FALSE]))
expr_avg <- as.matrix(expr_avg)
rownames(expr_avg) <- rownames(expr)
colnames(expr_avg) <- paste0("T", times)

# ---- Correlation transform (shape-based) ----
X <- log1p(expr_avg)
X_corr <- t(scale(t(X), center = TRUE, scale = FALSE))
X_corr[!is.finite(X_corr)] <- 0

# ---- kmeans ----
k <- 3
km_corr <- kmeans(X_corr, centers = k, nstart = 50, iter.max = 1000)

clusters_corr <- data.frame(gene = rownames(X_corr), cluster_corr = km_corr$cluster)
clusters_corr <- clusters_corr[order(clusters_corr$gene), ]
write.csv(clusters_corr, "clusters_correlation_k3.csv", row.names = FALSE)

# ---- cluster mean trajectories ----
means_corr <- t(sapply(1:k, function(cl) colMeans(X_corr[km_corr$cluster == cl, , drop = FALSE])))
rownames(means_corr) <- paste0("Cluster", 1:k)
colnames(means_corr) <- colnames(X_corr)

# ---- PCA of cluster means ----
pca_corr <- prcomp(means_corr, center = TRUE, scale. = FALSE)
pc <- pca_corr$x[, 1:2, drop = FALSE]

png("PCA_cluster_means_correlation.png", width = 900, height = 650)
plot(pc[,1], pc[,2], pch = 16, cex = 2, xlab = "PC1", ylab = "PC2",
     main = "PCA of Cluster Means (Correlation-based k-means)")
text(pc[,1], pc[,2], labels = rownames(means_corr), pos = 3)
dev.off()

# ---- Plot cluster mean trajectories ----
png("Cluster_mean_trajectories_correlation.png", width = 1100, height = 450)
par(mfrow = c(1, k), mar = c(4,4,3,1))
for (i in 1:k) {
  cl <- paste0("Cluster", i)
  plot(times, means_corr[cl, ], type = "l", lwd = 2,
       xlab = "Timepoint", ylab = "Mean (correlation space)", main = cl)
}
dev.off()
par(mfrow = c(1,1))