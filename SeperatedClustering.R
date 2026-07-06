library(plotly)
library(WGCNA)

set.seed(123)
load("data/ImmuneCounts.RData")

# Clean and prep (same as before)
df <- as.data.frame(counts, stringsAsFactors = FALSE)
rownames(df) <- df[[1]]
df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df))]

expr <- as.matrix(sapply(df, as.numeric))
rownames(expr) <- rownames(df)
expr[is.na(expr)] <- 0
expr <- expr[, order(as.numeric(sub("\\D+", "", colnames(expr))), sub("\\d+", "", colnames(expr)))]

X <- t(scale(t(log1p(expr))))
X[!is.finite(X)] <- 0
rownames(X) <- rownames(expr)

# Full clustering on all data (this is the reference)
km_full <- kmeans(X, centers = 3, nstart = 50, iter.max = 1000)

# Randomly split columns (samples) into two halves
cols <- 1:ncol(X)
set1 <- sample(cols, floor(ncol(X) / 2))
set2 <- setdiff(cols, set1)

X1 <- X[, set1]
X2 <- X[, set2]

# Cluster each half separately
km1 <- kmeans(X1, centers = 3, nstart = 50, iter.max = 1000)
km2 <- kmeans(X2, centers = 3, nstart = 50, iter.max = 1000)

# matchLabels makes sure cluster "1" in set2 corresponds to cluster "1" in set1
# so we're comparing apples to apples
matched_labels <- matchLabels(km2$cluster, km1$cluster)

# Fraction of genes that ended up in the same cluster in both halves
same_cluster <- mean(km1$cluster == matched_labels)
cat(sprintf("Fraction of genes in the same cluster across both splits: %.3f (%.1f%%)\n",
            same_cluster, same_cluster * 100))