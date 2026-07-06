############################################################
# HIERARCHICAL CLUSTERING DENDROGRAM (ONLY)
# Output:
#   hclust_dendrogram_k3.png
############################################################

set.seed(123)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ----------------------------------------------------------
# 1) Clean data
# ----------------------------------------------------------
df <- as.data.frame(counts, stringsAsFactors = FALSE)

# Handle gene names in first column if needed
first_col_numeric_test <- suppressWarnings(as.numeric(as.character(df[[1]])))
if (mean(is.na(first_col_numeric_test)) > 0.5) {
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
}

# Keep only timepoint columns (1A, 2B, etc.)
df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df)), drop = FALSE]

# Convert to numeric matrix
expr <- as.matrix(sapply(df, function(x) as.numeric(as.character(x))))
rownames(expr) <- rownames(df)
colnames(expr) <- colnames(df)
expr[is.na(expr)] <- 0

# ----------------------------------------------------------
# 2) Average A/B replicates
# ----------------------------------------------------------
time_num <- as.numeric(sub("^([0-9]+).*", "\\1", colnames(expr)))
times <- sort(unique(time_num))

expr_avg <- sapply(times, function(t) rowMeans(expr[, time_num == t, drop = FALSE]))
expr_avg <- as.matrix(expr_avg)
rownames(expr_avg) <- rownames(expr)

# ----------------------------------------------------------
# 3) Scale data (same as your clustering)
# ----------------------------------------------------------
X <- t(scale(t(log1p(expr_avg))))
X[!is.finite(X)] <- 0

# Remove genes with zero variance (prevents errors)
row_sd <- apply(X, 1, sd)
X <- X[row_sd > 0, ]

# ----------------------------------------------------------
# 4) Hierarchical clustering (correlation + average linkage)
# ----------------------------------------------------------
hc <- hclust(as.dist(1 - cor(t(X))), method = "average")

# ----------------------------------------------------------
# 5) Plot dendrogram
# ----------------------------------------------------------
png("hclust_dendrogram_k3.png", width = 1000, height = 600)

plot(hc,
     labels = FALSE,   # hide gene names (too many)
     main = "Gene Expression Dendrogram",
     xlab = "",
     sub = "")

# Draw cluster boxes (k = 3)
rect.hclust(hc, k = 3, border = c("red", "green", "blue"))

dev.off()

############################################################