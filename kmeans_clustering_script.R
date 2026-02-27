library(plotly)

set.seed(123)
setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")

# Gene names as row names
df <- as.data.frame(counts, stringsAsFactors = FALSE)
rownames(df) <- df[[1]]
df <- df[, grepl("^[0-9]+[A-Za-z]+$", colnames(df))]

# Convert to numeric, log scale the expression
expr <- as.matrix(sapply(df, as.numeric))
rownames(expr) <- rownames(df)
expr[is.na(expr)] <- 0
expr <- expr[, order(as.numeric(sub("\\D+", "", colnames(expr))), sub("\\d+", "", colnames(expr)))]

X <- t(scale(t(log1p(expr))))
X[!is.finite(X)] <- 0
rownames(X) <- rownames(expr)

# Cluster
km <- kmeans(X, centers = 3, nstart = 50, iter.max = 1000)

# PCA + plot
pca <- prcomp(X[, apply(X, 2, var) > 0], center = TRUE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

fig <- plot_ly(
  x = pca$x[,1], y = pca$x[,2],
  color = as.factor(km$cluster),
  text = rownames(X),
  hovertemplate = "<b>%{text}</b><br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra>Cluster %{color}</extra>",
  type = "scatter", mode = "markers",
  marker = list(size = 6, opacity = 0.8)
) %>% layout(
  title = "K-means (k=3) PCA",
  xaxis = list(title = paste0("PC1 (", pct[1], "%)")),
  yaxis = list(title = paste0("PC2 (", pct[2], "%)"))
)
annotate_figure(print(fig), top = "test")

