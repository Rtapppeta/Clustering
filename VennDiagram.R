############################################################
# VENN DIAGRAMS FOR OVERLAPPING GENES ACROSS SCHEMES
#
# For each polynomial group, this script:
#   1) rebuilds all k=4 clustering/grouping schemes
#   2) rebuilds the combined PCA of cluster means
#   3) finds the nearest cluster from each other scheme
#   4) makes a 5-set Venn diagram of overlapping genes
#
# Output:
#   venn_match_summary.csv
#   venn_POLY_<group>.png   for each polynomial group present
############################################################

set.seed(123)

# install once if needed
# install.packages("ggVennDiagram")
# install.packages("ggplot2")

library(ggplot2)
library(ggVennDiagram)

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

match_means <- function(raw_means, smooth_means) {
  perms <- all_perms(seq_len(nrow(raw_means)))
  scores <- sapply(perms, function(p) {
    sum(sapply(seq_len(nrow(raw_means)), function(i) {
      sqrt(sum((raw_means[i, ] - smooth_means[p[i], ])^2, na.rm = TRUE))
    }))
  })
  smooth_means[perms[[which.min(scores)]], , drop = FALSE]
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
# 5) Build schemes
# ----------------------------------------------------------
# kmeans
km_raw <- kmeans(X_raw, centers = k, nstart = 50, iter.max = 1000)
km_smooth <- kmeans(X_smooth, centers = k, nstart = 50, iter.max = 1000)

means_km_raw <- get_means(X_raw, km_raw$cluster, paste0("KM_RAW_", 1:k))
means_km_smooth <- get_means(X_smooth, km_smooth$cluster, paste0("KM_SMOOTH_", 1:k))
means_km_smooth <- match_means(means_km_raw, means_km_smooth)
rownames(means_km_smooth) <- paste0("KM_SMOOTH_", 1:k)

# hclust
cl_hc_raw <- safe_cor_hclust(X_raw, k = k)
cl_hc_smooth <- safe_cor_hclust(X_smooth, k = k)

means_hc_raw <- get_means(X_raw, cl_hc_raw, paste0("HC_RAW_", 1:k))
means_hc_smooth <- get_means(X_smooth, cl_hc_smooth, paste0("HC_SMOOTH_", 1:k))
means_hc_smooth <- match_means(means_hc_raw, means_hc_smooth)
rownames(means_hc_smooth) <- paste0("HC_SMOOTH_", 1:k)

# polynomial groups
poly_group <- sapply(seq_len(nrow(expr_avg)), function(i) {
  classify_polynomial_gene(expr_avg[i, ], times)
})

poly_assignments <- data.frame(
  gene = rownames(expr_avg),
  poly_group = poly_group,
  stringsAsFactors = FALSE
)

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
# 6) Combined PCA
# ----------------------------------------------------------
all_means <- rbind(
  means_km_raw,
  means_km_smooth,
  means_hc_raw,
  means_hc_smooth,
  means_poly
)

all_means <- all_means[rowSums(is.na(all_means)) < ncol(all_means), , drop = FALSE]

pca_all <- prcomp(all_means, center = TRUE, scale. = FALSE)
pc <- pca_all$x[, 1:2, drop = FALSE]

# ----------------------------------------------------------
# 7) Gene sets for each scheme
# ----------------------------------------------------------
gene_sets <- list()

# kmeans raw
for (i in 1:k) {
  gene_sets[[paste0("KM_RAW_", i)]] <- rownames(X_raw)[km_raw$cluster == i]
}

# kmeans smooth (original labels are fine for gene sets; PCA matching only affected means)
for (i in 1:k) {
  gene_sets[[paste0("KM_SMOOTH_", i)]] <- rownames(X_smooth)[km_smooth$cluster == i]
}

# hclust raw
for (i in 1:k) {
  gene_sets[[paste0("HC_RAW_", i)]] <- rownames(X_raw)[cl_hc_raw == i]
}

# hclust smooth
for (i in 1:k) {
  gene_sets[[paste0("HC_SMOOTH_", i)]] <- rownames(X_smooth)[cl_hc_smooth == i]
}

# polynomial groups
for (g in poly_present) {
  gene_sets[[paste0("POLY_", g)]] <- poly_assignments$gene[poly_assignments$poly_group == g]
}

# ----------------------------------------------------------
# 8) For each polynomial group, find nearest cluster in each other scheme
# ----------------------------------------------------------
poly_rows <- grep("^POLY_", rownames(all_means), value = TRUE)

match_summary <- data.frame()

for (poly_name in poly_rows) {
  
  poly_idx <- which(rownames(all_means) == poly_name)
  poly_pt <- pc[poly_idx, ]
  
  get_nearest <- function(prefix) {
    idx <- grep(paste0("^", prefix, "_"), rownames(all_means))
    d <- apply(pc[idx, , drop = FALSE], 1, function(z) sqrt(sum((z - poly_pt)^2)))
    rownames(all_means)[idx][which.min(d)]
  }
  
  nearest_km_raw <- get_nearest("KM_RAW")
  nearest_km_smooth <- get_nearest("KM_SMOOTH")
  nearest_hc_raw <- get_nearest("HC_RAW")
  nearest_hc_smooth <- get_nearest("HC_SMOOTH")
  
  match_summary <- rbind(
    match_summary,
    data.frame(
      poly_group = poly_name,
      nearest_km_raw = nearest_km_raw,
      nearest_km_smooth = nearest_km_smooth,
      nearest_hc_raw = nearest_hc_raw,
      nearest_hc_smooth = nearest_hc_smooth,
      stringsAsFactors = FALSE
    )
  )
  
  # 5-set venn
  sets_for_venn <- list(
    POLY = gene_sets[[poly_name]],
    KM_RAW = gene_sets[[nearest_km_raw]],
    KM_SMOOTH = gene_sets[[nearest_km_smooth]],
    HC_RAW = gene_sets[[nearest_hc_raw]],
    HC_SMOOTH = gene_sets[[nearest_hc_smooth]]
  )
  
  p <- ggVennDiagram(sets_for_venn, label_alpha = 0) +
    ggtitle(paste("Gene overlap for", poly_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    filename = paste0("venn_", poly_name, ".png"),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
}

write.csv(match_summary, "venn_match_summary.csv", row.names = FALSE)