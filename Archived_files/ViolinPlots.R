# ==========================================================
# Violin plots of MASE by cluster within each scheme
# ==========================================================

library(readr)
library(dplyr)
library(ggplot2)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/ImmuneCounts.RData")
stopifnot(exists("counts"))

# ----------------------------------------------------------
# 1. Read the spreadsheet
# ----------------------------------------------------------
# If your file is a CSV:
mase_df <- read_csv("gene_scheme_mase.csv")

# If your file is an Excel file instead, use this instead:
# library(readxl)
# mase_df <- read_xlsx("gene_scheme_mase.xlsx")

# ----------------------------------------------------------
# 2. Clean up columns
# ----------------------------------------------------------
plot_df <- mase_df %>%
  filter(!is.na(MASE), !is.na(Cluster), !is.na(Scheme)) %>%
  mutate(
    Scheme = factor(Scheme, levels = c("KM_RAW", "KM_SMOOTH", "HC_RAW", "HC_SMOOTH", "POLY")),
    Cluster = factor(Cluster)
  ) %>%
  filter(is.finite(MASE))

# ----------------------------------------------------------
# 3. Make violin plot
# ----------------------------------------------------------
p <- ggplot(plot_df, aes(x = Cluster, y = MASE, fill = Cluster)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  facet_wrap(~ Scheme, scales = "free_y") +
  labs(
    title = "MASE by Cluster Across Clustering Schemes",
    x = "Cluster",
    y = "MASE"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print(p)

violin_means_cluster <- aggregate(MASE ~ Scheme + Cluster, data = mase_df, FUN = function(x) mean(x, na.rm = TRUE))

print(violin_means_cluster)

write.csv(violin_means_cluster, "violin_plot_means_by_scheme_cluster.csv", row.names = FALSE)

# ----------------------------------------------------------
# 4. Save plot
# ----------------------------------------------------------
ggsave("MASE_violin_by_scheme.png", p, width = 10, height = 6, dpi = 300)
