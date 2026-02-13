#Oog Counts Kmeans Clustering 
library(tidymodels)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")

# Load real data
load("data/OogCounts.RData")

# Quick check - run these in console to see what you're working with
# glimpse(OogCounts)
# str(OogCounts)

set.seed(27)

#Initial known cluster reference: 
centers <- tibble(
  cluster = factor(1:3), 
  num_points = c(100, 150, 50),  # number points in each cluster
  x1 = c(5, 0, -3),              # x1 coordinate of cluster center
  x2 = c(-1, 1, -2)              # x2 coordinate of cluster center
)

counts <- 
  centers %>%
  mutate(
    x1 = map2(num_points, x1, rnorm),
    x2 = map2(num_points, x2, rnorm)
  ) %>% 
  select(-num_points) %>% 
  unnest(cols = c(x1, x2))

ggplot(counts, aes(x1, x2, color = cluster)) +
  geom_point(alpha = 0.3)

# Use OogCounts for real analysis
# Select only numeric columns, drop gene_name
points <- 
  OogCounts %>% 
  select(where(is.numeric))

# Test with k = 3
kclust <- kmeans(points, centers = 3)
kclust

# The effect of different choices of k, from 1 to 9, on this clustering
kclusts <- 
  tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(points, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, points)
  )

# Separate 3 different sets of data using tidy, augment, glance
clusters <- 
  kclusts %>%
  unnest(cols = c(tidied))

assignments <- 
  kclusts %>% 
  unnest(cols = c(augmented))

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

# Create elbow plot -> tot.withinss shows ideal number of k's
ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point()+
  labs(title = "Oogenesis Optimal K")


