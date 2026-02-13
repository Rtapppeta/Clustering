#Oogenesis Kmeans Elbow Plot


rm(list = ls())
library(tidymodels)
library(tidyverse)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/OogCounts.RData")
oog_counts <- counts

set.seed(27)
oog_points <- oog_counts %>% select(where(is.numeric))

kclusts <- 
  tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(oog_points, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, oog_points)
  )

clusters    <- kclusts %>% unnest(cols = c(tidied))
assignments <- kclusts %>% unnest(cols = c(augmented))
clusterings <- kclusts %>% unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  labs(title = "Oogenesis Optimal K")



