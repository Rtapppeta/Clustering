#Developmental Kmeans Elbow Plot
#final version 

rm(list = ls())
library(tidymodels)
library(tidyverse)

setwd("C:/Users/ritik/Downloads/Wunderlich Lab/Pleiotropy")
load("data/DevCounts.RData")
dev_counts <- counts

set.seed(27)
dev_points <- dev_counts %>% select(where(is.numeric))

kclusts <- 
  tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(dev_points, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, dev_points)
  )

clusters    <- kclusts %>% unnest(cols = c(tidied))
assignments <- kclusts %>% unnest(cols = c(augmented))
clusterings <- kclusts %>% unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  labs(title = "Developmental Optimal K")

