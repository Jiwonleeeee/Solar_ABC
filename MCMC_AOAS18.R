rm(list=ls())
load("~/Downloads/solar.rdata")
colnames(train.df)
set.seed(123)  # For reproducibility

library(dplyr)
train.df <- train.df %>% select(DateTime, time, GHI_Meas, NAM_GHI_inline, SREF212_GHI_inline,
                    site, lat, long)

start_date <- as.POSIXct("2014-12-01")
end_date <- as.POSIXct("2014-12-15")

# Filter data within the period
train.df <- train.df %>%
  filter(DateTime >= start_date & DateTime <= end_date)

## NA랑 0 일단 제거해볼까?
colnames(train.df)
train.df <- train.df %>%
  filter(if_all(c(GHI_Meas, NAM_GHI_inline, SREF212_GHI_inline), ~ !is.na(.) & . != 0))

locations <- train.df %>% 
  select(lat, long)

scaled_locations <- scale(locations)

# Perform k-means clustering with 50 groups

k <- 50  # Number of clusters
kmeans_result <- kmeans(scaled_locations, centers = k, nstart = 25)

# Add the cluster assignment to the original data
locations$cluster <- kmeans_result$cluster
table(locations$cluster)

train.df$group <- kmeans_result$cluster


## H = 50 / J = 3






