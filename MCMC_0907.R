load("~/Downloads/6month.rdata")

library(dplyr)

start_date <- as.POSIXct("2014-12-01")
end_date <- as.POSIXct("2014-12-15")

# Filter data within the period
filtered_data <- train.df %>%
  filter(DateTime >= start_date & DateTime <= end_date)

# length(unique(filtered_data$site))
# filtered_data$DateTime %in% as.POSIXct("2014-12-01 15:15:00")

# remove Na and 0
datause <- filtered_data %>% 
  filter(!is.na(GHI_Meas) & GHI_Meas!=0)


locations <- datause %>% 
  select(lat, long)

scaled_locations <- scale(locations)

# Perform k-means clustering with 50 groups
set.seed(123)  # For reproducibility
k <- 50  # Number of clusters
kmeans_result <- kmeans(scaled_locations, centers = k, nstart = 25)

# Add the cluster assignment to the original data
locations$cluster <- kmeans_result$cluster
table(locations$cluster)

datause$group <- kmeans_result$cluster

# datause %>%
#   filter(group==50) %>% 
#   select(site) %>% 
#   unique %>% 
#   nrow

colnames(datause)
datause <- datause %>% 
  select(GHI_Meas, NAM_GHI_inline, SREF212_GHI_inline, group)


subset(datause, group==1)$GHI_Meas %>% hist


library(mvtnorm)
ex <- rmvnorm(1, mean = rep(0,50))
ex %>% as.numeric %>% dmvnorm(mean=rep(0,50))
