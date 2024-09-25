cor(c(1,2,3,NA),c(2,2,3,4),use = "complete.obs")

## gather each y
min_length <- min(table(train.df$group))
H <- 50
y_mat <- matrix(0, H, min_length)
for(h in 1:50){
  y_mat[h,] <- train.df %>% 
    filter(group==h) %>% 
    pull(GHI_Meas) %>%
    .[1:min_length]
}
ty <- t(y_mat)

y_mat %>% 
  t %>% 
  cor %>% 
  { .[upper.tri(.) | lower.tri(.)] } %>% 
  range ## correlation is not strong


exp(-D) %>% 
  { .[upper.tri(.) | lower.tri(.)] } %>% 
  range ## correlation is not strong
