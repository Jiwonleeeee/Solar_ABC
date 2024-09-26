rm(list=ls())

library(Matrix)
library(plgp)
library(dplyr)
library(mvtnorm)
load("~/Downloads/solar.rdata")
colnames(train.df)
set.seed(123)  # For reproducibility

dIG <- function(x, shape, scale) {
  if (x <= 0) return(0)
  (scale^shape / gamma(shape)) * x^(-shape-1) * exp(-scale / x)
}


start_date <- as.POSIXct("2014-12-01")
end_date <- as.POSIXct("2014-12-15")

# Filter data within the period
train.df <- train.df %>%
  filter(DateTime >= start_date & DateTime <= end_date)

## Remove NA and 0 
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

train.df <- train.df %>% select(DateTime, GHI_Meas, NAM_GHI_inline, SREF212_GHI_inline,site, group)

## H = 50 / J = 3
## change the scale to log
colnames(train.df)
train.df[,c(2,3,4)] <- apply(train.df[,c(2,3,4)], 2, function(x){log(x)})
train.df$site <- as.numeric(train.df$site)

test_proportion <- 0.1
split_data <- train.df %>%
  mutate(row_id = row_number()) %>%  # Create a unique row ID within each group
  group_by(site, group) %>%          # Group again for sampling within each combination
  mutate(is_test = row_id %in% sample(row_id, size = ceiling(n() * test_proportion))) %>% 
  ungroup()

# Separate into test and train sets
test.df <- split_data %>% filter(is_test) %>% select(-is_test, -row_id)  # Test set
train.df <- split_data %>% filter(!is_test) %>% select(-is_test, -row_id)  # Train set


## y:array: H - I - J, group - site - observations
H_list <- unique(train.df$group)
H <- length(unique(train.df$group))
y <- list()
X <- list()
n_h <- vector()
sens_list <- list()
for(h in 1:H){
  
  y[[h]] <- list()
  X[[h]] <- list()
  set_h <- train.df %>% filter(group==H_list[h])
  n_sensor <- set_h$site %>% unique %>% length # number of sensor in each h group
  n_h[h] <- n_sensor
  sensor_list <- set_h$site %>% unique
  sens_list[[h]] <- sensor_list
  for(i in 1:n_sensor){
    y[[h]][[i]] <- set_h %>%
      filter(site == sensor_list[i]) %>%
      select(GHI_Meas, site, group) %>% 
      as.matrix()
    X[[h]][[i]] <- matrix(c(rep(1, nrow(y[[h]][[i]]) ), 
                            set_h %>% filter(site == sensor_list[i]) %>% pull(NAM_GHI_inline),
                            set_h %>% filter(site == sensor_list[i]) %>% pull(SREF212_GHI_inline)), ncol=3, byrow=F)
  } 
  # print(h)
}

# use
# y[[1]][[1]][,1] : just to know what site and group are assigned to each y and X

iter_total <- 50000
iter <- seq(1, iter_total, by=10) %>% length
# level 3 model parameters
J <- 3
mu_j <- tau_j <- rho_j <- numeric(J)
mu_j_store <- tau_j_store <- rho_j_store <- matrix(0, iter, J)

rho_j_hyp_mu <- c(-1.897,  -2.995, -2.995)
rho_j_hyp_sig <- c(1.196, 1.561, 1.561)

beta <- matrix(0, H, J)
beta_store <- array(NA, dim = c(iter, H, J))

sigma2_h <- numeric(H)
sigma2_h_store <- matrix(0, iter, H)

theta <- sigma2_hi <- list()
for(h in 1:H){
  theta[[h]] <- matrix(0, n_h[h], 3)
  sigma2_hi[[h]] <- numeric(n_h[h])
}
theta_store <- sigma2_hi_store <- list()

##### initial setting #####
## theta_hi : regression 
# for h=1,...,H and i=1,...,nh 
for(h in 1:H){
  for(i in 1:n_h[h]){
    theta[[h]][i,] <- lm(y[[h]][[i]][,1] ~ X[[h]][[i]][,2]+X[[h]][[i]][,3]) %>% coefficients()
    sigma2_hi[[h]][i] <- lm(y[[h]][[i]][,1] ~ X[[h]][[i]][,2]+X[[h]][[i]][,3])%>%
      summary() %>%
      {.$sigma^2}
  }
}

# beta -> mean of theta로 하자 
# sigma_h -> sample variance
for(h in 1:H){
  beta[h,] <- apply(theta[[h]], 2, mean)
  sigma2_h[h] <- as.vector(theta[[h]]) %>% var
}

## mu_j, tau_j, rho_j
mu_j <- apply(beta, 2, mean)
tau_j <- apply(beta, 2, var) ## tau^2
rho_j <- (apply(beta, 2, range)[2,]-apply(beta, 2, range)[1,])/10

##### initial setting done #####
center <- kmeans_result$centers
D <- distance(center) %>% sqrt # distance matrix
eps <- sqrt(.Machine$double.eps)


current_time <- Sys.time()
cnt <- 0
##### MCMC iteration starts? #####
for(t in 1:iter_total){
  
  Sigma_j_inv_store <- list()
  # mu_j, tau_j, rho_j
  for(j in 1:3){
    
    ### mu_j
    
    Sigma_j <- tau_j[j] * (exp(-D/rho_j[j]) + diag(eps, k)) # k = # of clusters (or H)
    one_vec <- matrix(1, nrow = 50, ncol = 1)
    Sigma_j_inv <- solve(Sigma_j)
    Sigma_j_inv_store[[j]] <- Sigma_j_inv
    
  mean <- c((t(one_vec) %*% Sigma_j_inv %*% beta[,j])/(t(one_vec) %*% Sigma_j_inv %*% one_vec))
    var <- c(1/(t(one_vec) %*% Sigma_j_inv %*% one_vec))
    
    mu_j[j] <- rnorm(1, mean, sqrt(var))
    
    ### tau_j
    alpha_0 <- 0.001
    beta_0 <- 0.001
    
    g_shape <- alpha_0 + H/2
    g_rate <- c((2*beta_0 + t(beta[,j]-mu_j[j]*one_vec)%*% (tau_j[j]*Sigma_j_inv) %*% (beta[,j]-mu_j[j]*one_vec))/2 )
    tau_j[j] <- 1/rgamma(1, shape = g_shape, rate = g_rate)
    
    # rho_j : no close form
    log_rho_j_star <- rnorm(1, log(rho_j[j]), 1)
    rho_j_star <- exp(log_rho_j_star)
    
    num <- dmvnorm(beta[,j], mu_j[j]*one_vec, sigma= tau_j[j] * (exp(-D/rho_j_star)), log=T) + dlnorm(rho_j_star, rho_j_hyp_mu[j], rho_j_hyp_sig[j], log=T) + log(rho_j_star)
    
    den <- dmvnorm(beta[,j], mu_j[j]*one_vec, sigma= tau_j[j] * (exp(-D/rho_j[j])), log=T) + dlnorm(rho_j[j], rho_j_hyp_mu[j], rho_j_hyp_sig[j], log=T) + log(rho_j[j])
    
    if(log(runif(1)) < (num-den)) rho_j[j] <- rho_j_star
    
    Sigma_j <- tau_j[j] * (exp(-D/rho_j[j]) + diag(eps, k))
    Sigma_j_inv <- solve(Sigma_j)
    Sigma_j_inv_store[[j]] <- Sigma_j_inv
  } # j=1,2,3 (level 3 model) done
  
  ## beta_h
  
  mu_beta <- matrix(rep(mu_j, each=H)) # HJ X 1 matrix
  Sigma_beta_stack_inv <- bdiag(Sigma_j_inv_store)
  tp <- solve(diag(sigma2_h))
  Sigma_theta_stack_inv <- bdiag(tp,tp,tp)
  
  tp <- lapply(theta, function(x){apply(x,2,mean)})
 theta_stack <-  matrix(matrix(unlist(tp), ncol=3, byrow=T), ncol=1)
 
 mu_beta_post <- solve((Sigma_beta_stack_inv + Sigma_theta_stack_inv)) %*% (Sigma_beta_stack_inv %*% mu_beta + Sigma_theta_stack_inv %*% theta_stack)
 
 Sigma_beta_post <- solve(Sigma_beta_stack_inv + Sigma_theta_stack_inv)
 
 Sigma_beta_post <- as.matrix(Sigma_beta_post)
 beta_stack <- rmvnorm(1, mu_beta_post, Sigma_beta_post)
 
 beta <- matrix(beta_stack, nrow= H, ncol=J, byrow=F)
 
  
  ## sigma2_h
 
  alpha_h <- beta_h <- 0.001
  
  for(h in 1:H){
    
    alpha_sigma2_h <- alpha_h + 3*n_h[h]/2
    

    tp <- theta[[h]] - matrix(rep(beta[h,], n_h[h]), ncol=3, byrow=T)
    s <- apply(tp, 1, function(u){sum(u^2)}) %>% sum
    
    # sum <- 0
    # for(i in 1:n_h[h]){
    #   sum <- sum +  c(t(theta[[h]][i,]-  beta[h,]) %*% (theta[[h]][i,]-  beta[h,]) )
    # }
    beta_sigma2_h <- 2*beta_h + s
    
    sigma2_h[h] <- 1/rgamma(1, shape=alpha_sigma2_h, rate=beta_sigma2_h)
    
  }

  
  ### level 2 done ###
  
  ## theta_hi and sigma_hi
  alpha_hi <- beta_hi <- 0.001
  for(h in 1:H){
    for(i in 1:n_h[h]){
      
      Sigma_theta <- solve( t(X[[h]][[i]])%*%(X[[h]][[i]])/sigma2_hi[[h]][i] + diag(3)/sigma2_h[h])
      
      mean_theta <- Sigma_theta %*% ( t(X[[h]][[i]])%*%y[[h]][[i]][,1]/sigma2_hi[[h]][i] + beta[h,]/sigma2_h[h]  )
      
      
      theta[[h]][i,] <- rmvnorm(1, mean_theta, Sigma_theta)
      
      alpha_sigma <- alpha_hi + length(y[[h]][[i]][,1])/2
      
      beta_sigma <- (2 * beta_hi + t(y[[h]][[i]][,1]-X[[h]][[i]]%*%theta[[h]][i,]) %*% (y[[h]][[i]][,1]-X[[h]][[i]]%*%theta[[h]][i,]) )/2 
      
      sigma2_hi[[h]][i] <- 1/rgamma(1, shape = alpha_sigma, rate=beta_sigma)
      
    }
  }
  
  ## save
  if(t%%10 == 1){
    
    cnt <- cnt + 1
 
    mu_j_store[cnt, ] <- mu_j
    tau_j_store[cnt, ] <- tau_j
    rho_j_store[cnt, ] <- rho_j
    
    beta_store[cnt, , ] <- beta
    sigma2_h_store[cnt, ] <- sigma2_h
    
    theta_store[[cnt]] <- theta
    sigma2_hi_store[[cnt]] <- sigma2_hi
  }
  
  print(t)
  
  
  
  
} ## MCMC end

Sys.time() - current_time


