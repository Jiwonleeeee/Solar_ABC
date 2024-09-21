rm(list=ls())

library(Matrix)
library(plgp)
library(dplyr)
library(mvtnorm)
load("~/Downloads/solar.rdata")
colnames(train.df)
set.seed(123)  # For reproducibility

dIG <- function(x, shape, scale) {
  if (x <= 0) return(0)  # Density is zero for non-positive x
  return((scale^shape / gamma(shape)) * (x^(-shape - 1)))
}


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

train.df <- train.df %>% select(DateTime, GHI_Meas, NAM_GHI_inline, SREF212_GHI_inline,site, group)

## H = 50 / J = 3
## change the scale to log
colnames(train.df)
train.df[,c(2,3,4)] <- apply(train.df[,c(2,3,4)], 2, function(x){log(x)})

## test data for 10000?
train.df$site <- as.numeric(train.df$site)
table(train.df$site) %>% range # 12 535
test.id <- sample(1:nrow(train.df), 10000) ## random selection of 10000 y_hij
test.df <- train.df[test.id,]
train.df <- train.df[-test.id,]

## y:array: H - I - J 순, group - site - observations
H_list <- unique(train.df$group)
H <- length(unique(train.df$group))
y <- list()
X <- list()
n_h <- vector()
for(h in 1:H){
  
  y[[h]] <- list()
  X[[h]] <- list()
  set_h <- train.df %>% filter(group==H_list[h])
  n_sensor <- set_h$site %>% unique %>% length # number of sensor in each h group
  n_h[h] <- n_sensor
  sensor_list <- set_h$site %>% unique
  
  for(i in 1:n_sensor){
    y[[h]][[i]] <- set_h[set_h$site==sensor_list[i], "GHI_Meas"]
    X[[h]][[i]] <- matrix(c(rep(1, length(y[[h]][[i]])), set_h[set_h$site==sensor_list[i], "NAM_GHI_inline"], set_h[set_h$site==sensor_list[i], "SREF212_GHI_inline"]), ncol=3, byrow=F)
  } 
  # print(h)
}

## MCMC save 할 거 지정하고 그냥 식따라서 코드 쓰면 됨 이제...
iter <- 100
# level 3 model parameters
J <- 3
mu_j <- tau_j <- rho_j <- numeric(J)
mu_j_store <- tau_j_store <- rho_j_store <- matrix(0, iter, J)

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
## theta_hi 들은 그냥 regression 하면 될 것 같고..
# for h=1,...,H and i=1,...,nh 마다 regression 해보자
for(h in 1:H){
  for(i in 1:n_h[h]){
    theta[[h]][i,] <- lm(y[[h]][[i]] ~ X[[h]][[i]][,2]+X[[h]][[i]][,3]) %>% coefficients()
    sigma2_hi[[h]][i] <- lm(y[[h]][[i]] ~ X[[h]][[i]][,2]+X[[h]][[i]][,3])%>%
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




##### MCMC iteration starts? #####
for(t in 1:iter){
  
  Sigma_j_inv_store <- list()
  # mu_j, tau_j, rho_j
  for(j in 1:3){
    
    ### mu_j
    
    Sigma_j <- tau_j[j] * (exp(-rho_j[j]*D) + diag(eps, k)) # k = # of clusters (or H)
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
    log_rho_j_star <- rnorm(1, log(rho_j[j]), 1.5)
    rho_j_star <- exp(log_rho_j_star)
    
    num <- dmvnorm(beta[,j], mu_j[j]*one_vec, sigma= tau_j[j] * (exp(-rho_j_star*D)), log=T) + log(dIG(rho_j_star, 0.001, 0.001)) + log(rho_j_star)
    
    den <- dmvnorm(beta[,j], mu_j[j]*one_vec, sigma= tau_j[j] * (exp(-rho_j[j]*D)), log=T) + log(dIG(rho_j[j], 0.001, 0.001)) + log(rho_j[j])
    
    if(log(runif(1)) < (num-den)) rho_j[j] <- rho_j_star
    
    
  } # j=1,2,3 (level 3 model) done
  
  ## beta_h
  
  mu_beta <- matrix(rep(mu_j, each=H)) # HJ X 1 matrix
  Sigma_beta_stack_inv <- bdiag(Sigma_j_inv_store)
  tp <- solve(diag(sigma2_h))
  Sigma_theta_stack_inv <- bdiag(tp,tp,tp)
  
  tp <- lapply(theta, function(x){apply(x,2,mean)})
 theta_stack <-  matrix(matrix(unlist(tp), ncol=3, byrow=T), ncol=1)
 
 mu_beta_post <- (Sigma_beta_stack_inv + Sigma_theta_stack_inv) %*% (Sigma_beta_stack_inv %*% mu_beta + Sigma_theta_stack_inv %*% theta_stack)
 
 Sigma_beta_post <- solve(Sigma_beta_stack_inv + Sigma_theta_stack_inv)
 
 beta_stack <- rmvnorm(1, mu_beta_post, Sigma_beta_post)
 
  
  ## sigma2_h
  
  
  
  
  
} ## MCMC end




