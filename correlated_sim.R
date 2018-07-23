library(MASS)
library(Matrix)
library(rstan)
source("pca_functions.R")
options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)
das_model <- stan_model("simple_multivar.stan")

gen_cor_cluster <- function(n_genes, n_experiments = 160, mean0_vec, prop_affec = 1, noise = F){
  
  naive_std_dv <- function(n, n0) {
    (1/n + 1/n0)/log(2)**2
  }
  
  log_fit <- function(n, n0) {
    log2(1 + n) - log2(1 + n0)
  }
  
  gen_fit_mixture <- function() {
    assign <- rbinom(n_genes, 1, prop_affec)
    return_vec <- rep(0, n_genes)
    for(i in 1:n_genes) {
      return_vec[i] <- ifelse(assign[i] == 1, runif(1, .1, .9) , 1)
    }
    return(return_vec)
  }
  
  gen_experiment <- function(means) {
    gen_mean <- log(round(means))
    
    cov_matrix <- matrix(0, ncol = ncol(cor_matrix), nrow = nrow(cor_matrix))
    for(i in 1:nrow(cor_matrix)) {
      for(j in 1:ncol(cor_matrix)){
        if(i == j) {
          cov_matrix[i,j] <- 1
          cov_matrix[i,j] <- log((cor_matrix[i,j]*sqrt(means[i]*means[j]))/(means[i]*means[j]) + 1)
        }
      }
    }
    return_vec <- sapply(exp(mvrnorm(1, gen_mean, Sigma = cov_matrix)), function(x) rpois(1, x))
    return(return_vec)
  } 
  
  gen_column <- function(m0, m1 = rep(1, n_genes), shared_percent = rep(0, n_genes)) {
    
    control_samples <- gen_experiment(m0)
    non_control <- gen_experiment(m1)
    
    return_vec <- rep(0, n_genes)
    indice <- sapply(shared_percent, function(x) {rbinom(1, 1, x)})
    
    for(i in 1:n_genes) {
      if(m0[i] == m1[i]) {
        indice[i] <- 0
      }
      return_vec[i] <- ifelse(indice[i] == 0, control_samples[i], non_control[i])
    }
    
    return(list(values = return_vec, labels = indice))
  }
  
  
  cor_matrix <- matrix(runif(n_genes*n_genes,.8,1), nrow = n_genes, ncol = n_genes)
  cor_matrix <- as.matrix(forceSymmetric(cor_matrix))
  
  
  
  #t0 <- rpois(n_genes , mean0)
  scale <- (1/(mean(mean0_vec)*length(mean0_vec)))*sum((mean0_vec - mean(mean0_vec))^2)
  shape <- mean(mean0_vec)/scale
  
  t0 <- rgamma(n_genes, scale = scale, shape = shape) + 1
  fit_kill <- gen_fit_mixture()
  result <- t0*fit_kill
  
  
  mixing <- runif(1 ,.5, .95)
  affected_cols <-  sample(1:n_experiments, n_experiments*mixing)
  shared_percent <- runif(n_genes, .5, .95) 
  control_time <- sapply(t0, function(x) {rpois(1,x)})
  
  count_matrix <- matrix(0, nrow = n_genes, ncol = n_experiments)
  label_matrix <-  matrix(0, nrow = n_genes, ncol = n_experiments)
  
  for(i in 1:n_experiments) {
    if(i %in% affected_cols) {
      temp <- gen_column(t0, result, shared_percent)
    } else {
      temp <- gen_column(t0)
    }
    count_matrix[,i] <- temp$values
    label_matrix[,i] <- temp$labels
  }
  
  sim_matrix <- matrix(0, nrow = n_genes, ncol = n_experiments)
  log_sim_matrix <- matrix(0, nrow = n_genes, ncol = n_experiments)
  
  for(i in 1:nrow(sim_matrix)) {
    logfit_sim <- sapply(count_matrix[i,], log_fit, n0 =  control_time[i])
    logfit_dev <- sapply(count_matrix[i,], naive_std_dv, n0 = control_time[i])
    
    log_sim_matrix[i,] <- logfit_sim
    sim_matrix[i,] <- logfit_sim/sqrt(.1**2 + logfit_dev)
  }
  return(
    list(
      t_matrix = sim_matrix,
      log_fit  = log_sim_matrix,
      label_matrix = label_matrix,
      count_matrix  = count_matrix,
      fit = matrix(fit_kill, ncol = 1),
      mixing = matrix(mixing * shared_percent, ncol = 1), 
      init = matrix(control_time, ncol = 1)
    )
  )
}

cor_wrapper <- function(num_clusters, n_experiments, max_group, prop_affec, mean0_vec) {
  collapse_func <- function(init, new_list) {
    for(i in 1:length(init)) {
      init[[i]] <- rbind(init[[i]], new_list[[i]])
    }
    return(init)
  }
  cluster_num <- sample(1:max_group, num_clusters, replace = T)
  
  list_o_lists <- lapply(cluster_num, gen_cor_cluster, n_experiments = n_experiments, prop_affec = prop_affec, mean0_vec = mean0_vec)
  
  return_list <- Reduce(collapse_func, list_o_lists)
  return(return_list)
}
