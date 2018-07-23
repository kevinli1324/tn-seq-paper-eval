source("pca_functions.R")

gen_t_data <- function(n_genes, n_experiments, prop_affec, mean0_vec = real_counts$set1IT002.Time0) {

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
      return_vec[i] <- ifelse(assign[i] == 1, runif(1,.15, .8), 1)
    }
    return(return_vec)
  }
  
 #t0 <- rpois(n_genes , mean0)
  
  scale <- (1/(mean(mean0_vec)*length(mean0_vec)))*sum((mean0_vec - mean(mean0_vec))^2)
  shape <- mean(mean0_vec)/scale
  
  t0 <- rgamma(n_genes, scale = scale, shape = shape) + 1
  
  
  fit_kill <- gen_fit_mixture()
  result <- t0*fit_kill
  
  gen_mixture <- function(mixing, noise_sd = 1,n, n0, num) {
    assignments <- rbinom(num, 1, 1 - mixing)
    
    return_vec <- rep(0, num)
    
    for(i in 1:num){
      return_vec[i] <- ifelse(assignments[[i]] == 0, rpois(1, n0), rpois(1, n))
    }
    
    if(n == n0) {
      assignments = rep(0, num)
    }
    
    return(list(samples = return_vec, labels = assignments))
  }
  
  mix_vec <- runif(n_genes, .3, .95)
  experiment_object <- Map(gen_mixture, mixing = mix_vec, n = result, n0 = t0, num = n_experiments)
  experiments <- lapply(experiment_object, function(x) {x[[1]]})
  experiment_labels <- lapply(experiment_object, function(x) {x[[2]]})
  control_time <- sapply(t0, function(x) {rpois(1,x)})
  
  sim_matrix <- matrix(0, nrow = n_genes, ncol = n_experiments)
  label_matrix <- matrix(0, nrow = n_genes, ncol = n_experiments)
  log_sim_matrix <- matrix(0, nrow = n_genes, ncol = n_experiments)
  
  for(i in 1:nrow(sim_matrix)) {
    logfit_sim <- sapply(experiments[[i]], log_fit, n0 =  control_time[[i]])
    logfit_dev <- sapply(experiments[[i]], naive_std_dv, n0 = control_time[[i]])
    
    log_sim_matrix[i,] <- logfit_sim
    label_matrix[i,] <- experiment_labels[[i]]
    sim_matrix[i,] <- logfit_sim/sqrt(.1**2 + logfit_dev)
  }
  
  return(list(log_data = log_sim_matrix, labels = label_matrix, t_data = sim_matrix, mix_vec = mix_vec))
  
  
}
