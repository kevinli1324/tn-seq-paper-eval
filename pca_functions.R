source("required_functions.R")
library(magrittr)
library(parallel)
das_model <- stan_model("simple_multivar.stan")

assign_groups <- function(data_matrix, princomps, single = F, method, clusters) {
  eval_angle <- function(vec1, vec2) {
    acos(sum(vec1*vec2)/(norm(vec1, type = "2") * norm(vec2, type = "2")))
  }
  
  svd_object <- svd(data_matrix)
  v <- svd_object$v
  u_temp <- svd_object$u

  if(nrow(data_matrix) < ncol(data_matrix)) {
    u <- matrix(0, nrow = nrow(v), ncol = nrow(v))
    u[1:nrow(u_temp), 1:ncol(u_temp)] <- u_temp
    z <-  t(u) %*% diag(svd_object$d, nrow = nrow(v))
  } else {
    u <- u_temp
    z <-  t(u) %*% diag(svd_object$d, nrow = ncol(t(u)))
  }

  if(method == "kmeans") {
    normalize_matrix <- t(apply(data_matrix, MARGIN = 1, function(x) {x/norm(x, type = "2")}))
    return_labels <- kmeans(normalize_matrix, princomps)
    return(return_labels$cluster)
  } else if(method == "r_kmeans") {
    normalize_matrix <-
      t(apply(u[, 1:princomps], MARGIN = 1, function(x) {
        if (norm(x, type = "2") == 0) {
          return(x)
        } else {
          x / norm(x, type = "2")
        }}))
    return_labels <- kmeans(normalize_matrix, clusters)
  	return(return_labels$cluster)
  }
  
  if(single) {
    return_labels <- 1:nrow(data_matrix)
    return(return_labels)
  }
  

  return_labels <- rep(0, nrow(data_matrix))
  princomp_store <- rep(0, princomps)
  for(row in 1:nrow(data_matrix)) {
      for(p in 1:princomps) {
        princomp_store[p] <- eval_angle(v[,p], u[row,])
      }
    return_labels[row] <- which(princomp_store == min(princomp_store))
  }
  return(return_labels)
}


gen_matrix <- function( mu0, mu1, theta0, sd, n = 162) {
  param_matrix <- matrix(c(mu0,mu1,theta0,sd),nrow=length(mu0))
  n_mix <- rep(1, length(mu0))#sample(1:5, nrow(param_matrix), replace=T)
  
  stor_sample <- matrix(nrow=sum(n_mix), ncol=n)
  stor_labels <- matrix(nrow=sum(n_mix), ncol=n)
  current_row = 1
  return_param_matrix <- matrix(0, nrow = sum(n_mix), ncol = 4)
  
  for (i in 1:nrow(param_matrix)){
    for(m in 1:n_mix[i]) {
      temp <- generate_sample(n,mu0[i],mu1[i],sd[i],theta0[i])
      stor_sample[current_row, ] <- 
        temp[[1]]
      stor_labels[current_row,] <- 
        temp[[2]]
      return_param_matrix[current_row, ] <- param_matrix[i,]
      current_row = current_row+1
    }
  }
  
  return(list(sample = stor_sample, labels = stor_labels, params = return_param_matrix))
}

return_partitions <- function(sample_matrix, label_matrix, princomps, params, single = F, method) {
  labels <- assign_groups(sample_matrix, princomps = princomps, single = single, method = method)
  
  return_list <- as.list(rep(NA, princomps))
  labeled_list <- as.list(rep(NA, princomps))
  param_list <- as.list(rep(NA, princomps))
  for(i in 1:princomps) {
    return_list[[i]] <- sample_matrix[which(labels == i),, drop = FALSE]
    labeled_list[[i]] <- label_matrix[which(labels == i),, drop = FALSE]
    param_list[[i]] <- params[which(labels == i),, drop = FALSE]
  }
  
  return(list(data = return_list, label = labels, label_data = labeled_list, params = params))
  
}

partition_data <- function(sample_matrix, princomps, method, single = F, clusters) {
  labels <- assign_groups(sample_matrix, princomps = princomps, single = single, method = method, clusters = clusters)
  return_list <- as.list(rep(NA, princomps))
  
  for(i in 1:princomps) {
    return_list[[i]] <- sample_matrix[which(labels == i),, drop = FALSE]
  }
  return(list(data = return_list, label = labels))
}

run_stan <- function(data, file, mean0 = 1) {
  if(nrow(data) == 0) {
    return(NA)
  } else if(is.null(data)) {
    return(NA)
  }
  stanFeed <- list(N = nrow(data ),  J = ncol(data),  y= data, t0 = mean0)
  fit <- stan(file = file, iter = 500, chains = 4, control = list(max_treedepth = 15), data = stanFeed)
  return(fit)
}



eval_stan_model <- function(stor_sample, stor_labels, stan_object) {
  
  if(class(stan_object) == "logical")  {
    return_frame <- data.frame(
      rank_entropy = numeric(), 
      rank_perf = numeric(),
      rank_false = numeric(),
      rank_neg_class = numeric(),
      rank_pos_class = numeric()
    )
    m_m_aij <- NULL
    x <-
      list(
        metrics = return_frame,
        label = m_m_aij
      )
    
    return(x)
  }
  
  
  models <- apply(X = as.data.frame(stan_object),MARGIN = 2, mean)
  
  extract_fix_param <- function(model_vec, row) {
    
    mu0 <- models[paste0("mu0[", row, "]")]
    mu1 <- models[paste0("mu1[", row, "]")]
    mixes <- models[paste0("theta[",row,"]")]
    sigma <- models[paste0("sigma[",row,"]")]
    
    return(list(mu0 = mu0, mu1 = mu1, sigma = sigma, mix = mixes))
  }
  
  m_m_perf <- sapply(1:nrow(stor_labels), function(i) {eval_accuracy(stor_sample[i,], stor_labels[i,], extract_fix_param(models,i), additive = T)})

  m_m_false <- sapply(1:nrow(stor_labels), function(i) {eval_accuracy(stor_sample[i,], stor_labels[i,], extract_fix_param(models,i), additive = T, falsePos = T)})

  m_m_entropy <- sapply(1:nrow(stor_labels), function(i) {eval_entropy(stor_sample[i,], stor_labels[i,], extract_fix_param(models,i), additive = T)})

  m_m_aij <- matrix(0, nrow = nrow(stor_labels), ncol = ncol(stor_labels))
  
  for(i in 1:nrow(stor_labels)) {
    m_m_aij[i,] <- return_aij(stor_sample[i,], extract_fix_param(models, i), additive = T) 
  }
  
  class_object <- class_metrics(stor_labels, m_m_aij)
  
  return_frame <- data.frame(
    mm_entropy = m_m_entropy, 
    mm_perf = m_m_perf,
    mm_false = m_m_false,
    mm_neg_class = class_object$neg_class,
    mm_pos_class = class_object$pos_class
  )
  
  x <-
    list(
      metrics = return_frame,
      label = m_m_aij
    )
  return(x)
  
}

eval_mclust_model <- function(data, labels, mclust_store) {
  if(class(mclust_store) == "logical")  {
    return_frame <- data.frame(
      rank_entropy = numeric(), 
      rank_perf = numeric(),
      rank_false = numeric(),
      rank_neg_class = numeric(),
      rank_pos_class = numeric()
    )
    m_m_aij <- NULL
    x <-
      list(
        metrics = return_frame,
        label = m_m_aij
      )
    
    return(x)
  }
  mc_perf <- sapply(1:nrow(data), function(i) {eval_accuracy(data[i,], labels[i,], get_mclust(mclust_store[[i]]))})
  
  mc_false <-  sapply(1:nrow(data), function(i) {eval_accuracy(data[i,], labels[i,], get_mclust(mclust_store[[i]]), falsePos = T)})
  
  mc_entropy <-  sapply(1:nrow(data), function(i) {eval_entropy(data[i,], labels[i,], get_mclust(mclust_store[[i]]))})
  
  mc_aij <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  
  for(i in 1:nrow(data)) {
    mc_aij[i,] <- return_aij(data[i,], get_mclust(mclust_store[[i]]))
  }
  
  class_object <- class_metrics(labels, mc_aij)
  return_frame <- data.frame(
    mc_entropy = mc_entropy, 
    mc_perf = mc_perf,
    mc_false = mc_false,
    mc_neg_class = class_object$neg_class,
    mc_pos_class = class_object$pos_class
  )
  
  
  x <-
    list(
      metrics = return_frame,
      mc_label = mc_aij
    )
  return(x)
  
}

eval_rank <- function(data, labels) {
  if(nrow(data) ==0)  {
    return_frame <- data.frame(
      rank_entropy = numeric(), 
      rank_perf = numeric(),
      rank_false = numeric(),
      rank_neg_class = numeric(),
      rank_pos_class = numeric()
    )
    m_m_aij <- NULL
    x <-
      list(
        metrics = return_frame,
        label = m_m_aij
      )
    
    return(x)
  }
  
  rank_aij <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  
  for(i in 1:nrow(data)) {
    rank_aij[i,] <- rank_method(data[i,])
  }
  
  rank_perf <- sapply(1:nrow(data), function(i){mean(labels[i,] == ifelse(rank_aij[i,] >.5, 1, 0))})
  
  rank_false <-  sapply(1:nrow(data), function(i) {sum(labels[i, ] != ifelse(rank_aij[i,] > .5, 1,0) & ifelse(rank_aij[i,] > .5, 1,0) == 1)/sum(labels[i,])})
  
  rank_entropy <-  sapply(1:nrow(data), function(i) {sum(-labels[i,]*log(rank_aij[i,]) - (1-labels[i,])*log(1 - rank_aij[i,]))})
  
  class_object <- class_metrics(labels, rank_aij)
  return_frame <- data.frame(
    rank_entropy = rank_entropy, 
    rank_perf = rank_perf,
    rank_false = rank_false,
    rank_neg_class = class_object$neg_class,
    rank_pos_class = class_object$pos_class
  )
  
  
  x <-
    list(
      metrics = return_frame,
      mc_label = rank_aij
    )
  return(x)
  
}


run_pca_simulation <- function(gen_data, gen_labels, princomps, param_matrix = matrix(0, nrow = nrow(gen_data), ncol = 4), single = F, multicore = F, method = "princomp") {
  das_model <- stan_model("simple_multivar.stan")
  
  groups <- return_partitions(sample_matrix = gen_data,gen_labels , princomps =  princomps, param_matrix, single = single, method = method)
  data_list <- groups$data
  label_list <- groups$label_data
  
  if(multicore){
    stan_models <- mclapply(data_list, run_stan)
    
  } else {
    stan_models <- lapply(data_list, run_stan)
  }
  
  mclust_model_list <- as.list(rep(NA, length(data_list)))
  for(i in 1:length(data_list)) {
    if(nrow(data_list[[i]]) == 0) {
      mclust_model_list[[i]] <- NA
    } else {
      mclust_model_list[[i]] <- apply(data_list[[i]], MARGIN = 1, Mclust, modelNames = "E", G = 2)
    }
  }
  
  if(multicore) {
    stan_perf <- mcMap(eval_stan_model, data_list , label_list, stan_models)
    
  } else {
    stan_perf <- Map(eval_stan_model, data_list , label_list, stan_models)
  }
  
  mclust_perf <- Map(eval_mclust_model, data_list, label_list, mclust_model_list)
  
  rank_perf <- Map(eval_rank, data_list, label_list)
  
  mm_frame <-
    data.frame(
      "mm_entropy" = numeric(),
      "mm_perf" = numeric(),
      "mm_false" = numeric(),
      "mm_neg_class" = numeric(),
      "mm_pos_class" = numeric()
    )
  
  for(i in 1:length(stan_perf))  {
    mm_frame <- rbind(mm_frame, stan_perf[[i]]$metrics)
  }
  
  mc_frame <- data.frame("mc_entropy" = numeric(), "mc_perf" = numeric(), "mc_false" = numeric(), "mc_neg_class" = numeric(), "mc_pos_class" = numeric())
  
  for(i in 1:length(mclust_perf)) {
    mc_frame <- rbind(mc_frame, mclust_perf[[i]]$metrics)
  }
  
  rank_frame <- data.frame("rank_entropy" = numeric(), "rank_perf" = numeric(), "rank_false" = numeric(), "rank_neg_class" = numeric(), "rank_pos_class" = numeric())
  
  for(i in 1:length(rank_perf)) {
    rank_frame <- rbind(rank_frame, rank_perf[[i]]$metrics)
  }
  
  return_frame <- cbind(mm_frame,mc_frame, rank_frame)
  
  return(list(metrics = return_frame, labels = list(stan = stan_perf,mclust =  mclust_perf), assign = groups$label, params = groups$params, data= groups))
}


pca_wrapper <- function(data, princomps, multicore = F, method = "princomp", file, mean0 = rep(1, nrow(data)), clusters = 1) {
  groups <- partition_data(sample_matrix = data, princomps = princomps, method = method, clusters = clusters)
  assign <- groups$label
  data_list <- groups$data
  if(multicore){
    stan_models <- mclapply(data_list, run_stan, file = file, mean0 = mean0)
    
  } else {
    stan_models <- lapply(data_list, run_stan, file = file, mean0 = mean0)
  }
  
  param_list <- lapply(stan_models, function(x) {
    if(class(x)[1] == "stanfit") {
      return(apply(as.data.frame(x), MARGIN = 2, mean))
    } else {
      return(NA)
    }
  })
  
  param_vec <- c()
  if (!any(grepl("mu1", names(param_list[[1]])))) {
    for (i in 1:length(param_list)) {
      vector_rename <- param_list[[i]]
      for (param in c("theta" , "sigma")) {
        position_index <- grepl(param, names(vector_rename))
        for (j in 1:length(which(assign == i))) {
          names(vector_rename)[position_index][j] <-
            paste(param, "[", which(assign == i)[j], "]", sep = "")
        }
      }
      
      position_index_1 <-
        grepl(",1]", names(vector_rename))
      position_index_2 <-
        grepl(",2]", names(vector_rename))
      
      for (j in 1:length(which(assign == i))) {
          names(vector_rename)[position_index_1][j] <-
            paste0("mu[", which(assign == i)[j], ",", 1, "]")
          
          names(vector_rename)[position_index_2][j] <-
            paste0("mu[", which(assign == i)[j], ",", 2, "]")
      }
      param_vec <- append(param_vec, vector_rename)
    }
  } else {
    for (i in 1:length(param_list)) {
      vector_rename <- param_list[[i]]
      for (param in c("theta", "mu1", "mu0", "sigma")) {
        position_index <- grepl(param, names(vector_rename))
        for (j in 1:length(which(assign == i))) {
          names(vector_rename)[position_index][j] <-
            paste(param, "[", which(assign == i)[j], "]", sep = "")
        }
      }
      param_vec <- append(param_vec, vector_rename)
    }
    
  }
  
  

  return(list(means = param_vec, models = stan_models, assign = assign))
}