source("correlated_sim.R")
source("pca_functions.R")
source("t_sim_functions.R")
################################# GENERATE DATA ###################################################
#generate simulations list
simulation_object <- cor_wrapper(num_clusters = 30, n_experiments = 162, max_group = 20, prop_affec = .5, mean0_vec = real_counts$set1IT002.Time0)

#t_matrix
t_matrix <- simulation_object$t_matrix
#logfit matirxx
logfit_matrix <- simulation_object$log_fit
#count_matrix
count_matrix <- simulation_object$count_matrix
#labels
label_matrix <- simulation_object$label_matrix
#fitness ratios
fit_vec <- as.vector(simulation_object$fit)
#mixing probabilities 
mix_vec <- as.vector(simulation_object$mixing)
#controls 
mean0 <- as.vector(simulation_object$init)

#pca mess
pca_object <- prcomp(logfit_matrix)
princomps <-  which(pca_object$sdev^2/sum(pca_object$sdev^2) < .01)[1]

stan_object <- pca_wrapper(data = logfit_matrix, princomps = princomps, method = "kmeans", file = "simple_multivar.stan")
#stan_count <- pca_wrapper(data = count_matrix, princomps = princomps, method = "kmeans", file = "pois_multivar.stan", mean0 = mean0)
mclust_list <- apply(logfit_matrix, MARGIN = 1, Mclust, G = 2, modelNames = "E") 

#calculate aij matrices 
stan_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol  = ncol(logfit_matrix))
mclust_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol = ncol(logfit_matrix))
rank_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol = ncol(logfit_matrix))
t_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol = ncol(logfit_matrix))

for(i in 1:nrow(logfit_matrix)) {
  stan_aij[i,] <- alpha_wrap( extract_fix_param(stan_object$means,i), logfit_matrix[i,], additive = F)
  mclust_aij[i,] <- alpha_wrap(get_mclust(mclust_list[[i]]), logfit_matrix[i,], additive = F)
  rank_aij[i,] <- rank_method(logfit_matrix[i,])
  t_aij[i,] <- mapply(paper_class,logfit_matrix[i,], t_matrix[i,])
}
#calculate classification metrics
stan_metrics <- as.data.frame(class_metrics(label_matrix, stan_aij))
stan_metrics$type <- "stan"

mclust_metrics <- as.data.frame(class_metrics(label_matrix, mclust_aij))
mclust_metrics$type <- "em"


ggplot(data = stan_metrics, )


rank_metrics <- as.data.frame(class_metrics(label_matrix, rank_aij))
rank_metrics$type <- "rank"


t_metrics <- as.data.frame(class_metrics(label_matrix, t_aij))
t_metrics$type <- "t"


#calculate cross entropies
stan_entropy <- data.frame(entropy = matrix_entropy(stan_aij, label_matrix), type = "stan")
mclust_entropy <-  data.frame(entropy = matrix_entropy(mclust_aij, label_matrix), type = "em")
rank_entropy <- data.frame(entropy = matrix_entropy(rank_aij, label_matrix), type = "rank")

master_entropy <- rbind(stan_entropy, mclust_entropy, rank_entropy)
master_entropy <- mutate(master_entropy, fit_ratio = fit_vec, mix = mix_vec, control = mean0)
plot_entropy <- melt(master_entropy, id.vars = c("type", "fit_ratio", "mix", "control"))

#class 
master_class <- rbind(stan_metrics, mclust_metrics, rank_metrics, t_metrics)
master_class <- mutate(master_class, fit_ratio = rep(fit_vec,4), mix = rep(mix_vec,4), control = rep(mean0,4))
plot_class <- melt(master_class, id.vars = c("type", "fit_ratio", "mix", "control"))

summary_class <- summarise(group_by(master_class, type), mean(pos_class, na.rm = T)) 
summary_class <- summarise(group_by(master_class, type), mean(raw)) 

#entropy metrics 

#####################################################################PLOTS###########################################

#classification distributione
class_plot <- ggplot(data = filter(plot_class, variable == "raw"), aes(value, fill = type)) + geom_histogram() + facet_grid(.~ type)

class_fit <- ggplot(data = filter(plot_class, variable == "raw"), aes(x = fit_ratio, y = value, color = type)) + geom_point() + facet_grid(.~ type)

class_mix <- ggplot(data = filter(plot_class, variable == "raw"), aes(x = mix, y = value, color = type)) + geom_point() + facet_grid(.~ type)

entropy_plot <- ggplot(data = filter(plot_entropy), aes(value, fill = type)) + geom_histogram() + facet_grid(.~ type)
entropy_fit <- ggplot(data = filter(plot_entropy), aes(x = fit_ratio, y = value, color = type)) + geom_point() + facet_grid(.~ type)
ggplot(data = filter(plot_entropy), aes(x = mix, y = value, color = type)) + geom_point() + facet_grid(.~ type)
ggplot(data = filter(plot_entropy), aes(x = mean0, y = value, color = type)) + geom_point() + facet_grid(.~ type)

## view
extract_plots_testparam()
extract_plots_testparam(193, stan_object$means,data = logfit_matrix)