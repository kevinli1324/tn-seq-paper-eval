source("correlated_sim.R")
source("pca_functions.R")
source("t_sim_functions.R")
################################# GENERATE DATA ###################################################
#generate simulations list6
simulation_object <- cor_wrapper(num_clusters = 30, n_experiments = 162, max_group = 5, prop_affec = .5, mean0_vec = real_counts$set1IT002.Time0)
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

stan_object <- pca_wrapper(data = logfit_matrix, princomps = 5, clusters = princomps,  method = "kmeans", file = "simple_multivar.stan", optimize = F)
#stan_count <- pca_wrapper(data = count_matrix, princomps = princomps, mexthod = "kmeans", file = "pois_multivar.stan", mean0 = mean0)
mclust_list <- apply(logfit_matrix, MARGIN = 1, Mclust, G = 2, modelNames = "E") 


# aij range calculation/ Bayesian confidence intervals
combined_frame <- combine_princomps(stan_object$models, stan_object$assign)
aij_dist <- matrix(0, nrow = nrow(logfit_matrix), ncol = 1000) 
store <- get_aij(logfit_matrix[1,1], combined_frame, 1 )

hist(get_aij(logfit_matrix[3,10], combined_frame, 1))

#calculate aij matrices 
stan_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol  = ncol(logfit_matrix))
mclust_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol = ncol(logfit_matrix))
rank_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol = ncol(logfit_matrix))
t_aij <- matrix(0, nrow = nrow(logfit_matrix), ncol = ncol(logfit_matrix))

for(i in 1:nrow(logfit_matrix)) {
  stan_aij[i,] <- alpha_wrap( extract_fix_param(stan_object$means,i), logfit_matrix[i,], additive = T)
  
  if(as.numeric(Mclust(logfit_matrix[i,], G = 1, modelNames = "E")$BIC) - 12 > as.numeric(mclust_list[[i]]$BIC)) {
    mclust_aij[i,] <- rep(1, ncol(logfit_matrix))
  } else {
    mclust_aij[i,] <- alpha_wrap(get_mclust(mclust_list[[i]]), logfit_matrix[i,], additive = F)
  }
  
  rank_aij[i,] <- rank_method(logfit_matrix[i,])
  t_aij[i,] <- mapply(paper_class,logfit_matrix[i,], t_matrix[i,])
}
#calculate classification metrics
stan_metrics <- as.data.frame(class_metrics(label_matrix, stan_aij))
stan_metrics$type <- "stan"

mclust_metrics <- as.data.frame(class_metrics(label_matrix, mclust_aij))
mclust_metrics$type <- "em"




rank_metrics <- as.data.frame(class_metrics(label_matrix, rank_aij))
rank_metrics$type <- "rank"


t_metrics <- as.data.frame(class_metrics(label_matrix, t_aij))
t_metrics$type <- "t"


#calculate cross entropies
stan_entropy <- data.frame(entropy = matrix_entropy(stan_aij, label_matrix), type = "stan")
mclust_entropy <-  data.frame(entropy = matrix_entropy(mclust_aij, label_matrix), type = "em")
rank_entropy <- data.frame(entropy = matrix_entropy(rank_aij, label_matrix), type = "rank")


master_entropy <- rbind(stan_entropy, mclust_entropy, rank_entropy)
master_entropy <- mutate(master_entropy, fit_ratio = rep(fit_vec,3), mix = rep(mix_vec,3), control = rep(mean0,3))
plot_entropy <- melt(master_entropy, id.vars = c("type", "fit_ratio", "mix", "control"))
summary_entropy <- summarise(group_by(master_entropy, type), mean(entropy, na.rm = T)) 



#class 
master_class <- rbind(stan_metrics, mclust_metrics, rank_metrics, t_metrics)
master_class <- mutate(master_class, fit_ratio = rep(fit_vec,4), mix = rep(mix_vec,4), control = rep(mean0,4))
plot_class <- melt(master_class, id.vars = c("type", "fit_ratio", "mix", "control"))

summary_class <- summarise(group_by(master_class, type), mean(pos_class, na.rm = T)) 
summary_class <- summarise(group_by(master_class, type), mean(raw)) 


summary_false <- summarise(group_by(master_class, type), mean(false_positive)) 

#entropy metrics 

#####################################################################PLOTS###########################################

#classification distributione
class_plot <- ggplot(data = filter(plot_class, variable == "raw"), aes(value, fill = type)) + geom_histogram() + facet_grid(.~ type)

ggplot(data = filter(plot_class, variable == "raw"), aes(y = value, x = type, fill= type)) + geom_boxplot() + xlab("Classifier") + ylab("Classification Rate") + ggtitle("Classification Rate for Classifiers")


class_fit <- ggplot(data = filter(plot_class, variable == "raw"), aes(x = fit_ratio, y = value, color = type)) + geom_point() + facet_grid(.~ type)

class_mix <- ggplot(data = filter(plot_class, variable == "raw"), aes(x = mix, y = value, color = type)) + geom_point() + facet_grid(.~ type)

entropy_plot <- ggplot(data = filter(plot_entropy), aes(value, fill = type)) + geom_histogram() + facet_grid(.~ type)
entropy_fit <- ggplot(data = filter(plot_entropy), aes(x = fit_ratio, y = value, color = type)) + geom_point() + facet_grid(.~ type)
ggplot(data = filter(plot_entropy), aes(x = mix, y = value, color = type)) + geom_point() + facet_grid(.~ type)
ggplot(data = filter(plot_entropy), aes(x = mean0, y = value, color = type)) + geom_point() + facet_grid(.~ type)

ggplot(data = plot_entropy, aes(x = type, y = value)) + geom_boxplot() + xlab("Classifier") + ylab("Cross Entropy") + ggtitle("Cross Entropy of Classifications")

ggplot(data = filter(plot_class, variable == "false_positive"), aes(y = value, x = type, fill = type)) + geom_boxplot() + xlab("Classifier") + ylab("FDR") + ggtitle("FDR for Classifiers")

ggplot(data = filter(plot_class, variable == "pos_class"), aes(y = value, x = type, fill = type)) + geom_boxplot() + xlab("Classifier") + ylab("FDR") + ggtitle("Positive Classification Rate")
ggplot(data = filter(plot_entropy), aes(y = value, x = type, fill = type)) + geom_boxplot() + xlab("Classifier") + ylab("Entropy") + ggtitle("Cross Entropy for Classifiers")

## view
which(stan_metrics$raw < .2)
#look ast posterior
extract_plots_testparam(54, stan_object$means, logfit_matrix)
draws <-  as.data.frame(stan_object$models[[5]])

plot(draws$`mu0[35]`, draws$`mu1[35]`, xlab = "mu0 posterior draws", ylab = "mu1 posterior draws", main = "mu1 and mu0 Posterior Draws")
hist(draws$`theta[35]`)
hist(draws$'theta[42]')
which(stan_)