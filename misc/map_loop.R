library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(future)
library(future.apply)
library(loo)
library(dplyr)
document()

# Set the parallel plan to use all local cores (currently not used as future package doesn't recognize functions in prior.R)
plan(multisession)
# Defining the random seed
set.seed(123)

# Defines the upper bounds for the quantiles
rho0 <- 1 # Controls the size of kappa
a0 <- 2 # Controls the size of v
sigma_u0 <- 2 # controls standard deviation of field
sigma_epsilon0 <- 0.1 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

# Setting mean of the field

m_u <- 0
# Calculates the log prior density function of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(
  sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)

log_not_pc_prior <- log_gaussian_prior_quantile(
  sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)

# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(10, 0.01), c(10, 10), c(0, 10), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 1))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

number_of_loops <- 100 # number of iterations
maxit_MAP <- 600
number_of_weights <- 500
confidence_level <- 0.05
results <- vector("list", number_of_loops) # Pre-allocates a list for m iterations
start_time <- Sys.time()
for (i in 1:number_of_loops) {
  tryCatch(
    {
      # Simulate parameters from PC prior
      true_params <- sim_theta_pc_quantile(
        alpha = alpha, sigma_u0 = sigma_u0,
        sigma_epsilon0 = sigma_epsilon0,
        a0 = a0, rho0 = rho0, m = 1
      )
      # Extract true parameters
      log_kappa <- true_params$log_kappa
      kappa <- exp(log_kappa)
      v <- true_params$v
      log_sigma_u <- true_params$log_sigma_u
      log_sigma_epsilon <- true_params$log_sigma_epsilon
      aniso <- list(rep(kappa, n), matrix(v, n, 2))

      # Sample from noisy data
      x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso)
      A <- Matrix::Diagonal(n, 1)
      y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)

      # Calculating the MAP under each prior knowing simulated data
      # Takes around 20s with mesh size 1
      map_pc <- MAP_prior(
        log_prior = log_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
        theta0 = unlist(true_params)
      )

      map_not_pc <- MAP_prior(
        log_prior = log_not_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
        theta0 = unlist(true_params)
      )

      # Laplace approximation
      mu_Laplace_pc <- map_pc$par
      mu_Laplace_not_pc <- map_not_pc$par
      Q_Laplace_pc <- -map_pc$hessian
      Q_Laplace_not_pc <- -map_not_pc$hessian
      Covariance_Laplace_pc <- solve(Q_Laplace_pc)
      Covariance_Laplace_not_pc <- solve(Q_Laplace_not_pc)
      std_dev_estimates_Laplace_pc <- sqrt(diag(Covariance_Laplace_pc))
      std_dev_estimates_Laplace_not_pc <- sqrt(diag(Covariance_Laplace_not_pc))
      m_u <- 0
      log_posterior_pc <- function(theta) {
        log_kappa <- theta[1]
        v <- theta[2:3]
        log_sigma_u <- theta[4]
        log_sigma_epsilon <- theta[5]
        log_posterior_prior(
          log_prior = log_pc_prior,
          mesh = mesh, log_kappa = log_kappa, v = v,
          log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
          y = y, A = A, m_u = m_u
        )
      }
      log_posterior_not_pc <- function(theta) {
        log_kappa <- theta[1]
        v <- theta[2:3]
        log_sigma_u <- theta[4]
        log_sigma_epsilon <- theta[5]
        log_posterior_prior(
          log_prior = log_not_pc_prior,
          mesh = mesh, log_kappa = log_kappa, v = v,
          log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
          y = y, A = A, m_u = m_u
        )
      }
      # Takes about 4s each with mesh size (1,1) and with 100 weights. Scales linearly in number of weights
      importance_pc <- log_unnormalized_importance_weights_and_integrals(
        log_posterior_density = log_posterior_pc,
        mu_Laplace = mu_Laplace_pc, Q_Laplace = Q_Laplace_pc,
        n_weights = number_of_weights, q = confidence_level
      )
      importance_not_pc <- log_unnormalized_importance_weights_and_integrals(
        log_posterior_density = log_posterior_not_pc,
        mu_Laplace = mu_Laplace_not_pc, Q_Laplace = Q_Laplace_not_pc,
        n_weights = number_of_weights, q = confidence_level
      )
      KL <- data.frame(
        KL_unsmoothed_pc_not_pc = KL_discrete_log_unnormalized_weights(
          importance_pc$log_unnormalized_weights,
          importance_not_pc$log_unnormalized_weights
        ),
        KL_smoothed_pc_not_pc = KL_discrete_log_unnormalized_weights(
          importance_pc$log_unnormalized_weights_smoothed,
          importance_not_pc$log_unnormalized_weights_smoothed
        ),
        KL_unsmoothed_pc_smoothed_pc = KL_discrete_log_unnormalized_weights(
          importance_pc$log_unnormalized_weights,
          importance_pc$log_unnormalized_weights_smoothed
        ),
        KL_unsmoothed_not_pc_smoothed_not_pc = KL_discrete_log_unnormalized_weights(
          importance_not_pc$log_unnormalized_weights,
          importance_not_pc$log_unnormalized_weights_smoothed
        )
      )

      # CIs
      confidence_intervals_Laplace_pc <- importance_pc$confidence_intervals_Laplace
      confidence_intervals_importance_pc <- importance_pc$confidence_intervals_importance
      confidence_intervals_importance_smoothed_pc <- importance_pc$confidence_intervals_importance_smoothed
      confidence_intervals_Laplace_not_pc <- importance_not_pc$confidence_intervals_Laplace
      confidence_intervals_importance_not_pc <- importance_not_pc$confidence_intervals_importance
      confidence_intervals_importance_smoothed_not_pc <- importance_not_pc$confidence_intervals_importance_smoothed
      parameter_within_confidence_intervals_Laplace_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_Laplace_pc)
      parameter_within_confidence_intervals_importance_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_importance_pc)
      parameter_within_confidence_intervals_importance_smoothed_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_importance_smoothed_pc)
      parameter_within_confidence_intervals_Laplace_not_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_Laplace_not_pc)
      parameter_within_confidence_intervals_importance_not_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_importance_not_pc)
      parameter_within_confidence_intervals_importance_smoothed_not_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_importance_smoothed_not_pc)


      # Accumulate results
      pc_results <- list(
        MAP_estimate = map_pc$par, # a 5d vector
        MAP_value = map_pc$value,
        convergence = map_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_pc$par - unlist(true_params)), # a 5d vector
        covariance_estimate = Covariance_Laplace_pc, # a 5x5 matrix
        std_dev_estimates_Laplace = std_dev_estimates_Laplace_pc, # a 5d vector
        confidence_intervals_Laplace = confidence_intervals_Laplace_pc,
        confidence_intervals_importance = confidence_intervals_importance_pc,
        confidence_intervals_importance_smoothed = confidence_intervals_importance_smoothed_pc,
        true_parameter_within_c_interval_Laplace = parameter_within_confidence_intervals_Laplace_pc, # a 5d vector
        true_parameter_within_c_interval_importance = parameter_within_confidence_intervals_importance_pc,
        true_parameter_within_c_interval_importance_smoothed = parameter_within_confidence_intervals_importance_smoothed_pc,
        importance_pc = importance_pc
      )
      # Not-PC results
      not_pc_results <- list(
        MAP_estimate = map_not_pc$par, # a 5d vector
        MAP_value = map_not_pc$value,
        convergence = map_not_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_not_pc$par - unlist(true_params)), # a 5d vector
        covariance_estimate = Covariance_Laplace_not_pc, # a 5x5 matrix
        std_dev_estimates_Laplace = sqrt(diag(Covariance_Laplace_not_pc)), # a 5d vector
        confidence_intervals_Laplace = confidence_intervals_Laplace_not_pc,
        confidence_intervals_importance = confidence_intervals_importance_not_pc,
        confidence_intervals_importance_smoothed = confidence_intervals_importance_smoothed_not_pc,
        true_parameter_within_c_interval_Laplace = parameter_within_confidence_intervals_Laplace_not_pc, # a 5d vector
        true_parameter_within_c_interval_importance = parameter_within_confidence_intervals_importance_not_pc,
        true_parameter_within_c_interval_importance_smoothed = parameter_within_confidence_intervals_importance_smoothed_not_pc,
        importance_not_pc = importance_not_pc
        # prob_MAP_greater = not_pc_prob_map                                      # a 5d vector
      )

      # Store results
      results[[i]] <- list(
        true_params = true_params,
        pc = pc_results,
        not_pc = not_pc_results,
        KL = KL
      )
    },
    error = function(e) {
      e
    }
  )
}

results
names <- rownames(results[[1]]$pc$confidence_intervals_Laplace)


# Histogram of distances for PC and not PC and mean.
par(mfrow = c(ceiling(length(results[[1]]$pc$distance_vector) / 2), 2))
mean_distances_pc <- c()
for (i in seq_len(length(results[[1]]$pc$distance_vector))) {
  all_distances <- c()
  for (j in seq_along(results)) {
    distance <- results[[j]]$pc$distance_vector[i]
    all_distances <- c(all_distances, distance)
  }
  mean_distance_pc <- mean(all_distances)
  mean_distances_pc <- c(mean_distances_pc, mean_distance_pc)
  hist(all_distances, main = paste("Distance for", names[[i]]), xlab = "Distance")
}
par(mfrow = c(ceiling(length(results[[1]]$not_pc$distance_vector) / 2), 2))

mean_distances_not_pc <- c()
for (i in seq_len(length(results[[1]]$not_pc$distance_vector))) {
  all_distances <- c()
  for (j in seq_along(results)) {
    distance <- results[[j]]$not_pc$distance_vector[i]
    all_distances <- c(all_distances, distance)
  }
  mean_distance_not_pc <- mean(all_distances)
  mean_distances_not_pc <- c(mean_distances_not_pc, mean_distance_not_pc)
  hist(all_distances, main = paste("Distance for", names[[i]]), xlab = "Distance")
}

# We also get the mean standard deviation estimates for the Laplace approximation
std_dev_estimates_Laplace_pc <- sapply(results, function(x) x$pc$std_dev_estimates_Laplace)
std_dev_estimates_Laplace_not_pc <- sapply(results, function(x) x$not_pc$std_dev_estimates_Laplace)
mean_std_dev_pc <- rowMeans(std_dev_estimates_Laplace_pc)
mean_std_dev_not_pc <- rowMeans(std_dev_estimates_Laplace_not_pc)

names(mean_distances_pc) <- names(results[[1]]$pc$distance_vector)
names(mean_distances_not_pc) <- names(results[[1]]$not_pc$distance_vector)
print(mean_distances_pc)
print(mean_std_dev_pc)
print(mean_distances_not_pc)
print(mean_std_dev_not_pc)



# Confidence Laplace intervals. Histogram of lengths for PC and not PC and mean.
par(mfrow = c(ceiling(length(results[[1]]$pc$confidence_intervals_Laplace) / 2), 2))
mean_lengths_pc <- c()
for (i in seq_len(nrow(results[[1]]$pc$confidence_intervals_Laplace))) {
  all_lengths <- c()

  # For each run
  for (j in seq_along(results)) {
    length <- diff(results[[j]]$pc$confidence_intervals_Laplace[i, ])
    all_lengths <- c(all_lengths, length)
  }
  mean_length_pc <- mean(all_lengths)
  mean_lengths_pc <- c(mean_lengths_pc, mean_length_pc)
  hist(all_distances, main = paste("Length of Laplace_pc CI for", names[[i]]), xlab = "Length")
}
par(mfrow = c(ceiling(length(results[[1]]$not_pc$confidence_intervals_Laplace) / 2), 2))

mean_lengths_not_pc <- c()

for (i in seq_len(nrow(results[[1]]$not_pc$confidence_intervals_Laplace))) {
  all_lengths <- c()
  for (j in seq_along(results)) {
    length <- diff(results[[j]]$not_pc$confidence_intervals_Laplace[i, ])
    all_lengths <- c(all_lengths, length)
  }
  mean_length_not_pc <- mean(all_lengths)
  mean_lengths_not_pc <- c(mean_lengths_not_pc, mean_length_not_pc)
  hist(all_lengths, main = paste("Length of Laplace_not_pc CI for", names[[i]]), xlab = "Distance")
}

print(mean_lengths_pc)
print(mean_lengths_not_pc)


## KL divergence between pc and not pc. Histograms
KL <- data.frame(
  KL_unsmoothed_pc_not_pc = sapply(results, function(x) x$KL$KL_unsmoothed_pc_not_pc),
  KL_smoothed_pc_not_pc = sapply(results, function(x) x$KL$KL_smoothed_pc_not_pc),
  KL_unsmoothed_pc_smoothed_pc = sapply(results, function(x) x$KL$KL_unsmoothed_pc_smoothed_pc),
  KL_unsmoothed_not_pc_smoothed_not_pc = sapply(results, function(x) x$KL$KL_unsmoothed_not_pc_smoothed_not_pc)
)

# Histogram of the KL divergences
par(mfrow = c(2, 2))
hist(KL$KL_unsmoothed_pc_not_pc, main = "KL(unsmoothed PC, not PC)", xlab = "KL divergence")
hist(KL$KL_smoothed_pc_not_pc, main = "KL(smoothed PC, not PC)", xlab = "KL divergence")
hist(KL$KL_unsmoothed_pc_smoothed_pc, main = "KL(unsmoothed PC, smoothed PC)", xlab = "KL divergence")
hist(KL$KL_unsmoothed_not_pc_smoothed_not_pc, main = "KL(unsmoothed not PC, smoothed not PC)", xlab = "KL divergence")

# We also store and print the mean of the KL divergences

KL_mean <- data.frame(
  KL_unsmoothed_pc_not_pc = mean(KL$KL_unsmoothed_pc_not_pc),
  KL_smoothed_pc_not_pc = mean(KL$KL_smoothed_pc_not_pc),
  KL_unsmoothed_pc_smoothed_pc = mean(KL$KL_unsmoothed_pc_smoothed_pc),
  KL_unsmoothed_not_pc_smoothed_not_pc = mean(KL$KL_unsmoothed_not_pc_smoothed_not_pc)
)
print(KL_mean)

# Percentage of times the true parameter is within the confidence interval

within_ci <- data.frame()

unlist(results[[j]][[method]][paste0("true_parameter_within_c_interval_", type)])[1]
within <- results[[j]][[method]][paste0("true_parameter_within_c_interval_", type)][i] * 1
within_ci <- data.frame(matrix(ncol = 0, nrow = 5))
for (method in c("pc", "not_pc")) {
  for (type in c("Laplace", "importance_smoothed", "importance")) {
    method_results <- c()
    for (i in 1:5) {
      all_results <- c()
      for (j in seq_along(results)) {
        within <- unlist(results[[j]][[method]][paste0("true_parameter_within_c_interval_", type)])[i] * 1
        all_results <- c(all_results, within)
      }
      proportion_within <- mean(all_results)
      method_results <- c(method_results, proportion_within)
    }
    # Add the method results as a new column in the data frame
    within_ci[paste(method, type, sep = "_")] <- method_results
  }
}
# Set the row names of the data frame to the parameter names
rownames(within_ci) <- rownames(results[[1]]$pc$confidence_intervals_Laplace)

# Create a bar plot for each parameter

for (i in seq_len(nrow(within_ci))) {
  dev.new()
  midpoints <- barplot(unlist(within_ci[i, ]), main = rownames(within_ci)[i], ylab = "Proportion within CI", ylim = c(0, 1))
  text(x = midpoints, y = unlist(within_ci[i, ]) + 0.02, labels = round(unlist(within_ci[i, ]), 2), pos = 3, cex = 0.8)
}
KL_Laplace <- data.frame(
  KL_unsmoothed_pc_Laplace = sort(sapply(results, function(x) x$pc$importance_pc$KL_divergence_importance_Laplace)),
  KL_unsmoothed_not_pc_Laplace = sort(sapply(results, function(x) x$not_pc$importance_not_pc$KL_divergence_importance_Laplace)),
  KL_smoothed_pc_Laplace = sort(sapply(results, function(x) x$pc$importance_pc$KL_divergence_smoothed_importance_Laplace)),
  KL_smoothed_not_pc_Laplace = sort(sapply(results, function(x) x$not_pc$importance_not_pc$KL_divergence_smoothed_importance_Laplace))
)

# Histogram of the KL divergences between the Laplace and importance samples
par(mfrow = c(2, 2))
hist(KL_Laplace$KL_unsmoothed_pc_Laplace, main = "KL(unsmoothed PC, Laplace)", xlab = "KL divergence")
hist(KL_Laplace$KL_unsmoothed_not_pc_Laplace, main = "KL(unsmoothed not PC, Laplace)", xlab = "KL divergence")
hist(KL_Laplace$KL_smoothed_pc_Laplace, main = "KL(smoothed PC, Laplace)", xlab = "KL divergence")
hist(KL_Laplace$KL_smoothed_not_pc_Laplace, main = "KL(smoothed not PC, Laplace)", xlab = "KL divergence")

# We also store the mean of the KL divergences
KL_Laplace_mean <- data.frame(
  KL_unsmoothed_pc_Laplace = mean(KL_Laplace$KL_unsmoothed_pc_Laplace),
  KL_unsmoothed_not_pc_Laplace = mean(KL_Laplace$KL_unsmoothed_not_pc_Laplace),
  KL_smoothed_pc_Laplace = mean(KL_Laplace$KL_smoothed_pc_Laplace),
  KL_smoothed_not_pc_Laplace = mean(KL_Laplace$KL_smoothed_not_pc_Laplace)
)
print(KL_Laplace_mean)


end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
