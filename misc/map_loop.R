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
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 2 # control standard deviation of noise
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
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(5, 5))
nodes <- mesh$loc
n <- mesh$n
par(mfrow = c(1, 1))
plot(mesh)

number_of_loops <- 2 # number of iterations
maxit_MAP <- 600
number_of_weights <- 100
confidence_level <- 0.05
results <- vector("list", number_of_loops) # Pre-allocates a list for m iterations

for (i in 1:number_of_loops) {
  start_time <- Sys.time()
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
      x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
      # To only observe the field at some of the nodes we set A to be rectangualr of size mxn
      m <- round(n / 4)
      A <- matrix(0, m, n)
      A[1:m, 1:m] <- diag(m)
      y <- A %*% x + rnorm(m, 0, exp(log_sigma_epsilon))



      # Calculating the MAP under each prior knowing simulated data
      # Takes around 20s with mesh size 1 (554 degrees of freedom) and scales linearly in degrees of freedom
      delta <- rnorm(5, 0, 1) # Used to randomize starting point of MAP
      map_pc <- MAP_prior(
        log_prior = log_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
        theta0 = unlist(true_params) + delta
      )

      map_not_pc <- MAP_prior(
        log_prior = log_not_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
        theta0 = unlist(true_params) + delta
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
      # Takes about 4s each with mesh size (1,1) and with 100 weights. Scales linearly in number of weights and
      # degrees of freedom
      importance_pc <- log_unnormalized_importance_weights_and_integrals(
        log_posterior_density = log_posterior_pc,
        mu_Laplace = mu_Laplace_pc, Q_Laplace = Q_Laplace_pc,
        n_weights = number_of_weights, q = confidence_level, true_params = unlist(true_params)
      )
      importance_not_pc <- log_unnormalized_importance_weights_and_integrals(
        log_posterior_density = log_posterior_not_pc,
        mu_Laplace = mu_Laplace_not_pc, Q_Laplace = Q_Laplace_not_pc,
        n_weights = number_of_weights, q = confidence_level, true_params = unlist(true_params)
      )
      # KL <- data.frame(
      #   KL_unsmoothed_pc_not_pc = KL_discrete_log_unnormalized_weights(
      #     importance_pc$log_unnormalized_weights,
      #     importance_not_pc$log_unnormalized_weights
      #   ),
      #   KL_smoothed_pc_not_pc = KL_discrete_log_unnormalized_weights(
      #     importance_pc$log_unnormalized_weights_smoothed,
      #     importance_not_pc$log_unnormalized_weights_smoothed
      #   ),
      #   KL_unsmoothed_pc_smoothed_pc = KL_discrete_log_unnormalized_weights(
      #     importance_pc$log_unnormalized_weights,
      #     importance_pc$log_unnormalized_weights_smoothed
      #   ),
      #   KL_unsmoothed_not_pc_smoothed_not_pc = KL_discrete_log_unnormalized_weights(
      #     importance_not_pc$log_unnormalized_weights,
      #     importance_not_pc$log_unnormalized_weights_smoothed
      #   )
      # )

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
        not_pc = not_pc_results
        # ,KL = KL
      )
    },
    error = function(e) {
      e
    }
  )
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  formatted_time <- seconds_to_hours_minutes_seconds(execution_time)
  print(paste("Iteration", i, "took", formatted_time))
  time_left <- as.numeric(difftime(end_time, start_time, units = "secs")) * (number_of_loops - i)
  print(paste("Estimated time left:", seconds_to_hours_minutes_seconds(time_left)))
}
# Eliminates NULL results
not_null_indices <- sapply(results, function(x) !is.null(x$pc$importance_pc$log_unnormalized_weights))
results <- results[not_null_indices]
prior_types <- c("pc", "not_pc")
approximation_types <- c("Laplace", "importance", "importance_smoothed")
parameter_names <- rownames(results[[1]]$pc$confidence_intervals_Laplace)

mean_distances <- list()
mean_std_dev <- list()

for (prior_type in prior_types) {
  par(mfrow = c(3, 2))
  mean_distances[[prior_type]] <- c()
  for (i in seq_len(length(results[[1]][[prior_type]]$distance_vector))) {
    all_distances <- c()
    for (j in seq_along(results)) {
      distance <- results[[j]][[prior_type]]$distance_vector[i]
      all_distances <- c(all_distances, distance)
    }
    mean_distance <- mean(all_distances)
    mean_distances[[prior_type]] <- c(mean_distances[[prior_type]], mean_distance)
    hist(all_distances, main = paste("Distance to MAP for", parameter_names[[i]], prior_type), xlab = "Distance")
  }

  std_dev_estimates_Laplace <- do.call(rbind, lapply(results, function(x) x[[prior_type]]$std_dev_estimates_Laplace))
  mean_std_dev[[prior_type]] <- colMeans(std_dev_estimates_Laplace)

  names(mean_distances[[prior_type]]) <- names(results[[1]][[prior_type]]$distance_vector)
  print(paste("Mean distances for", prior_type))
  print(mean_distances[[prior_type]])
  print(paste("Mean standard deviations for", prior_type))
  print(mean_std_dev[[prior_type]])
}



# Credible intervals. Histogram of lengths.


calculate_lengths <- function(results, prior_type, approximation_type) {
  ci_type <- paste0("confidence_intervals_", approximation_type)
  mean_lengths <- c()

  for (i in seq_along(parameter_names)) {
    all_lengths <- c()

    for (j in seq_along(results)) {
      length <- diff(results[[j]][[prior_type]][[ci_type]][i, ])
      all_lengths <- c(all_lengths, length)
    }

    mean_length <- mean(all_lengths)
    mean_lengths <- c(mean_lengths, mean_length)
    hist(all_lengths, main = paste("Length of", ci_type, "of", parameter_names[[i]], "for", prior_type), xlab = "Length")
  }

  print(paste("Mean lengths for", prior_type, ci_type))
  print(mean_lengths)
}


# Set layout
layout_matrix <- rbind(c(1, 2, 3), c(4, 5, 0))
layout(layout_matrix)

for (prior_type in prior_types) {
  for (approximation_type in approximation_types) {
    calculate_lengths(results, prior_type, approximation_type)
  }
}

# Percentage of times the true parameter is within the confidence interval

within_ci <- data.frame()

within_ci <- data.frame(matrix(ncol = 0, nrow = 5))
for (prior_type in prior_types) {
  for (approximation_type in approximation_types) {
    prior_type_results <- c()
    for (i in 1:5) {
      all_results <- c()
      for (j in seq_along(results)) {
        within <- unlist(results[[j]][[prior_type]][paste0("true_parameter_within_c_interval_", approximation_type)])[i] * 1
        all_results <- c(all_results, within)
      }
      proportion_within <- mean(all_results)
      prior_type_results <- c(prior_type_results, proportion_within)
    }
    # Add the prior_type results as a new column in the data frame
    within_ci[paste(prior_type, approximation_type, sep = "_")] <- prior_type_results
  }
}
# Set the row names of the data frame to the parameter names
rownames(within_ci) <- rownames(results[[1]]$pc$confidence_intervals_Laplace)

# Create a bar plot for each parameter

layout(layout_matrix)

for (i in seq_len(nrow(within_ci))) {
  midpoints <- barplot(unlist(within_ci[i, ]), main = rownames(within_ci)[i], ylab = "Proportion within CI", ylim = c(0, 1))
  text(x = midpoints, y = unlist(within_ci[i, ]) + 0.02, labels = round(unlist(within_ci[i, ]), 2), pos = 3, cex = 0.8)
}
KL_Laplace <- data.frame(
  KL_unsmoothed_pc_Laplace = sort(unlist(sapply(results, function(x) x$pc$importance_pc$KL_divergence_importance_Laplace))),
  KL_unsmoothed_not_pc_Laplace = sort(unlist(sapply(results, function(x) x$not_pc$importance_not_pc$KL_divergence_importance_Laplace))),
  KL_smoothed_pc_Laplace = sort(unlist(sapply(results, function(x) x$pc$importance_pc$KL_divergence_smoothed_importance_Laplace))),
  KL_smoothed_not_pc_Laplace = sort(unlist(sapply(results, function(x) x$not_pc$importance_not_pc$KL_divergence_smoothed_importance_Laplace)))
)

# Histogram of the KL divergences between the Laplace and importance samples
par(mfrow = c(2, 2))
hist(KL_Laplace$KL_unsmoothed_pc_Laplace, main = "KL(unsmoothed PC, Laplace)", xlab = "KL divergence")
hist(KL_Laplace$KL_unsmoothed_not_pc_Laplace, main = "KL(unsmoothed not PC, Laplace)", xlab = "KL divergence")
hist(KL_Laplace$KL_smoothed_pc_Laplace, main = "KL(smoothed PC, Laplace)", xlab = "KL divergence")
hist(KL_Laplace$KL_smoothed_not_pc_Laplace, main = "KL(smoothed not PC, Laplace)", xlab = "KL divergence")
ggplot(rbind(
  data.frame(model = "PC", generator = "Laplace", IS = "unsmoothed", KL = KL_Laplace$KL_unsmoothed_pc_Laplace),
  data.frame(model = "PC", generator = "Laplace", IS = "smoothed", KL = KL_Laplace$KL_smoothed_pc_Laplace),
  data.frame(model = "NOT_PC", generator = "Laplace", IS = "unsmoothed", KL = KL_Laplace$KL_unsmoothed_not_pc_Laplace),
  data.frame(model = "NOT_PC", generator = "Laplace", IS = "smoothed", KL = KL_Laplace$KL_smoothed_not_pc_Laplace)
)) +
  stat_ecdf(aes(KL, col = model, linetype = IS))
# +  facet_wrap(vars(IS, model))
# We also store the mean of the KL divergences
KL_Laplace_mean <- data.frame(
  KL_unsmoothed_pc_Laplace = mean(KL_Laplace$KL_unsmoothed_pc_Laplace),
  KL_unsmoothed_not_pc_Laplace = mean(KL_Laplace$KL_unsmoothed_not_pc_Laplace),
  KL_smoothed_pc_Laplace = mean(KL_Laplace$KL_smoothed_pc_Laplace),
  KL_smoothed_not_pc_Laplace = mean(KL_Laplace$KL_smoothed_not_pc_Laplace)
)
print(KL_Laplace_mean)

# We check whether the quantiles P[theta<=theta^*] are uniformly distributed
par(mfrow = c(2, 3))
for (prior_type in prior_types) {
  for (approximation_type in approximation_types) {
    probabilities <- sapply(results, function(x) x[[prior_type]][[paste0("importance_", prior_type)]][[paste0("probabilities_", approximation_type)]])
    for (i in length(parameter_names)) {
      hist(probabilities[i, ], main = paste("Histogram of P[theta<=MAP] for", parameter_names[[i]], "for", prior_type, approximation_type), xlab = "P[theta<=theta^*]")
      # We also get the mean of the probabilities
      print(paste("Mean of P[theta<=theta^*] for", parameter_names[[i]], "for", prior_type, approximation_type))
      print(mean(probabilities[i, ]))
    }
  }
}



# # Average variance of log unnormalized weights
# var_weights <- sapply(results, function(x) var(normalize_log_weights(x$pc$importance_pc$log_unnormalized_weights)))
# hist(var_weights, main = "Histogram of variances for PC", xlab = "Variance")
#
# var_weights_smoothed <- sapply(results, function(x) var(normalize_log_weights(x$pc$importance_pc$log_unnormalized_weights_smoothed)))
# hist(var_weights_smoothed, main = "Histogram of variances for PC smoothed", xlab = "Variance")


# We get the vector of k diagnostics and show the log weights are similar
par(mfrow = c(1, 1))
k_diagnostics <- sapply(results, function(x) x$pc$importance_pc$k_diagnostic)
hist(k_diagnostics, main = "Histogram of k diagnostics", xlab = "k diagnostic")
par(mfrow = c(1, 2))
j <- 1
hist(results[[j]]$pc$importance_pc$log_unnormalized_weights_smoothed, main = "Log weights", xlab = "Log weight")
hist(results[[j]]$pc$importance_pc$log_unnormalized_weights, main = "Log weights", xlab = "Log weight")
