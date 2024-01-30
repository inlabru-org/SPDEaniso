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

log_priors <- list(
  pc = log_pc_prior,
  not_pc = log_not_pc_prior,
  uniform = log_uniform_prior
)
prior_types <- setNames(as.list(names(log_priors)), names(log_priors))
approximation_types <- list("Gaussian_median", "importance", "importance_smoothed")
approximation_types <- setNames(approximation_types, approximation_types)

log_not_pc_prior <- log_gaussian_prior_quantile(
  sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)
L <- 10 # Length of domain
log_uniform_prior <- log_prior_uniform(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L)

# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(2, 2))
nodes <- mesh$loc
n <- mesh$n
par(mfrow = c(1, 1))
plot(mesh)

number_of_loops <- 500 # number of iterations
maxit_MAP <- 600
number_of_weights <- 1000
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
      # To only observe the field at some of the nodes we set A to be rectangular of size m x n
      m <- round(n / 10)
      A <- matrix(0, m, n)
      A[1:m, 1:m] <- diag(m)
      y <- A %*% x + rnorm(m, 0, exp(log_sigma_epsilon))



      # Calculating the MAP under each prior knowing simulated data
      # Takes around 20s with mesh size 1 (554 degrees of freedom) and scales linearly in degrees of freedom
      delta <- rnorm(5, 0, 1) # Used to randomize starting point of MAP
      maps <- lapply(log_priors, function(log_prior) {
        MAP_prior(
          log_prior = log_prior, mesh = mesh,
          y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
          theta0 = unlist(true_params) + delta
        )
      })


      # Gaussian_median approximations
      mus_Gaussian_median <- lapply(maps, function(map) map$par)
      Qs_Gaussian_median <- lapply(maps, function(map) -map$hessian)
      Covariances_Gaussian_median <- lapply(Qs_Gaussian_median, function(Q) solve(Q))
      std_dev_estimates_Gaussian_median <- lapply(Covariances_Gaussian_median, function(Covariance) sqrt(diag(Covariance)))

      m_u <- 0

      log_posteriors <- lapply(log_priors, function(log_prior) {
        function(theta) {
          log_kappa <- theta[1]
          v <- theta[2:3]
          log_sigma_u <- theta[4]
          log_sigma_epsilon <- theta[5]
          log_posterior_prior(
            log_prior = log_prior,
            mesh = mesh, log_kappa = log_kappa, v = v,
            log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
            y = y, A = A, m_u = m_u
          )
        }
      })

      # Importance sampling
      importances <- lapply(prior_types, function(prior_type) {
        log_posterior <- log_posteriors[[prior_type]]
        mu_Gaussian_median <- mus_Gaussian_median[[prior_type]]
        Q_Gaussian_median <- Qs_Gaussian_median[[prior_type]]
        log_unnormalized_importance_weights_and_integrals(
          log_posterior_density = log_posteriors[[prior_type]],
          mu_Gaussian_median = mu_Gaussian_median, Q_Gaussian_median = Q_Gaussian_median,
          n_weights = number_of_weights, q = confidence_level, true_params = unlist(true_params)
        )
      })

      # CIs
      confidence_intervals <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
          importances[[prior_type]][[paste0("confidence_intervals_", approximation_type)]]
        })
      })

      true_parameter_is_within_CI <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
          parameter_within_confidence_intervals(true_params, confidence_intervals[[prior_type]][[approximation_type]])
        })
      })


      # Accumulate results
      results_accumulator <- function(prior_type) {
        # Return a list of the calculated values
        list(
          MAP_estimate = maps[[prior_type]]$par,
          MAP_value = maps[[prior_type]]$value,
          convergence = maps[[prior_type]]$convergence,
          distance_vector = abs(maps[[prior_type]]$par - unlist(true_params)),
          covariance_estimate = Covariances_Gaussian_median[[prior_type]],
          std_dev_estimates_Gaussian_median = std_dev_estimates_Gaussian_median[[prior_type]],
          confidence_intervals = lapply(approximation_types, function(approximation_type) {
            confidence_intervals[[prior_type]][[approximation_type]]
          }),
          true_parameter_within_c_interval = lapply(approximation_types, function(approximation_type) {
            true_parameter_is_within_CI[[prior_type]][[approximation_type]]
          }),
          importance = importances[[prior_type]]
        )
      }


      # Store results
      partial_results <- lapply(prior_types, results_accumulator)
      results[[i]] <- c(
        list(true_params = true_params),
        partial_results, prior_types
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
not_null_indices <- sapply(results, function(x) !is.null(x$pc$importance$log_unnormalized_weights))
results <- results[not_null_indices]
parameter_names <- rownames(results[[1]]$pc$confidence_intervals_Gaussian_median)

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

  std_dev_estimates_Gaussian_median <- do.call(rbind, lapply(results, function(x) x[[prior_type]]$std_dev_estimates_Gaussian_median))
  mean_std_dev[[prior_type]] <- colMeans(std_dev_estimates_Gaussian_median)

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

within_ci <- data.flog_priors

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
rownames(within_ci) <- rownames(results[[1]]$pc$confidence_intervals_Gaussian_median)

# Create a bar plot for each parameter

layout(layout_matrix)

for (i in seq_len(nrow(within_ci))) {
  midpoints <- barplot(unlist(within_ci[i, ]), main = rownames(within_ci)[i], ylab = "Proportion within CI", ylim = c(0, 1))
  text(x = midpoints, y = unlist(within_ci[i, ]) + 0.02, labels = round(unlist(within_ci[i, ]), 2), pos = 3, cex = 0.8)
}
KL <- data.frame(
  KL_unsmoothed_pc_Gaussian_median = sort(unlist(sapply(results, function(x) x$pc$importance$KL_divergence_importance_Gaussian_median))),
  KL_unsmoothed_not_pc_Gaussian_median = sort(unlist(sapply(results, function(x) x$not_pc$importance$KL_divergence_importance_Gaussian_median))),
  KL_smoothed_pc_Gaussian_median = sort(unlist(sapply(results, function(x) x$pc$importance$KL_divergence_smoothed_importance_Gaussian_median))),
  KL_smoothed_not_pc_Gaussian_median = sort(unlist(sapply(results, function(x) x$not_pc$importance$KL_divergence_smoothed_importance_Gaussian_median)))
)

ggplot(rbind(
  data.frame(model = "PC", generator = "Gaussian_median", IS = "unsmoothed", KL = KL$KL_unsmoothed_pc_Gaussian_median),
  data.frame(model = "PC", generator = "Gaussian_median", IS = "smoothed", KL = KL$KL_smoothed_pc_Gaussian_median),
  data.frame(model = "NOT_PC", generator = "Gaussian_median", IS = "unsmoothed", KL = KL$KL_unsmoothed_not_pc_Gaussian_median),
  data.frame(model = "NOT_PC", generator = "Gaussian_median", IS = "smoothed", KL = KL$KL_smoothed_not_pc_Gaussian_median)
)) +
  stat_ecdf(aes(KL, col = model, linetype = IS))
# We also store the mean of the KL divergences
KL_Gaussian_median_mean <- data.frame(
  KL_unsmoothed_pc_Gaussian_median = mean(KL$KL_unsmoothed_pc_Gaussian_median),
  KL_unsmoothed_not_pc_Gaussian_median = mean(KL$KL_unsmoothed_not_pc_Gaussian_median),
  KL_smoothed_pc_Gaussian_median = mean(KL$KL_smoothed_pc_Gaussian_median),
  KL_smoothed_not_pc_Gaussian_median = mean(KL$KL_smoothed_not_pc_Gaussian_median)
)
print(KL_Gaussian_median_mean)

log_priors
all_probabilities <- data.frame()

for (prior_type in prior_types) {
  for (approximation_type in approximation_types) {
    probabilities <- sapply(results, function(x) x[[prior_type]][["importance"]][[paste0("probabilities_", approximation_type)]])
    for (i in seq_along(parameter_names)) {
      # # If any element in probabilities[i, ] is numeric(0), replace it with 0
      # probabilities[i, ] <- sapply(probabilities[i, ], function(x) if (length(x) == 0) 0 else x)
      # Create a data frame for the current probabilities and add it to all_probabilities
      df <- data.frame(prob = unlist(probabilities[i, ]), parameter = parameter_names[[i]], prior = prior_type, approximation = approximation_type)
      all_probabilities <- rbind(all_probabilities, df)
    }
  }
}





# Plot the ECDF for all probabilities
ggplot(all_probabilities) +
  stat_ecdf(aes(prob, col = prior, linetype = approximation)) +
  geom_abline(slope = 1, intercept = 0, color = "red") + # Add this line
  facet_wrap(~parameter)

# Now we calculate the Kolmogorov-Smirnov statistic for each parameter and represent it in a point plot
KS_results <- data.frame()

for (i in seq_along(parameter_names)) {
  # Get the probabililog_priorse current parameter
  probabilities <- all_probabilities[all_probabilities$parameter == parameter_names[[i]], ]
  # Calculate the KS statistic for each prior and approximation type
  for (prior_type in prior_types) {
    for (approximation_type in approximation_types) {
      KS_result <- ks.test(probabilities[probabilities$prior == prior_type & probabilities$approximation == approximation_type, ]$prob, "punif")
      # Add the KS statistic and p-value for the current parameter to the data frame
      KS_results <- rbind(KS_results, data.frame(parameter = parameter_names[[i]], prior = prior_type, approximation = approximation_type, statistic = KS_result$statistic, p_value = KS_result$p.value))
    }
  }
}

ggplot(KS_results) +
  geom_point(aes(x = parameter, y = statistic, color = prior, shape = approximation)) +
  facet_wrap(~parameter)

ggplot(KS_results) +
  geom_point(aes(x = parameter, y = p_value, color = prior, shape = approximation)) +
  facet_wrap(~parameter)




# # Average variance of log unnormalized weights
# var_weights <- sapply(results, function(x) var(normalize_log_weights(x$pc$importance_pc$log_unnormalized_weights)))
# hist(var_weights, main = "Histogram of variances for PC", xlab = "Variance")
#
# var_weights_smoothed <- sapply(results, function(x) var(normalize_log_weights(x$pc$importance_pc$log_unnormalized_weights_smoothed)))
# hist(var_weights_smoothed, main = "Histogram of variances for PC smoothed", xlab = "Variance")


# We get the vector of k diagnostics and show the log weights are similar
par(mfrow = c(1, 1))
k_diagnostics <- sapply(results, function(x) x$pc$importance$k_diagnostic)
hist(k_diagnostics, main = "Histogram of k diagnostics", xlab = "k diagnostic")
par(mfrow = c(1, 2))
j <- 1
hist(results[[j]]$pc$importance$log_unnormalized_weights_smoothed, main = "Log weights", xlab = "Log weight")
hist(results[[j]]$pc$importance$log_unnormalized_weights, main = "Log weights", xlab = "Log weight")
