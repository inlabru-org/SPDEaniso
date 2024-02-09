library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(future)
library(future.apply)
library(dplyr)
library(loo)
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

L <- 10
width_uniform <- Inf
log_uniform_prior <- log_prior_uniform(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, width_support_factor = width_uniform)
shape <- 1.1
width_beta <- 20
# width_beta <- 2 # Changing to not get extreme values in simulations
log_beta_prior <- log_prior_beta(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, shape = shape, width_support_factor = width_beta)
log_priors <- list(
  pc = log_pc_prior,
  not_pc = log_not_pc_prior,
  uniform = log_uniform_prior,
  beta = log_beta_prior
)
prior_types <- setNames(as.list(names(log_priors)), names(log_priors))
approximation_types <- list("Gaussian_median", "importance", "importance_smoothed")
approximation_types <- setNames(approximation_types, approximation_types)
# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 3))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)
n_observations <- 15
observations <- L * matrix(runif(n_observations * 2), ncol = 2)
A <- fm_basis(mesh, loc = observations)

number_of_loops <- 200 # number of iterations
maxit_MAP <- 600
number_of_weights <- 5000
credible_level <- 0.05
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
      # #Simulate parameters from beta
      # true_params <- sim_theta_beta(
      #   sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
      #   a0 = a0, rho0 = rho0, L = L, shape = shape,
      #   width_support_factor = width_beta)
      # Simulate parameters from not_pc
      # true_params <- sim_not_pc(
      #   alpha = alpha, sigma_u0 = sigma_u0,
      #   sigma_epsilon0 = sigma_epsilon0,
      #   a0 = a0, rho0 = rho0
      # )

      log_kappa <- true_params$log_kappa
      kappa <- exp(log_kappa)
      v <- true_params$v
      log_sigma_u <- true_params$log_sigma_u
      log_sigma_epsilon <- true_params$log_sigma_epsilon
      aniso <- list(rep(kappa, n), matrix(v, n, 2))

      # Sample from noisy data
      x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
      y <- A %*% x + rnorm(n_observations, 0, exp(log_sigma_epsilon))

      # delta <- rnorm(5, 0, 1) # Used to randomize starting point of MAP
      delta <- 0
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
          log_posterior_density = log_posterior,
          mu_Gaussian_median = mu_Gaussian_median, Q_Gaussian_median = Q_Gaussian_median,
          n_weights = number_of_weights, q = credible_level, true_params = unlist(true_params)
        )
      })

      # CIs
      credible_intervals <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
          importances[[prior_type]][[paste0("credible_intervals_", approximation_type)]]
        })
      })

      true_parameter_is_within_CI <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
          parameter_within_credible_intervals(true_params, credible_intervals[[prior_type]][[approximation_type]])
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
          credible_intervals = lapply(approximation_types, function(approximation_type) {
            credible_intervals[[prior_type]][[approximation_type]]
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
# Results obtained simulating parameters from PC priors and using a mesh size of 1, 15 observations, 200 iterations, 5000 weights, a credible level of 0.05 a width of uniform =inf and for beta a multiplier of 20.
# saveRDS(results, "results_pc_1_15_200_5000_005_wu_inf_wb_20.rds")
# results <- readRDS("Simulation_results/results_pc_1_15_200_5000_005_wu_inf_wb_20.rds")
parameter_names <- rownames(results[[1]]$pc$credible_intervals$Gaussian_median)
# Plots ecdf of distances to MAP using ggplot
plot_distances_to_MAP <- function(results, prior_types) {
  all_distances <- data.frame()
  for (prior_type in prior_types) {
    distances_to_MAP <- lapply(results, function(x) x[[prior_type]]$distance_vector)
    distances_to_MAP <- do.call(rbind, distances_to_MAP)
    distances_to_MAP <- as.data.frame(distances_to_MAP)
    distances_to_MAP$iteration <- seq_len(nrow(distances_to_MAP))
    distances_to_MAP <- reshape2::melt(distances_to_MAP, id.vars = "iteration")
    distances_to_MAP$prior_type <- prior_type
    all_distances <- rbind(all_distances, distances_to_MAP)
  }

  ggplot(all_distances) +
    stat_ecdf(aes(value, color = prior_type)) +
    facet_wrap(~variable)
}

plot_distances_to_MAP(results, prior_types)

# Mean distances and standard deviation estimates
mean_distances <- lapply(prior_types, function(prior_type) {
  mean_distances <- sapply(1:5, function(i) {
    all_distances <- sapply(seq_along(results), function(j) {
      results[[j]][[prior_type]]$distance_vector[i]
    })
    mean(all_distances)
  })

  std_dev_estimates_Gaussian_median <- do.call(rbind, lapply(results, function(x) x[[prior_type]]$std_dev_estimates_Gaussian_median))
  mean_std_dev <- colMeans(std_dev_estimates_Gaussian_median)

  names(mean_distances) <- parameter_names
  print(paste("Mean distances for", prior_type))
  print(mean_distances)
  print(paste("Mean standard deviations for", prior_type))
  print(mean_std_dev)
})


# THere we test how to get the dataframe of credible interval lengths using lapplyand for each approximation type

lengths_df <- lapply(prior_types, function(prior_type) {
  lapply(approximation_types, function(approximation_type) {
    lengths <- lapply(parameter_names, function(parameter_name) {
      all_lengths <- sapply(seq_along(results), function(j) {
        length <- diff(results[[j]][[prior_type]]$credible_intervals[[approximation_type]][parameter_name, ])
        length
      })
      all_lengths
    })
    lengths <- do.call(cbind, lengths)
    colnames(lengths) <- parameter_names
    lengths
  })
})

mean_lengths_df <- lapply(lengths_df, function(prior_type) {
  lapply(prior_type, function(approximation_type) {
    colMeans(approximation_type)
  })
})
print(mean_lengths_df)

# We show a cdf of lengths_df using ggplot
plot_CI_lengths <- function(lengths_df, prior_types, approximation_types) {
  all_lengths <- data.frame()

  for (prior_type in prior_types) {
    for (approximation_type in approximation_types) {
      lengths <- lengths_df[[prior_type]][[approximation_type]]
      lengths <- as.data.frame(lengths)
      lengths$iteration <- seq_len(nrow(lengths))
      lengths <- reshape2::melt(lengths, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
      lengths$prior_type <- prior_type
      lengths$approximation_type <- approximation_type
      all_lengths <- rbind(all_lengths, lengths)
    }
  }

  ggplot(all_lengths) +
    stat_ecdf(aes(value, color = prior_type, linetype = approximation_type)) +
    facet_wrap(~variable)
}

plot_CI_lengths(lengths_df, prior_types, approximation_types)



# Percentage of times the true parameter is within the credible interval

within_ci <- lapply(prior_types, function(prior_type) {
  lapply(approximation_types, function(approximation_type) {
    rowMeans(sapply(results, function(x) x[[prior_type]][["true_parameter_within_c_interval"]][[approximation_type]]))
  })
})
# Convert the list to a data frame
within_ci <- do.call(rbind, lapply(names(within_ci), function(prior_type) {
  do.call(rbind, lapply(names(within_ci[[prior_type]]), function(approximation_type) {
    data.frame(
      prior_type = prior_type,
      approximation_type = approximation_type,
      parameter = factor(names(within_ci[[prior_type]][[approximation_type]]), levels = c("log_kappa", "v1", "v2", "log_sigma_u", "log_sigma_epsilon")),
      value = c(within_ci[[prior_type]][[approximation_type]])
    )
  }))
}))

# Plot the points using ggplot
ggplot(within_ci, aes(x = approximation_type, y = value, color = prior_type)) +
  geom_point() +
  geom_text(aes(label = round(value, 2)), vjust = -0.5) +
  facet_wrap(~parameter) +
  labs(x = "Approximation Type", y = "Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


KL_approx_types <- list(importance = "importance", smoothed_importance = "smoothed_importance")
KL <- lapply(prior_types, function(prior_type) {
  lapply(KL_approx_types, function(approximation_type) {
    sapply(results, function(x) {
      kl_values <- x[[prior_type]]$importance[[paste0("KL_divergence_", approximation_type, "_Gaussian_median")]]
      kl_values <- replace(kl_values, is.na(kl_values), Inf) # Replace NA values with Inf
      kl_values
    })
  })
})

# ggplot of ecdf of KL divergences
plot_KL_divergences <- function(KL, prior_types, approximation_types) {
  all_KL <- data.frame()
  for (prior_type in prior_types) {
    for (approximation_type in KL_approx_types) {
      KL_divergences <- KL[[prior_type]][[approximation_type]]
      KL_divergences <- as.data.frame(KL_divergences)
      KL_divergences$iteration <- seq_len(nrow(KL_divergences))
      KL_divergences <- reshape2::melt(KL_divergences, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
      KL_divergences$prior_type <- prior_type
      KL_divergences$approximation_type <- approximation_type
      all_KL <- rbind(all_KL, KL_divergences)
    }
  }

  ggplot(all_KL) +
    stat_ecdf(aes(value, color = prior_type, linetype = approximation_type)) +
    facet_wrap(~variable)
}

plot_KL_divergences(KL, prior_types, approximation_types)

# We calculate the mean of the KL divergences
KL_Gaussian_median_mean <- lapply(prior_types, function(prior_type) {
  lapply(KL_approx_types, function(approximation_type) {
    mean(KL[[prior_type]][[approximation_type]])
  })
})
print(KL_Gaussian_median_mean)


# Probabilities that marginal posterior is smaller than true parameter
all_probabilities <- data.frame()

for (prior_type in prior_types) {
  for (approximation_type in approximation_types) {
    probabilities <- sapply(results, function(x) x[[prior_type]][["importance"]][[paste0("probabilities_", approximation_type)]])
    for (i in seq_along(parameter_names)) {
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
  # Get the probabilities for the current parameter
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

# ks.test(x<-runif(5000),"punif")

# Now we do the CDF of the complexity of the model
complexity <- lapply(prior_types, function(prior_type) {
  lapply(approximation_types[2:3], function(approximation_type) {
    sapply(results, function(x) {
      complexity <- x[[prior_type]]$importance[[paste0("complexity_", approximation_type)]]
      complexity
    })
  })
})

# ggplot of ecdf of complexity
plot_complexity <- function(complexity, prior_types, approximation_types) {
  all_complexity <- data.frame()
  for (prior_type in prior_types) {
    for (approximation_type in approximation_types) {
      complexity_values <- complexity[[prior_type]][[approximation_type]]
      complexity_values <- as.data.frame(complexity_values)
      complexity_values$iteration <- seq_len(nrow(complexity_values))
      complexity_values <- reshape2::melt(complexity_values, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
      complexity_values$prior_type <- prior_type
      complexity_values$approximation_type <- approximation_type
      all_complexity <- rbind(all_complexity, complexity_values)
    }
  }

  ggplot(all_complexity) +
    stat_ecdf(aes(value, color = prior_type, linetype = approximation_type)) +
    facet_wrap(~variable)
}

plot_complexity(complexity, prior_types[2:3], approximation_types[2:3])

# We calculate the mean of the complexity
complexity_mean <- lapply(prior_types, function(prior_type) {
  lapply(approximation_types[2:3], function(approximation_type) {
    mean(complexity[[prior_type]][[approximation_type]])
  })
})
print(complexity_mean)


# We get the vector of k diagnostics and show the log weights are similar
par(mfrow = c(1, 1))
k_diagnostics <- sapply(results, function(x) x$pc$importance$k_diagnostic)
hist(k_diagnostics, main = "Histogram of k diagnostics", xlab = "k diagnostic")
par(mfrow = c(1, 2))
j <- 1
hist(results[[j]]$pc$importance$log_unnormalized_weights_smoothed, main = "Log weights", xlab = "Log weight")
hist(results[[j]]$pc$importance$log_unnormalized_weights, main = "Log weights", xlab = "Log weight")
