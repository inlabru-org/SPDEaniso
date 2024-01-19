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

number_of_loops <- 1 # number of iterations
maxit_MAP <- 600
number_of_weights <- 100
confidence_level <- 0.001
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

      # Confidence intervals
      confidence_intervals_Laplace_pc <- importance_pc$confidence_intervals_Laplace
      confidences_intervals_importance_pc <- importance_pc$confidence_intervals_importance
      confidence_intervals_importance_smoothed_pc <- importance_pc$confidence_intervals_importance_smoothed
      confidence_intervals_Laplace_not_pc <- importance_not_pc$confidence_intervals_Laplace
      confidences_intervals_importance_not_pc <- importance_not_pc$confidence_intervals_importance
      confidence_intervals_importance_smoothed_not_pc <- importance_not_pc$confidence_intervals_importance_smoothed
      parameter_within_confidence_intervals_Laplace_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_Laplace_pc)
      parameter_within_confidence_intervals_importance_pc <- parameter_within_confidence_intervals(true_params, confidences_intervals_importance_pc)
      parameter_within_confidence_intervals_importance_smoothed_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_importance_smoothed_pc)
      parameter_within_confidence_intervals_Laplace_not_pc <- parameter_within_confidence_intervals(true_params, confidence_intervals_Laplace_not_pc)
      parameter_within_confidence_intervals_importance_not_pc <- parameter_within_confidence_intervals(true_params, confidences_intervals_importance_not_pc)
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
        confidences_intervals_importance = confidences_intervals_importance_pc,
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
        covariance_estimate = Q_Laplace_not_pc, # a 5x5 matrix
        std_dev_estimates_Laplace = sqrt(diag(Q_Laplace_pc)), # a 5d vector
        true_parameter_within_c_interval_Laplace = parameter_within_confidence_intervals_Laplace_not_pc, # a 5d vector
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




# This gives the 5d mean vector of the error in each parameter
# distance_vectors <- sapply(results, function(x) x$pc$distance_vector)
# distance_vectors
#
#
# intervals <- sapply(results, function(x) x$pc$true_parameter_within_c_interval_Laplace)
# intervals

# results <- future_lapply(1:m, function(i) {
#   true_params <- sim_theta_pc_quantile(alpha = alpha, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = 1)
#   # Extract true parameters
#   log_kappa <- true_params$log_kappa
#   kappa <- exp(log_kappa)
#   v <- true_params$v
#   log_sigma_u <- true_params$log_sigma_u
#   log_sigma_epsilon <- true_params$log_sigma_epsilon
#   aniso <- list(rep(kappa, n), matrix(v, n, 2))
#
#   # Sample from noisy data
#   x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso)
#   A <- Matrix::Diagonal(n, 1)
#   y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)
#
#   # Calculating the MAP under each prior knowing simulated data
#   map_pc <- MAP_prior(
#     log_prior = log_pc_prior, mesh = mesh,
#     y = y, A = A, m_u = m_u, max_iterations = maxit
#   )
#
#   map_not_pc <- MAP_prior(
#     log_prior = log_not_pc_prior, mesh = mesh,
#     y = y, A = A, m_u = m_u, max_iterations = maxit
#   )
#
#   # PC results
#   pc_results <- list(
#     MAP_estimate = map_pc$par,
#     convergence = map_pc$convergence,
#     distance_vector = abs(map_pc$par - unlist(true_params)),
#     covariance_estimate = solve(-map_pc$hessian),
#     std_dev_estimates_Laplace = sqrt(diag(solve(-map_pc$hessian)))
#   )
#
#   # Not-PC results
#   not_pc_results <- list(
#     MAP_estimate = map_not_pc$par,
#     convergence = map_not_pc$convergence,
#     distance_vector = abs(map_not_pc$par - unlist(true_params)),
#     covariance_estimate = solve(-map_not_pc$hessian),
#     std_dev_estimates_Laplace = sqrt(diag(solve(-map_not_pc$hessian)))
#   )
#
#   # Return the results for this iteration
#   list(
#     true_params = true_params,
#     pc = pc_results,
#     not_pc = not_pc_results
#   )
# })
