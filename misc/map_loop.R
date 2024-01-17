library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(future)
library(future.apply)
document()

# Set the parallel plan to use all local cores (currently not used as future package doesn't recognize functions in prior.R)
plan(multisession)
# Defining the random seed
set.seed(123)

# Defines the upper bounds for the quantiles
a0 <- 2 # Controls the size of v
rho0 <- 0.1 # Controls the size of kappa
sigma0 <- 1.5 # Controls the size of v in non PC priors
sigma_u0 <- 2 # controls standard deviation of field
sigma_epsilon0 <- 0.1 # control standard deviation of noise

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
number_of_weights <- 1000
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




      # PC results
      pc_results <- list(
        MAP_estimate = map_pc$par, # a 5d vector
        MAP_value = map_pc$value,
        convergence = map_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_pc$par - unlist(true_params)), # a 5d vector
        covariance_estimate = solve(-map_pc$hessian), # a 5x5 matrix
        std_dev_estimates = sqrt(diag(solve(-map_pc$hessian))), # a 5d vector
        within_c_interval = c(1, 1, 1, 1, 1) * (abs(map_pc$par - unlist(true_params))
        < 1.96 * sqrt(diag(solve(-map_pc$hessian)))) # a 5d vector
        # prob_MAP_greater = pc_prob_map                                      # a 5d vector
      )
      # Not-PC results
      not_pc_results <- list(
        MAP_estimate = map_not_pc$par, # a 5d vector
        MAP_value = map_not_pc$value,
        convergence = map_not_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_not_pc$par - unlist(true_params)), # a 5d vector
        covariance_estimate = solve(-map_not_pc$hessian), # a 5x5 matrix
        std_dev_estimates = sqrt(diag(solve(-map_pc$hessian))), # a 5d vector
        within_c_interval = c(1, 1, 1, 1, 1) * (abs(map_not_pc$par - unlist(true_params))
        < 1.96 * sqrt(diag(solve(-map_not_pc$hessian)))) # a 5d vector
        # prob_MAP_greater = not_pc_prob_map                                      # a 5d vector
      )

      # Store results
      results[[i]] <- list(
        true_params = true_params,
        pc = pc_results,
        not_pc = not_pc_results
      )
    },
    error = function(e) {}
  )
}

# sapply(results, function(x) x$pc$std_dev_estimates)
results




# This gives the 5d mean vector of the error in each parameter
# distance_vectors <- sapply(results, function(x) x$pc$distance_vector)
# distance_vectors
#
#
# intervals <- sapply(results, function(x) x$pc$within_c_interval)
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
#     std_dev_estimates = sqrt(diag(solve(-map_pc$hessian)))
#   )
#
#   # Not-PC results
#   not_pc_results <- list(
#     MAP_estimate = map_not_pc$par,
#     convergence = map_not_pc$convergence,
#     distance_vector = abs(map_not_pc$par - unlist(true_params)),
#     covariance_estimate = solve(-map_not_pc$hessian),
#     std_dev_estimates = sqrt(diag(solve(-map_not_pc$hessian)))
#   )
#
#   # Return the results for this iteration
#   list(
#     true_params = true_params,
#     pc = pc_results,
#     not_pc = not_pc_results
#   )
# })
