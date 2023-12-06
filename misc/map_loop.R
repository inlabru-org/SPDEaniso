library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(future)
library(future.apply)

# Set the parallel plan to use all local cores (currently not used as future package doesn't recognize functions in prior.R)
plan(multisession)
# Defining the random seed
set.seed(123)

# Defines the upper bounds for the quantiles
a0 <- 2 # Controls the size of v
rho0 <- 0.1 # Controls the size of kappa
sigma0 <- 1.5 # Controls the size of v in non PC priors
sigmau0 <- 2 # controls standard deviation of field
sigmaepsilon0 <- 0.1 # control standard deviation of noise

# Defines the quantile
alpha <- 0.01
# Calculates the hyperparameters of PC by hand
lambda1 <- lambda1_quantile(a0) # Controls the size of v
lambda <- lambda_quantile(rho0 = rho0, lambda1 = lambda1) # Controls the size of kappa
lambda_u <- lambda_variance_quantile(sigma0 = sigmau0) # Controls the variance of the field
lambda_epsilon <- lambda_variance_quantile(sigma0 = sigmaepsilon0) # Controls the variance of the noise
hyper_pc <- list(
  lambda = lambda, lambda1 = lambda1, lambda_u = lambda_u,
  lambda_epsilon = lambda_epsilon
)
true_params <- sim_theta_pc_quantile(
  alpha = alpha, sigmau0 = sigmau0,
  sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1
)
print(hyper_pc)
print(true_params)
# Calculates the hyperparameters of not PC by hand
sigma2_quantile_v(a0)
lambda_quantile_kappa(rho0)


# Setting mean of the field
m_u <- 0
# Calculates the log prior density function of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(
  sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)

log_not_pc_prior <- log_gaussian_prior_quantile(
  sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)
# Testing value of priors
true_params <- sim_theta_pc_quantile(alpha = alpha, sigmau0 = sigmau0,
                                     sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1)
log_kappa <- true_params$log_kappa
v <- true_params$v
log_sigma_u <- true_params$log_sigma_u
log_sigma_epsilon <- true_params$log_sigma_epsilon
# Testing that log_pc_prior is the same when using parameters and quantiles
log_pc_prior(
  log_kappa = log_kappa, v = v,
  log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
)
log_pc_prior_theta(
  lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon,
  lambda_u = lambda_u, log_kappa = log_kappa, v = v,
  log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
)
# Comparing against value with arbitraty hyperparameters. With arbitrary parameters should be smaller in general.
log_pc_prior_theta(
  lambda = 1, lambda1 = 1, lambda_epsilon = 1, lambda_u = 1,
  log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u,
  log_sigma_epsilon = log_sigma_epsilon
)
# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(10, 0.01), c(10, 10), c(0, 10), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(0.25, 0.25))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)


m <- 1 # number of iterations
maxit <- 10
m <- 10 # number of iterations
maxit <- 1200
results <- vector("list", m) # Pre-allocates a list for m iterations
start_time <- Sys.time()
for (i in 1:m) {
  tryCatch(
    {
      # Simulate parameters from PC prior
      true_params <- sim_theta_pc_quantile(
        alpha = alpha, sigmau0 = sigmau0,
        sigmaepsilon0 = sigmaepsilon0,
        a0 = a0, rho0 = rho0, m = 1
      )
      # true_params <- list(
      #   log_kappa = -0.5,
      #   v = c(0, 0),
      #   log_sigma_u = 0,
      #   log_sigma_epsilon = -3)
      # Extract true parameters
      log_kappa <- true_params$log_kappa
      kappa <- exp(log_kappa)
      v <- true_params$v
      log_sigma_u <- true_params$log_sigma_u
      log_sigma_epsilon <- true_params$log_sigma_epsilon
      aniso <- list(rep(kappa, n), matrix(v, n, 2))

      # Sample from noisyy data
      x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso)
      A <- Matrix::Diagonal(n, 1)
      y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)

      # Calculating the MAP under each prior knowing simulated data
      map_pc <- MAP_prior(
        logprior = log_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, maxiterations = maxit,
        theta0 = unlist(true_params)
      )

      map_not_pc <- MAP_prior(
        logprior = log_not_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, maxiterations = maxit,
        theta0 = unlist(true_params)
      )

      # PC results
      pc_results <- list(
        MAP_estimate = map_pc$par, # a 5d vector
        convergence = map_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_pc$par - unlist(true_params)), # a 5d vector
        covariance_estimate = solve(-map_pc$hessian), # a 5x5 matrix
        std_dev_estimates = sqrt(diag(solve(-map_pc$hessian))), # a 5d vector
        within_cinterval = c(1, 1, 1, 1, 1) * (abs(map_pc$par - unlist(true_params))
        < 1.96 * sqrt(diag(solve(-map_pc$hessian)))) # a 5d vector
        # prob_MAP_greater = pc_prob_map                                      # a 5d vector
      )
      # Not-PC results
      not_pc_results <- list(
        MAP_estimate = map_not_pc$par, # a 5d vector
        convergence = map_not_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_not_pc$par - unlist(true_params)), # a 5d vector
        covariance_estimate = solve(-map_not_pc$hessian), # a 5x5 matrix
        std_dev_estimates = sqrt(diag(solve(-map_pc$hessian))), # a 5d vector
        within_cinterval = c(1, 1, 1, 1, 1) * (abs(map_not_pc$par - unlist(true_params))
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
end_time <- Sys.time()

# Calculate the duration
duration <- end_time - start_time
print(paste("Total time taken: ", duration))
sapply(results, function(x) x$pc$std_dev_estimates)
# This gives the 5d mean vector of the error in each parameter
# distance_vectors <- sapply(results, function(x) x$pc$distance_vector)
# distance_vectors
#
#
# intervals <- sapply(results, function(x) x$pc$within_cinterval)
# intervals

# results <- future_lapply(1:m, function(i) {
#   true_params <- sim_theta_pc_quantile(alpha = alpha, sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1)
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
#     logprior = log_pc_prior, mesh = mesh,
#     y = y, A = A, m_u = m_u, maxiterations = maxit
#   )
#
#   map_not_pc <- MAP_prior(
#     logprior = log_not_pc_prior, mesh = mesh,
#     y = y, A = A, m_u = m_u, maxiterations = maxit
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





results[[10]][[1]]
true_params
log_posterior()
