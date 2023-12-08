#In this file, parameters are simulated from the PC prior,
#the posterior distribution given a noisy observation is calculated both for pc and non-pc priors
#and is maximized while keeping the noise parameter log_sigma_epsilon fixed.
library(devtools)
document()
library(SPDEaniso)
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

# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(10, 0.01), c(10, 10), c(0, 10), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 1))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

m <- 100 # number of iterations
maxit <- 600
results <- vector("list", m) # Pre-allocates a list for m iterations
for (i in 1:m) {
  tryCatch(
    {
      # Simulate parameters from PC prior
      true_params <- sim_theta_pc_quantile(
        alpha = alpha, sigmau0 = sigmau0,
        sigmaepsilon0 = sigmaepsilon0,
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
        logprior = log_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, maxiterations = maxit,
        theta0 = unlist(true_params), log_sigma_epsilon = log_sigma_epsilon
      )

      map_not_pc <- MAP_prior(
        logprior = log_not_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, maxiterations = maxit,
        theta0 = unlist(true_params), log_sigma_epsilon = log_sigma_epsilon
      )

      # PC results
      pc_results <- list(
        MAP_estimate = map_pc$par, # a 5d vector
        MAP_value <- map_pc$value,
        convergence = map_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_pc$par - unlist(true_params)[1:4]), # a 5d vector
        covariance_estimate = solve(-map_pc$hessian), # a 5x5 matrix
        std_dev_estimates = sqrt(diag(solve(-map_pc$hessian))), # a 5d vector
        within_cinterval = c(1, 1, 1, 1) * (abs(map_pc$par - unlist(true_params)[1:4])
                                               < 1.96 * sqrt(diag(solve(-map_pc$hessian)))), # a 5d vector
        # To calculate the probability that the MAP estimate is greater than the true parameter we estimate the
        # marginal posterior of the posterior using the marginal standard deviations of the posterior. So this probability is calculated
        # using the CDF of the normal distribution with mean equal to the MAP estimate and standard deviation equal to the
        # marginal standard deviation of the posterior
         prob_MAP_greater = pnorm(map_pc$par, mean = unlist(true_params)[1:4], sd = sqrt(diag(solve(-map_pc$hessian))))

      )
    #  Not-PC results
      not_pc_results <- list(
        MAP_estimate = map_not_pc$par, # a 5d vector
        MAP_value <- map_not_pc$value,
        convergence = map_not_pc$convergence, # convergence = 0 no convergence =1
        distance_vector = abs(map_not_pc$par - unlist(true_params)[1:4]), # a 5d vector
        covariance_estimate = solve(-map_not_pc$hessian), # a 5x5 matrix
        std_dev_estimates = sqrt(diag(solve(-map_pc$hessian))), # a 5d vector
        within_cinterval = c(1, 1, 1, 1) * (abs(map_not_pc$par - unlist(true_params)[1:4])
                                               < 1.96 * sqrt(diag(solve(-map_not_pc$hessian)))), # a 5d vector
        prob_MAP_greater = pnorm(map_not_pc$par, mean = unlist(true_params)[1:4], sd = sqrt(diag(solve(-map_not_pc$hessian))))
      )

      # Store results
      results[[i]] <- list(
        true_params = true_params,
        pc = pc_results
         , not_pc = not_pc_results
      )
    },
    error = function(e) {}
  )
}

# This calculates the proportion of times the MAP estimate is in the 95% confidence intervalfor each parameter ignoring
# results with any NANs
pc_within_cinterval <- lapply(results, function(x) x$pc$within_cinterval)
#This transforms the list into a matrix
pc_within_cinterval <- do.call(rbind, pc_within_cinterval)
#This calculates the mean of each column ignoring NA values
pc_within_cinterval <- colMeans(pc_within_cinterval, na.rm = TRUE)
pc_within_cinterval

# This calculates the proportion of times the MAP estimate is in the 95% confidence intervalfor each parameter ignoring
# results with any NANs
not_pc_within_cinterval <- lapply(results, function(x) x$not_pc$within_cinterval)
#This transforms the list into a matrix
not_pc_within_cinterval <- do.call(rbind, not_pc_within_cinterval)
#This calculates the mean of each column ignoring NA values
not_pc_within_cinterval <- colMeans(not_pc_within_cinterval, na.rm = TRUE)
not_pc_within_cinterval

#This plots a histogram of the probability that the MAP estimate is greater than the true parameter
# It should be uniform if the MAP estimate is unbiased
pc_prob_MAP_greater <- lapply(results, function(x) x$pc$prob_MAP_greater)
pc_prob_MAP_greater <- do.call(rbind, pc_prob_MAP_greater)
pc_prob_MAP_greater <- colMeans(pc_prob_MAP_greater, na.rm = TRUE)
hist(pc_prob_MAP_greater)

