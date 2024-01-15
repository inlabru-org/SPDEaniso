# Tests that functions to calculate hyperparameters given quantiles work
# Tests that parameters samples with PC-priors verify quantiles
# Tests that PC_prior defined through quantiles is the same as PC_prior defined through hyperparameters
# Tests that MAP_prior(log_pc_prior) works and gives a larger posterior density than the true parameters
# Tests that MAP_prior with PC priors returns the same thing as MAP using PC priors defined through hyperparameters
library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
document()
# Defining the random seed
set.seed(123)

# Defines the upper bounds for the quantiles
rho0 <- 1 # Controls the size of kappa in PC and non PC priors
a0 <- 2 # Controls the size of v in PC and non PC priors
sigmau0 <- 2 # controls standard deviation of field
sigmaepsilon0 <- 0.1 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
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

# Simulates the "true parameters" from the pc_prior.
true_params <- sim_theta_pc_quantile(
  alpha = alpha, sigmau0 = sigmau0,
  sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1
)
# Prints the hyperparameters and the true parameters
print(hyper_pc)
print(true_params)

# Extracting simulated parameters
log_kappa <- true_params$log_kappa
kappa <- exp(log_kappa)
v <- true_params$v
log_sigma_u <- true_params$log_sigma_u
log_sigma_epsilon <- true_params$log_sigma_epsilon

print("Testing that the true parameters are correct. They should verify the quantiles with high probability.")
print("The correlation l sqrt(8)/exp(log_kappa) is greater than rho0")
print(sqrt(8) / exp(log_kappa) > rho0)
print("The anisotropy ratio exp(||v||) is smaller than a0")
print(exp(sqrt(v[1]^2 + v[2]^2)) < a0)
print("The standard deviation of the field exp(log_sigma_u) is smaller than sigmau0")
print(exp(log_sigma_u) < sigmau0)
print("The standard deviation of the noise exp(log_sigma_epsilon) is smaller than sigmaepsilon0")
print(exp(log_sigma_epsilon) < sigmaepsilon0)


# Calculates the hyperparameters of not PC by hand
sigma2_quantile_v(a0)
lambda_quantile_kappa(rho0)

# Setting the mean of the field
m_u <- 0

# Defines the log prior density function of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(
  sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)

log_not_pc_prior <- log_gaussian_prior_quantile(
  sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)

### Testing value of priors ###
print("Testing that log_pc_prior is the same when using quantiles and hyper_parameters calculated by hand.")
abs(log_pc_prior(
  log_kappa = log_kappa, v = v,
  log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
) - log_pc_prior_theta(
  lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon,
  lambda_u = lambda_u, log_kappa = log_kappa, v = v,
  log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
)) < 1e-10


print("Comparing against value with arbitrary hyperparameters. With arbitrary hyperparameters
      the prior density should be smaller in general.")
log_pc_prior_theta(
  lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon,
  lambda_u = lambda_u, log_kappa = log_kappa, v = v,
  log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
) > log_pc_prior_theta(
  lambda = 1, lambda1 = 1, lambda_epsilon = 1, lambda_u = 1,
  log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u,
  log_sigma_epsilon = log_sigma_epsilon
)


# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(10, 0.01), c(10, 10), c(0, 10), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 1))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

# Sample from noisy data
aniso <- list(rep(kappa, n), matrix(v, n, 2))
x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso)
A <- Matrix::Diagonal(n, 1)
y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)

# Testing the MAP_prior function for log_prior = log_pc_prior
maxit <- 600
tryCatch({
  # Calculating the MAP under each prior knowing simulated data
  map_pc <- MAP_prior(
    logprior = log_pc_prior, mesh = mesh,
    y = y, A = A, m_u = m_u, maxiterations = maxit,
    theta0 = unlist(true_params)
    # ,log_sigma_epsilon = log_sigma_epsilon
  )

  error <- function(e) {}
})
#### TESTING MAP_prior ####

# Defining the log posterior density as a function of the parameters using the log_pc_prior
log_posterior_pc <- function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
  log_posterior_prior(
    logprior = log_pc_prior,
    mesh = mesh, log_kappa = log_kappa, v = v,
    log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
    y = y, A = A, m_u = m_u
  )
}

print("Checking if map is larger than posterior value at true parameters. ")
map_pc$value > log_posterior_pc(
  log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u,
  log_sigma_epsilon = log_sigma_epsilon
)

print("Checking if value for map is the same as posterior value at map parameters. ")
map_pc$value == log_posterior_pc(
  log_kappa = map_pc$par[1], v = map_pc$par[2:3],
  log_sigma_u = map_pc$par[4], log_sigma_epsilon = map_pc$par[5]
)

print("Checking if MAP_prior with PC priors returns the same thing as MAP using PC priors defined through hyperparameters")

maphyper <- MAP_pc(mesh, lambda, lambda1, lambda_epsilon, lambda_u,
  y, A, m_u,
  maxiterations = maxit, theta0 = unlist(true_params)
)
maphyper$value == map_pc$value & maphyper$par == map_pc$par

print("Checking if Hessian of MAP_prior is invertible and calculating marginal standard deviation")
hessian <- map_pc$hessian
hessian_inv <- solve(hessian)
sqrt(-diag(hessian_inv))

########## UNNORMALIZED LAPLACE APPROXIMATION TO THE POSTERIOR############
log_laplace <- function(theta) {
    logGdensity(
        x = theta, mu = map_pc$par, Q = -hessian
    )+map_pc$value-logGdensity(
        x = map_pc$par, mu = map_pc$par, Q = -hessian
    )
}

print("Plotting the posterior distribution of the parameters")

plotter <- function(map) {
  # Defines the log_pc_posterior density as a function of log_kappa keeping the other parameters fixed at the MAP
  log_posterior_pc_log_kappa <- function(log_kappa) {
    log_posterior_pc(
      log_kappa = log_kappa, v = map$par[2:3],
      log_sigma_u = map$par[4], log_sigma_epsilon = map$par[5]
    )
  }
  # Defines the log_pc_posterior density as a function of v1 keeping the other parameters fixed at the MAP
  log_posterior_pc_v1 <- function(v1) {
    log_posterior_pc(
      log_kappa = map$par[1], v = c(v1, map$par[3]),
      log_sigma_u = map$par[4], log_sigma_epsilon = map$par[5]
    )
  }
  # Defines the log_pc_posterior density as a function of v2 keeping the other parameters fixed at the MAP
  log_posterior_pc_v2 <- function(v2) {
    log_posterior_pc(
      log_kappa = map$par[1], v = c(map$par[2], v2),
      log_sigma_u = map$par[4], log_sigma_epsilon = map$par[5]
    )
  }
  # Defines the log_pc_posterior density as a function of log_sigma_u keeping the other parameters fixed at the MAP
  log_posterior_pc_log_sigma_u <- function(log_sigma_u) {
    log_posterior_pc(
      log_kappa = map$par[1], v = map$par[2:3],
      log_sigma_u = log_sigma_u, log_sigma_epsilon = map$par[5]
    )
  }
  # Defines the log_pc_posterior density as a function of log_sigma_epsilon keeping the other parameters fixed at the MAP
  log_posterior_pc_log_sigma_epsilon <- function(log_sigma_epsilon) {
    log_posterior_pc(
      log_kappa = map$par[1], v = map$par[2:3],
      log_sigma_u = map$par[4], log_sigma_epsilon = log_sigma_epsilon
    )
  }
log_laplace_log_kappa <- function(log_kappa) {
    log_laplace(
      theta = c(log_kappa, map$par[2:5])
    )
  }

  log_laplace_v1 <- function(v1) {
    log_laplace(
      theta = c(map$par[1], v1, map$par[3:5])
    )
  }

  log_laplace_v2 <- function(v2) {
    log_laplace(
      theta = c(map$par[1:2], v2, map$par[4:5])
    )
  }

  log_laplace_log_sigma_u <- function(log_sigma_u) {
    log_laplace(
      theta = c(map$par[1:3], log_sigma_u, map$par[5])
    )
  }

  log_laplace_log_sigma_epsilon <- function(log_sigma_epsilon) {
    log_laplace(
      theta = c(map$par[1:4], log_sigma_epsilon)
    )
  }
    l <- 2 # l of the partition
  # Defines the partitions for plotting
  n_points <- 51 # Number of points in the partition centered at MAP_prior value of kappa

  partition_log_kappa <- seq(map$par[1] - l, map$par[1] + l, length.out = n_points)
  partition_v1 <- seq(map$par[2] - l, map$par[2] + l, length.out = n_points)
  partition_v2 <- seq(map$par[3] - l, map$par[3] + l, length.out = n_points)
  partition_log_sigma_u <- seq(map$par[4] - l, map$par[4] + l, length.out = n_points)
  partition_log_sigma_epsilon <- seq(map$par[5] - l, map$par[5] + l, length.out = n_points)

  # Apply the function to each point in the partition
  posterior_values_log_kappa <- sapply(partition_log_kappa, log_posterior_pc_log_kappa)
  posterior_values_v1 <- sapply(partition_v1, log_posterior_pc_v1)
  posterior_values_v2 <- sapply(partition_v2, log_posterior_pc_v2)
  posterior_values_log_sigma_u <- sapply(partition_log_sigma_u, log_posterior_pc_log_sigma_u)
  posterior_values_log_sigma_epsilon <- sapply(partition_log_sigma_epsilon, log_posterior_pc_log_sigma_epsilon)
  laplace_values_log_kappa <- sapply(partition_log_kappa, log_laplace_log_kappa)
  laplace_values_v1 <- sapply(partition_v1, log_laplace_v1)
  laplace_values_v2 <- sapply(partition_v2, log_laplace_v2)
  laplace_values_log_sigma_u <- sapply(partition_log_sigma_u, log_laplace_log_sigma_u)
  laplace_values_log_sigma_epsilon <- sapply(partition_log_sigma_epsilon, log_laplace_log_sigma_epsilon)


  # Plot the results with a vertical line at the MAP_prior value of kappa
  plot(partition_log_kappa, posterior_values_log_kappa, type = "l", xlab = "log_kappa", ylab = "log density")
  points(partition_log_kappa, laplace_values_log_kappa, type = "l", col = "blue")
  abline(v = map$par[1], col = "red")
  legend("bottomleft", legend = c("posterior", "Laplace"), col = c("black", "blue"), pch=1)

  # Plot the results with a vertical line at the MAP_prior value of v1
  plot(partition_v1, posterior_values_v1, type = "l", xlab = "v1", ylab = "log density")
  points(partition_v1, laplace_values_v1, type = "l", col = "blue")
  abline(v = map$par[2], col = "red")
  legend("bottomleft", legend = c("posterior", "Laplace"), col = c("black", "blue"), pch=1)


  # Plot the results with a vertical line at the MAP_prior value of v2
  plot(partition_v2, posterior_values_v2, type = "l", xlab = "v2", ylab = "log density")
  points(partition_v2, laplace_values_v2, type = "l", col = "blue")
  abline(v = map$par[3], col = "red")
  legend("bottomleft", legend = c("posterior", "Laplace"), col = c("black", "blue"), pch=1)


  # Plot the results with a vertical line at the MAP_prior value of log_sigma_u
  plot(partition_log_sigma_u, posterior_values_log_sigma_u, type = "l", xlab = "log_sigma_u", ylab = "log density")
  points(partition_log_sigma_u, laplace_values_log_sigma_u, type = "l", col = "blue")
  abline(v = map$par[4], col = "red")
  legend("bottomleft", legend = c("posterior", "Laplace"), col = c("black", "blue"), pch = 1)

  # Plot the results with a vertical line at the MAP_prior value of log_sigma_epsilon
  plot(partition_log_sigma_epsilon, posterior_values_log_sigma_epsilon, type = "l", xlab = "log_sigma_epsilon", ylab = "log density")
  points(partition_log_sigma_epsilon, laplace_values_log_sigma_epsilon, type = "l", col = "blue")
  abline(v = map$par[5], col = "red")
  legend("bottomleft", legend = c("posterior", "Laplace"), col = c("black", "blue"), pch = 1)
}

plotter(map = map_pc)

# Showing true parameters, MAP, and standard deviation of the parameters
unlist(true_params)
map_pc$par
sqrt(-diag(hessian_inv))
