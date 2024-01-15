# Implements Pareto smoothed importance sampling for the PC prior
library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(loo)
document()
# Defining the random seed
set.seed(123)

# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(10, 0.01), c(10, 10), c(0, 10), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 1))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

# ---------------------- SIMULATING DATA Y ----------------------

# Defines the upper bounds for the quantiles
rho0 <- 1 # Controls the size of kappa in PC and non PC priors
a0 <- 2 # Controls the size of v in PC and non PC priors
sigmau0 <- 2 # controls standard deviation of field
sigmaepsilon0 <- 0.1 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

# Simulates the "true parameters" from the pc_prior.
true_params <- sim_theta_pc_quantile(
    alpha = alpha, sigmau0 = sigmau0,
    sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1
)

log_kappa <- true_params$log_kappa
kappa <- exp(log_kappa)
v <- true_params$v
log_sigma_u <- true_params$log_sigma_u
log_sigma_epsilon <- true_params$log_sigma_epsilon


# Sample from noisy data
aniso <- list(rep(kappa, n), matrix(v, n, 2))
x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso)
A <- Matrix::Diagonal(n, 1)
y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)

# Defines the log prior density function of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(
    sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
    a0 = a0, rho0 = rho0, alpha = alpha
)
# Defining the log posterior density as a function of the parameters using the log_pc_prior
m_u <- 0
log_posterior_pc <- function(theta) {
    log_kappa <- theta[1]
    v <- theta[2:3]
    log_sigma_u <- theta[4]
    log_sigma_epsilon <- theta[5]
    log_posterior_prior(
        logprior = log_pc_prior,
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
        y = y, A = A, m_u = m_u
    )
}

########## LAPLACE APPROXIMATION TO THE POSTERIOR############
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

hessian <- map_pc$hessian
sigma <- -solve(hessian)

log_laplace <- function(theta) {
    logGdensity(
        x = theta, mu = map_pc$par, Q = -hessian
    )
}
# ---------------------- IMPORTANCE SAMPLING ----------------------
# Simulate from Laplace approximation
n_sim <- 100
theta_sim <- MASS::mvrnorm(n_sim, map_pc$par, sigma)
log_ratio_function <- function(theta) {
    log_posterior_pc(theta)  - log_laplace(theta)
}

# Calculate the importance weights
log_importance_ratios <- apply(theta_sim, 1, log_ratio_function)
psis_result <- psis(-log_importance_ratios, r_eff = NA)
weights_smoothed <- exp(psis_result$log_weights)
weights_smoothed_normalized <- weights_smoothed / sum(weights_smoothed)
print(paste("The value of k is:", psis_result$diagnostics$pareto_k))
print(paste("The variance of the importance weights is:", var(weights_smoothed)))

#Compare the mean and covariance of the importance weights with the MAP and -hessian^{-1}
mean_importance_smoothed <- apply(theta_sim, 2, weighted.mean, w = weights_smoothed_normalized)
print(mean_importance_smoothed)
apply(theta_sim, 2, weighted.mean, w = weights_smoothed_normalized)
(mean_importance_smoothed - map_pc$par) / map_pc$par
cov_importance_smoothed <- stats::cov.wt(theta_sim, wt = as.vector(weights_smoothed))$cov
sweep(cov_importance_smoothed - sigma, 2, diag(cov_importance_smoothed), "/")

#No smoothing
weights <- exp(log_importance_ratios)
weights_normalized <- weights / sum(weights)
print(paste("The variance of the importance weights is:", var(weights)))
mean_importance <- apply(theta_sim, 2, weighted.mean, w = weights_normalized)
print(mean_importance)
(mean_importance - map_pc$par) / map_pc$par

cov_importance <- stats::cov.wt(theta_sim, wt = weights)$cov
sweep(cov_importance - sigma, 2, diag(cov_importance), "/")

dim(weights)
dim(weights_smoothed)
sigma
cov_importance_smoothed
cov_importance
