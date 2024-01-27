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
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(2, 2))
nodes <- mesh$loc
n <- mesh$n
# plot(mesh)

# ---------------------- SIMULATING DATA Y ----------------------

# Defines the upper bounds for the quantiles
rho0 <- 1 # Controls the size of kappa in PC and non PC priors
a0 <- 2 # Controls the size of v in PC and non PC priors
sigma_u0 <- 2 # controls standard deviation of field
sigma_epsilon0 <- 0.1 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

# Simulates the "true parameters" from the pc_prior.
true_params <- sim_theta_pc_quantile(
    alpha = alpha, sigma_u0 = sigma_u0,
    sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = 1
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
    sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
    a0 = a0, rho0 = rho0, alpha = alpha
)
# testing
# Defining the log posterior density as a function of the parameters using the log_pc_prior
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

########## Gaussian_median APPROXIMATION TO THE POSTERIOR############
maxit <- 600
tryCatch({
    # Calculating the MAP under each prior knowing simulated data
    map_pc <- MAP_prior(
        log_prior = log_pc_prior, mesh = mesh,
        y = y, A = A, m_u = m_u, max_iterations = maxit,
        theta0 = unlist(true_params)
        # ,log_sigma_epsilon = log_sigma_epsilon
    )

    error <- function(e) {}
})

Q_Gaussian_median <- -map_pc$hessian
covariance_Gaussian_median <- solve(Q_Gaussian_median)



log_Gaussian_median <- function(theta, Q) {
    logGdensity(
        x = theta, mu = map_pc$par, Q = Q
    )
}
# ---------------------- IMPORTANCE SAMPLING ----------------------
# Simulate from Gaussian_median approximation
n_sim <- 1000
precision_multiplier <- 1.2
theta_sim <- MASS::mvrnorm(n_sim, map_pc$par, covariance_Gaussian_median / precision_multiplier)
log_ratio_function <- function(theta) {
    log_posterior_pc(theta) - log_Gaussian_median(theta, Q_Gaussian_median * precision_multiplier)
}

# Calculate the importance weights
log_importance_ratios <- apply(theta_sim, 1, log_ratio_function)

# Subtract the mean to avoid numerical issues as it shouldn't change the result
log_importance_ratios_2 <- log_importance_ratios - max(log_importance_ratios)


psis_result <- psis(log_importance_ratios_2, r_eff = NA)
# subsample_ratios <- sample(log_importance_ratios_2, 1000, replace = FALSE)
# psis_result <- psis(subsample_ratios, r_eff = NA)
psis_result$diagnostics

# Plots the log_weights against the log posterior to see when the weights are small and large depending on the posterior. Smoothed and unsmoothed to see how the smoothing works
log_posterior_pc_values <- apply(theta_sim, 1, log_posterior_pc)
ggplot(data.frame(lppv = log_posterior_pc_values, psis_lw = psis_result$log_weights, lw = log_importance_ratios_2), aes(x = lppv, y = psis_lw)) +
    geom_point(aes(color = "psis")) +
    geom_point(aes(y = lw, color = "original")) +
    ggtitle("Pareto smoothed importance sampling") +
    xlab("log posterior") +
    ylab("log smoothed weights")



# Comparison of weights vs smoothed weights
plot(exp(log_importance_ratios_2), exp(psis_result$log_weights), main = "smoothed weights", pch = 20)
plot.ecdf(exp(log_importance_ratios_2) / mean(exp(log_importance_ratios_2)), col = "blue")
plot.ecdf(exp(psis_result$log_weights) / mean(exp(psis_result$log_weights)), add = TRUE, col = "red")



## Histogram of log_weights and comparison with smoothed weights
# hist(psis_result$log_weights)
# hist(log_importance_ratios_2)
# hist(log_importance_ratios_2 - psis_result$log_weights)
# log_importance_ratios_2 - psis_result$log_weights
