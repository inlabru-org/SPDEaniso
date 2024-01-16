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

# Subtract the mean to avoid numerical issues as it shouldn't change the result
log_importance_ratios_2 <- log_importance_ratios - mean(log_importance_ratios)

psis_result <- psis(-log_importance_ratios_2, r_eff = NA)
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
weights <- exp(log_importance_ratios_2)
weights_normalized <- weights / sum(weights)
print(paste("The variance of the importance weights is:", var(weights)))
mean_importance <- apply(theta_sim, 2, weighted.mean, w = weights_normalized)
print(mean_importance)
(mean_importance - map_pc$par) / map_pc$par

cov_importance <- stats::cov.wt(theta_sim, wt = weights)$cov
sweep(cov_importance - sigma, 2, diag(cov_importance), "/")


## We repeat the same but now only let log_kappa vary
log_laplace_log_kappa <- function(log_kappa) {
    log_laplace(c(log_kappa, map_pc$par[2:5]))
}
log_posterior_pc_log_kappa <- function(log_kappa) {
    log_posterior_pc(c(log_kappa, map_pc$par[2:5]))
}
log_ratio_function_log_kappa <- function(log_kappa) {
    log_posterior_pc_log_kappa(log_kappa)  - log_laplace_log_kappa(log_kappa)
}

# Simulate from Laplace approximation
n_sim <- 1000
log_kappa_sim <- rnorm(n_sim, map_pc$par[1], sqrt(sigma[1, 1]))
log_importance_ratios_log_kappa <- sapply(log_kappa_sim, log_ratio_function_log_kappa)

# Subtract the mean to avoid numerical issues as it shouldn't change the result
log_importance_ratios_2_log_kappa <- log_importance_ratios_log_kappa - mean(log_importance_ratios_log_kappa)

psis_result_log_kappa <- psis(-log_importance_ratios_2_log_kappa, r_eff = NA)
weights_smoothed_log_kappa <- exp(psis_result_log_kappa$log_weights)
weights_smoothed_normalized_log_kappa <- weights_smoothed_log_kappa / sum(weights_smoothed_log_kappa)
print(paste("The value of k is:", psis_result_log_kappa$diagnostics$pareto_k))
print(paste("The variance of the importance weights is:", var(weights_smoothed_log_kappa)))

##We repeat the same but now only let v vary
log_laplace_v <- function(v) {
    log_laplace(c(map_pc$par[1], v, map_pc$par[4:5]))
}
log_posterior_pc_v <- function(v) {
    log_posterior_pc(c(map_pc$par[1], v, map_pc$par[4:5]))
}
log_ratio_function_v <- function(v) {
    log_posterior_pc_v(v)  - log_laplace_v(v)
}

# Simulate from Laplace approximation
n_sim <- 1000
v_sim <- MASS::mvrnorm(n_sim, map_pc$par[2:3], sigma[2:3, 2:3])
log_importance_ratios_v <- apply(v_sim, 1, log_ratio_function_v)

# Subtract the mean to avoid numerical issues as it shouldn't change the result
log_importance_ratios_2_v <- log_importance_ratios_v - mean(log_importance_ratios_v)

psis_result_v <- psis(-log_importance_ratios_2_v, r_eff = NA)
weights_smoothed_v <- exp(psis_result_v$log_weights)
weights_smoothed_normalized_v <- weights_smoothed_v / sum(weights_smoothed_v)
print(paste("The value of k is:", psis_result_v$diagnostics$pareto_k))
print(paste("The variance of the importance weights is:", var(weights_smoothed_v)))

##We repeat the same but now only let log_sigma_u vary
log_laplace_log_sigma_u <- function(log_sigma_u) {
    log_laplace(c(map_pc$par[1:3], log_sigma_u, map_pc$par[5]))
}
log_posterior_pc_log_sigma_u <- function(log_sigma_u) {
    log_posterior_pc(c(map_pc$par[1:3], log_sigma_u, map_pc$par[5]))
}
log_ratio_function_log_sigma_u <- function(log_sigma_u) {
    log_posterior_pc_log_sigma_u(log_sigma_u)  - log_laplace_log_sigma_u(log_sigma_u)
}

# Simulate from Laplace approximation
n_sim <- 1000
log_sigma_u_sim <- rnorm(n_sim, map_pc$par[4], sqrt(sigma[4, 4]))
log_importance_ratios_log_sigma_u <- sapply(log_sigma_u_sim, log_ratio_function_log_sigma_u)

# Subtract the mean to avoid numerical issues as it shouldn't change the result
log_importance_ratios_2_log_sigma_u <- log_importance_ratios_log_sigma_u - mean(log_importance_ratios_log_sigma_u)

psis_result_log_sigma_u <- psis(-log_importance_ratios_2_log_sigma_u, r_eff = NA)
weights_smoothed_log_sigma_u <- exp(psis_result_log_sigma_u$log_weights)
weights_smoothed_normalized_log_sigma_u <- weights_smoothed_log_sigma_u / sum(weights_smoothed_log_sigma_u)
print(paste("The value of k is:", psis_result_log_sigma_u$diagnostics$pareto_k))
print(paste("The variance of the importance weights is:", var(weights_smoothed_log_sigma_u)))

##We repeat the same but now only let log_sigma_epsilon vary
log_laplace_log_sigma_epsilon <- function(log_sigma_epsilon) {
    log_laplace(c(map_pc$par[1:4], log_sigma_epsilon))
}
log_posterior_pc_log_sigma_epsilon <- function(log_sigma_epsilon) {
    log_posterior_pc(c(map_pc$par[1:4], log_sigma_epsilon))
}
log_ratio_function_log_sigma_epsilon <- function(log_sigma_epsilon) {
    log_posterior_pc_log_sigma_epsilon(log_sigma_epsilon)  - log_laplace_log_sigma_epsilon(log_sigma_epsilon)
}

# Simulate from Laplace approximation
n_sim <- 1000
log_sigma_epsilon_sim <- rnorm(n_sim, map_pc$par[5], sqrt(sigma[5, 5]))
log_importance_ratios_log_sigma_epsilon <- sapply(log_sigma_epsilon_sim, log_ratio_function_log_sigma_epsilon)

# Subtract the mean to avoid numerical issues as it shouldn't change the result
log_importance_ratios_2_log_sigma_epsilon <- log_importance_ratios_log_sigma_epsilon - mean(log_importance_ratios_log_sigma_epsilon)

psis_result_log_sigma_epsilon <- psis(-log_importance_ratios_2_log_sigma_epsilon, r_eff = NA)
weights_smoothed_log_sigma_epsilon <- exp(psis_result_log_sigma_epsilon$log_weights)
weights_smoothed_normalized_log_sigma_epsilon <- weights_smoothed_log_sigma_epsilon / sum(weights_smoothed_log_sigma_epsilon)
print(paste("The value of k is:", psis_result_log_sigma_epsilon$diagnostics$pareto_k))
print(paste("The variance of the importance weights is:", var(weights_smoothed_log_sigma_epsilon)))


#now we repeat the same but let all the parameters vary except one, we define a general function for this
pareto_smoothing <- function(i) {
    # Simulate from Laplace approximation in four variables and keep i fixes
    n_sim <- 500
    sigma <- sigma[-i, -i]
    non_fixed_params <- MASS::mvrnorm(n_sim, map_pc$par[-i], sigma)
    if (i == 1) {
        theta_sim <- cbind(rep(map_pc$par[i], n_sim), non_fixed_params)
    } else if(i == 5) {
        theta_sim <- cbind(non_fixed_params, rep(map_pc$par[i], n_sim))
    } else {
        theta_sim <- cbind(non_fixed_params[, 1:(i - 1)], rep(map_pc$par[i], n_sim), non_fixed_params[, i:4])
    }
 log_importance_ratios <- apply(theta_sim, 1, log_ratio_function)

    # Subtract the mean to avoid numerical issues as it shouldn't change the result
    log_importance_ratios_2 <- log_importance_ratios - mean(log_importance_ratios)

    psis_result <- psis(-log_importance_ratios_2, r_eff = NA)
    weights_smoothed <- exp(psis_result$log_weights)
    weights_smoothed_normalized <- weights_smoothed / sum(weights_smoothed)
    print(paste("The value of k is:", psis_result$diagnostics$pareto_k))
    print(paste("The variance of the importance weights is:", var(weights_smoothed)))
    # return(weights_smoothed_normalized)
}
pareto_smoothing(1)

