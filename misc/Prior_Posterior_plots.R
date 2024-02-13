library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(cowplot)
document()
# Defining the random seed
set.seed(223)

# Defines the upper bounds for the quantiles
# rho0 <- 1 # Controls the size of kappa in PC and non PC priors
rho0 <- 1
# a0 <- 2 # Controls the size of v in PC and non PC priors
a0 <- 2
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 0.5 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

m_u <- 0

# Defines the log prior density function of theta for PC and non-PC priors
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
log_beta_prior <- log_prior_beta(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, shape = shape, width_support_factor = width_beta)
log_priors <- list(
  pc = log_pc_prior,
  not_pc = log_not_pc_prior,
  uniform = log_uniform_prior,
  beta = log_beta_prior
)
prior_types <- setNames(as.list(names(log_priors)), names(log_priors))



# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 3))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)
n_observations <- 15
observations <- 10 * matrix(runif(n_observations * 2), ncol = 2)
A <- fm_basis(mesh, loc = observations)

# Defining the log posterior
log_posteriors <- lapply(log_priors, function(log_prior) {
  function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
    log_posterior_prior(
      log_prior = log_prior,
      mesh = mesh, log_kappa = log_kappa, v = v,
      log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
      y = y, A = A, m_u = m_u
    )
  }
})
maxit <- 600
tryCatch({
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
  x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
  y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(nrow(A))
  # true_params <- sim_theta_beta(
  #   sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  #   a0 = a0, rho0 = rho0, L = L, shape = shape,
  #   width_support_factor = width_beta)

  delta <- rnorm(5, 0, 1) # Used to randomize starting point of MAP
  # Calculating the MAP under each prior knowing simulated data
  map_pc <- MAP_prior(
    log_prior = log_pc_prior, mesh = mesh,
    y = y, A = A, m_u = m_u, max_iterations = maxit,
    theta0 = unlist(true_params) + delta
    # ,log_sigma_epsilon = log_sigma_epsilon
  )

  error <- function(e) {}
})

# The plots are different when rho0 and a0 are large. eg (2,4) more for (4,10). But are still quite similar.
prior_posterior_plotter(theta_fixed = map_pc$par, log_priors = log_priors, log_posteriors = log_posteriors, l = 3, n_points = 200, n_parameters_to_plot = 2, path = "Simulation_images/prior_posterior_plots.pdf")
