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
#rho0 <- 1 # Controls the size of kappa in PC and non PC priors
rho0 <- 1
#a0 <- 2 # Controls the size of v in PC and non PC priors
a0 <- 2
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 2 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

# Simulates the "true parameters" from the pc_prior.
true_params <- sim_theta_pc_quantile(
  alpha = alpha, sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = 1
)




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
a0_inf <- 1.01
log_uniform_prior <- log_prior_uniform(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, a0_inf = a0_inf, rho0 = rho0, L = L, width_support_factor = width_uniform)
shape <- 1.1
width_beta <- 20
log_beta_prior <- log_prior_beta(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, a0_inf = a0_inf, rho0 = rho0, L = L, shape = shape, width_support_factor = width_beta)

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
observations <- 10*matrix(runif(n_observations * 2), ncol = 2)
A <- fm_basis(mesh, loc = observations)

# Sample from noisy data
aniso <- list(rep(kappa, n), matrix(v, n, 2))
x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(nrow(A))

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
plotter <- function(theta_fixed = map_pc$par, log_priors = log_priors, log_posteriors = log_posteriors, l = 2, n_points = 10, together = TRUE, n_parameters_to_plot = 3) {
  function_types <- list(prior = "prior", posterior = "posterior")
  ########## UNNORMALIZED Gaussian_median APPROXIMATION TO THE POSTERIOR############

  ### UNNORMALIZED LOG FUNCTION SO THEY ALL START AT 0###
  unnormalize_prior_and_posterior <- function(log_prior) {
    function(log_kappa, v1, v2, log_sigma_u, log_sigma_epsilon) {
      log_prior(
        log_kappa = log_kappa, v = c(v1, v2),
        log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
      ) - log_prior(
        log_kappa = theta_fixed[1], v = theta_fixed[2:3],
        log_sigma_u = theta_fixed[4], log_sigma_epsilon = theta_fixed[5]
      )
    }
  }

  unnormalized_priors <- lapply(log_priors, unnormalize_prior_and_posterior)
  unnormalized_posteriors <- lapply(log_posteriors, unnormalize_prior_and_posterior)
  unnormalized_priors_and_posteriors <- list(prior = unnormalized_priors, posterior = unnormalized_posteriors)

  # Restricting the functions to one parameter
  restricting_function_to_one_parameter <- function(f, x0) {
    f_list <- lapply(seq_along(x0), function(i) {
      function(x) {
        x0_copy <- x0
        x0_copy[i] <- x
        unname(do.call(f, as.list(x0_copy)))
      }
    })
    names(f_list) <- names(x0)
    f_list
  }
  restricted_priors_and_posteriors <- lapply(unnormalized_priors_and_posteriors, function(f) {
    lapply(f, restricting_function_to_one_parameter, theta_fixed)
  })


  # Getting data for plotting
  partitions <- lapply(seq_along(theta_fixed), function(i) {
    seq(theta_fixed[i] - l, theta_fixed[i] + l, length.out = n_points)
  })
  names(partitions) <- names(theta_fixed)

  plot_data <- do.call(rbind, lapply(names(restricted_priors_and_posteriors), function(function_type) {
    do.call(rbind, lapply(names(restricted_priors_and_posteriors[[function_type]]), function(prior_type) {
      do.call(rbind, lapply(names(restricted_priors_and_posteriors[[function_type]][[prior_type]]), function(parameter_name) {
        # Calculate the function values
        x_values <- partitions[[parameter_name]]
        y_values <- sapply(x_values, restricted_priors_and_posteriors[[function_type]][[prior_type]][[parameter_name]])

        data.frame(
          x = x_values,
          Value = y_values,
          Parameter = parameter_name,
          FunctionType = function_type,
          PriorType = prior_type,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))

  plots <- list()
  for (parameter in unique(plot_data$Parameter)) {
    parameter_data <- subset(plot_data, Parameter == parameter)

    p <- ggplot(parameter_data, aes(x = x, y = exp(Value), color = PriorType, linetype = FunctionType)) +
      geom_line() +
      geom_vline(xintercept = unlist(true_params)[parameter], color = "purple") + # Move xintercept outside of aes()
      labs(title = paste("Unnormalized Priors and Posteriors for", parameter), x = parameter, y = "Log density") +
      theme_minimal()

    # Add the plot to the list
    plots[[parameter]] <- p
  }
  if (together) {
    require(gridExtra)
    do.call(grid.arrange, c(plots, ncol = 3))
  } else {
    # Print the plots
    plots[1:n_parameters_to_plot]
  }
}
alpha <- 0.05

#The plots are different when rho0 and a0 are large. eg (2,4) more for (4,10). But are still quite similar.


plotter(theta_fixed = map_pc$par, log_priors = log_priors, log_posteriors = log_posteriors, l = 2, n_points = 50, together = FALSE, n_parameters_to_plot = 2)
