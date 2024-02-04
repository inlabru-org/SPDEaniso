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
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 2 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

# Calculates the hyperparameters of PC by hand
lambda1 <- lambda1_quantile(a0) # Controls the size of v
lambda <- lambda_quantile(rho0 = rho0, lambda1 = lambda1) # Controls the size of kappa
lambda_u <- lambda_variance_quantile(sigma0 = sigma_u0) # Controls the variance of the field
lambda_epsilon <- lambda_variance_quantile(sigma0 = sigma_epsilon0) # Controls the variance of the noise
hyper_pc <- list(
  lambda = lambda, lambda1 = lambda1, lambda_u = lambda_u,
  lambda_epsilon = lambda_epsilon
)

# Simulates the "true parameters" from the pc_prior.
true_params <- sim_theta_pc_quantile(
  alpha = alpha, sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = 1
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
print("The standard deviation of the field exp(log_sigma_u) is smaller than sigma_u0")
print(exp(log_sigma_u) < sigma_u0)
print("The standard deviation of the noise exp(log_sigma_epsilon) is smaller than sigma_epsilon0")
print(exp(log_sigma_epsilon) < sigma_epsilon0)


# Calculates the hyperparameters of not PC by hand
sigma2_quantile_v(a0)
lambda_quantile_kappa(rho0)

# Setting the mean of the field
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

# Uniform prior
L <- 10
log_uniform_prior <- log_prior_uniform(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L)

log_priors <- list(
  pc = log_pc_prior,
  not_pc = log_not_pc_prior,
  uniform = log_uniform_prior
)
prior_types <- setNames(as.list(names(log_priors)), names(log_priors))

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
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(3, 3))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

# Sample from noisy data
aniso <- list(rep(kappa, n), matrix(v, n, 2))
x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
A <- Matrix::Diagonal(n, 1)
y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)

# Testing the MAP_prior function for log_prior = log_pc_prior
maxit <- 600
tryCatch({
  delta <- rnorm(5, 0, 1) # Used to randomize starting point of MAP
  # Calculating the MAP under each prior knowing simulated data
  map_pc <- MAP_prior(
    log_prior = log_uniform_prior, mesh = mesh,
    y = y, A = A, m_u = m_u, max_iterations = maxit,
    theta0 = unlist(true_params),
    do_u_want_hessian = FALSE
    # ,log_sigma_epsilon = log_sigma_epsilon
  )

  error <- function(e) {}
})
#### TESTING MAP_prior ####

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

print("Checking if map is larger than posterior value at true parameters. ")
map_pc$value > log_posteriors$pc(log_kappa, v, log_sigma_u, log_sigma_epsilon)

print("Checking if value for map is the same as posterior value at map parameters. ")
map_pc$value == log_posteriors$pc(map_pc$par[1], map_pc$par[2:3], map_pc$par[4], map_pc$par[5])

# print("Checking if MAP_prior with PC priors returns the same thing as MAP using PC priors defined through hyperparameters")

# map_hyper <- MAP_pc(mesh, lambda, lambda1, lambda_epsilon, lambda_u,
#   y, A, m_u,
#   max_iterations = maxit, theta0 = unlist(true_params) + delta
# )
# map_hyper$value == map_pc$value & map_hyper$par == map_pc$par



plotter <- function(map = map_pc, log_priors = log_priors, log_posteriors = log_posteriors, l = 2, n_points = 10, together = TRUE) {
  function_types <- list(prior = "prior", posterior = "posterior")
  ########## UNNORMALIZED Gaussian_median APPROXIMATION TO THE POSTERIOR############
  log_Gaussian_median <- function(theta) {
    logGdensity(
      x = theta, mu = map_pc$par, Q = -map_pc$hessian
    ) - logGdensity(
      x = map_pc$par, mu = map_pc$par, Q = -map_pc$hessian
    )
  }
  ### UNNORMALIZED LOG FUNCTION SO THEY ALL START AT 0###
  unnormalize_prior_and_posterior <- function(log_prior) {
    function(log_kappa, v1, v2, log_sigma_u, log_sigma_epsilon) {
      log_prior(
        log_kappa = log_kappa, v = c(v1, v2),
        log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
      ) - log_prior(
        log_kappa = map_pc$par[1], v = map_pc$par[2:3],
        log_sigma_u = map_pc$par[4], log_sigma_epsilon = map_pc$par[5]
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
    lapply(f, restricting_function_to_one_parameter, map_pc$par)
  })


  # Getting data for plotting
  partitions <- lapply(seq_along(map_pc$par), function(i) {
    seq(map_pc$par[i] - l, map_pc$par[i] + l, length.out = n_points)
  })
  names(partitions) <- names(map_pc$par)

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

  # Create a list to store the ggplot objects
  plots <- list()

  # For each parameter, create a ggplot of the functions
  for (parameter in unique(plot_data$Parameter)) {
    # Subset the data for the current parameter
    parameter_data <- subset(plot_data, Parameter == parameter)

    # Create a ggplot
    # Create a ggplot
    p <- ggplot(parameter_data, aes(x = x, y = Value, color = PriorType, linetype = FunctionType)) +
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
    plots
  }
}

plotter(map = map_pc, log_priors = log_priors, log_posteriors = log_posteriors, l = 4, n_points = 20, together = FALSE)
