library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)

#Defining the random seed
set.seed(123)

# Defines the upper bounds for the quantiles
a0 <- 10
rho0 <- 0.1
sigma0 <- 10
sigmau0 <- 10
sigmaepsilon0 <- 1

# Defines the quantile
alpha <- 0.01

#Mesh definition
library(sf)
boundary_sf = st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10),c(0,0)))))
boundary = fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(  boundary = boundary,  max.edge = c(0.5, 0.5))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)




m=3 #number of iterations
results <- vector("list", m)  # Pre-allocates a list for s iterations

for (i in 1:m) {
    # Simulate parameters from PC prior
  true_params <- sim_theta_pc(lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u, m = 1)
    # Extract true parameters
  log_kappa <- true_params$log_kappa
    v <- true_params$v
    log_sigma_u <- true_params$log_sigma_u
    log_sigma_epsilon <- true_params$log_sigma_epsilon

  #Sample from noisyy data
  x = fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
A = Matrix::Diagonal(n, 1)
y = A %*% x + exp(log_sigma_epsilon) * stats::rnorm(nrow(Q))

#Calculating the MAP under each prior knowing simulated data
map_pc <- MAP_prior(logprior = log_pc_prior ,mesh = mesh,
           y= y, A = A, m_u =m_u, maxiterations = 10)
map_not_pc <- MAP_prior(logprior = log_not_pc_prior ,mesh = mesh,
                    y= y, A = A, m_u =m_u, maxiterations = 10)

  # PC results
  pc_results <- list(
    MAP_estimate = map_pc$par,                                          # a 5d vector
    distance_vector = abs(MAP_estimate - true_params),                  # a 5d vector
    covariance_estimate = solve(-map_pc$hessian),                     # a 5x5 matrix
    std_dev_estimates = sqrt(diag(covariance_estimate)),                # a 5d vector
    #prob_MAP_greater = pc_prob_map                                      # a 5d vector
  )
  # Not-PC results
  not_pc_results <- list(
    MAP_estimate = map_not_pc$par,                                          # a 5d vector
    distance_vector = abs(MAP_estimate - true_params),                  # a 5d vector
    covariance_estimate = solve(-map_not_pchessian),                     # a 5x5 matrix
    std_dev_estimates = sqrt(diag(covariance_estimate)),                # a 5d vector
    #prob_MAP_greater = not_pc_prob_map                                      # a 5d vector
  )

  # Store results
  results[[i]] <- list(
    true_params = true_params,
    pc = pc_results,
    not_pc = not_pc_results
  )
}

