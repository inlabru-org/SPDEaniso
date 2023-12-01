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
a0 <- 2
rho0 <- 0.1
sigma0 <- 1.5
sigmau0 <- 2
sigmaepsilon0 <- 0.1

# Defines the quantile
alpha <- 0.01
#Setting mean of the field
m_u <-0
#Calculates the log prior of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
                                      a0 = a0, rho0 = rho0, alpha = alpha)

log_not_pc_prior <- log_gaussian_prior_quantile(sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
                                                a0 = a0, rho0 = rho0, alpha = alpha)

#Mesh definition
library(sf)
boundary_sf = st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10),c(0,0)))))
boundary = fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(  boundary = boundary,  max.edge = c(0.75, 0.75))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)
# set.seed(123)
# debugonce(sim_theta_pc)
# true_params <- sim_theta_pc_quantile(alpha = alpha, sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1)
# lambda1<- lambda1_quantile(a0,alpha)
# lambda_quantile(rho0 = rho0,lambda1 = lambda1,alpha=alpha)

m=1 #number of iterations
maxit=10
results <- vector("list", m)  # Pre-allocates a list for s iterations

for (i in 1:m) {
  tryCatch({
    # Simulate parameters from PC prior
  true_params <- sim_theta_pc_quantile(alpha = alpha, sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, m = 1)
    # Extract true parameters
  log_kappa <- true_params$log_kappa
  kappa <- exp(log_kappa)
  v <- true_params$v
  log_sigma_u <- true_params$log_sigma_u
  log_sigma_epsilon <- true_params$log_sigma_epsilon
  aniso <- list(rep(kappa, n), matrix(v, n, 2))

  #Sample from noisyy data
  x = fm_aniso_basis_weights_sample(x = mesh, aniso = aniso)
  A = Matrix::Diagonal(n, 1)
  y = A %*% x + exp(log_sigma_epsilon) * stats::rnorm(n)

  #Calculating the MAP under each prior knowing simulated data
  map_pc <- MAP_prior(logprior = log_pc_prior ,mesh = mesh,
           y= y, A = A, m_u =m_u, maxiterations = maxit)

  map_not_pc <- MAP_prior(logprior = log_not_pc_prior ,mesh = mesh,
                    y= y, A = A, m_u =m_u, maxiterations = maxit)

  # PC results
  pc_results <- list(
    MAP_estimate = map_pc$par,                                          # a 5d vector
    distance_vector = abs(map_pc$par - unlist(true_params)),                  # a 5d vector
    covariance_estimate = solve(-map_pc$hessian),                     # a 5x5 matrix
    std_dev_estimates = sqrt(diag(solve(-map_pc$hessian)))                # a 5d vector
    #prob_MAP_greater = pc_prob_map                                      # a 5d vector
  )
  # Not-PC results
  not_pc_results <- list(
    MAP_estimate = map_not_pc$par,                                          # a 5d vector
    distance_vector = abs(map_not_pc$par - unlist(true_params)),                  # a 5d vector
    covariance_estimate = solve(-map_not_pchessian),                     # a 5x5 matrix
    std_dev_estimates = sqrt(diag(solve(-map_pc$hessian)))                # a 5d vector
    #prob_MAP_greater = not_pc_prob_map                                      # a 5d vector
  )

  # Store results
  results[[i]] <- list(
    true_params = true_params,
    pc = pc_results,
    not_pc = not_pc_results
  )
  },
  error= function(e){})
}

