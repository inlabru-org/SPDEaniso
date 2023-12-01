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

#Anisotropy parameters
kappa <- exp(-0.5); log_kappa <- log(kappa);
v <- c(0,0.01)

#Correlation range calculation
rho <- sqrt(8)/kappa/sqrt(exp(sqrt(v[1]^2 + v[2]^2)))

#Field and noise parameters
sigma_u <- 0.5; log_sigma_u <- log(sigma_u)
sigma_epsilon <- exp(-3) ; log_sigma_epsilon <- log(sigma_epsilon)

#Mesh definition
library(sf)
boundary_sf = st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10),c(0,0)))))
boundary = fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(  boundary = boundary,  max.edge = c(0.5, 0.5))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)

#Defining anisotropy
kappa_values <- rep(kappa, n)
vec_values <- matrix(v, n, 2, byrow = TRUE)
aniso <- list(kappa = kappa_values, vec = vec_values)
#Calculating precision and defining mean
Q <-fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
m_u = as.vector(rep(0,n))

#Sampling from noisy data
x = fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
A = Matrix::Diagonal(n, 1)
y = A %*% x + exp(log_sigma_epsilon) * stats::rnorm(nrow(Q))

#Calculates the log prior of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
                                      a0 = a0, rho0 = rho0, alpha = alpha)

log_not_pc_prior <- log_gaussian_prior_quantile(sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0,
                                       a0 = a0, rho0 = rho0, alpha = alpha)

#Calculates the log posterior for PC priors and non PC priors
log_posterior_pc_value <- log_posterior_prior(logprior = log_pc_prior, mesh = mesh,
log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon, y = y, A = A, m_u = m_u)

log_posterior_notpc_value <- log_posterior_prior(logprior = log_not_pc_prior,
mesh = mesh, log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon, y = y, A = A, m_u = m_u)




#Optimizing over log(kappa), v, log(sigma_u) log(sigma_noise)
map_pc <- MAP_prior(logprior = log_pc_prior ,mesh = mesh,
           y= y, A = A, m_u =m_u, maxiterations = 2000)
map_not_pc <- MAP_prior(logprior = log_not_pc_prior ,mesh = mesh,
                    y= y, A = A, m_u =m_u, maxiterations = 1000)
print(map_pc)
cov2cor(solve(-map_pc$hessian))
par_full <- map_pc$par
real_par_full <- c(log_kappa,v,log_sigma_u, log_sigma_epsilon)
print(map_pc$par-real_par_full)
sqrt(diag(solve(-map_pc$hessian)))

print(map_not_pc)
cov2cor(solve(-map_not_pc$hessian))
par_full <- map_not_pc$par
real_par_full <- c(log_kappa,v,log_sigma_u, log_sigma_epsilon)
print(map_not_pc$par-real_par_full)
sqrt(diag(solve(-map_not_pc$hessian)))


m=3 #number of iterations
results <- vector("list", m)  # Pre-allocates a list for s iterations

for (i in 1:m) {
    # Simulate parameters from PC prior
  true_params <- sim_theta_pc(lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u, m = 1)
  #



  

  # PC results
  pc_results <- list(
    MAP_estimates = pc_map_est,            # a 5d vector
    covariance_estimates = pc_cov_est,     # a 5x5 matrix
    std_dev_estimates = pc_std_dev,        # a 5d vector
    prob_MAP_greater = pc_prob_map         # a single value or 5d vector
  )

  # Not-PC results
  not_pc_results <- list(
    MAP_estimates = not_pc_map_est,        # a 5d vector
    covariance_estimates = not_pc_cov_est, # a 5x5 matrix
    std_dev_estimates = not_pc_std_dev,    # a 5d vector
    prob_MAP_greater = not_pc_prob_map     # a single value or 5d vector
  )

  # Store results
  results[[i]] <- list(
    true_params = true_params,
    pc = pc_results,
    not_pc = not_pc_results
  )
}

