library(devtools)
# Defines the upper bounds for the quantiles
a0 <- 10
rho0 <- 0.1
sigma0 <- 10
sigmau0 <- 10
sigmaepsilon0 <- 1
# Defines the quantile
alpha <- 0.01
#Defines the true value of the Parameters
kappa <- exp(-0.5); log_kappa <- log(kappa);
v <- c(0,0.0001)
log_sigma_epsilon <- 0.4
log_sigma_u <- 0.1



##PC PRIORS##
# Finds the hyperparameters for the PC priors on the anisotropy given the quantiles
lambda1 <- lambda1_quantile(a0 =a0,alpha_a = alpha)
lambda <- lambda_quantile(alpha=alpha, rho0=rho0, lambda1=lambda1)
#Calulate the hyperparameters of the PC priors for variances
lambda_u <- lambda_variance_quantile(alpha_sigma = alpha, sigma0 = sigmau0)
lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha, sigma0 = sigmaepsilon0)

#GAUSSIAN PRIORS#
# Finds the hyperparameters for the anisotropy given the quantiles
sigma2_v <- sigma2_quantile_v(alpha_v = alpha, a0 = a0)
lambda_k <- lambda_quantile_kappa(alpha_k = alpha, rho0 = rho0)
# Finds the hyperparameters for the variances given the quantiles
lambda_u <- lambda_variance_quantile(alpha_sigma = alpha, sigma0 = sigmau0)
lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha, sigma0 = sigmaepsilon0)
# Calculates the log prior for the PC priors
log_pc_prior <- log_pc_prior_quantile(sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, alpha = alpha)

# Calculates the log prior for the Gaussian priors
log_gaussian_prior <- log_gaussian_prior_quantile(sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0, alpha = alpha)

# Calculates the log prior at the true value of the parameters for each prior

log_pc_prior(log_kappa, v, log_sigma_u, log_sigma_epsilon)
log_gaussian_prior(log_kappa, v, log_sigma_u, log_sigma_epsilon)

