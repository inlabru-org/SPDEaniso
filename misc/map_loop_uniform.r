library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(future)
library(future.apply)
library(dplyr)
library(tidyr)
library(loo)
document()

# Set the parallel plan to use all local cores (currently not used as future package doesn't recognize functions in prior.R)
plan(multisession)
# Defining the random seed
set.seed(123)

# Defines the upper bounds for the quantiles
rho0 <- 1 # Controls the size of kappa
a0 <- 2 # Controls the size of v
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 2 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors
# Defines the quantile
alpha <- 0.01

# Setting mean of the field

m_u <- 0
# Calculates the log prior density function of theta for PC and non-PC priors
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
# width_beta <- 2 # Changing to not get extreme values in simulations
log_beta_prior <- log_prior_beta(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, shape = shape, width_support_factor = width_beta)
log_priors <- list(
    pc = log_pc_prior,
    not_pc = log_not_pc_prior,
    uniform = log_uniform_prior,
    beta = log_beta_prior
)
prior_types <- setNames(as.list(names(log_priors)), names(log_priors))
approximation_types <- list("Gaussian_median", "importance", "importance_smoothed")
approximation_types <- setNames(approximation_types, approximation_types)
# Mesh definition
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)))))
boundary <- fm_as_segm(boundary_sf)
mesh <- fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 3))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)
n_observations <- 15
observations <- L * matrix(runif(n_observations * 2), ncol = 2)
A <- fm_basis(mesh, loc = observations)

number_of_loops <- 200 # number of iterations
maxit_MAP <- 600
number_of_weights <- 5000
credible_level <- 0.05
results_uniform <- vector("list", number_of_loops) # Pre-allocates a list for m iterations

for (i in 1:number_of_loops) {
    start_time <- Sys.time()
    tryCatch(
        {
            # Simulate parameters from beta
            true_params <- sim_theta_uniform(
                sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
                a0 = a0, rho0 = rho0, L = L,
                width_support_factor = 2
            )


            log_kappa <- true_params$log_kappa
            kappa <- exp(log_kappa)
            v <- true_params$v
            log_sigma_u <- true_params$log_sigma_u
            log_sigma_epsilon <- true_params$log_sigma_epsilon
            aniso <- list(rep(kappa, n), matrix(v, n, 2))

            # Sample from noisy data
            x <- fm_aniso_basis_weights_sample(x = mesh, aniso = aniso, log_sigma = log_sigma_u)
            y <- A %*% x + rnorm(n_observations, 0, exp(log_sigma_epsilon))

            # delta <- rnorm(5, 0, 1) # Used to randomize starting point of MAP
            maps <- lapply(log_priors, function(log_prior) {
                MAP_prior(
                    log_prior = log_prior, mesh = mesh,
                    y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
                    theta0 = unlist(true_params) + delta
                )
            })

            # Gaussian_median approximations
            mus_Gaussian_median <- lapply(maps, function(map) map$par)
            Qs_Gaussian_median <- lapply(maps, function(map) -map$hessian)
            Covariances_Gaussian_median <- lapply(Qs_Gaussian_median, function(Q) solve(Q))
            std_dev_estimates_Gaussian_median <- lapply(Covariances_Gaussian_median, function(Covariance) sqrt(diag(Covariance)))

            m_u <- 0

            log_posteriors <- lapply(log_priors, function(log_prior) {
                function(theta) {
                    log_kappa <- theta[1]
                    v <- theta[2:3]
                    log_sigma_u <- theta[4]
                    log_sigma_epsilon <- theta[5]
                    log_posterior_prior(
                        log_prior = log_prior,
                        mesh = mesh, log_kappa = log_kappa, v = v,
                        log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
                        y = y, A = A, m_u = m_u
                    )
                }
            })

            # Importance sampling
            importances <- lapply(prior_types, function(prior_type) {
                log_posterior <- log_posteriors[[prior_type]]
                mu_Gaussian_median <- mus_Gaussian_median[[prior_type]]
                Q_Gaussian_median <- Qs_Gaussian_median[[prior_type]]
                log_unnormalized_importance_weights_and_integrals(
                    log_posterior_density = log_posterior,
                    mu_Gaussian_median = mu_Gaussian_median, Q_Gaussian_median = Q_Gaussian_median,
                    n_weights = number_of_weights, q = credible_level, true_params = unlist(true_params)
                )
            })

            # CIs
            credible_intervals <- lapply(prior_types, function(prior_type) {
                lapply(approximation_types, function(approximation_type) {
                    importances[[prior_type]][[paste0("credible_intervals_", approximation_type)]]
                })
            })

            true_parameter_is_within_CI <- lapply(prior_types, function(prior_type) {
                lapply(approximation_types, function(approximation_type) {
                    parameter_within_credible_intervals(true_params, credible_intervals[[prior_type]][[approximation_type]])
                })
            })


            # Accumulate results
            results_accumulator <- function(prior_type) {
                # Return a list of the calculated values
                list(
                    MAP_estimate = maps[[prior_type]]$par,
                    MAP_value = maps[[prior_type]]$value,
                    convergence = maps[[prior_type]]$convergence,
                    distance_vector = abs(maps[[prior_type]]$par - unlist(true_params)),
                    covariance_estimate = Covariances_Gaussian_median[[prior_type]],
                    std_dev_estimates_Gaussian_median = std_dev_estimates_Gaussian_median[[prior_type]],
                    credible_intervals = lapply(approximation_types, function(approximation_type) {
                        credible_intervals[[prior_type]][[approximation_type]]
                    }),
                    true_parameter_within_c_interval = lapply(approximation_types, function(approximation_type) {
                        true_parameter_is_within_CI[[prior_type]][[approximation_type]]
                    }),
                    importance = importances[[prior_type]]
                )
            }


            # Store results
            partial_results <- lapply(prior_types, results_accumulator)
            results_uniform[[i]] <- c(
                list(true_params = true_params),
                partial_results, prior_types
            )
        },
        error = function(e) {
            e
        }
    )
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    formatted_time <- seconds_to_hours_minutes_seconds(execution_time)
    print(paste("Iteration", i, "took", formatted_time))
    time_left <- as.numeric(difftime(end_time, start_time, units = "secs")) * (number_of_loops - i)
    print(paste("Estimated time left:", seconds_to_hours_minutes_seconds(time_left)))
}
# Eliminates NULL results
not_null_indices <- sapply(results_uniform, function(x) !is.null(x$pc$importance$log_unnormalized_weights))
results_uniform <- results_uniform[not_null_indices]
# Results obtained simulating parameters from uniform priors and using a mesh size of 1, 15 observations, 200 iterations, 5000 weights, a credible level of 0.05 a width of uniform =2 (to prevent extreme results) and for beta a multiplier of 20 (to prevent errors in optimization).
# saveRDS(results_uniform, "results_uniform_1_15_200_5000_005_wu_2_wb_20.rds")
#results_uniform <- readRDS("Simulation_results/results_uniform_1_15_200_5000_005_wu_2_wb_20.rds")
parameter_names <- rownames(results_uniform[[1]]$pc$credible_intervals$Gaussian_median)

# DISTANCES true parameter to MAP
plt_distances_to_MAP(results = results_uniform, prior_types = prior_types, path = "Simulation_images/Distances_to_map_uniform.png")
mean_distance_and_std_dev <- mean_distance_to_MAP_and_std_dev_of_Gaussian_approximation(results = results_uniform, prior_types = prior_types)

# CI: Credible intervals
plt_CI_lengths_and_get_mean_lengths(results = results_uniform, prior_types = prior_types, approximation_types = approximation_types, parameter_names = parameter_names, path = "Simulation_images/CI_lengths_uniform.png")
plt_frequency_true_parameter_in_CI(results = results_uniform, prior_types = prior_types, approximation_types = approximation_types, parameter_names = parameter_names, path = "Simulation_images/within_CI_uniform.png")


# KL divergences
KL_approx_types <- list(importance = "importance", smoothed_importance = "smoothed_importance")
plt_KL_and_get_mean_KL(results = results_uniform, prior_types = prior_types, approximation_types = KL_approx_types, path = "Simulation_images/KL_uniform.png")

# PROBABILITIES that marginal posterior is smaller than true parameter
plt_probabilities(results = results_uniform, prior_types = prior_types, approximation_types = approximation_types, parameter_names = parameter_names, path = "Simulation_images/probabilities_uniform.png")

# KS TEST FOR EACH PARAMETER

plt_KS(results = results_uniform, prior_types = prior_types, approximation_types = approximation_types, parameter_names = parameter_names, path1 = "Simulation_images/KS_distance_uniform.png", path2 = "Simulation_images/KS_pvalue_uniform.png")

# ks.test(x<-runif(5000),"punif")

# COMPLEXITY
plt_complexity_and_get_mean_complexity(results = results_uniform, prior_types = prior_types, approximation_types = approximation_types, path = "Simulation_images/complexity_uniform.png")
plt_complexity_and_get_mean_complexity(results = results_uniform, prior_types = prior_types[1:2], approximation_types = approximation_types, path = "Simulation_images/Complexity_12_uniform.png")


# K diagnostics and checking weights are similar
plt_k_diagnostics(results = results_uniform, prior_types = prior_types, path = "Simulation_images/k_diagnostics_uniform.png")
plt_weights_cdf(results = results_uniform, prior_types = prior_types, path = "Simulation_images/weights_uniform.png")

