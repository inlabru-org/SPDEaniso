#' @importFrom Matrix t solve

#' @title fdist
#' @description Calculates the component of the distance from the flexible to the base model depending on |v|
#'
#' @param r A positive number representing the norm of the anisotropy vector
#' @return The calculated component of the distance from the flexible to the base model depending on |v|
#' @export
#' @examples
#' r <- 0.5
#' fdist(r)
fdist <- function(r) {
  return(sqrt((1 / (48 * pi)) * (3 * cosh(2 * r) + 1)))
}


#' @title fdist_prime
#' @description Calculates the derivative of the component of the distance from the flexible to the base model depending on |v|
#'
#' @param r A positive number representing the norm of the anisotropy vector
#' @return The calculated derivative of the component of the distance from the flexible to the base model depending on |v|
#' @export
#' @examples
#' r <- 0.5
#' fdist_prime(r)
fdist_prime <- function(r) {
  return(sinh(2 * r) / (4 * sqrt(pi * cosh(2 * r) + pi / 3)))
}


#' @title Marginal PC prior density of log(kappa)
#' @description Calculates  the marginal PC prior pi_{log(kappa)} of the inverse correlation range kappa
#'
#' @param lambda A hyperparameter controlling the penalization of the distance from the base model as a function of kappa
#' @param lambda1 A hyperparameter controlling the penalization of the distance from the base model as a function of v
#' @param log_kappa Logarithm of the inverse correlation range kappa.
#'
#' @return The calculated marginal prior density of log(kappa) .
#' @export
#' @examples
#' log_kappa <- -0.3
#' lambda <- 0.5
#' lambda1 <- 1
#' result <- PC_prior_log_kappa(log_kappa, lambda, lambda1)
PC_prior_log_kappa <- function(log_kappa, lambda, lambda1) {
  kappa <- exp(log_kappa)
  c <- 2 * sqrt(3 * pi)
  fact <- exp(-((kappa * lambda) / c)) * lambda * lambda1 / (kappa * lambda + lambda1)^2
  change_variable <- kappa
  return(change_variable * fact * (1 + (kappa * lambda + lambda1) / c))
}

#' @title Marginal PC prior density of r = |v|
#' @description Calculates  the marginal PC prior, pi_r, of the norm of the anisotropy vector
#'
#' @param lambda1 A hyperparameter controlling the penalization of the distance from the base model as a function of r.
#' @param log_r The logarithm of |v|, which controls the magnitude of the anisotropy
#'
#' @return The calculated marginal prior density of log(r).
#' @export
#' @examples
#' log_r <- -0.5
#' lambda1 <- 1
#' result <- PC_prior_r(log_r, lambda1)
PC_prior_r <- function(log_r, lambda1) {
  r <- exp(log_r)
  numerator <- sqrt(3) * exp(-((lambda1 * (-2 + sqrt(1 + 3 * cosh(2 * r)))) / (4 * sqrt(3 * pi)))) * lambda1 * sinh(2 * r)
  denominator <- 4 * sqrt(pi + 3 * pi * cosh(2 * r))
  change_variable <- r
  result <- change_variable * numerator / denominator
  return(result)
}

#' @title Marginal PC prior density of v
#' @description Calculates  the marginal PC prior, pi_v, of the anisotropy vector
#'
#' @param lambda1 A hyperparameter controlling the penalization of the distance from the base model as a function of v
#' @param v A two dimensional vector that controls the direction and magnitude of anisotropy
#'
#' @return The calculated marginal prior density of v.
#' @export
#' @examples
#' v <- c(1, 2)
#' lambda1 <- 1
#' result <- PC_prior_v(v, lambda1)
PC_prior_v <- function(v, lambda1) {
  log_norm_v <- 0.5 * log(sum(v^2))
  return(PC_prior_r(log_norm_v, lambda1) / (2 * pi * exp(log_norm_v)))
  return(result)
}




#' @title PC Prior on log kappa, v
#' @description Calculates  the PC prior for (log(kappa), v) based on given hyperparameters.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param log_kappa Logarithm of inverse correlation range.
#' @param v Vector that controls anisotropy.
#'
#' @return The calculated PC prior density for (log(kappa),v).
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' log_kappa <- 0.5
#' v <- c(1, 2)
#' pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
pc_prior_aniso <- function(lambda, lambda1, log_kappa, v) {
  kappa <- exp(log_kappa)
  change_variable <- kappa

  v_norm <- sqrt(sum(v^2))
  f_val <- fdist(v_norm)
  f_prime_val <- fdist_prime(v_norm)

  term1 <- (lambda * lambda1 * abs(f_prime_val) * f_val) / (2 * pi * v_norm)
  term2 <- exp(-lambda1 * (f_val - fdist(0)) - lambda * f_val * kappa)

  return(term1 * term2 * change_variable)
}

#' @title Log PC Prior Calculation for anisotropy parameters.
#' @description Calculates  the log of the PC prior for (kappa, v) based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param log_kappa Logarithm of inverse correlation range.
#' @param v Vector that controls anisotropy.
#'
#' @return The calculated log PC prior value of (log(kappa),v).
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' log_kappa <- 0.5
#' v <- c(1, 2)
#' log_pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
log_pc_prior_aniso <- function(lambda, lambda1, log_kappa, v) {
  kappa <- exp(log_kappa)
  change_variable <- log_kappa
  v_norm <- sqrt(sum(v^2))

  f_val <- fdist(v_norm)
  f_prime_val <- fdist_prime(v_norm)

  term1 <- log(lambda) + log(lambda1) + log(f_prime_val) + log(f_val) - log(2 * pi * v_norm)
  term2 <- -lambda1 * (f_val - fdist(0)) - lambda * f_val * kappa

  return(term1 + term2 + change_variable)
}

#' @title Log PC Prior calculation for log variance of noise (and u)
#' @description Calculates  the log of the PC prior for sigma_epsilon based on given hyperparameters and vectors.
#'
#' @param lambda_epsilon A hyperparameter controlling the size of epsilon.
#' @param log_sigma_epsilon The logarithm of the variance of the additive noise.
#'
#' @return The calculated log PC prior value of log(sigma_epsilon).
#' @export
#' @examples
#' lambda_epsilon <- 1
#' log_sigma_epsilon <- 0
#' log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
log_pc_prior_noise_variance <- function(lambda_epsilon, log_sigma_epsilon) {
  sigma_epsilon <- exp(log_sigma_epsilon)
  change_variable <- log_sigma_epsilon
  # Calculates  the logarithm of exponential density.
  return(log(lambda_epsilon) - lambda_epsilon * sigma_epsilon + change_variable)
}

#' @title Log PC Prior calculation for theta=(log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @description Calculates  the log of the PC prior for theta=(log(kappa), v, log(sigma_u), log(sigma_epsilon)) based on given hyperparameters and vectors.
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of epsilon.
#' @param lambda_u A hyperparameter controlling the size of sigma_u.
#' @param log_kappa Logarithm of inverse correlation range.
#' @param v Vector that controls anisotropy.
#' @param log_sigma_u Logarithm of the variance of the field.
#' @param log_sigma_epsilon The logarithm of the variance of the additive noise.
#'
#' @return The calculated log PC prior value of theta=(log(kappa), v, log(sigma_u), log(sigma_epsilon)).
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' lambda_epsilon <- 1
#' lambda_u <- 1
#' log_kappa <- 0.5
#' v <- c(1, 2)
#' log_sigma_u <- 0
#' log_sigma_epsilon <- 0
#' log_pc_prior_theta(lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u, log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon)
log_pc_prior_theta <- function(lambda, lambda1, lambda_epsilon, lambda_u, log_kappa, v, log_sigma_u, log_sigma_epsilon) {
  # Calculates  the log PC prior for (log(kappa), v)
  log_pc_prior_aniso_value <- log_pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)

  # Calculates  the log PC prior for log(sigma_u)
  log_pc_prior_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)

  # Calculates  the log PC prior for log(sigma_epsilon)
  log_pc_prior_sigma_epsilon_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)

  # Calculates  the log PC prior for theta
  log_pc_value <- log_pc_prior_aniso_value + log_pc_prior_sigma_u_value + log_pc_prior_sigma_epsilon_value
  return(log_pc_value)
}

#' @title Sparse Matrix Determinant using Cholesky Factorization
#' @description Calculates the determinant of a sparse matrix using Cholesky factorization.
#'
#' @param Q A sparse matrix.
#'
#' @return The determinant of the sparse matrix.
#' @export
#' @examples
#' library(Matrix)
#' Q <- Matrix(c(4, -1, -1, 4), nrow = 2, sparse = TRUE)
#' sparse_log_determinant_chol(Q)
sparse_log_determinant_chol <- function(Q) {
  # Perform Cholesky factorization
  chol_fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)

  # Extract square diagonal elements of the Cholesky factor L
  diag_L <- Matrix::diag(chol_fact)

  # Calculates  determinant using Cholesky factor
  det_val <- sum(log(diag_L))

  return(det_val)
}



#' @title Calculates  the norm ||x||_Q:= <Qx,x>
#'
#' @description
#' Calculates  the term <Qx, x> appearing in Gaussian density
#'
#' @param Q A sparse matrix representing the precision matrix
#' @param x A numeric vector representing the point x
#'
#' @return The calculated norm term
#' @export
norm_Q <- function(Q, x) {
  norm <- sum(x * (Q %*% x))
  return(norm)
}

#' @title Calculates  the log Gaussian density
#'
#' @description
#' Calculates  the log Gaussian density using the precision matrix Q and mean mu
#'
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#' @param Q A sparse matrix representing the precision matrix
#'
#' @return The calculated log Gaussian density
#' @export
logGdensity <- function(x, mu, Q) {
  # Calculates  determinant using Cholesky factorization
  log_det_val <- sparse_log_determinant_chol(Q)

  # Calculates  the norm term
  norm_term <- norm_Q(Q, x - mu)

  # Calculates  the log Gaussian density
  0.5 * (log_det_val - norm_term - nrow(Q) * log(2 * pi))
}

#' @title Calculates  the log-posterior density for parameters (log_kappa,v, log_sigma_u, log_sigma_epsilon) with PC prior.
#'
#' @description
#' Calculates  the log-posterior density of parameters (log(kappa),v, log(sigma_u), log(epsilon)))
#' given a linear noisy observation y= A*u + epsilon
#' Uses based on the prior density and the likelihood.
#' Only stationary parameters are accepted.
#' Value is up to an additive constant depending only on y
#'
#' @param mesh The mesh
#' @param log_kappa Logarithm of inverse correlation range
#' @param v 2D vector that controls anisotropy
#' @param log_sigma_u Variance of field u. If unspecified, it is assumed to be 0.
#' @param log_sigma_epsilon Variance of noise
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param lambda_u A hyperparameter controlling the size of sigma_u.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given, it is transformed into (m_u, m_u,..., m_u)
#'
#' @return The calculated log-posterior density
#' @export


log_pc_posterior <- function(mesh, log_kappa, v, log_sigma_u = 0, log_sigma_epsilon, lambda, lambda1, lambda_epsilon, lambda_u = 1, y, A, m_u) {
  kappa <- exp(log_kappa)
  sigma_epsilon <- exp(log_sigma_epsilon)

  # Calculates log-prior density
  log_pc_value <- log_pc_prior_theta(lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u, log_kappa = log_kappa, v = v, log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon)

  # Calculates anisotropy
  n <- nrow(mesh$loc)
  kappa_values <- rep(kappa, n)
  vec_values <- matrix(v, n, 2, byrow = TRUE)
  aniso <- list(kappa = kappa_values, vec = vec_values)

  # Calculates log-density of the distribution of u at m_u knowing (kappa, v)
  Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
  if (length(m_u) == 1) {
    m_u <- rep(m_u, n)
  }
  u <- m_u
  logGdty_prior <- logGdensity(x = u, mu = m_u, Q = Q_u)

  # Calculates Q_epsilon,  Q_{u|y,theta} and m_{u|y,theta}
  Q_epsilon <- Matrix::Diagonal(n, 1 / sigma_epsilon^2)
  Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
  m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

  # Calculates  log-density of the posterior distribution of u given y and theta
  logGdty_posterior <- logGdensity(x = u, mu = m_uy_theta, Q = Q_uy_theta)

  # Calculates  log-density of the observation of y given u, theta
  logGdty_observation <- logGdensity(x = y, mu = A %*% u, Q = Q_epsilon)

  # Calculates  log-posterior density
  log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}

#' @title Calculates the MAP estimate for linear noisy observation of the field using PC priors on all parameters.
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters are accepted.
#'
#' @param mesh The mesh
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#' @param max_iterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL, it is estimated by the MAP
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export


MAP_pc <- function(mesh, lambda, lambda1, lambda_epsilon, lambda_u, y, A, m_u, max_iterations = 300, log_sigma_epsilon = NULL, theta0 = c(-0.5, c(0.1, 0.1), 0, -3)) {
  if (missing(log_sigma_epsilon)) {
    # Maximizes the log-posterior density over (log_kappa, v, log_sigma_u, log_sigma_epsilon)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      log_sigma_epsilon <- theta[5]
      return(log_pc_posterior(
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
        lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
        y = y, A = A, m_u = m_u
      ))
    }
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE, method = "BFGS"))

    # Maximizes the log-posterior density over (log_kappa, v, log_sigma_u)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      return(log_pc_posterior(
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
        lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
        y = y, A = A, m_u = m_u
      ))
    }
    theta0 <- theta0[1:4]
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE, method = "BFGS"))
  }
}

#' @title Log Gaussian density 1D
#' @description Calculates  the log of the Gaussian density with mean mu and variance sigma^2
#'
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#' @param log_sigma A numeric vector representing the variance log sigma
#'
#' @return The calculated log Gaussian density
#' @export
#' @examples
#' x <- 1
#' mu <- 0
#' log_sigma <- 1
#' log_gaussian_density(x, mu, log_sigma)
log_gaussian_density <- function(x, mu, log_sigma) {
  sigma <- exp(log_sigma)
  -0.5 * log(2 * pi) - log_sigma - 0.5 * (x - mu)^2 / sigma^2
}

#' @title Calculates  the log-posterior density for parameters (log_kappa,v, log_sigma_u, log_sigma_epsilon)
#' with a general prior on the anisotropy and PC priors on log(sigma_u) and log(sigma_epsilon).
#'
#' @description
#' Calculates  the log-posterior density of parameters (log(kappa),v, log(sigma_u), log(epsilon)))
#' given a linear noisy observation y= A*u + epsilon
#' Uses based on the prior density and the likelihood.
#' Only stationary parameters are accepted.
#' Value is up to an additive constant depending only on y
#'
#' @param log_prior_aniso A function that calculates the log prior of (kappa, v)
#' @param mesh The mesh
#' @param log_kappa Logarithm of inverse correlation range
#' @param v 2D vector that controls anisotropy
#' @param log_sigma_u Variance of field u. If unspecified, it is assumed to be 0.
#' @param log_sigma_epsilon Variance of noise
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param lambda_u A hyperparameter controlling the size of sigma_u.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given, it is transformed into (m_u, m_u,..., m_u)
#'
#' @return The calculated log-posterior density
#' @export
log_posterior_prior_on_aniso <- function(log_prior_aniso, mesh, log_kappa, v, log_sigma_u = 0, log_sigma_epsilon, lambda_epsilon, lambda_u = 1, y, A, m_u) {
  kappa <- exp(log_kappa)
  sigma_epsilon <- exp(log_sigma_epsilon)

  # Calculates log-prior
  log_prior_aniso_value <- log_prior_aniso(log_kappa = log_kappa, v = v)
  log_pc_noise_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
  log_pc_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
  log_pc_value <- log_prior_aniso_value + log_pc_noise_value + log_pc_sigma_u_value

  # Calculates anisotropy
  n <- nrow(mesh$loc)
  kappa_values <- rep(kappa, n)
  vec_values <- matrix(v, n, 2, byrow = TRUE)
  aniso <- list(kappa = kappa_values, vec = vec_values)

  # Calculates log-density of the distribution of u at m_u knowing (kappa, v)
  Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
  if (length(m_u) == 1) {
    m_u <- rep(m_u, n)
  }
  u <- m_u
  logGdty_prior <- logGdensity(x = u, mu = m_u, Q = Q_u)

  # Calculates Q_epsilon,  Q_{u|y,theta} and m_{u|y,theta}
  Q_epsilon <- Matrix::Diagonal(n, 1 / sigma_epsilon^2)
  Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
  m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

  # Calculates  log-density of the posterior distribution of u given y and theta
  logGdty_posterior <- logGdensity(x = u, mu = m_uy_theta, Q = Q_uy_theta)

  # Calculates  log-density of the observation of y given u, theta
  logGdty_observation <- logGdensity(x = y, mu = A %*% u, Q = Q_epsilon)

  # Calculates  log-posterior density
  log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}

#' @title Calculates  the log-posterior density for parameters (log_kappa,v, log_sigma_u, log_sigma_epsilon) with a general prior.
#'
#' @description
#' Calculates  the log-posterior density of parameters (log(kappa),v, log(sigma_u), log(epsilon)))
#' given a linear noisy observation y= A*u + epsilon
#' Uses based on the prior density and the likelihood.
#' Only stationary parameters are accepted.
#' Value is up to an additive constant depending only on y
#'
#' @param log_prior A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param mesh The mesh
#' @param log_kappa Logarithm of inverse correlation range
#' @param v 2D vector that controls anisotropy
#' @param log_sigma_u Variance of field u. If unspecified, it is assumed to be 0.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given, it is transformed into (m_u, m_u,..., m_u)

#' @return The calculated log-posterior
#' @export
log_posterior_prior <- function(log_prior, mesh, log_kappa, v, log_sigma_u = 0, log_sigma_epsilon, y, A, m_u) {
  kappa <- exp(log_kappa)
  sigma_epsilon <- exp(log_sigma_epsilon)

  # Calculates log-prior density
  log_prior_value <- log_prior(log_kappa, v, log_sigma_u, log_sigma_epsilon)

  # Calculates anisotropy
  n <- nrow(mesh$loc)
  if (length(y) > n) {
    stop("y cannot be a vector with length larger than the size of the mesh")
  }
  if (length(y) != nrow(A)) {
    stop("y and A must have the same number of rows")
  }
  kappa_values <- rep(kappa, n)
  vec_values <- matrix(v, n, 2, byrow = TRUE)
  aniso <- list(kappa = kappa_values, vec = vec_values)

  # Calculates log-density of the distribution of u at m_u knowing (kappa, v)
  Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
  if (length(m_u) == 1) {
    m_u <- rep(m_u, n)
  }
  u <- m_u
  logGdty_prior <- logGdensity(x = u, mu = m_u, Q = Q_u)

  # Calculates Q_epsilon,  Q_{u|y,theta} and m_{u|y,theta}
  Q_epsilon <- Matrix::Diagonal(length(y), 1 / sigma_epsilon^2)
  Q_uy_theta <- Q_u + t(A) %*% Q_epsilon %*% A
  m_uy_theta <- solve(Q_uy_theta, Q_u %*% m_u + t(A) %*% Q_epsilon %*% y)

  # Calculates  log-density of the posterior distribution of u given y and theta
  logGdty_posterior <- logGdensity(x = u, mu = m_uy_theta, Q = Q_uy_theta)

  # Calculates  log-density of the observation of y given u, theta
  logGdty_observation <- logGdensity(x = y, mu = A %*% u, Q = Q_epsilon)

  # Calculates  log-posterior density
  unname(log_prior_value + logGdty_prior + logGdty_observation - logGdty_posterior)
}
#' @title Calculates the MAP estimate for linear noisy observation of field with a general prior on theta = (log(kappa),v, log(sigma_u), log(sigma_epsilon)).
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters are accepted.
#'
#' @param log_prior A function that calculates the log prior of (log(kappa), v, log(sigma_u), log(sigma_epsilon)).
#' If not specified, it is assumed to be the function that returns 0.
#' @param mesh The mesh
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#' @param max_iterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL, it is estimated by the MAP
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (log(0.5), 1, 2, 1, 1)
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export
MAP_prior <- function(log_prior = function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
                        return(0)
                      }, mesh, y, A, m_u, log_sigma_epsilon = NULL, max_iterations = 300, theta0 = c(-0.5, c(0.1, 0.1), 0, -3)) {
  if (missing(log_sigma_epsilon) || is.null(log_sigma_epsilon)) {
    # Maximizes the log-posterior density over (log_kappa, v, log_sigma_u, log_sigma_epsilon)
    # First uses Nelder-Mead to find a good starting point for the optimization
    # and then uses BFGS to maximize the log-posterior density
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      log_sigma_epsilon <- theta[5]
      tryCatch(
        {
          log_posterior_prior(
            log_prior = log_prior,
            mesh = mesh, log_kappa = log_kappa, v = v,
            log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
            y = y, A = A, m_u = m_u
          )
        },
        error = function(e) {
          -Inf
        }
      )
    }
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE, method = "BFGS"))
  } else {
    # Maximizes the log-posterior density over (log_kappa, v, log_sigma_u)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      tryCatch(
        {
          log_posterior_prior(
            log_prior = log_prior,
            mesh = mesh, log_kappa = log_kappa, v = v,
            log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
            y = y, A = A, m_u = m_u
          )
        },
        error = function(e) {
          -Inf
        }
      )
    }
    theta0 <- theta0[1:4]
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = TRUE, method = "BFGS"))
  }
}

#' @title Simulation of anisotropy parameters (log(kappa), v) from the PC prior
#'
#' @description Simulates anisotropy parameters (log(kappa), v) from the PC prior
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param m Number of samples, by default 1.
#'
#' @return A list with two elements: log(kappa) and v
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' m <- 1
#' result <- sim_aniso_pc(lambda = lambda, lambda1 = lambda1, m = m)
sim_aniso_pc <- function(lambda, lambda1, m = 1) {
  # Initialize a list to store results
  results <- vector("list", m)

  # Calculate the CDF of Rayleigh distribution for a given value
  R <- function(x) 1 - exp(-x^2 / 2)

  # Inverse function f^{-1}
  f_inv <- function(x) 0.5 * acosh(16 * pi * x^2 - 1 / 3)

  for (i in 1:m) {
    # Generate Y vector from standard normal distribution
    Y <- rnorm(3)

    # Calculations for v1, v2, and kappa
    radius <- sqrt(Y[1]^2 + Y[2]^2)
    common_term <- fdist(0) - (log(1 - R(radius)) / lambda1)

    v1 <- f_inv(common_term) * Y[1] / radius
    v2 <- f_inv(common_term) * Y[2] / radius
    kappa <- -(common_term)^-1 * log(1 - pnorm(Y[3])) / lambda
    log_kappa <- log(kappa)
    v <- c(v1, v2)

    # Store the results as a vector in the list
    results[[i]] <- list(log_kappa = log_kappa, v = v)
  }
  if (m == 1) {
    return(results[[1]])
  }
  return(results)
}

#' @title Simulation of theta under PC priors
#' @description Simulates theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon)) from the PC prior
#'
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param lambda_u A hyperparameter controlling the size of sigma_u.
#' @param m Number of samples, by default 1.
#'
#' @return A list with four elements: log_kappa, v, log_sigma_u, log_sigma_epsilon
#' @export
#' @examples
#' lambda <- 1
#' lambda1 <- 1
#' lambda_epsilon <- 1
#' lambda_u <- 1
#' m <- 10
#' result <- sim_theta_pc(lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u, m = m)
sim_theta_pc <- function(lambda, lambda1, lambda_u, lambda_epsilon, m = 1) {
  # Initialize a list to store results
  results <- vector("list", m)

  for (i in 1:m) {
    # Simulate kappa and v
    kappa_v <- sim_aniso_pc(lambda = lambda, lambda1 = lambda1, m = 1)

    # Simulate sigma_u and sigma_epsilon from exponential distributions
    sigma_u <- rexp(1, rate = lambda_u)
    sigma_epsilon <- rexp(1, rate = lambda_epsilon)

    # Store all parameters
    results[[i]] <- list(log_kappa = kappa_v$log_kappa, v = kappa_v$v, log_sigma_u = log(sigma_u), log_sigma_epsilon = log(sigma_epsilon))
  }
  if (m == 1) {
    return(results[[1]])
  }
  return(results)
}

#' @title Normalized importance weights
#' @description Given log unnormalized weights, calculates the normalized importance weights
#' @param log_unnormalized_weights A vector of log unnormalized weights
#' @return A vector of normalized importance weights
#' @export
#' @examples
#' log_unnormalized_weights <- log(c(1, 2, 3))
#' normalize_log_weights(log_unnormalized_weights)
normalize_log_weights <- function(log_unnormalized_weights) {
  max_log_unnormalized_weight <- max(log_unnormalized_weights)
  # Subtract the max to avoid numerical issues as it shouldn't change the result
  weights <- exp(log_unnormalized_weights - max_log_unnormalized_weight)
  as.vector(weights / sum(weights))
}



#' @title Calculate confidence intervals for each entry of a weighted sample
#' @description Calculates the confidence intervals for each entry of a weighted sample.
#' @param theta A matrix of size n x d representing the weighted sample
#' @param w A vector of length n representing the weights of the unnormalized weighted sample.
#' @param q The significance level
#' @return A matrix of size d x 2 representing the confidence intervals for each entry of the weighted sample
#' @export
calculate_confidence_intervals_importance <- function(theta, w, q) {
  # Normalize the weights
  w <- w / sum(w)

  # Function to calculate the confidence interval for one dimension
  confidence_interval <- function(theta, w, q) {
    df <- data.frame(theta = theta, w = w)

    # Calculate cumulative weights
    df <- df %>%
      arrange(theta) %>% # Sort by theta
      mutate(cum_w = cumsum(w)) # Calculate cumulative weights and add to the data frame as a new column called cum_w

    # Find lower and upper bounds
    lower_bound <- df$theta[which(df$cum_w >= q / 2)[1]]
    upper_bound <- df$theta[which(df$cum_w >= 1 - q / 2)[1]]

    return(c(lower_bound, upper_bound))
  }

  # Apply the function to each dimension of theta
  t(apply(theta, 2, confidence_interval, w, q))
}
#' @title Calculate the probabilities P[theta_i <= theta^*_i]
#' @description Calculates the probabilities P[theta_i <= theta^*_i] using the weighted empirical CDF.
#' @param theta_0 A vector representing the values theta^*_i.
#' @param theta A vector representing the sample.
#' @param w A vector of the same length as theta representing the weights.
#' @return A vector of probabilities P[theta_i <= theta^*_i].
calculate_probabilities <- function(theta_0, theta, w) {
  w <- w / sum(w)
  # Function to calculate the probabilities for one dimension
  probability <- function(theta, w, theta_0) {
    df <- data.frame(theta = theta, w = w)

    # Calculate cumulative weights
    df <- df %>%
      arrange(theta) %>% # Sort by theta
      mutate(cum_w = cumsum(w)) # Calculate cumulative weights and add to the data frame as a new column called cum_w

    # Finds probabilities
    probabilities <- df$cum_w[which(df$theta <= theta_0)]
    return(probabilities[length(probabilities)])
  }

  # Apply the function to each dimension of theta
  sapply(1:length(theta_0), function(i) probability(theta[, i], w, theta_0[i]))
}
#' @title Calculate confidence intervals for Gaussian
#' @description Calculates the confidence intervals for a Gaussian distribution.
#' @param mu A vector of length d representing the mean of the Gaussian distribution
#' @param standard_deviation A vector of length d representing the standard deviation of the Gaussian distribution
#' @param q The significance level
#' @return A matrix of size d x 2 representing the confidence intervals for each entry of the Gaussian distribution
#' @export
#' @examples
#' mu <- c(0, 0)
#' standard_deviation <- c(1, 1)
#' q <- 0.05
#' calculate_confidence_intervals_gaussian(mu, standard_deviation, q)
calculate_confidence_intervals_gaussian <- function(mu, standard_deviation, q) {
  # Calculate the confidence intervals for each dimension
  ci <- cbind(mu - qnorm(1 - q / 2) * standard_deviation, mu + qnorm(1 - q / 2) * standard_deviation)
  return(ci)
}

#' @title Parameter within confidence intervals
#' @description Checks if the parameter is within the confidence intervals
#' @param parameter A vector of length d representing the parameter
#' @param confidence_intervals A matrix of size d x 2 representing the confidence intervals
#' @return A vector of length d representing whether the parameter is within the confidence intervals
#' @export
#' @examples
#' mu <- c(0, 0)
#' standard_deviation <- c(1, 1)
#' q <- 0.05
#' confidence_intervals <- calculate_confidence_intervals_gaussian(mu, standard_deviation, q)
#' theta <- c(0.5, 0.5)
#' parameter_within_confidence_intervals(theta, confidence_intervals)
#'
parameter_within_confidence_intervals <- function(parameter, confidence_intervals) {
  # Check if the parameter is within the confidence intervals
  parameter <- unlist(parameter)
  parameter_within_confidence_intervals <- (parameter >= confidence_intervals[, 1]) & (parameter <= confidence_intervals[, 2])
  return(parameter_within_confidence_intervals)
}


#' @title Calculate unnormalized log importance weights, mean, variance, and KL divergences
#' @description Calculates the unnormalized log importance weights, mean, variance, and KL divergences between the importance approximation and the Laplace approximation
#' @param log_posterior_density A function that calculates the unnormalized log posterior of theta
#' @param mu_Laplace A 5D vector representing the mean of the Laplace approximation
#' @param Q_Laplace A 5x5 matrix representing the precision matrix of the Laplace approximation
#' @param n_weights Number of weights to be calculated
#' @param q The significance level. By default, set to 0.05.
#'
#' @return A list with three elements: log_unnormalized_weights, log_unnormalized_weights_smoothed and the k-diagnostic for PSIS
#'
#' @export
#' @examples
#' # Defines the log prior density function of theta for PC and non-PC priors
#' log_pc_prior <- log_pc_prior_quantile(
#'   sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
#'   a0 = a0, rho0 = rho0, alpha = alpha
#' )
#' # Defining the log posterior density as a function of the parameters using the log_pc_prior
#' m_u <- 0
#' log_posterior_pc <- function(theta) {
#'   log_kappa <- theta[1]
#'   v <- theta[2:3]
#'   log_sigma_u <- theta[4]
#'   log_sigma_epsilon <- theta[5]
#'   log_posterior_prior(
#'     log_prior = log_pc_prior,
#'     mesh = mesh, log_kappa = log_kappa, v = v,
#'     log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
#'     y = y, A = A, m_u = m_u
#'   )
#' }
#' # Calculates the importance weights, mean, variance, and KL divergences
#' mu_Laplace <- c(0, 0, 0, 0, 0)
#' Q_Laplace <- diag(5)
#' log_unnormalized_importance_weights_and_integrals(log_posterior_pc, mu_Laplace, Q_Laplace, n_weights = 1000)
log_unnormalized_importance_weights_and_integrals <- function(log_posterior_density, mu_Laplace, Q_Laplace, n_weights, q = 0.05, true_params) {
  # Calculate the importance weights
  log_Laplace_density <- function(theta) {
    logGdensity(
      x = theta, mu = mu_Laplace, Q = Q_Laplace
    )
  }

  log_ratio_function <- function(theta) {
    log_posterior_density(theta) - log_Laplace_density(theta)
  }
  covariance_Laplace <- solve(Q_Laplace)
  theta_sim_importance <- MASS::mvrnorm(n_weights, mu_Laplace, covariance_Laplace)
  # Subtract the max to avoid numerical issues as it shouldn't change the result
  log_importance_ratios <- apply(theta_sim_importance, 1, log_ratio_function)
  log_importance_ratios <- log_importance_ratios - max(log_importance_ratios)

  psis_result <- psis(log_importance_ratios, r_eff = NA)
  log_weights_smoothed <- psis_result$log_weights

  # Calculate the mean, variance and KL divergences with unsmoothed weights
  weights_normalized <- normalize_log_weights(log_importance_ratios)
  n_eff <- 1 / sum(weights_normalized^2)
  mean_importance <- apply(theta_sim_importance, 2, weighted.mean, w = weights_normalized)
  cov_importance <- stats::cov.wt(theta_sim_importance, wt = weights_normalized)$cov
  marginal_variance_importance <- diag(cov_importance)

  # Calculate mean and variance using smoothed weights
  weights_smoothed_normalized <- normalize_log_weights(log_weights_smoothed)
  mean_smoothed_importance <- apply(theta_sim_importance, 2, weighted.mean, w = weights_smoothed_normalized)
  cov_smoothed_importance <- stats::cov.wt(theta_sim_importance, wt = weights_smoothed_normalized)$cov
  marginal_variance_smoothed_importance <- diag(cov_smoothed_importance)

  # Calculate a 95% confidence interval for each component of the parameters around its mean
  std_dev_Laplace <- sqrt(diag(covariance_Laplace))
  confidence_intervals_Laplace <- calculate_confidence_intervals_gaussian(mu_Laplace, std_dev_Laplace, q)
  confidence_intervals_importance <- calculate_confidence_intervals_importance(theta_sim_importance, weights_normalized, q)
  confidence_intervals_importance_smoothed <- calculate_confidence_intervals_importance(theta_sim_importance, weights_smoothed_normalized, q)
  # Calculates the probabilities P[theta_i <= true_params_i]
  probabilities_Laplace <- pnorm(true_params, mu_Laplace, std_dev_Laplace)
  probabilities_importance <- calculate_probabilities(true_params, theta_sim_importance, weights_normalized)
  probabilities_importance_smoothed <- calculate_probabilities(true_params, theta_sim_importance, weights_smoothed_normalized)

  # A discrete distribution (importance approximation to posterior) is not absolutely continuous with respect to a continuous one (Laplace approximation to posterior) so the KL divergence is not defined. As a result, we replace the Laplace approximation with its importance approximation to calculate the KL divergence.
  log_weights_Laplace <- apply(theta_sim_importance, 1, log_Laplace_density)
  weights_Laplace_normalized <- normalize_log_weights(log_weights_Laplace)
  KL_divergence_importance_Laplace <- sum(weights_normalized * log(weights_normalized)) + log(n_weights)
  KL_divergence_smoothed_importance_Laplace <- sum(weights_smoothed_normalized * log(weights_smoothed_normalized)) + log(n_weights)

  list(
    log_unnormalized_weights = log_importance_ratios,
    log_unnormalized_weights_smoothed = log_weights_smoothed,
    psis_result = psis_result,
    k_diagnostic = psis_result$diagnostics$pareto_k,
    n_eff = n_eff,
    n_eff_smoothed = psis_result$diagnostics$n_eff,
    mean_importance = mean_importance,
    marginal_variance_importance = marginal_variance_importance,
    mean_smoothed_importance = mean_smoothed_importance,
    marginal_variance_smoothed_importance = marginal_variance_smoothed_importance,
    confidence_intervals_Laplace = confidence_intervals_Laplace,
    confidence_intervals_importance = confidence_intervals_importance,
    confidence_intervals_importance_smoothed = confidence_intervals_importance_smoothed,
    probabilities_Laplace = probabilities_Laplace,
    probabilities_importance = probabilities_importance,
    probabilities_importance_smoothed = probabilities_importance_smoothed,
    KL_divergence_importance_Laplace = KL_divergence_importance_Laplace,
    KL_divergence_smoothed_importance_Laplace = KL_divergence_smoothed_importance_Laplace
  )
}

#' @title KL discrete log unnormalized weights
#' @description Calculates the KL divergence between two discrete distributions with the same support using the log unnormalized weights
#' @param log_unnormalized_weights_1 A vector of log unnormalized weights for the first distribution
#' @param log_unnormalized_weights_2 A vector of log unnormalized weights for the second distribution
#'
#' @return The KL divergence between the two distributions
#'
#' @export
#' @examples
#' log_unnormalized_weights_1 <- log(c(1, 2, 3))
#' log_unnormalized_weights_2 <- log(c(3, 2, 1))
#' KL_discrete_log_unnormalized_weights(log_unnormalized_weights_1, log_unnormalized_weights_2)
KL_discrete_log_unnormalized_weights <- function(log_unnormalized_weights_1, log_unnormalized_weights_2) {
  weights_1_normalized <- normalize_log_weights(log_unnormalized_weights_1)
  weights_2_normalized <- normalize_log_weights(log_unnormalized_weights_2)
  return(sum(weights_1_normalized * (log(weights_1_normalized) - log(weights_2_normalized))))
}

#' @title Convert seconds to hours, minutes, and seconds
#' @description Converts seconds to hours, minutes, and seconds
#' @param seconds A number representing the number of seconds
#' @return A string representing the number of hours, minutes, and seconds
#' @export
#' @examples
#' seconds <- 100
#' seconds_to_hours_minutes_seconds(seconds)
#' seconds <- 10000
#' seconds_to_hours_minutes_seconds(seconds)
seconds_to_hours_minutes_seconds <- function(seconds) {
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  seconds <- round(seconds %% 60)
  if (hours == 0 && minutes == 0) {
    return(paste0(seconds, "s"))
  }
  if (hours == 0) {
    return(paste0(minutes, "m ", seconds, "s"))
  }
  return(paste0(hours, "h ", minutes, "m ", seconds, "s"))
}
