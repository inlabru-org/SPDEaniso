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
  logGdty <- 0.5 * (log_det_val - norm_term - nrow(Q) * log(2 * pi))

  return(logGdty)
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
#' @param maxiterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL, it is estimated by the MAP
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export


MAP_pc <- function(mesh, lambda, lambda1, lambda_epsilon, lambda_u, y, A, m_u, maxiterations = 300, log_sigma_epsilon = NULL, theta0 = c(-0.5, c(0.1, 0.1), 0, -3)) {
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
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE, method = "BFGS"))

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
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE, method = "BFGS"))
  }
}

#' @title Log Gaussian density 1D
#' @description Calculates  the log of the Gaussian density with mean mu and variance sigma^2
#'
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#' @param logsigma A numeric vector representing the variance log sigma
#'
#' @return The calculated log Gaussian density
#' @export
#' @examples
#' x <- 1
#' mu <- 0
#' logsigma <- 1
#' log_gaussian_density(x, mu, logsigma)
log_gaussian_density <- function(x, mu, logsigma) {
  sigma <- exp(logsigma)
  return(-0.5 * log(2 * pi) - logsigma - 0.5 * (x - mu)^2 / sigma^2)
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
#' @param logprior_aniso A function that calculates the log prior of (kappa, v)
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
log_posterior_general <- function(logprior_aniso, mesh, log_kappa, v, log_sigma_u = 0, log_sigma_epsilon, lambda, lambda1, lambda_epsilon, lambda_u = 1, y, A, m_u) {
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
#' @param logprior A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
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
log_posterior_prior <- function(logprior, mesh, log_kappa, v, log_sigma_u = 0, log_sigma_epsilon, y, A, m_u) {
  kappa <- exp(log_kappa)
  sigma_epsilon <- exp(log_sigma_epsilon)

  # Calculates log-prior density
  log_prior_value <- logprior(log_kappa, v, log_sigma_u, log_sigma_epsilon)

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
  log_posterior_val <- log_prior_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}
#' @title Calculates the MAP estimate for linear noisy observation of field with a general prior on theta = (log(kappa),v, log(sigma_u), log(sigma_epsilon)).
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters are accepted.
#'
#' @param logprior A function that calculates the log prior of (log(kappa), v, log(sigma_u), log(sigma_epsilon)).
#' If not specified, it is assumed to be the function that returns 0.
#' @param mesh The mesh
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#' @param maxiterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL, it is estimated by the MAP
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (log(0.5), 1, 2, 1, 1)
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export
MAP_prior <- function(logprior = function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
                        return(0)
                      }, mesh, y, A, m_u, log_sigma_epsilon = NULL, maxiterations = 300, theta0 = c(-0.5, c(0.1, 0.1), 0, -3)) {
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
            logprior = logprior,
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
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE, method = "BFGS"))
  } else {
    # Maximizes the log-posterior density over (log_kappa, v, log_sigma_u)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      tryCatch(
        {
          log_posterior_prior(
            logprior = logprior,
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
    theta0 <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE)$par
    return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations / 2), hessian = TRUE, method = "BFGS"))
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
