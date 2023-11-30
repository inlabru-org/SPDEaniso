#' @importFrom Matrix t solve

#' @title fdist
#' @description Calculates the compenent of the distance from the flexible to the base model depending on |v|
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
#' @description Calculates the derivative of the compenent of the distance from the flexible to the base model depending on |v|
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
#' @param lambda A hyperparameter controlling the penalization of the distance from the base moedel as a function of kappa
#' @param lambda1 A hyperparameter controlling the penalization of the distance from the base moedel as a function of v
#' @param log_kappa Logarithm of the inverse correlation range kappa.
#'
#' @return The calculated marginal prior density of log(kappa) .
#' @export
#' @examples
#' log_kappa <- -0.3
#' lambda <- 0.5
#' lambda1 <- 1
#' result <- PC_prior_kappa(log_kappa, lambda, lambda1)
PC_prior_kappa <- function(log_kappa, lambda, lambda1) {
  kappa <- exp(log_kappa)
  c <- 2 * sqrt(3 * pi)
  fact <- exp(-((kappa * lambda) / c)) * lambda * lambda1 / (kappa * lambda + lambda1)^2
  change_variable <- kappa
  return(change_variable * fact * (1 + (kappa * lambda + lambda1) / c))
}

#' @title Marginal PC prior density of r = |v|
#' @description Calculates  the marginal PC prior, pi_r, of the norm of the anisotropy vector
#'
#' @param lambda1 A hyperparameter controlling the penalization of the distance from the base moedel as a function of r.
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
#' @param lambda1 A hyperparameter controlling the penalization of the distance from the base moedel as a function of v
#' @param v A two dimensional vector which controls the direction and magintude of anisotropy
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




#' @title PC Prior on kappa, v
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
#' pc_prior(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
pc_prior <- function(lambda, lambda1, log_kappa, v) {
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
  change_variable <- kappa
  v_norm <- sqrt(sum(v^2))
  
  f_val <- fdist(v_norm)
  f_prime_val <- fdist_prime(v_norm)

  term1 <- log(lambda) + log(lambda1) + log(f_prime_val) + log(f_val) - log(2 * pi * v_norm)
  term2 <- -lambda1 * (f_val - fdist(0)) - lambda * f_val * kappa

  return(term1 + term2 + change_variable)
}

#' @title Log PC Prior calcualtion for log variance of noise (and u)
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
  change_variable <- sigma_epsilon
  # Calculates  the logarithm of exponential density.
  return(log(lambda_epsilon) - lambda_epsilon * sigma_epsilon + sigma_epsilon)
}

#' @title Sparse Matrix Determinant using Cholesky Factorization
#' @description Calculates  the determinant of a sparse matrix using Cholesky factorization.
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
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given is transformed into (m_u, m_u,..., m_u)
#'
#' @return The calculated log-posterior
#' @export


log_posterior <- function(mesh, log_kappa, v, log_sigma_u = 0, log_sigma_epsilon, lambda, lambda1, lambda_epsilon, lambda_u = 1, y, A, m_u) {
  kappa <- exp(log_kappa)
  sigma_epsilon <- exp(log_sigma_epsilon)

  # Calculates log-prior
  log_pc_aniso_value <- log_pc_prior_aniso(lambda = lambda, lambda1 = lambda1, log_kappa = log_kappa, v = v)
  log_pc_noise_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
  log_pc_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
  log_pc_value <- log_pc_aniso_value + log_pc_noise_value + log_pc_sigma_u_value

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

  # Calculates  log-posterior
  log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}

#' @title Calculates  the MAP estimate for linear noisy observation of field.
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters accepted.
#'
#' @param mesh The mesh
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#' @param maxiterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL it is estimated by the MAP
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export


MAP <- function(mesh, lambda, lambda1, lambda_epsilon, lambda_u, y, A, m_u, maxiterations = 300, log_sigma_epsilon = NULL) {
  if (missing(log_sigma_epsilon)) {
    # Optimizes the log-posterior over (log_kappa, v, log_sigma_u, log_sigma_epsilon)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      log_sigma_epsilon <- theta[5]
      return(log_posterior(
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
        lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
        y = y, A = A, m_u = m_u
      ))
    }
    aniso0 <- c(log(0.5), c(1, 2), 1, 1)
    # To do: calculate the gradient of log posterior
    # gradient= grad_log_posterior(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u)
    return(optim(par = aniso0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations), hessian = TRUE))
  } else {
    # Optimizes the log-posterior over (log_kappa, v, log_sigma_u)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      return(log_posterior(
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
        lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
        y = y, A = A, m_u = m_u
      ))
    }
    aniso0 <- c(1, c(0.1, 0.1), 1)
    return(optim(par = aniso0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations), hessian = TRUE))
  }
}

#' @title Calculates  the gradient of the log posterior of a linear observation y = A u + noise
#'
#' @description
#' Calculates  the gradient of the log posterior of a linear observation y = A u + noise
#'
#' @param mesh The mesh
#' @param kappa Inverse correlation range
#' @param v Vector that controls anisotropy
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param Q_epsilon A sparse matrix of size nxn representing the noise precision matrix
#' @param m_u A vector with length n representing the prior mean m_u
#'
#' @return The calculated log-posterior
#' @export


grad_log_posterior <- function(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u) {
  # Calculation
  grad_log_prior <- 0
  # d|Q|= |Q| Tr(Q^{-1}dQ)
  # Need to differentiate Q = C_kappa + 2 G_v + G_v C_kp^{-1} G_v
  # Equivalent to differentiating C, and G component by component.
  # Should we add the calculation of gradient over each triangle to the calcCaniso, calcGaniso, calcQaniso functions
  # dC^{-1} = - dC * M * dC   (don't forget lumped mass)
  # Should I calDo we need to differentitate manifold case?
  grad_logGdty_prior <- 0
  grad_logGdty_posterior <- 0
  return(grad_log_prior + grad_logGdty_prior + grad_logGdty_posterior)
}

#' @title Calculates  the MAP estimate for linear noisy observation of field.
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters accepted.
#'
#' @param mesh The mesh
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#'
#' @return The calculated log-posterior
#' @export


MAPold <- function(mesh, lambda, lambda1, lambda_epsilon, y, A, m_u, maxiiterations = 300) {
  # Writes the log-posterior as a function of (log_kappa, v, sigma)
  log_post <- function(theta) {
    log_kappa <- theta[1]
    v <- theta[2:3]
    log_sigma_epsilon <- theta[4]
    return(log_posterior(mesh = mesh, log_kappa = log_kappa, v = v, log_sigma_epsilon = log_sigma_epsilon, lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, y = y, A = A, m_u = m_u))
  }
  aniso0 <- c(log(0.5), c(1, 2), 1)
  # To do: calculate the gradient of log posterior
  # gradient= grad_log_posterior(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u)
  return(optim(par = aniso0, fn = log_post, control = list(fnscale = -1, maxit = maxiiterations)))
}

#' @title Log gaussian 1D density
#' @description Calculates  the log of the gaussian density with mean mu and variance sigma^2
#'
#' @param x A numeric vector representing the point x
#' @param mu A numeric vector representing the mean mu
#' @param logsigma A numeric vector representing the variance logsigma
#'
#' @return The calculated log gaussian density
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
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given is transformed into (m_u, m_u,..., m_u)
#'
#' @return The calculated log-posterior
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

  # Calculates  log-posterior
  log_posterior_val <- log_pc_value + logGdty_prior + logGdty_observation - logGdty_posterior

  return(log_posterior_val)
}
#' @title Calculates  the MAP estimate for linear noisy observation of field with a general prior on (kappa,v).
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters accepted.
#'
#' @param logprior_aniso A function that calculates the log prior of (kappa, v)
#' @param mesh The mesh
#' @param lambda A hyperparameter controlling the size of kappa.
#' @param lambda1 A hyperparameter controlling the size of |v|.
#' @param lambda_epsilon A hyperparameter controlling the size of sigma_epsilon.
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param m_u A vector with length n representing the prior mean m_u
#' @param maxiterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL it is estimated by the MAP
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export
MAPgeneral <- function(logprior_aniso, mesh, lambda, lambda1, lambda_epsilon, lambda_u, y, A, m_u, maxiterations = 300, log_sigma_epsilon = NULL) {
  if (missing(log_sigma_epsilon)) {
    # Optimizes the log-posterior over (log_kappa, v, log_sigma_u, log_sigma_epsilon)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      log_sigma_epsilon <- theta[5]
      return(log_posterior_general(logprior_aniso = logprior_aniso, 
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
        lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
        y = y, A = A, m_u = m_u
      ))
    }
    aniso0 <- c(log(0.5), c(1, 2), 1, 1)
    # To do: calculate the gradient of log posterior
    # gradient= grad_log_posterior(mesh, kappa, v, lambda, lambda1, y, A, Q_epsilon, m_u)
    return(optim(par = aniso0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations), hessian = TRUE))
  } else {
    # Optimizes the log-posterior over (log_kappa, v, log_sigma_u)
    log_post <- function(theta) {
      log_kappa <- theta[1]
      v <- theta[2:3]
      log_sigma_u <- theta[4]
      return(log_posterior_general(logprior_aniso = logprior_aniso,
        mesh = mesh, log_kappa = log_kappa, v = v,
        log_sigma_epsilon = log_sigma_epsilon, log_sigma_u = log_sigma_u,
        lambda = lambda, lambda1 = lambda1, lambda_epsilon = lambda_epsilon, lambda_u = lambda_u,
        y = y, A = A, m_u = m_u
      ))
    }
    aniso0 <- c(1, c(0.1, 0.1), 1)
    return(optim(par = aniso0, fn = log_post, control = list(fnscale = -1, maxit = maxiterations), hessian = TRUE))
  }
}


#' @title Gaussian prior on anisotropy given quantiles
#' @description Calculates  the log of the gaussian prior on the (log(kappa),v) supposing log(|v|)~N(1,sigma_v^2) and log(kappa)~N(1,sigma_k^2)
#' such that the anisotropy ratio exp(2|v|) is larger than a0 with probability alpha and the correlation range sqrt{16}/kappa is smaller than rho0 with probability alpha
#'
#' @param alpha_k A quantile in (0,1) for kappa
#' @param alpha_v A quantile in (0,1) for |v|
#' @param rho0 A surprisingly small correlation range
#' @param a0 A surprisingly high ratio of anisotropy
#'
#' @return The calculated log of the gaussian prior on the (log(kappa),v)
#' @export
#' @examples
#' alpha_k <- 0.01
#' alpha_v <- 0.01
#' a0 <- 10
#' rho0 <- 1
#' log_prior_aniso_gaussian <- log_prior_aniso_quantile(alpha_k = alpha_k, beta_v = beta_v, rho0 = rho0, a0 = a0)
log_prior_aniso_quantile_gaussian <- function(alpha_k, alpha_v, rho0, a0) {
  #This warns the user if alpha is not in (0,1)
  if (alpha_k <= 0 | alpha_k >= 1) {
    warning("alpha should be in (0,1)")
  }
  #This warns the user if beta is not in (0,1)
  if (alpha_v <= 0 | alpha_v >= 1) {
    warning("beta should be in (0,1)")
  }
  #This warns the user if rho0 is not greater than 0
  if (rho0 > 0) {
    warning("rho0 should be greater than 0")
  }
  #This warns the user if a0 is not greater than 1
  if (a0 <= 1) {
    warning("a0 should be greater than 1")
  }
  sigma_v <- sigma_quantile_v(alpha_v = alpha_v, a0 = a0)
  sigma_k <- sigma_quantile_kappa(alpha_k = alpha_k, rho0 = rho0)
  log_prior_aniso <- function(log_kappa, v) {
    v_norm <- sqrt(sum(v^2))
    log_prior_aniso <- log_gaussian_density(x = log_kappa, mu = 1, logsigma = log(sigma_k)) + log_gaussian_density(x = log(v_norm), mu = 1, logsigma = log(sigma_v)) - 2*log(v_norm)
    return(log_prior_aniso)
  }
  return(log_prior_aniso)
}

#' @title Hyperparameter for variance of gaussian field given quantiles
#' @description Calculates  the hyperparameter lambda such that if sigma~Exp(lambda) then sigma>sigma0 with probability alpha
#'
#' @param alpha_sigma A quantile in (0,1)
#' @param sigma0 A surprisingly high number
#'
#' @return The calculated hyperparameter lambda
#' @export
#' @examples
#' alpha_sigma <- 0.01
#' sigma0 <- 10
#' lambda_sigma <- lambda_variance_quantile(alpha_sigma = alpha_sigma, sigma0 = sigma0)
lambda_variance_quantile <- function(alpha, sigma0) {
  #This warns the user if alpha is not in (0,1)
  if (alpha_sigma <= 0 | alpha_sigma >= 1) {
    warning("alpha_sigma should be in (0,1)")
  }
  if (sigma0 <= 0) {
    warning("sigma0 should be greater than 0")
  }
  lambda_sigma <- -log(alpha_sigma) / sigma0
  return(lambda_sigma)
}

#' @title Log prior on anisotropy, noise and variance of field supposing log(|v|)~N(1,sigma_v^2) and log(kappa)~N(1,sigma_k^2)
#' and with PC priors on noise and variance of field, given certain quantiles.
#'
#' @description
#' Calculates  the log of the prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' supposing log(|v|)~N(1,sigma_v^2) and log(kappa)~N(1,sigma_k^2)
#' and with PC priors on noise and variance of field, given certain quantiles.
#' such that the anisotropy ratio exp(2|v|) is larger than a0 with probability alpha and the correlation range sqrt{16}/kappa is smaller than rho0 with probability alpha
#'
#' @param alpha_u A quantile in (0,1) for the variance of the field
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise
#' @param alpha_k A quantile in (0,1) for kappa
#' @param alpha_v A quantile in (0,1) for |v|
#' @param sigmau0 A surprisingly high variance of field, in (0,infinity)
#' @param sigmaepsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#'
#' @return The calculated log of the prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' alpha_k <- 0.01
#' alpha_v <- 0.01
#' sigmau0 <- 10
#' sigmaepsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' log_prior_gaussian <- log_prior_aniso_quantile(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, alpha_k = alpha_k, alpha_v = alpha_v, sigmau0 = sigmau0, sigmaepsilon0 = sigmaepsilon0, a0 = a0, rho0 = rho0)
log_prior_quantile_gaussian <- function(alpha_u, alpha_epsilon, alpha_k, alpha_v, sigmau0, sigmaepsilon0, a0, rho0) {
  #This warns the user if alpha is not in (0,1)
  if (alpha_u <= 0 | alpha_u >= 1) {
    warning("alpha_u should be in (0,1)")
  }
  if (alpha_epsilon <= 0 | alpha_epsilon >= 1) {
    warning("alpha_epsilon should be in (0,1)")
  }
  if (alpha_k <= 0 | alpha_k >= 1) {
    warning("alpha_k should be in (0,1)")
  }
  if (alpha_v <= 0 | alpha_v >= 1) {
    warning("alpha_v should be in (0,1)")
  }
  #This warns the user if a0 is not greater than 1
  if (a0 <= 1) {
    warning("a0 should be greater than 1")
  }
  #This warns the user if rho0 is not greater than 0
  if (rho0 > 0) {
    warning("rho0 should be greater than 0")
  }
  #This warns the user if sigmau0 is not greater than 0
  if (sigmau0 <= 0) {
    warning("sigmau0 should be greater than 0")
  }
  #This warns the user if sigmaepsilon0 is not greater than 0
  if (sigmaepsilon0 <= 0) {
    warning("sigmaepsilon0 should be greater than 0")
  }
  sigma_v <- sigma_quantile_v(alpha_v = alpha_v, a0 = a0)
  sigma_k <- sigma_quantile_kappa(alpha_k = alpha_k, rho0 = rho0)
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigmau0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigmaepsilon0)
  log_prior <- function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
    variance_term <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon) + log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
    aniso_term <- log_prior_aniso_quantile_gaussian(alpha_k = alpha_k, alpha_v = alpha_v, rho0 = rho0, a0 = a0)(log_kappa = log_kappa, v = v)
    return(variance_term + aniso_term)
  }
  return(log_prior)
}
