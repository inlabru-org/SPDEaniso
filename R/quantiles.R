#' @title Hyperparameter lambda1 for given quantile of the PC prior of |v|
#' @description Calculates the hyperparameter lambda1 for a given quantile alpha_a and large anisotropy a0 of the PC prior of |v|
#'
#' @param alpha_a A quantile in [0,1] of the PC prior of |v|
#' @param a0 A surprisingly high ratio of anisotropy
#' @export
#' @examples
#' alpha_a <- 0.01
#' a0 <- 10
#' lambda1 <- lambda1_quantile(a0 = a0, alpha_a = alpha_a)
lambda1_quantile <- function(a0, alpha_a = 0.01) {
    if (alpha_a <= 0 | alpha_a >= 1) {
        warning("alpha_a should be in (0,1)")
    }
    if (a0 <= 1) {
        warning("a0 should be greater than 1")
    }
    r0 <- log(a0) / 2
    lambda1 <- -log(alpha_a) / (fdist(r0) - fdist(0))
    return(lambda1)
}
#' @description Calculates the hyperparameter lambda for a given quantile alpha, tolerance kappa0, and lambda1 of the PC prior of kappa
#'
#' @param alpha A quantile in (0,1) of the PC prior of kappa
#' @param rho0 A surprisingly large correlation range
#' @export
#' @examples
#' alpha <- 0.01
#' rho0 <- 0.1
#' lambda1 <- 1
#' lambda <- lambda_quantile(alpha = alpha, rho0 = rho0, lambda1 = lambda1)
lambda_quantile <- function(rho0, lambda1, alpha = 0.01) {
  #This warns the user if alpha is not in (0,1)
  if (alpha <= 0 | alpha >= 1) {
    warning("alpha should be in (0,1)")
  }
  kappa0 <- sqrt(16)/ rho0
  product <- lambda1 * fdist(0)
  lambert <- lamW::lambertW0(product * exp(product) / (1-alpha))
  lambda <- (lambert/ fdist(0)- lambda1) / kappa0
  return(lambda)
}
#' @title This calculates sigma such that if x ~ N(0, sigma^2) then x >= x0 with probability alpha
#' @description Calculates the parameter sigma^2 such that if x ~ N(0, sigma^2) then x > x0 with probability alpha
#'
#' @param alpha A quantile in (0,1)
#' @param x0 A surprisingly high number
#'
#' @return The calculated parameter sigma^2
#' @export
#' @examples
#' alpha <- 0.01
#' x0 <- 10
#' sigma <- sigma_quantile(alpha = alpha, x0 = x0)
sigma_quantile <- function(alpha, x0) {
    if (alpha <= 0 | alpha >= 1) {
        warning("alpha should be in (0,1)")
    }
    sigma <- x0 / qnorm(1 - alpha)
    return(sigma)
}

#' @title Parameter sigma such that if v1,v2 ~ N(0, sigma^2) the anisotropy ratio exp(2|v|) is larger than a0 with probability alpha
#' @description Calculates the parameter sigma^2 such that if v1,v2 are iid N(0, sigma^2), then the anisotropy ratio exp(2|v|) is larger than a0 with probability alpha
#'
#' @param alpha_v A quantile in (0,1)
#' @param a0 A surprisingly high ratio of anisotropy
#'
#' @return The calculated parameter sigma^2
#' @export
#' @examples
#' alpha_v <- 0.01
#' a0 <- 10
#' sigma_v <- sigma_quantile_v(alpha_v = alpha_v, a0 = a0)
sigma_quantile_v <- function(alpha_v, a0) {
    if (alpha_v <= 0 | alpha_v >= 1) {
        warning("alpha_v should be in (0,1)")
    }
    if (a0 <= 1) {
        warning("a0 should be greater than 1")
    }
    r02 <- (log(a0) / 2)^2
    # r^2 follows a exponential distribution with rate 1/2sigma^2
    sigma_v <- sqrt(-r02 / (2 * log(alfa_v)))
    return(sigma_v)
}
