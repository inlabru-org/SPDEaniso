% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantiles.R
\name{log_prior_uniform}
\alias{log_prior_uniform}
\title{Log uniform prior on anisotropy, noise and variance of field}
\usage{
log_prior_uniform(
  sigma_u0,
  sigma_epsilon0,
  a0 = NULL,
  rho0 = NULL,
  L = NULL,
  width_support_factor = 2,
  alpha = 0.01,
  alpha_u = NULL,
  alpha_epsilon = NULL
)
}
\arguments{
\item{sigma_u0}{A surprisingly high variance of field, in (0,infinity)}

\item{sigma_epsilon0}{A surprisingly high variance of noise, in (0,infinity)}

\item{a0}{A surprisingly high ratio of anisotropy, in (1,infinity), anisotropy ratio is not allowed to go below a0/2 or above 2a0.}

\item{rho0}{A surprisingly small correlation range, correlation is not allowed to go below rho0/2}

\item{L}{A surprisingly large correlation range (e.g. size of domain), correlation is not allowed to go above L}

\item{width_support_factor}{A factor to multiply the width of the support of the uniform distribution. By default 2. If Inf,an improper uniform prior is used}

\item{alpha}{A quantile in (0,1) for all the parameters. By default, alpha = 0.01}

\item{alpha_u}{A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha}

\item{alpha_epsilon}{A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha}
}
\value{
The calculated log of the prior on  theta=(log(kappa),v, log(sigma_u), log(sigma_epsilon))
}
\description{
Calculates  the log of the prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon)) supposing:
kappa ~ Uniform(-log(w\emph{L)+1/2}log(8),-log(rho0/w)+1/2\emph{log(8)), so that rho= sqrt{8}/kappa has support \link{rho0/w,w*L}
|v1|,|v2| ~ Uniform(0,log(w}a0)/sqrt(2)), so that exp(|v|) has support \link{1,a0*w}
and with PC priors on noise and variance of field, given certain quantiles.
}
\examples{
alpha_u <- 0.01
alpha_epsilon <- 0.01
sigma_u0 <- 10
sigma_epsilon0 <- 1
a0 <- 10
rho0 <- 0.1
L <- 100

log_prior <- log_prior_uniform(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, width_support_factor = 2)
}
