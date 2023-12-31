% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCPrior.R
\name{sim_aniso_pc}
\alias{sim_aniso_pc}
\title{Simulation of anisotropy parameters (log(kappa), v) from the PC prior}
\usage{
sim_aniso_pc(lambda, lambda1, m = 1)
}
\arguments{
\item{lambda}{A hyperparameter controlling the size of kappa.}

\item{lambda1}{A hyperparameter controlling the size of |v|.}

\item{m}{Number of samples, by default 1.}
}
\value{
A list with two elements: log(kappa) and v
}
\description{
Simulates anisotropy parameters (log(kappa), v) from the PC prior
}
\examples{
lambda <- 1
lambda1 <- 1
m <- 1
result <- sim_aniso_pc(lambda = lambda, lambda1 = lambda1, m = m)
}
