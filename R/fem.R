#' @title Anisotropic FEM
#' @description Calculates FEM matrices for anisotropic field.
#' @param mesh An `fm_mesh_2d` object
#' @param aniso A `list(kappa, v)`.
#' Calculates anisotropic structure matrices for
#'  an anisotropic operator \eqn{\kappa^2-\nabla\cdot H \nabla}{kappa^2-div H grad}.
#' Here \eqn{\kappa>0,v=(v_1,v_2)\in\mathbb{R}^2}{kappa>0,v=(v1,v2) in R^2} and
#'  \eqn{H=e^{|v|}\tilde{v}\tilde{v}^T+e^{-|v|}\tilde{v}_\perp\tilde{v}^T_\perp}{H = H = exp(|v|) * v_tilde * v_tilde' + exp(-|v|) * v_tilde_perp * v_tilde_perp'}., where
#' and \eqn{\tilde{v}=|v| e^{i \alpha /2 }, \alpha := \arctan(v_2 /v_1)}{v_tilde=|v|exp(i alpha /2),arctan(v2 /v1)}.
#' @examples
# mesh <- fm_rcdt_2d_inla(globe = 1)
#' fem3 <- fm_fem_aniso(mesh, aniso = list(kappa = rep(1, mesh$n), v = matrix(0, mesh$n, 3)))
#' @return `fm_fem_aniso`: A list with elements `c0`, `c1`, `g1`, `g2` `va`, `ta`,
#'
#' @export
fm_fem_aniso <- function(mesh, aniso, ...) {
  if (!is.list(aniso) || length(aniso) != 2) {
    stop("'aniso' must be NULL or a list of length 2.")
  }
  result <- fmesher_fem_aniso(
    mesh_loc = mesh$loc,
    mesh_tv = mesh$graph$tv - 1L,
    aniso = aniso,
    options = list()
  )
  result
}

