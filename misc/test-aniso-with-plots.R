#Tests that functions fm_aniso_sample and fm_aniso precision for solving nonstationary SPDE and sampling work
#Plots field and covariance and variancee
library(devtools)
library(ggplot2)

# Parameter values stationary, anisotropic field
kp <- 1
v1 <- 0
v2 <- 1
lambda <- exp(sqrt(v1^2 + v2^2))
stretch <- sqrt(lambda)


# Parameter values Matérn field
nu <- 2 - 2 / 2
rh <- sqrt(8 * nu) / kp


# Kappa and vector field
kappa <- function(x) {
  return(kp)
}

vec <- function(x) {
  #v <- c(x[1],x[2])
   v <- c(max(-4,min(4,x[1])),max(-4,min(4,x[2])))
  #v <- c(max(-4,min(4,-x[2])),max(-4,min(4,x[1])))
  return(v)
}
vec <- function(x) {
  alpha <- 2*atan2(x[2],x[1])
  norm <- sqrt(x[1]^2+x[2]^2)
  v <- norm * c(-sin(alpha),cos(alpha))
  return(v)
}
# Square mesh for field
l <- 4
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(c(-l, -l), c(l, -l), c(l, l), c(-l, l), c(-l, -l)))))
#boundary_sf_2 <- st_sfc(st_polygon(list(rbind(c(-l, -l), c(l+1e-6, -l-1e-5), c(l, l), c(-l, l), c(-l, -l)))))
boundary_sf_zoom <- st_sfc(st_polygon(list(rbind(c(-1, 1), c(1, 1), c(1, 3), c(-1, 3), c(-1, 1)))))
boundary <- fm_as_segm(boundary_sf)

mesh <- fm_mesh_2d_inla(
  boundary = boundary,
  max.edge = c(0.2, 0.4)
)
mesh <- fm_mesh_2d_inla(
  loc = fm_lattice_2d(x=seq(-l,l,length.out=65),y=seq(-l,l,length.out=65))$loc,
  max.edge = c(0.4, 0.4)
)
nodes1 <- mesh$loc
#plot(mesh)

# Defining anisotropy
kappa_values <- apply(nodes1, 1, kappa)
vec_values <- t(apply(nodes1, 1, vec))
aniso <- list(kappa_values, vec_values)

# Sample of Matérn field u' with correlation range kp
sample_matern <- fm_matern_sample(mesh, alpha = 2, rho = rh, sigma = 1)

# Sample of anisotropic field u, should be equal to  u'(H^{-1/2}x)
sample_aniso <- fm_aniso_sample(mesh, aniso)
sample_matern <- fm_matern_sample(mesh, alpha = 2, rho = rh, sigma = 1)

# Data for plotting
field_matern <- data.frame(
  x = nodes1[, 1],
  y = nodes1[, 2],
  u = sample_matern
)
field_aniso <- data.frame(
  x = nodes1[, 1],
  y = nodes1[, 2],
  u = sample_aniso
)

# Definining pixels for plotting
pxl <- fm_pixels(mesh,dims=c(150,150), mask = boundary_sf)
# pxl$uisotropic <- fm_evaluate(mesh,
#   loc = pxl,
#   field = field_matern$u
# )
pxl$uaniso <- fm_evaluate(mesh,
  loc = pxl,
  field = field_aniso$u
)

# # Plotting Matérn field
# ggplot(pxl) +
#   geom_tile(aes(geometry = geometry, fill = uisotropic),
#     stat = "sf_coordinates", alpha = 1
#   ) +
#   scale_fill_gradientn(colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits = c(min(field_matern$u), max(field_matern$u))) +
#   coord_equal() +
#   xlab("X Coordinate") +
#   ylab("Y Coordinate") +
#   labs(fill = "Field u")

# Plotting anisotropic field
ggplot(pxl) +
  geom_tile(aes(geometry = geometry, fill = uaniso),
    stat = "sf_coordinates", alpha = 1
  ) +
  scale_fill_gradientn(colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits = c(-max(abs(field_aniso$u)), max(abs(field_aniso$u)))) +
  coord_equal() +
  xlab("X Coordinate") +
  ylab("Y Coordinate") +
  labs(fill = "Field u")

#Covariance plot
Q <-fm_aniso_precision(mesh, aniso)
A1 <- Matrix::Diagonal(nrow(Q))
A2 <- fm_basis(mesh,cbind(2,2))
K <- fm_covariance(Q,A1,A2)
pxl$K <- fm_evaluate(mesh,
                          loc = pxl,
                          field = as.vector(K)
)

ggplot(pxl) +
  geom_tile(aes(geometry = geometry, fill = K),
            stat = "sf_coordinates", alpha = 1
  ) +
  scale_fill_gradientn(colours = c( "#FFFFFFFF", "#FF0000FF"), limits = c(min(K), max(abs(K)))) +
  coord_equal() +
  xlab("X Coordinate") +
  ylab("Y Coordinate") +
  labs(fill = "Covariance")
#+   ggtitle("Covariance with (0,0)")

#Variance plot
V <- diag(INLA::inla.qinv(Q+Matrix::Diagonal(nrow(Q),1e-6)))
pxl$V <- fm_evaluate(mesh,
                     loc = pxl,
                     field = as.vector(sqrt(V))
)^2

ggplot(pxl) +
  geom_tile(aes(geometry = geometry, fill = V),
            stat = "sf_coordinates", alpha = 1
  ) +
  scale_fill_gradientn(colours = c( "#0000FFFF", "#FF0000FF"), limits = c(0, max(pxl$V))) +
  coord_equal() +
  xlab("X Coordinate") +
  ylab("Y Coordinate") +
  labs(fill = "Variance")
#+  gg(mesh,alpha=0)


