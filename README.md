# SPDEaniso: Tool for working with anisotropic fields using SPDEs
This package can be used to simulate from and perform Bayesian inference on anisotropic Gaussian fields using stochastic partial differential equations (SPDEs). The package is based on the [`fmesher`](https://github.com/inlabru-org/fmesher) package. Anisotropy is implemented through a paretrization of the anisotropic SPDE

```math
(\kappa^2-\nabla\cdot \mathbf{H}({\mathbf{v}})\nabla)u=\kappa\sigma\mathcal{W}
```

The parameters are $\kappa >0$, a two dimensional vector $v$ and $\sigma>0$. These parameters control the length scale and anisotropy, respectively.

Spatially varying parameters $\kappa(x),v(x)$ are supported for simulation. For Bayesian simulation we consider a linear, noisy observation process

```math
\mathbf{y} = \mathbf{A}\mathbf{u} + \mathbf{\epsilon},
```

where $\mathbf{A}$ is the observation matrix and $\mathbf{\epsilon}\sim\mathcal{N}(0,\sigma_{\mathbf{\epsilon}}^2\mathbf{I})$ is a vector of independent Gaussian noise. The package supports Bayesian inference on 
```math
\theta:=(\kappa, \mathbf{v}, \sigma,\sigma_{\mathbf{\epsilon}}) ,
```
for spatially _constant_ parameters. Here penalized complexity priors are used both for the anisotropy parameters (see...) and the [noise parameters](https://arxiv.org/abs/1403.4630#:~:text=Proper%20priors%20are%20defined%20to%20penalise%20the%20complexity,both%20in%20the%20univariate%20and%20the%20multivariate%20case.).
