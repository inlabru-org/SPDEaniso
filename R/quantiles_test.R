a0 <- 10
alpha <- 0.01
lambda1 <- lambda1_quantile(a0 =a0,alpha_a = alpha)

lambda <- lambda_quantile(alpha=alpha, rho0=rho0, lambda1=lambda1)

x0=1000
sigma_quantile(alpha = alpha, x0=x0)
sigma_v <- sigma_quantile_v(alpha_v = alpha, a0 = a0)
