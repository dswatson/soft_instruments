# Load libraries, register cores
library(data.table)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Simulate data (inspired by Hartwig et al., 2017)
sim_dat <- function(n, d_z, rho, alpha, r2_x, r2_y) {
  # SNPs in biallelic Hardy-Weinberg equilibrium
  pr_z <- runif(d_z, min = 0.1, max = 0.9)
  z <- sapply(pr_z, function(p) rbinom(n, size = 1, prob = p))
  colnames(z) <- paste0('z', seq_len(d_z))
  # Extract covariance matrix
  var_z <- sapply(pr_z, function(p) p * (1 - p))
  Sigma_z <- diag(var_z)
  # Simulate standardized residual vectors
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
  eps <- matrix(rnorm(n * 2), ncol = 2)
  eps <- eps %*% chol(Sigma_eps)
  # Draw random weights for Z (adapted from Hartford et al., 2020)
  # Then calculate eta_x from r2_x
  nu1 <- runif(d_z, min = 0.01, max = 0.2)
  sigma_zx <- sd(sqrt(0.1) * as.numeric(z %*% nu1))
  beta <- (sqrt(0.1) / sigma_zx) * nu1
  var_x <- as.numeric(t(beta) %*% Sigma_z %*% beta) / r2_x 
  eta_x <- sqrt(var_x * (1 - r2_x))
  x <- as.numeric(z %*% beta) + eps[, 1] * eta_x
  # And again for Y, although we need to solve a quadratic equation for eta_y
  nu2 <- runif(d_z, min = 0.01, max = 0.2)
  sigma_zy <- sd(sqrt(0.1) * as.numeric(z %*% nu2))
  gamma <- (sqrt(0.1) / sigma_zy) * nu2 # Not imposing any sparsity here...
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% gamma) + alpha^2 * var_x + 
    2 * alpha * as.numeric(beta %*% Sigma_z %*% gamma)
  var_y <- var_mu / r2_y
  b <- 2 * alpha * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * alpha + eps[, 2] * eta_y
  dat <- data.table(z, 'x' = x, 'y' = y)
  params <- list(
    'alpha' = alpha, 'beta' = beta, 'gamma' = gamma, 'eps_x' = eps_x, 
    'eps_y' = eps_y, 'rho' = rho, 'r2_x' = r2_x, 'r2_y' = r2_y
  )
  out <- list('dat' = dat, 'params' = params)
  return(out)
}
# For example:
sim <- sim_dat(n = 1e4, d_z = 5, rho = 0.25, alpha = 1, r2_x = 0.5, r2_y = 0.5)






