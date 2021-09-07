# Load libraries, register cores
library(data.table)
library(MultiRNG)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param z_cnt Logical. If \code{TRUE}, Z is multivariate normal; else, 
#'   multivariate binomial.
#' @param z_rho Autocorrelation of the Toeplitz matrix for Z.
#' @param rho Correlation between residuals for X and Y.
#' @param alpha Average treatment effect of X on Y.
#' @param r2_x Proportion of variance explained for X.
#' @param r2_y Proportion of variance explained for Y.
#' @param prop_valid Proportion of candidate instruments that are valid, 
#'   i.e. have no direct effect on Y.

# Simulate data (inspired by Hartwig et al., 2017)
sim_dat <- function(n, d_z, z_cnt, z_rho,
                    rho, alpha, r2_x, r2_y, prop_valid) {
  # What proportion are valid?
  valid_cnt <- round(prop_valid * d_z)
  prop_valid <- valid_cnt / d_z
  # Draw Z's
  if (z_rho != 0) {
    Sigma_z <- toeplitz(z_rho^(0:(d_z - 1)))
  }
  if (z_cnt) {
    z <- matrix(rnorm(n * d_z), ncol = d_z)
    if (z_rho == 0) {
      Sigma_z < diag(rep(1, d_z))
    } else {
      z <- z %*% chol(Sigma_z)
    }
  } else {
    pr_z <- runif(d_z, min = 0.1, max = 0.9)
    var_z <- pr_z * (1 - pr_z)
    if (z_rho == 0) {
      z <- sapply(pr_z, function(p) rbinom(n, 1, p))
      Sigma_z <- diag(var_z)
    } else {
      z <- draw.correlated.binary(no.row = n, d = d_z, prop.vec = pr_z, 
                                  corr.mat = Sigma_z)
      Sigma_z <- diag(sqrt(var_z)) %*% Sigma_z %*% diag(sqrt(var_z))
    }
  }
  colnames(z) <- paste0('z', seq_len(d_z))
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
  gamma <- (sqrt(0.1) / sigma_zy) * nu2 * (1 - prop_valid)
  if (prop_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% gamma) + alpha^2 * var_x + 
    2 * alpha * as.numeric(beta %*% Sigma_z %*% gamma)
  var_y <- var_mu / r2_y
  b <- 2 * alpha * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * alpha + eps[, 2] * eta_y
  dat <- data.table(z, 'x' = x, 'y' = y)
  params <- list(
    'alpha' = alpha, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'r2_x' = r2_x, 'r2_y' = r2_y
  )
  out <- list('dat' = dat, 'params' = params)
  return(out)
}
# For example:
sim <- sim_dat(n = 1e4, d_z = 5, z_cnt = 6, z_rho = 0.25, 
               rho = 0.25, alpha = 1, r2_x = 0.5, r2_y = 0.5, prop_valid = 0.5)
head(sim$dat)
sim$params




