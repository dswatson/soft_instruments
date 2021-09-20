# Set working directory
setwd('~/Documents/UCL/soft_instruments')

# Load libraries, register cores
library(data.table)
library(MultiRNG)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Define simulation grid
sim_idx <- expand_grid(
  z_rho = c(0, 0.5),
  rho = c(0.25, 0.75, -0.25, -0.75),
  pr_valid = c(0, 0.25, 0.5)
) %>% 
  mutate(idx = row_number()) %>%
  select(idx, z_rho, rho, pr_valid)

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param z_cnt Logical. If \code{TRUE}, Z is multivariate normal; else, 
#'   multivariate binomial.
#' @param z_rho Autocorrelation of the Toeplitz matrix for Z.
#' @param rho Correlation between residuals for X and Y.
#' @param alpha Average treatment effect of X on Y.
#' @param r2_x Proportion of variance explained for X.
#' @param r2_y Proportion of variance explained for Y.
#' @param pr_valid Proportion of candidate instruments that are valid, 
#'   i.e. have no direct effect on Y.
#' @param s_idx Simulation index

# Simulate data (inspired by Hartwig et al., 2017)
sim_dat <- function(n, d_z, z_cnt, z_rho, rho, alpha, r2_x, r2_y, pr_valid, s_idx) {
  # What proportion are valid?
  valid_cnt <- round(pr_valid * d_z)
  pr_valid <- valid_cnt / d_z
  # Draw Z's
  if (z_rho != 0) {
    Sigma_z <- toeplitz(z_rho^(0:(d_z - 1)))
  }
  if (z_cnt) {
    z <- matrix(rnorm(n * d_z), ncol = d_z)
    if (z_rho == 0) {
      Sigma_z <- diag(rep(1, d_z))
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
  gamma <- (sqrt(0.1) / sigma_zy) * nu2 * (1 - pr_valid)
  if (pr_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% gamma) + alpha^2 * var_x + 
    2 * alpha * as.numeric(beta %*% Sigma_z %*% gamma)
  var_y <- var_mu / r2_y
  b <- 2 * alpha * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * alpha + eps[, 2] * eta_y
  # Export results
  dat <- data.table(z, 'x' = x, 'y' = y)
  fwrite(dat, paste0('./simulations/sim', s_idx, '.csv'))
  params <- list(
    'alpha' = alpha, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'r2_x' = r2_x, 'r2_y' = r2_y
  )
  saveRDS(params, paste0('./simulations/params', s_idx, '.rds'))
}
# Execute in parallel
foreach(aa = sim_idx$idx) %dopar%
  sim_dat(n = 1e4, d_z = 8, z_cnt = TRUE, z_rho = sim_idx$z_rho[aa],
          rho = sim_idx$rho[aa], alpha = 1, r2_x = 0.5, r2_y = 0.5, 
          pr_valid = sim_idx$pr_valid[aa], aa)

# And also
fwrite(sim_idx, './simulations/sim_idx.csv')


# Problem: pr_z imposes some constraints on Sigma_z in the binomial case,
# as correlations cannot exceed upper limits imposed by expectations

# Do we really want to bound beta and gamma away from 0?

# Simulate data 
sim_dat <- function(n, d_z, z_rho, rho, alpha, r2_x, r2_y, pr_valid, s_idx) {
  # What proportion are valid?
  valid_cnt <- round(pr_valid * d_z)
  pr_valid <- valid_cnt / d_z
  # Draw Z's
  var_z <- 1 / d_z
  z <- matrix(rnorm(n * d_z, sd = sqrt(var_z)), ncol = d_z)
  if (z_rho == 0) {
    Sigma_z <- diag(var_z, d_z)
  } else {
    Sigma_z <- toeplitz(z_rho^(0:(d_z - 1)))
    z <- z %*% chol(Sigma_z)
  }
  colnames(z) <- paste0('z', seq_len(d_z))
  # Simulate standardized residual vectors
  rho <- rho * sample(c(-1, 1), 1)
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
  eps <- matrix(rnorm(n * 2), ncol = 2)
  eps <- eps %*% chol(Sigma_eps)
  # Draw random beta, set var_mu = r2_x
  beta <- runif(d_z, min = -1, 1)
  fn <- function(lambda) {
    var_mu <- as.numeric(t(beta * lambda) %*% Sigma_z %*% (beta * lambda))
    (var_mu - r2_x)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = 0, upper = 10)$par
  beta <- beta * lambda
  var_x <- 1
  eta_x <- sqrt(var_x * (1 - r2_x))
  x <- as.numeric(z %*% beta) + eps[, 1] * eta_x
  # Set alpha to account for a third, a half or two thirds of signal variance
  if (alpha == 'low') {
    alpha <- sqrt(r2_y / 3)
  } else if (alpha == 'med') {
    alpha <- sqrt(r2_y / 2) 
  } else if (alpha == 'high') {
    alpha <- sqrt(2 * r2_y / 3) 
  }
  alpha <- alpha * sample(c(-1, 1), 1)
  # Draw random gamma
  gamma <- runif(d_z, min = -1, 1) / (d_z * (1 - pr_valid))
  if (pr_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  fn <- function(lambda) {
    var_mu <- as.numeric(t(gamma * lambda) %*% Sigma_z %*% (gamma * lambda)) +
      alpha^2 * var_x + 2 * alpha * as.numeric(beta %*% Sigma_z %*% (gamma * lambda))
    (var_mu - r2_y)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = 0, upper = 10)$par
  gamma <- gamma * lambda
  # Solve a quadratic equation for eta_y
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% gamma) + alpha^2 * var_x + 
    2 * alpha * as.numeric(beta %*% Sigma_z %*% gamma)
  var_y <- 1
  b <- 2 * alpha * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * alpha + eps[, 2] * eta_y
  # Export results
  dat <- data.table(z, 'x' = x, 'y' = y)
  params <- list(
    'alpha' = alpha, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'r2_x' = r2_x, 'r2_y' = r2_y
  )
  return(list('dat' = dat, 'params' = params))
}






