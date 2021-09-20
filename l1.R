### Experiment in L1 geometry ###
# Hyperparameters
n <- 1e4
rho <- 0.25
alpha <- 0.5
r2_x <- 0.75
r2_y <- 0.75
# Sample univariate Z
z <- rnorm(n)
# Simulate standardized residual vectors
Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
eps <- matrix(rnorm(n * 2), ncol = 2)
eps <- eps %*% chol(Sigma_eps)
# Set beta to enforce r2_x
beta <- sqrt(r2_x)
var_x <- 1
eta_x <- sqrt(var_x * (1 - r2_x))
x <- z * beta + eps[, 1] * eta_x
# Set gamma to enforce r2_y
bb <- 2 * alpha * beta 
cc <- alpha^2 * var_x - r2_y
gamma <- c(
  (-bb + sqrt(bb^2 - 4 * cc)) / 2,
  (-bb - sqrt(bb^2 - 4 * cc)) / 2
)
var_y <- 1
bb <- 2 * alpha * rho * eta_x
eta_y <- (-bb + sqrt(bb^2 - 4 * (r2_y - var_y))) / 2
y <- z * gamma[1] + x * alpha + eps[, 2] * eta_y  # COULD BE GAMMA[2] THO

# Compute over range of L1 bounds, plot results
fn <- function(tau) {
  gammas <- c(tau, -tau)
  alphas <- sapply(gammas, function(g) {
    rhs <- y - z * g - eps[, 2] * eta_y
    tmp <- data.frame('x' = x, 'rhs' = rhs)
    ols <- lm(rhs ~ 0 + x, tmp)
    tidy(ols)$estimate[1]
  })
  out <- data.frame('l1' = tau, 'lo' = alphas[1], 'hi' = alphas[2])
  return(out)
}
df <- foreach(taus = seq(0, 5, 0.01), .combine = rbind) %dopar% fn(taus)
ggplot(df) + 
  geom_ribbon(aes(l1, ymin = lo, ymax = hi), alpha = 0.5) + 
  theme_bw() + 
  labs(x = 'L1(gamma)', y = 'alpha')

# Slope of top and bottom lines appears to be given by beta?
tidy(lm(hi ~ l1, df))


# Compute rho-alpha for rho \in [-1, 1] and solve for gamma at each point
# Rejection sample: L1(gamma) <= tau?


# Is there a falsifiable range of gamma that constrains the range of rho?

# Not sure about the intercept tho...

# Wonder how this extends to the multidimensional case?
n <- 1e4
rho <- 0.5
alpha <- 0.5
r2_x <- 0.75
r2_y <- 0.75
d_z <- 4

loop_fn <- function(d_z, rho_z, rho, alpha, r2_x, r2_y) {
  ### STEP 1: SIMULATE PARAMETERS TO CONFORM TO STRUCTURAL ASSUMPTIONS ###
  # Covariance of Z
  if (rho_z == 0) {
    Sigma_z <- diag(1, d_z)
  } else {
    Sigma_z <- toeplitz(rho_z^(0:(d_z - 1)))
    z <- z %*% chol(Sigma_z)
  }
  colnames(Sigma_z) <- paste0('z', seq_len(d_z))
  # Simulate standardized residual vectors
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
  eps <- matrix(rnorm(n * 2), ncol = 2)
  eps <- eps %*% chol(Sigma_eps)
  # Set beta to enforce r2_x
  beta <- runif(d_z, min = -1, max = 1)
  fn <- function(lambda) {
    var_mu <- as.numeric(t(beta * lambda) %*% Sigma_z %*% (beta * lambda))
    (var_mu - r2_x)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = 0, upper = 10)$par
  beta <- beta * lambda
  var_x <- 1
  eta_x <- sqrt(var_x * (1 - r2_x))
  # Set gamma to enforce r2_y
  gamma <- runif(d_z, min = -1, max = 1)
  fn <- function(lambda) {
    var_mu <- as.numeric(t(gamma * lambda) %*% Sigma_z %*% (gamma * lambda)) +
      alpha^2 * var_x + 2 * alpha * as.numeric(beta %*% Sigma_z %*% (gamma * lambda))
    (var_mu - r2_y)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = 0, upper = 10)$par
  gamma <- gamma * lambda
  var_mu <- r2_y
  var_y <- 1
  b <- 2 * alpha * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  Theta_z <- solve(Sigma_z)
  Theta_z2 <- Theta_z %*% Theta_z
  Sigma_zx <- Sigma_z %*% beta
  Sigma_xz <- t(Sigma_zx)
  Sigma_zy <- Sigma_z %*% gamma + Sigma_zx %*% alpha
  Sigma_yz <- t(Sigma_zy)
  sigma_xy <- Sigma_xz %*% gamma + alpha * var_x + eta_x * eta_y * rho
  ### STEP 2: COMPUTE ALPHA AND GAMMA FROM RHO
  aa <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zx)
  bb <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zy)
  cc <- as.numeric(Sigma_yz %*% Theta_z2 %*% Sigma_zx)
  dd <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
  ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
  ff <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
  from_rho <- function(rho_in) {
    gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho_in^2)))
    hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho_in^2)))
    ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho_in^2)) + ff - var_y
    if (rho_in < 0) {
      alpha_hat <- as.numeric((-hh + sqrt(hh^2 - gg * ii)) / gg)
    } else {
      alpha_hat <- as.numeric((-hh - sqrt(hh^2 - gg * ii)) / gg)
    }
    gamma_hat <- as.numeric(Theta_z %*% (Sigma_zy - alpha_hat * Sigma_zx))
    out <- data.frame(
      'rho_in' = rho_in, 'alpha_hat' = alpha_hat,
      'l1' = sum(abs(gamma_hat)), 'l2' = sum(gamma_hat^2)
    )
    return(out)
  }
  # Loop over rhos
  rhos <- seq(-0.99, 0.99, length.out = 200)
  out <- foreach(r = rhos, .combine = rbind) %do% from_rho(r)
  # Record configurations
  out <- out %>%
    mutate('d_z' = d_z, 'rho_z' = rho_z, 'rho' = rho, 'alpha' = alpha, 
           'r2_x' = r2_x, 'r2_y' = r2_y)
  return(out)
}
df <- foreach(d_zs = c(2, 3, 4), .combine = rbind) %:%
  foreach(rho_zs = c(0, 0.5), .combine = rbind) %:%
  foreach(rhos = c(-0.75, -0.25, 0.25, 0.75), .combine = rbind) %:%
  foreach(alphas = c(-0.2, 0.2), .combine = rbind) %dopar%
  loop_fn(1e4, d_zs, rho_zs, rhos, alphas, 0.5, 0.5)

# Plot results
df %>%
  filter(rho_z == 0, alpha == 0.2) %>%
  ggplot(aes(l1, alpha_hat)) + 
  geom_path() + 
  geom_hline(aes(yintercept = alpha), color = 'red', linetype = 'dashed') +
  theme_bw() + 
  labs(x = 'L1(gamma)', y = 'alpha') + 
  facet_grid(d_z ~ rho)

df %>%
  filter(rho_z == 0, alpha == 0.2) %>%
  ggplot(aes(rho_in, l1)) + 
  geom_path() + 
  geom_vline(aes(xintercept = rho), color = 'red', linetype = 'dashed') +
  theme_bw() + 
  labs(x = 'rho', y = 'L1(gamma)') + 
  facet_grid(d_z ~ rho)






# When gamma = 0-vector
cc <- alpha^2 * var_x - var_y
eta_y <- (-b + sqrt(b^2 - 4 * cc)) / 2

# Compute over range of L1 bounds, plot results
# There are 2^d_z extrema of the hypercube within which all gammas 
# must exist. Let's visit them!

# NOOOO THIS DOESN'T WORK B/C THESE GAMMA VECTORS ARE NOT CONSISTENT WITH
# THE DATA COVARIANCE
fn <- function(tau) {
  gammas <- list(
    c(tau, 0), c(-tau, 0), c(0, tau), c(0, -tau)
  )
  alphas <- sapply(gammas, function(g) {
    rhs <- y - as.numeric(z %*% g) - eps[, 2] * eta_y
    tmp <- data.frame('x' = x, 'rhs' = rhs)
    ols <- lm(rhs ~ 0 + x, tmp)
    tidy(ols)$estimate[1]
  })
  out <- data.frame('l1' = tau, 'lo' = min(alphas), 'hi' = max(alphas))
  return(out)
}
df <- foreach(taus = seq(0, 5, 0.01), .combine = rbind) %dopar% fn(taus)
ggplot(df) + 
  geom_ribbon(aes(l1, ymin = lo, ymax = hi), alpha = 0.5) + 
  theme_bw() + 
  labs(x = 'L1(gamma)', y = 'alpha')
# Slope of top and bottom lines appears to be given by beta?
tidy(lm(hi ~ l1, df))


# gamma as function of rho?
gamma_fn <- function(rho_in) {
  # Data covariance
  Theta_z <- solve(Sigma_z)
  Theta_z2 <- Theta_z %*% Theta_z
  Sigma_zx <- Sigma_z %*% beta
  Sigma_xz <- t(Sigma_zx)
  Sigma_zy <- Sigma_z %*% gamma + Sigma_zx %*% alpha
  Sigma_yz <- t(Sigma_zy)
  sigma_xy <- Sigma_xz %*% gamma + alpha * var_x + eta_x * eta_y * rho
  # Define
  aa <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zx)
  bb <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zy)
  cc <- as.numeric(Sigma_yz %*% Theta_z2 %*% Sigma_zx)
  dd <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
  ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
  ff <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
  alpha_fn <- function(rho) {
    gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
    hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
    ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho^2)) + ff - var_y
    if (rho < 0) {
      alpha <- (-hh + sqrt(hh^2 - gg * ii)) / gg 
    } else {
      alpha <- (-hh - sqrt(hh^2 - gg * ii)) / gg 
    }
    return(alpha)
  }
  alpha_hat <- as.numeric(alpha_fn(rho_in))
  gamma_hat <- as.numeric(Theta_z %*% (Sigma_zy - alpha_hat * Sigma_zx))
  out <- data.frame('l1' = sum(abs(gamma_hat)), 
                    'l2' = sum(gamma_hat^2),
                    'rho' = rho_in, 
                    'alpha' = alpha_hat)
  return(out)
}
df <- foreach(rhos = seq(-0.99, 0.99, length.out = 200), .combine = rbind) %dopar%
  gamma_fn(rhos)
ggplot(df, aes(l1, rho)) + 
  geom_point() + 
  theme_bw()
ggplot(df, aes(l1, alpha)) + 
  geom_point() + 
  theme_bw()
ggplot(df, aes(l2, rho)) + 
  geom_point() + 
  theme_bw()
ggplot(df, aes(l2, alpha)) + 
  geom_point() + 
  theme_bw()


# Feasible region will take the form
j - k * tau <= alpha <= j + k * tau

# Set this to zero and solve for alpha
Theta_z %*% (Sigma_zy - alpha * Sigma_zx) = rep(0, d_z)
Theta_z %*% Sigma_zy - alpha * Theta_z %*% Sigma_zx = rep(0, d_z)


alpha * Theta_z %*% Sigma_zx = Theta_z %*% Sigma_zy


alpha * Sigma_zx = Sigma_zy

alpha = Sigma_zy / Sigma_zx 

alpha <- mean(Sigma_zy / Sigma_zx) # ???






