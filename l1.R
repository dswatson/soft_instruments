
### Experiment in L1 geometry ###
# Hyperparameters
n <- 1e4
rho <- 0.25
alpha <- 0.5
r2_x <- 0.75
r2_y <- 0.5
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

# Not sure about the intercept tho...

# Wonder how this extends to the multidimensional case?
n <- 1e4
rho <- 0.25
alpha <- 0.5
r2_x <- 0.75
r2_y <- 0.5
d_z <- 2
# Sample multivariate Z
z <- matrix(rnorm(n * d_z), ncol = d_z)
Sigma_z <- diag(1, d_z)
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
x <- z * beta + eps[, 1] * eta_x
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
y <- as.numeric(z %*% gamma) + x * alpha + eps[, 2] * eta_y


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