



soft_iv <- function(dat, rho) {
  # Estimate eta_x
  f1 <- lm(x ~ ., data = select(dat, -y))
  eta_x <- sqrt(mean((residuals(f1)^2)))
  # Estimate data covariance
  Sigma_z <- cov(select(dat, starts_with('z')))
  Theta_z <- solve(Sigma_z)
  #Theta_z2 <- Theta_z %*% Theta_z
  Sigma_zy <- cov(select(dat, starts_with('z')), dat$y)
  Sigma_zx <- cov(select(dat, starts_with('z')), dat$x)
  Sigma_yz <- cov(dat$y, select(dat, starts_with('z')))
  Sigma_xz <- cov(dat$x, select(dat, starts_with('z')))
  sigma_xy <- cov(dat$x, dat$y)
  var_x <- var(dat$x)
  var_y <- var(dat$y)
  # Simplify
  aa <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
  bb <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
  cc <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
  dd <- (bb - var_x) * (1 + ((bb - var_x) / (eta_x^2 * rho^2)))
  ee <- -(aa - sigma_xy) * (1 + ((bb - var_x) / (eta_x^2 * rho^2)))
  ff <- ((aa - sigma_xy)^2 / (eta_x^2 * rho^2)) + cc - var_y
  # Quadratics
  out <- data.frame(
    alpha_lo = (-ee - sqrt(ee^2 - dd * ff)) / dd,
    alpha_hi = (-ee + sqrt(ee^2 - dd * ff)) / dd 
  )
  return(out)
}
loop_fn <- function(AA, BB, CC, DD, EE, FF) {
  tmp <- sim_dat(
    n = 2000, d_z = 4, z_cnt = TRUE, z_rho = BB,
    rho = CC, alpha = 1, r2_x = 0.5, r2_y = 0.5, 
    prop_valid = DD, AA
  )
  alphas <- soft_iv(tmp, CC)
  out <- data.table(
    'alpha_lo' = alphas$alpha_lo, 'alpha_hi' = alphas$alpha_hi,
    'z_rho' = BB, 'rho' = CC, 'r2_x' = DD, 'r2_y' = EE, 'pr_valid' = FF
  )
}

df <- foreach(idxs = seq_len(2000), .combine = rbind) %:%
  foreach(z_rhos = c(0, 0.5), .combine = rbind) %:%
  foreach(rhos = c(0.25, 0.5, 0.75), .combine = rbind) %:%
  foreach(r2_x = c(0.25, 0.5, 0.75), .combine = rbind) %:%
  foreach(r2_y = c(0.25, 0.5, 0.75), .combine = rbind) %:%
  foreach(props = c(0, 0.5), .combine = rbind)  %dopar% 
  loop_fn(idxs, z_rhos, rhos, props)


df[pr_valid == 0] %>%
  ggplot(aes(alpha_lo)) + 
  geom_histogram(bins = 50) +  
  geom_vline(xintercept = 1, color = 'red', linetype = 'dashed') +
  theme_bw() + 
  facet_grid(z_rho ~ rho)

df %>%
  filter(z_rho == 0, props == 0) %>%
  pivot_longer(starts_with('alpha'), names_to = 'root', values_to = 'alpha') %>%
  ggplot(aes(rho, alpha, fill = root)) + 
  geom_boxplot() + 
  scale_fill_d3() + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  facet_grid(r2_x ~ r2_y)



soft_iv <- function(Sigma, eta_x, rho) {
  Theta_z <- solve(Sigma_z)
  Theta_z2 <- Theta_z %*% Theta_z
  Sigma_zx <- Sigma_z %*% beta
  Sigma_xz <- t(Sigma_zx)
  Sigma_zy <- Sigma_z %*% gamma + Sigma_zx %*% alpha
  #var_x <- t(beta) %*% Sigma_z %*% beta + eta_x^2
  sigma_xy <- Sigma_xz %*% gamma + alpha %*% var_x + eta_x * eta_y * rho
 # var_y <- t(gamma) %*% Sigma_z %*% gamma + 2 * alpha * Sigma_xz %*% gamma + 
 #   alpha^2 * var_x + 2 * alpha * eta_x * eta_y * rho + eta_y^2
  
}

# NEW IDEA: Don't care about r^2. eta_x = eta_y = 1.
z <- matrix(rnorm(n * d_z), ncol = d_z)
Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
eps <- matrix(rnorm(n * 2), ncol = 2)
eps <- eps %*% chol(Sigma_eps)
beta <- runif(d_z)
gamma <- runif(d_z)
x <- z %*% beta + eps[, 1]
y <- z %*% gamma + x * alpha + eps[, 2]

# Now...
eta_x <- eta_y <- alpha <- 1
Sigma_z <- diag(1, d_z)
Sigma_zx <- Sigma_z %*% beta
Sigma_xz <- t(Sigma_zx)
Sigma_zy <- Sigma_z %*% gamma + Sigma_zx %*% alpha
var_x <- t(beta) %*% Sigma_z %*% beta + eta_x^2

var_y <- var_x + 2 * (Sigma_xz %*% gamma + rho) + t(gamma) %*% Sigma_z %*% gamma + 1
var(y)

aa <- var_x
bb <- 2 * (Sigma_xz %*% gamma + eta_x * eta_y * rho)
cc <- t(gamma) %*% Sigma_z %*% gamma + eta_y^2 - var_y
alpha <- c(
  (-bb + sqrt(bb^2 - 4 * aa * cc)) / (2 * aa),
  (-bb - sqrt(bb^2 - 4 * aa * cc)) / (2 * aa)
)

# var of Y should be covar ZY, covar XY, covar XZ, plus noise










