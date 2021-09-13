# Solve for alpha given rho
soft_iv <- function(dat, rho) {
  # Estimate eta_x
  f1 <- lm(x ~ ., data = select(dat, -y))
  eta_x <- sqrt(mean((residuals(f1)^2)))
  # Estimate data covariance
  d_z <- sum(grepl('z', colnames(dat)))
  Sigma <- cov(dat)
  Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
  Theta_z <- solve(Sigma_z)
  Theta_z2 <- Theta_z %*% Theta_z
  Sigma_zy <- matrix(Sigma[seq_len(d_z), ncol(Sigma)], ncol = 1)
  Sigma_yz <- t(Sigma_zy)
  Sigma_zx <- matrix(Sigma[seq_len(d_z), ncol(Sigma) - 1], ncol = 1)
  Sigma_xz <- t(Sigma_zx)
  sigma_xy <- Sigma[ncol(Sigma), ncol(Sigma) - 1]
  var_x <- Sigma[ncol(Sigma) - 1, ncol(Sigma) - 1]
  var_y <- Sigma[ncol(Sigma), ncol(Sigma)]
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
loop_fn <- function(idx_b, alpha_b, z_rho_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 2000, d_z = 4, z_cnt = TRUE, z_rho = z_rho_b,
    rho = rho_b, alpha = alpha_b, r2_x = r2_xb, r2_y = r2_yb, 
    pr_valid = prop_b, idx_b
  )
  alphas <- soft_iv(tmp, rho_b)
  out <- data.table(
    'alpha_lo' = alphas$alpha_lo, 'alpha_hi' = alphas$alpha_hi,
    'ACE' = alpha_b, 'z_rho' = z_rho_b, 'rho' = rho_b,  
    'r2_x' = r2_xb, 'r2_y' = r2_yb, 'pr_valid' = prop_b
  )
}
# Execute in parallel
df <- foreach(idxs = seq_len(500), .combine = rbind) %:%
  foreach(alphas = c(-1, 1), .combine = rbind) %:%
  foreach(z_rhos = c(0, 0.5), .combine = rbind) %:%
  foreach(rhos = c(-0.5, 0.5), .combine = rbind) %:%
  foreach(props = c(0, 0.5), .combine = rbind)  %dopar% 
  loop_fn(idxs, alphas, z_rhos, rhos, 0.5, 0.5, props)

# All about varying the sign of rho!
df <- foreach(idxs = seq_len(500), .combine = rbind) %:%
  foreach(alphas = c(-1, 1), .combine = rbind) %:%
  foreach(rhos = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), .combine = rbind) %dopar% 
  loop_fn(idxs, alphas, 0, rhos, 0.5, 0.5, 0)

# In histogram form (no r2_x or r2_y)
df[pr_valid == 0 & ACE == 1] %>%
  ggplot(aes(alpha_lo)) + 
  geom_histogram(bins = 50) +  
  geom_vline(xintercept = 1, color = 'red', linetype = 'dashed') +
  theme_bw() + 
  facet_grid(z_rho ~ rho)
df[pr_valid == 1] %>%
  ggplot(aes(alpha_lo)) + 
  geom_histogram(bins = 50) +  
  geom_vline(xintercept = 1, color = 'red', linetype = 'dashed') +
  theme_bw() + 
  facet_grid(z_rho ~ rho)

# Boxplots
df %>%
  rename(lower = alpha_lo, upper = alpha_hi) %>%
  filter(z_rho == 0) %>%
  pivot_longer(c('lower', 'upper'), names_to = 'root', values_to = 'alpha') %>%
  ggplot(aes(as.factor(rho), alpha, fill = root)) + 
  geom_boxplot() + 
  scale_fill_d3() + 
  geom_hline(aes(yintercept = ACE), color = 'red', linetype = 'dashed') +
  theme_bw() +
  facet_grid(ace2 ~ r2_x) + 
  labs(x = 'Unobserved Confounding', y = 'Average Causal Effect')

df %>%
  rename(lower = alpha_lo, upper = alpha_hi) %>%
  mutate(ace_class = if_else(ACE > 0, 'ACE = 1', 'ACE = -1'),
         rho_class = if_else(rho > 0, 'Positive Confounding', 
                             'Negative Confounding')) %>%
  pivot_longer(c('lower', 'upper'), names_to = 'root', values_to = 'alpha') %>%
  ggplot(aes(as.factor(rho), alpha, color = root)) + 
  geom_boxplot() + 
  scale_color_d3() + 
  geom_hline(aes(yintercept = ACE), color = 'red', linetype = 'dashed') +
  theme_bw() +
  facet_grid(ace_class ~ rho_class) + 
  labs(x = 'Unobserved Confounding', y = 'Average Causal Effect')


# Bootstrapping for better results, standard errors, etc.
soft_iv <- function(dat, rho, n_boot) {
  d_z <- sum(grepl('z', colnames(dat)))
  boot_loop <- function(b) {
    # Draw bootstrap sample
    tmp <- dat[sample.int(nrow(dat), replace = TRUE)]
    # Estimate eta_x
    f1 <- lm(x ~ ., data = select(tmp, -y))
    eta_x <- sqrt(mean((residuals(f1)^2)))
    # Estimate data covariance
    Sigma <- cov(tmp)
    Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
    Theta_z <- solve(Sigma_z)
    Theta_z2 <- Theta_z %*% Theta_z
    Sigma_zy <- matrix(Sigma[seq_len(d_z), ncol(Sigma)], ncol = 1)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[seq_len(d_z), ncol(Sigma) - 1], ncol = 1)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[ncol(Sigma), ncol(Sigma) - 1]
    var_x <- Sigma[ncol(Sigma) - 1, ncol(Sigma) - 1]
    var_y <- Sigma[ncol(Sigma), ncol(Sigma)]
    # Simplify
    aa <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    bb <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    cc <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    dd <- (bb - var_x) * (1 + ((bb - var_x) / (eta_x^2 * rho^2)))
    ee <- -(aa - sigma_xy) * (1 + ((bb - var_x) / (eta_x^2 * rho^2)))
    ff <- ((aa - sigma_xy)^2 / (eta_x^2 * rho^2)) + cc - var_y
    # Return roots
    out <- data.frame(
      lo = (-ee - sqrt(ee^2 - dd * ff)) / dd,
      hi = (-ee + sqrt(ee^2 - dd * ff)) / dd 
    )
    return(out)
  }
  # Execute, export
  alphas <- foreach(i = seq_len(n_boot), .combine = rbind) %do% boot_loop(i)
  out <- data.frame(
      'root' = c('lower', 'upper'), 
     'alpha' = c(mean(alphas$lo), mean(alphas$hi)),
        'se' = c(sd(alphas$lo), sd(alphas$hi)),
    'rho_in' = rho
  )
  return(out)
}
# Loop over rhos in parallel
loop_fn <- function(idx_b, alpha_b, z_rho_b, rho_b, r2_xb, r2_yb, prop_b) {
  # Draw dataset
  tmp <- sim_dat(
    n = 1000, d_z = 4, z_cnt = TRUE, z_rho = z_rho_b,
    rho = rho_b, alpha = alpha_b, r2_x = r2_xb, r2_y = r2_yb, 
    pr_valid = prop_b, idx_b
  )
  # Loop over rhos
  out <- foreach(rhos = seq(-0.99, 0.99, 0.02), .combine = rbind) %do%
    soft_iv(tmp, rhos, n_boot = 100)
  out <- out %>%
    mutate('ACE' = alpha_b, 'z_rho' = z_rho_b, 'rho' = rho_b,  
           'r2_x' = r2_xb, 'r2_y' = r2_yb, 'pr_valid' = prop_b)
  return(out)
}
df <- foreach(alphas = c(-1, 1), .combine = rbind) %:%
  foreach(rhos = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), .combine = rbind) %dopar% 
  loop_fn(1, alphas, 0, rhos, 0.5, 0.5, 0)

# Plot
df %>%
  mutate(ace_class = if_else(ACE > 0, 'ACE = 1', 'ACE = -1'),
         rho_class = if_else(rho > 0, 'Positive Confounding', 
                             'Negative Confounding')) %>%
  filter(rho > 0) %>%
  ggplot(aes(rho_in, alpha, color = root)) + 
  geom_point(size = 0.25) + 
  geom_line() + 
  geom_errorbar(aes(ymin = alpha - se, ymax = alpha + se)) +
  geom_hline(aes(yintercept = ACE), color = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = rho), color = 'red', linetype = 'dashed') +
  scale_color_d3() +
  theme_bw() +
  facet_grid(ace_class ~ rho) + 
  labs(x = 'Unobserved Confounding', y = 'Average Causal Effect')

# LOWER ROOT IS RIGHT WHEN CONFOUNDING IS POSITIVE!!!
# UPPER ROOT IS RIGHT WHEN CONFOUNDING IS NEGATIVE!!!
# Upshot: there's always one alpha per rho. If rho < 0, return upper root;
# else, return lower root.
# Means that alpha is a monotonically decreasing function of rho...hm....






loop_fn <- function(idx_b, alpha_b, z_rho_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 2000, d_z = 4, z_cnt = TRUE, z_rho = z_rho_b,
    rho = rho_b, alpha = alpha_b, r2_x = r2_xb, r2_y = r2_yb, 
    pr_valid = prop_b, idx_b
  )
  f_x <- lm(x ~ ., data = select(tmp, -y))
  tmp <- tmp %>% mutate(x = fitted(f_x))
  f_y <- lm(y ~ ., data = tmp)
  out <- data.table(
    'rho' = rho_b, 'rho_hat' = cor(residuals(f_x), residuals(f_y)),
    'alpha' = alpha_b, 'z_rho' = z_rho_b, 
    'r2_x' = r2_xb, 'r2_y' = r2_yb, 'pr_valid' = prop_b
  )
  return(out)
}
# Execute in parallel
df <- foreach(idxs = seq_len(100), .combine = rbind) %:%
  foreach(alphas = c(-1, 1), .combine = rbind) %:%
  foreach(z_rhos = c(0, 0.5), .combine = rbind) %:%
  foreach(rhos = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), .combine = rbind) %:%
  foreach(r2_xs = c(0.25, 0.75), .combine = rbind) %:%
  foreach(r2_ys = c(0.25, 0.75), .combine = rbind) %:%
  foreach(props = c(0, 0.5), .combine = rbind)  %dopar% 
  loop_fn(idxs, alphas, z_rhos, rhos, r2_xs, r2_ys, props)

df %>%
  filter(z_rho == 0, pr_valid == 0) %>%
  ggplot(aes(rho, rho_hat)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = 'red') + 
  theme_bw() + 
  facet_grid(r2_x ~ r2_y)

df %>%
  filter(z_rho == 0.5, pr_valid == 0) %>%
  ggplot(aes(rho, rho_hat, color = as.factor(alpha))) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = 'red') + 
  scale_color_d3() +
  theme_bw() + 
  facet_grid(r2_x ~ r2_y)

df %>%
  filter(z_rho == 0, pr_valid == 0.5) %>%
  ggplot(aes(rho, rho_hat)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = 'red') + 
  theme_bw() + 
  facet_grid(r2_x ~ r2_y)

df %>%
  filter(z_rho == 0.5, pr_valid == 0.5) %>%
  ggplot(aes(rho, rho_hat)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = 'red') + 
  theme_bw() + 
  facet_grid(r2_x ~ r2_y)




# For selecting tau: 
# regress y ~ Z to get an upper bound (this will overestimate the Z-signal)
# regress y ~ Z + X to get a lower bound (this will underestimate the Z signal)
# True L2 norm of gamma should be somewhere in between
# You probably want to be conservative and underestimate tau...


ggplot(df, aes(rho_in, alpha)) + 
  geom_point(size = 0.25) + 
  geom_line() + 
  geom_errorbar(aes(ymin = alpha - se, ymax = alpha + se)) +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0.5, color = 'red', linetype = 'dashed') +
  theme_bw()




# Can we ever get ee - var_x to flip signs? Nope!
loop_fn <- function(idx_b, alpha_b, z_rho_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 1000, d_z = 4, z_cnt = TRUE, z_rho = z_rho_b,
    rho = rho_b, alpha = alpha_b, r2_x = r2_xb, r2_y = r2_yb, 
    pr_valid = prop_b, idx_b
  )
  Sigma <- cov(tmp)
  Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
  Theta_z <- solve(Sigma_z)
  Sigma_zx <- matrix(Sigma[seq_len(d_z), ncol(Sigma) - 1], ncol = 1)
  Sigma_xz <- t(Sigma_zx)
  var_x <- Sigma[ncol(Sigma) - 1, ncol(Sigma) - 1]
  ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
  out <- data.table(
    'ACE' = alpha_b, 'z_rho' = z_rho_b, 'rho' = rho_b,  
    'r2_x' = r2_xb, 'r2_y' = r2_yb, 'pr_valid' = prop_b, 
    'theta' = ee - var_x
  )
  return(out)
}
# Execute in parallel
df <- foreach(idxs = seq_len(100), .combine = rbind) %:%
  foreach(alphas = c(-1, 1), .combine = rbind) %:%
  foreach(z_rhos = c(-0.5, 0, 0.5), .combine = rbind) %:%
  foreach(rhos = c(-0.5, 0.5), .combine = rbind) %:%
  foreach(props = c(0, 0.5), .combine = rbind) %:%
  foreach(r2_xs = c(0.25, 0.75), .combine = rbind) %:% 
  foreach(r2_ys = c(0.25, 0.75), .combine = rbind) %dopar% 
  loop_fn(idxs, alphas, z_rhos, rhos, r2_xs, r2_ys, props)

# Do we believe that causal effects are a monotonically decreasing function of
# unobserved confounding?



