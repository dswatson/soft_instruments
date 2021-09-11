# Solve for alpha given rho
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
loop_fn <- function(idx_b, alpha_b, z_rho_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 2000, d_z = 4, z_cnt = TRUE, z_rho = z_rho_b,
    rho = rho_b, alpha = alpha_b, r2_x = r2_xb, r2_y = r2_yb, 
    prop_valid = prop_b, idx_b
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
  boot_loop <- function(b) {
    # Draw bootstrap sample
    tmp <- dat[sample.int(nrow(dat), replace = TRUE)]
    # Estimate eta_x
    f1 <- lm(x ~ ., data = select(tmp, -y))
    eta_x <- sqrt(mean((residuals(f1)^2)))
    # Estimate data covariance
    Sigma_z <- cov(select(tmp, starts_with('z')))
    Theta_z <- solve(Sigma_z)
    Sigma_zy <- cov(select(tmp, starts_with('z')), tmp$y)
    Sigma_zx <- cov(select(tmp, starts_with('z')), tmp$x)
    Sigma_yz <- cov(tmp$y, select(tmp, starts_with('z')))
    Sigma_xz <- cov(tmp$x, select(tmp, starts_with('z')))
    sigma_xy <- cov(tmp$x, tmp$y)
    var_x <- var(tmp$x)
    var_y <- var(tmp$y)
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
    prop_valid = prop_b, idx_b
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
  ggplot(aes(rho, alpha)) + 
  geom_point(size = 0.25) + 
  geom_line() + 
  geom_errorbar(aes(ymin = alpha - se, ymax = alpha + se)) +
  geom_hline(aes(yintercept = ACE), color = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = rho), color = 'red', linetype = 'dashed') +
  theme_bw() +
  facet_grid(ace_class ~ rho_class) + 
  labs(x = 'Unobserved Confounding', y = 'Average Causal Effect')

# LOWER ROOT IS RIGHT WHEN CONFOUNDING IS POSITIVE!!!
# UPPER ROOT IS RIGHT WHEN CONFOUNDING IS NEGATIVE!!!




