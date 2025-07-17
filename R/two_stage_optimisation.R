tune_parameters <- function(
    beta_psi = 1,
    beta_z_biased_fixed = 1.5,
    beta_z_occ = -1,
    beta_psi_bounds = c(-15, 0),
    alpha_occ_start = -1,
    beta_year_z_biased_start = 0.2,
    mean_change_z_start = 0.1,
    beta_year_start = 0.1,
    target_mean_occ_year1 = 0.3,
    target_occ_ratio = 0.9,
    target_corr_RY = -0.15,
    target_delta_corr = -0.1,
    target_sample_frac_year1 = 0.1,
    target_sample_frac_year2 = 0.15,
    nSites = 10000,
    nYears = 2,
    sd_change_z = 0.1,
    rho_XZ = 0.5,
    nVisits = 1
) {
  set.seed(1)
  z <- rnorm(nSites)
  
  find_alpha_psi <- function(beta_year, z_mat, beta_psi, target_frac, year_idx) {
    root_fn <- function(alpha) {
      eta <- alpha + beta_psi * z_mat[, year_idx] + beta_year * (year_idx - 1)
      p <- plogis(eta)
      mean(p) - target_frac
    }
    uniroot(root_fn, interval = beta_psi_bounds)$root
  }
  
  objective_fn <- function(par) {
    alpha_occ <- par[1]
    beta_year_z_biased <- par[2]
    mean_change_z <- par[3]
    beta_year <- par[4]
    
    # Create z matrix
    z_mat <- matrix(NA, nrow = nSites, ncol = nYears)
    z_mat[, 1] <- z
    for (t in 2:nYears) {
      z_mat[, t] <- z + rnorm(nSites, mean = mean_change_z, sd = sd_change_z)
    }
    
    # Find alpha_psi using uniroot to match year 1 target
    alpha_psi <- find_alpha_psi(beta_year, z_mat, beta_psi, target_sample_frac_year1, 1)
    
    # Check year 2 fraction
    eta2 <- alpha_psi + beta_psi * z_mat[, 2] + beta_year
    sample_frac2 <- mean(plogis(eta2))
    
    # Simulate occupancy
    occupancy <- simulateOccupancyStatic(
      nYears = nYears,
      nSites = nSites,
      z_mat = z_mat,
      beta0 = alpha_occ,
      beta_z = beta_z_occ,
      plot = FALSE
    )$simulated_data
    
    mean_occ_year1 <- mean(occupancy[, 1])
    occ_ratio <- mean(occupancy[, 2]) / mean(occupancy[, 1])
    
    # Simulate sample inclusion
    S <- simulateSample(
      Y = occupancy,
      nVisits = nVisits,
      z_mat = z_mat,
      alpha_psi = alpha_psi,
      beta_psi = beta_z_biased_fixed,
      beta_year = beta_year,
      beta_year_z = beta_year_z_biased,
      plot = FALSE
    )
    
    R <- S$sampled_once
    Y <- occupancy
    corr_RY1 <- cor(R[, 1], Y[, 1])
    corr_RY2 <- cor(R[, 2], Y[, 2])
    delta_corr <- corr_RY2 - corr_RY1
    
    cat(sprintf(
      "alpha_occ = %.3f, mean_occ = %.3f, occ_ratio = %.3f\n",
      alpha_occ, mean_occ_year1, occ_ratio
    ))
    cat(sprintf(
      "corr_RY1 = %.3f, corr_RY2 = %.3f, delta_corr = %.3f\n",
      corr_RY1, corr_RY2, delta_corr
    ))
    cat(sprintf(
      "sample fractions: year1 = %.3f, year2 = %.3f\n",
      target_sample_frac_year1, sample_frac2
    ))
    
    # Loss
    loss <- (mean_occ_year1 - target_mean_occ_year1)^2 +
      (occ_ratio - target_occ_ratio)^2 +
      2 * (corr_RY1 - target_corr_RY)^2 +
      5 * (delta_corr - target_delta_corr)^2 +
      (sample_frac2 - target_sample_frac_year2)^2 * 100  # extra weight
    
    return(loss)
  }
  
  result <- optim(
    par = c(alpha_occ_start, beta_year_z_biased_start, mean_change_z_start, beta_year_start),
    fn = objective_fn,
    method = "L-BFGS-B",
    lower = c(-5, 0, 0, 0),
    upper = c(4, 10, 1, 4),
    control = list(maxit = 30000, factr = 1e-14, pgtol = 1e-14)
  )
  
  names(result$par) <- c("alpha_occ", "beta_year_z_biased", "mean_change_z", "beta_year")
  return(result)
}


result <- tune_parameters(
  target_mean_occ_year1 = 0.05,
  target_occ_ratio = 0.9,
  target_corr_RY = -0.1,
  target_delta_corr = -0.05,
  target_sample_frac_year1 = 0.1,
  target_sample_frac_year2 = 0.15
)

print(result$par)
