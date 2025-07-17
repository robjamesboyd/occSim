tune_parameters <- function(
    alpha_start = -0.2,
    beta_year_z_biased_start = -0.05,
    mean_change_z_start = 0.15,
    alpha_psi_start = 0,
    beta_year_start = 1,
    
    target_mean_occ_year1 = 0.3, 
    target_occ_ratio = 0.9,
    target_corr_RY = -0.15,
    target_delta_corr = -0.1,
    target_sample_frac_yr1 = 0.10,
    target_sample_frac_yr2 = 0.15,
    
    nSites = 10000,
    nYears = 2,
    beta_z_occ = -1,
    beta_z_biased_fixed = 1.5,
    sd_change_z = 0.1,
    rho_XZ = 0.5,
    beta_psi = 1,
    nVisits = 1
) {
  z <- rnorm(nSites)  # baseline Z
  
  objective_fn <- function(par) {
    alpha_occ <- par[1]
    beta_year_z_biased <- par[2]
    mean_change_z <- par[3]
    alpha_psi <- par[4]
    beta_year <- par[5]
    
    # Generate Z matrix
    z_mat <- matrix(NA, nSites, nYears)
    z_mat[, 1] <- z
    for (t in 2:nYears) {
      z_mat[, t] <- z + rnorm(nSites, mean = mean_change_z, sd = sd_change_z)
    }
    
    # Simulate occupancy
    occupancy <- simulateOccupancyStatic(
      nYears = nYears, nSites = nSites,
      z_mat = z_mat,
      beta0 = alpha_occ,
      beta_z = beta_z_occ,
      plot = FALSE
    )$simulated_data
    
    mean_occ_year1 <- mean(occupancy[, 1])
    occ_ratio <- mean(occupancy[, 2]) / mean(occupancy[, 1])
    
    # Simulate sampling
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
    
    # Correlations
    corr_RY1 <- cor(R[, 1], Y[, 1])
    corr_RY2 <- cor(R[, 2], Y[, 2])
    delta_corr <- corr_RY2 - corr_RY1
    
    # Sample fractions
    sample_frac_yr1 <- mean(R[, 1])
    sample_frac_yr2 <- mean(R[, 2])
    
    cat("alpha_occ =", alpha_occ, "mean_occ_year1 =", round(mean_occ_year1, 3), "\n")
    cat("occ_ratio =", round(occ_ratio, 3), "delta_corr =", round(delta_corr, 3), "\n")
    cat("corr_RY1 =", corr_RY1, "corr_RY2 =", corr_RY2, "\n")
    cat("sample fractions =", round(sample_frac_yr1, 3), round(sample_frac_yr2, 3), "\n")
    
    # Loss function
    loss <- 1.5 * (mean_occ_year1 - target_mean_occ_year1)^2 +
      1.5 * (occ_ratio - target_occ_ratio)^2 +
      2.0 * (corr_RY1 - target_corr_RY)^2 +
      3.0 * (delta_corr - target_delta_corr)^2 +
      100.0 * (sample_frac_yr1 - target_sample_frac_yr1)^2 +
      100.0 * (sample_frac_yr2 - target_sample_frac_yr2)^2
    
    
    return(loss)
    
    
  }
  
  result <- optim(
    par = c(alpha_start, beta_year_z_biased_start, mean_change_z_start, alpha_psi_start, beta_year_start),
    fn = objective_fn,
    method = "L-BFGS-B",
    control = list(maxit = 25000, factr = 1e-14, pgtol = 1e-14),
    lower = c(-7, 0, 0, -10, 0),
    upper = c(-1, 5, 1, -2, 2)
  )
  
  names(result$par) <- c("alpha_occ", "beta_year_z_biased", "mean_change_z", "alpha_psi", "beta_year")
  return(result)
}



result <- tune_parameters(
  target_mean_occ_year1 = 0.3,
  target_occ_ratio = 0.9,
  target_corr_RY = -0.15,
  target_delta_corr = -0.1,
  target_sample_frac_yr1 = 0.10,
  target_sample_frac_yr2 = 0.15
)

print(result$par)


