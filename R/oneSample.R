oneSample <- function(sample,
                      nSites,
                      nYears,
                      z,
                      w,
                      alpha_initial,
                      beta_initial,
                      alpha_colonisation,
                      beta_colonisation,
                      alpha_extinction,
                      beta_extinction,
                      target_f_biased_sample,
                      beta_sample_inclusion,
                      beta_detection, 
                      nVisits_biased,
                      nVisits_SRS,
                      sample_inclusion_prob_SRS,
                      target_error_biased_sample) {
  
  if (sample==1) {
    occupancy_matrix <- simulateOccupancy(
      nYears = nYears,
      nSites = nSites,
      explanatoryVar = z,
      discreteExplanatoryVar = FALSE,
      beta0_extinction = alpha_extinction,
      beta0_colonization = alpha_colonisation,
      beta0_initial = alpha_initial,
      beta_extinction = beta_extinction,
      beta_colonization = beta_colonisation,
      beta_initial = beta_initial,
      randomInitial = FALSE,
      randomInitialPsi = 0.5,
      nBins = 5,
      plot = FALSE
    )$simulated_data
  }
  # Use simulateOccupancy() to get the true states
 
  
  
  # 2. -------------------------------
  # SIMULATE SAMPLING VISITS
  # --------------------------------
  
  # small simple random sample
  SRS <- simulateSample(
    Y = occupancy_matrix,
    nVisits = nVisits_SRS,
    z = z,
    alpha_psi = log(sample_inclusion_prob_SRS / (1 - sample_inclusion_prob_SRS)),
    beta_psi = 0,
    fixed_panel = TRUE,
    plot = FALSE
  )
  
  # large biased sample
  S <- simulateSample(
    Y = occupancy_matrix,
    nVisits = nVisits_biased,
    z = z,
    alpha_psi = find_sampling_intercept(target_f = target_f_biased_sample, Z=z, V=nVisits_biased,beta=beta_sample_inclusion),
    beta_psi = beta_sample_inclusion,
    plot = FALSE
  )
  
  # 3. -------------------------------
  # SIMULATE DETECTIONS CONDITIONAL ON OCCUPANCY + SAMPLING
  # --------------------------------

  # detection for srs
  D_SRS <- simulateDetection(
    Y = occupancy_matrix,
    S = SRS$S,
    w = w,
    alpha_p = log(0.99/(1-0.99)), # essentially perfect detection
    beta_p = 0,
    plot = FALSE
  )
  D_SRS
  
  # detection for biased sample
  D <- simulateDetection(
    Y = occupancy_matrix,
    S = S$S,
    w = w,
    alpha_p = find_detection_intercept(target_error=target_error_biased_sample, w=w, V=nVisits_biased,beta=beta_detection),
    beta_p = beta_detection,
    plot = FALSE
  )
  
  # Helper function: compute observed detection matrix from 3D array
  detection_matrix <- function(D_array) {
    apply(D_array, c(1, 2), function(x) as.integer(any(x == 1)))
  }
  
  # detected occupancy at SRS sites and biased sample sites
  obs_SRS <- detection_matrix(D_SRS$D)  # [nSites x nYears]
  obs_biased <- detection_matrix(D$D)
  
  sampled_SRS <- apply(SRS$S, c(1, 2), max)
  sampled_biased <- apply(S$S, c(1, 2), max)
  
  # Identify sites sampled in both years (i.e. sampled in all columns)
  sampled_both_years_SRS <- rowSums(sampled_SRS) == nYears
  sampled_both_years_biased <- rowSums(sampled_biased) == nYears
  
  combined_obs <- pmax(obs_SRS, obs_biased)  # at least one detected
  combined_sampled <- pmax(sampled_SRS, sampled_biased)
  
  sampled_both_years_combined <- rowSums(combined_sampled) == nYears
  
  mean_obs_SRS <- sapply(1:nYears, function(t) {
    idx <- which(sampled_both_years_SRS)
    mean(obs_SRS[idx, t])
  })
  
  mean_obs_biased <- sapply(1:nYears, function(t) {
    idx <- which(sampled_both_years_biased)
    mean(obs_biased[idx, t])
  })
  
  mean_obs_combined <- sapply(1:nYears, function(t) {
    idx <- which(sampled_both_years_combined)
    mean(combined_obs[idx, t])
  })
  
  ## get true mean occupancy at sites sampled in both years
  
  # Mean true occupancy for those sites, per year
  true_occ_SRS <- colMeans(occupancy_matrix[sampled_both_years_SRS, ])
  true_occ_biased <- colMeans(occupancy_matrix[sampled_both_years_biased, ])
  true_occ_combined <- colMeans(occupancy_matrix[sampled_both_years_combined, ])
  
  actual_error_SRS <- mean_obs_SRS - colMeans(occupancy_matrix)
  actual_error_biased <- mean_obs_biased - colMeans(occupancy_matrix)
  actual_error_combined <- mean_obs_combined - colMeans(occupancy_matrix)
  
  coverage_error_SRS <- true_occ_SRS - colMeans(occupancy_matrix)
  coverage_error_biased <- true_occ_biased - colMeans(occupancy_matrix)
  coverage_error_combined <- true_occ_combined - colMeans(occupancy_matrix)
  
  prediction_error_SRS <- mean_obs_SRS - true_occ_SRS
  prediction_error_biased <- mean_obs_biased - true_occ_biased
  prediction_error_combined <- mean_obs_combined - true_occ_combined
  
  return(data.frame(year=1:nYears,
    coverage_error_SRS=coverage_error_SRS,
    coverage_error_biased=coverage_error_biased,
    coverage_error_combined=coverage_error_combined,
    actual_error_SRS=actual_error_SRS,
    actual_error_biased=actual_error_biased,
    actual_error_combined=actual_error_combined,
    prediction_error_SRS=prediction_error_SRS,
    prediction_error_biased=prediction_error_biased,
    prediction_error_combined=prediction_error_combined,
    true_sampled_occ_SRS=true_occ_SRS,
    true_sampled_occ_biased=true_occ_biased,
    true_sampled_occ_combined=true_occ_combined,
    est_SRS=mean_obs_SRS,
    est_biased=mean_obs_biased,
    est_combined=mean_obs_combined,
    population_mean=colMeans(occupancy_matrix),
    sample = sample))
  
}
