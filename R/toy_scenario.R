# Load or source your three functions here:
# simulateOccupancy()
# simulateSample()
# simulateDetection()

# Set seed for reproducibility
set.seed(42)

# 1. -------------------------------
# SIMULATE TRUE OCCUPANCY STATES
# --------------------------------

oneSample <- function(nSites,
                      nYears,
                      z,
                      w,
                      alpha_initial_extinction,
                      beta_initial_extinction,
                      alpha_initial_colonisation,
                      beta_initial_colonisation,
                      alpha_colonisation,
                      beta_colonisation,
                      alpha_extinction,
                      beta_extinction,
                      alpha_sample_inclusion,
                      beta_sample_inclusion,
                      alpha_detection,
                      beta_detection, 
                      nVisits_biased,
                      nVisits_SRS,
                      per_visit_detection_prob) {
  
  
  # Use simulateOccupancy() to get the true states
  occupancy_matrix <- simulateOccupancy(
    nYears = nYears,
    nSites = nSites,
    explanatoryVar = z,
    discreteExplanatoryVar = FALSE,
    beta0_extinction = -3,
    beta0_colonization = -7,
    beta0_initial = NULL,
    beta_extinction = 2,
    beta_colonization = 2,
    beta_initial = NULL,
    randomInitial = TRUE,
    randomInitialPsi = 0.5,
    nBins = 5
  )$simulated_data
  
  
  # 2. -------------------------------
  # SIMULATE SAMPLING VISITS
  # --------------------------------
  nVisits_SRS <- 1
  
  per_visit_prob <- 0.01 / nVisits_SRS
  
  # small simple random sample
  SRS <- simulateSample(
    Y = occupancy_matrix,
    nVisits = nVisits_SRS,
    z = z,
    alpha_psi = log(per_visit_prob/(1-per_visit_prob)),
    beta_psi = 0
  )
  S$summary_df
  
  nVisits <- 3
  
  # large biased sample
  S <- simulateSample(
    Y = occupancy_matrix,
    nVisits = nVisits,
    z = z,
    alpha_psi = find_sampling_intercept(target_f=0.1, Z=z, V=nVisits,beta=3),
    beta_psi = 3
  )
  S$summary_df
  
  # 3. -------------------------------
  # SIMULATE DETECTIONS CONDITIONAL ON OCCUPANCY + SAMPLING
  # --------------------------------
  # Use a different covariate for detection probability, or reuse z
  w <- rnorm(nSites)  # could represent site accessibility or observer skill
  # detection for srs
  D_SRS <- simulateDetection(
    Y = occupancy_matrix,
    S = SRS$S,
    w = w,
    alpha_p = log(0.99/(1-0.99)),
    beta_p = 0
  )
  D_SRS
  
  # detection for biased sample
  D <- simulateDetection(
    Y = occupancy_matrix,
    S = S$S,
    w = w,
    alpha_p = find_detection_intercept(target_error=0.5, w=w, V=nVisits,beta=2),
    beta_p = 2
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
  
  mean_obs_SRS <- sapply(1:nYears, function(t) {
    idx <- which(sampled_SRS[, t] == 1)
    mean(obs_SRS[idx, t])
  })
  
  mean_obs_biased <- sapply(1:nYears, function(t) {
    idx <- which(sampled_biased[, t] == 1)
    mean(obs_biased[idx, t])
  })
  
  combined_obs <- pmax(obs_SRS, obs_biased)  # at least one detected
  combined_sampled <- pmax(sampled_SRS, sampled_biased)
  
  mean_obs_combined <- sapply(1:nYears, function(t) {
    idx <- which(combined_sampled[, t] == 1)
    mean(combined_obs[idx, t])
  })
  
  ## get true mean occupancy at sites sampled in both years
  
  # Identify sites sampled in both years (i.e. sampled in all columns)
  sampled_both_years_SRS <- rowSums(sampled_SRS) == nYears
  sampled_both_years_biased <- rowSums(sampled_biased) == nYears
  sampled_both_years_combined <- rowSums(combined_sampled) == nYears
  
  # Mean true occupancy for those sites, per year
  true_occ_SRS <- colMeans(occupancy_matrix[sampled_both_years_SRS, ])
  true_occ_biased <- colMeans(occupancy_matrix[sampled_both_years_biased, ])
  true_occ_combined <- colMeans(occupancy_matrix[sampled_both_years_combined, ])
  
  actual_error_SRS <- mean_obs_SRS - colMeans(occupancy_matrix)
  actual_error_biased <- mean_obs_biased - colMeans(occupancy_matrix)
  actual_error_combined <- mean_obs_combined - colMeans(occupancy_matrix)
  
  coverage_error_SRS <- true_occ_SRS - colMeans(occupancy_matrix)
  coverage_error_biased <- true_occ_SRS - colMeans(occupancy_matrix)
  coverage_error_combined <- true_occ_SRS - colMeans(occupancy_matrix)
  
  prediction_error_SRS <- mean_obs_SRS - true_occ_SRS
  prediction_error_biased <- mean_obs_biased - true_occ_SRS
  prediction_error_combined <- mean_obs_combined - true_occ_SRS
  
  return(data.frame(coverage_error_SRS=coverage_error_SRS,
                    coverage_error_biased=coverage_error_biased,
                    coverage_error_combined=coverage_error_combined,
                    actual_error_SRS=coverage_error_SRS,
                    actual_error_biased=coverage_error_biased,
                    actual_error_combined=coverage_error_combined,
                    prediction_error_SRS=coverage_error_SRS,
                    prediction_error_biased=coverage_error_biased,
                    prediction_error_combined=coverage_error_combined))
  
}
                      
nSites <- 10000
nYears <- 2

# Simulate a continuous site-level covariate (z)
z <- runif(nSites, min = 0, max = 1)

# Use simulateOccupancy() to get the true states
occupancy_matrix <- simulateOccupancy(
  nYears = nYears,
  nSites = nSites,
  explanatoryVar = z,
  discreteExplanatoryVar = FALSE,
  beta0_extinction = -3,
  beta0_colonization = -7,
  beta0_initial = -1.000696,
  beta_extinction = 2,
  beta_colonization = 2,
  beta_initial = 2,
  randomInitial = FALSE,
  randomInitialPsi = NA,
  nBins = 5
)$simulated_data


# 2. -------------------------------
# SIMULATE SAMPLING VISITS
# --------------------------------
nVisits_SRS <- 1

per_visit_prob <- 0.01 / nVisits_SRS

# small simple random sample
SRS <- simulateSample(
  Y = occupancy_matrix,
  nVisits = nVisits_SRS,
  z = z,
  alpha_psi = log(per_visit_prob/(1-per_visit_prob)),
  beta_psi = 0
)
S$summary_df

nVisits <- 3

# large biased sample
S <- simulateSample(
  Y = occupancy_matrix,
  nVisits = nVisits,
  z = z,
  alpha_psi = find_sampling_intercept(target_f=0.1, Z=z, V=nVisits,beta=3),
  beta_psi = 3
)
S$summary_df
S$S
# 3. -------------------------------
# SIMULATE DETECTIONS CONDITIONAL ON OCCUPANCY + SAMPLING
# --------------------------------
# Use a different covariate for detection probability, or reuse z
w <- rnorm(nSites)  # could represent site accessibility or observer skill
# detection for srs
D_SRS <- simulateDetection(
  Y = occupancy_matrix,
  S = SRS$S,
  w = w,
  alpha_p = log(0.99/(1-0.99)),
  beta_p = 0
)
D_SRS

# detection for biased sample
D <- simulateDetection(
  Y = occupancy_matrix,
  S = S$S,
  w = w,
  alpha_p = find_detection_intercept(target_error=0.5, w=w, V=nVisits,beta=2),
  beta_p = 2
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

mean_obs_SRS <- sapply(1:nYears, function(t) {
  idx <- which(sampled_SRS[, t] == 1)
  mean(obs_SRS[idx, t])
})

mean_obs_biased <- sapply(1:nYears, function(t) {
  idx <- which(sampled_biased[, t] == 1)
  mean(obs_biased[idx, t])
})

combined_obs <- pmax(obs_SRS, obs_biased)  # at least one detected
combined_sampled <- pmax(sampled_SRS, sampled_biased)

mean_obs_combined <- sapply(1:nYears, function(t) {
  idx <- which(combined_sampled[, t] == 1)
  mean(combined_obs[idx, t])
})

## get true mean occupancy at sites sampled in both years

# Identify sites sampled in both years (i.e. sampled in all columns)
sampled_both_years_SRS <- rowSums(sampled_SRS) == nYears
sampled_both_years_biased <- rowSums(sampled_biased) == nYears
sampled_both_years_combined <- rowSums(combined_sampled) == nYears

# Mean true occupancy for those sites, per year
true_occ_SRS <- colMeans(occupancy_matrix[sampled_both_years_SRS, ])
true_occ_biased <- colMeans(occupancy_matrix[sampled_both_years_biased, ])
true_occ_combined <- colMeans(occupancy_matrix[sampled_both_years_combined, ])

actual_error_SRS <- mean_obs_SRS - colMeans(occupancy_matrix)
actual_error_biased <- mean_obs_biased - colMeans(occupancy_matrix)
actual_error_combined <- mean_obs_combined - colMeans(occupancy_matrix)

coverage_error_SRS <- true_occ_SRS - colMeans(occupancy_matrix)
coverage_error_biased <- true_occ_SRS - colMeans(occupancy_matrix)
coverage_error_combined <- true_occ_SRS - colMeans(occupancy_matrix)

prediction_error_SRS <- mean_obs_SRS - true_occ_SRS
prediction_error_biased <- mean_obs_biased - true_occ_SRS
prediction_error_combined <- mean_obs_combined - true_occ_SRS

return(data.frame(coverage_error_SRS=coverage_error_SRS,
                  coverage_error_biased=coverage_error_biased,
                  coverage_error_combined=coverage_error_combined,
                  actual_error_SRS=coverage_error_SRS,
                  actual_error_biased=coverage_error_biased,
                  actual_error_combined=coverage_error_combined,
                  prediction_error_SRS=coverage_error_SRS,
                  prediction_error_biased=coverage_error_biased,
                  prediction_error_combined=coverage_error_combined))
# 4. -------------------------------
# BASIC CHECKS
# --------------------------------

# Check dimensions
print(dim(occupancy_matrix))  # Should be [nSites, nYears]
print(dim(SRS$S))                 # Should be [nSites, nYears, nVisits]
print(dim(D_SRS$D))                 # Should match S

# Pick one site and inspect
site <- 1
cat("True occupancy for site", site, ":\n")
print(occupancy_matrix[site, ])

cat("\nSampling matrix for site", site, ":\n")
print(S[site, , ])

cat("\nDetection matrix for site", site, ":\n")
print(D[site, , ])

# How many detections were made total?
cat("\nTotal detections across all sites and years:\n")
print(sum(D))

# Confirm that there are no detections where either occupancy or sampling = 0
violations <- 0
for (i in 1:nSites) {
  for (t in 1:nYears) {
    for (v in 1:nVisits) {
      if (D[i, t, v] == 1 && (occupancy_matrix[i, t] == 0 || S[i, t, v] == 0)) {
        violations <- violations + 1
      }
    }
  }
}
cat("\nDetection violations (should be 0):", violations, "\n")

occupancy_matrix


# Simulate fixed occupancy
Y <- matrix(1, nrow = 100, ncol = 5)  # all sites occupied

# Simulate all sites visited always
S <- array(1, dim = c(100, 5, 3))  # all visits occur

# Set detection to be perfect
w <- rep(0, 100)
detect_result <- simulateDetection(Y, S, w, alpha_p = 10, beta_p = 0, plot = TRUE)
detect_result$error_summary


manySamples <- lapply(1:1000, FUN=oneSample,
       nSites = 10000,
          nYears = 2,
          z = z,
          w = w,
          alpha_initial = -1.000696,
          beta_initial = 2,
          alpha_extinction = -3,
          beta_extinction = 2,
          alpha_colonisation = -7,
          beta_colonisation = 2,
          nVisits_SRS = 1,
          nVisits_biased = 3,
          target_f_biased_sample = 0.1,
          target_error_biased_sample = 0.5,
          sample_inclusion_prob_SRS = 0.01,
          beta_sample_inclusion = 3,
          beta_detection = 2)


manySamples <- do.call("rbind", manySamples)

df <- manySamples[manySamples$year == 1, ]

# Bias^2 component
bias_squared <- (mean(df$est_biased) - mean(df$population_mean))^2

# Estimator variance
estimator_variance <- var(df$est_biased)

# Variance of population slice (true mean at sampled sites)
slice_variance <- var(df$true_sampled_occ_biased)

# Covariance (error coupling)
error_covariance <- cov(df$est_biased, df$true_sampled_occ_biased)

# Full identity-based MSE decomposition
MSE_decomposed <- bias_squared + estimator_variance + slice_variance - 2 * error_covariance

# Direct empirical MSE
MSE_empirical <- mean((df$est_biased - df$population_mean)^2)

print(c(MSE_decomposed = MSE_decomposed, MSE_empirical = MSE_empirical))

library(dplyr)
library(tidyr)
library(ggplot2)

# Subset to year 1
df <- manySamples[manySamples$year == 1, ]

# Function to compute components for a given scenario
mse_components <- function(est, true_slice, pop_mean) {
  bias2 <- (mean(est) - mean(pop_mean))^2
  var_est <- var(est)
  var_slice <- var(true_slice)
  cov_term <- cov(est, true_slice)
  data.frame(
    Bias_squared = bias2,
    Estimator_variance = var_est,
    Coverage_variance = var_slice,
    Error_covariance = -2 * cov_term
  )
}

# Compute for each scenario
srs <- mse_components(df$est_SRS, df$true_sampled_occ_SRS, df$population_mean)
biased <- mse_components(df$est_biased, df$true_sampled_occ_biased, df$population_mean)
combined <- mse_components(df$est_combined, df$true_sampled_occ_combined, df$population_mean)

# Combine into one data frame
mse_df <- bind_rows(srs, biased, combined)
mse_df$Scenario <- c("Random_sample", "Unstructured", "Combined")

# Reshape for plotting
mse_long <- mse_df %>%
  pivot_longer(cols = -Scenario, names_to = "Component", values_to = "Value")

# Order components (optional aesthetic)
mse_long$Component <- factor(mse_long$Component,
                             levels = c("Bias_squared", "Estimator_variance", "Coverage_variance", "Error_covariance"))

# Plot
png("plot.png", width = 7, height = 5, units = "in", res = 500)
print(ggplot(mse_long, aes(x = Scenario, y = Value, fill = Component)) +
  geom_bar(stat = "identity") +
  labs(title = "MSE Decomposition by Scenario",
       y = "MSE", x = "Scenario", fill = "Component") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal())
dev.off()

