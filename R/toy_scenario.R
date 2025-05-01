# Load or source your three functions here:
# simulateOccupancy()
# simulateSample()
# simulateDetection()

# Set seed for reproducibility
set.seed(42)

# 1. -------------------------------
# SIMULATE TRUE OCCUPANCY STATES
# --------------------------------
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

per_visit_prob <- 0.01 / nVisits

# small simple random sample
SRS <- simulateSample(
  Y = occupancy_matrix,
  nVisits = nVisits_SRS,
  z = z,
  alpha_psi = log(per_visit_prob/(1-per_visit_prob)),
  beta_psi = 0
)
S$summary_df

nVists <- 3

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

