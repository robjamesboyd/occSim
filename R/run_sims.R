library(dplyr)

# Set simulation parameters
nSites <- 10000
nYears <- 2
nStrata <- 10
nVisits_biased <- 1
nVisits_SRS <- 1
sample_inclusion_prob_SRS <- 0.01
target_error_biased_sample <- 0.1
beta_year_biased <- 0
beta_year_z_biased <- 1.5
beta_detection <- 1
sd_change_z <- 0.1

# Simulation replicates per scenario
n_reps <- 100

# Create w and z (same across reps if desired)
w <- rnorm(nSites)
z <- rnorm(nSites)

# Function to wrap oneSample call with parameters from a row
run_scenario <- function(row, scenario_id) {
  replicate_out <- lapply(1:n_reps, function(rep) {
    oneSample(
      sample = rep,
      nSites = nSites,
      nYears = nYears,
      nStrata = nStrata,
      z = z,
      w = w,
      alpha_occ = row$alpha_occ,
      beta_z_occ = row$beta_z,
      alpha_biased = row$alpha_sample,
      beta_z_biased = row$beta_z_sample,
      beta_year_biased = beta_year_biased,
      beta_year_z_biased = beta_year_z_biased,
      beta_detection = beta_detection,
      nVisits_biased = nVisits_biased,
      nVisits_SRS = nVisits_SRS,
      sample_inclusion_prob_SRS = sample_inclusion_prob_SRS,
      target_error_biased_sample = target_error_biased_sample,
      rho = row$rho_XZ,
      mean_change_z = row$mean_change_z,
      sd_change_z = sd_change_z
    ) %>%
      mutate(scenario = scenario_id, rep = rep)
  })
  
  write.csv(do.call(rbind, replicate_out), file = paste0(scenario_id, "out.csv"))
  do.call(rbind, replicate_out)
}

# Run for all 64 parameter sets
all_results <- lapply(1:nrow(param_set_sample), function(i) {
  run_scenario(param_set_sample[i, ], scenario_id = i)
})

# Combine all into one data frame
final_results <- do.call(rbind, all_results)

# Save if desired
write.csv(final_results, "simulation_results.csv", row.names = FALSE)
