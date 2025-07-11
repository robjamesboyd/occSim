simulate_scenarios_from_priors <- function(nDraws = 100, nSites = 10000, nYears = 2) {
  # Storage for each draw
  results <- data.frame()
  
  for (i in 1:nDraws) {
    # Sample parameters from priors
    alpha_occ <- runif(1, -2.5, -0.5)
    beta_z <- runif(1, -2, 2)
    mean_change_z <- runif(1, -0.2, 0.2)
    alpha_sample <- runif(1, -4, -1)
    beta_z_sample <- runif(1, 0.5, 2)
    rho_XZ <- runif(1, 0.3, 0.9)
    
    # Generate Z and X
    z <- rnorm(nSites)
    z_mat <- matrix(NA, nSites, nYears)
    x_mat <- matrix(NA, nSites, nYears)
    z_mat[, 1] <- z
    x_mat[, 1] <- rho_XZ * z + sqrt(1 - rho_XZ^2) * rnorm(nSites)
    
    for (t in 2:nYears) {
      z_mat[, t] <- z + rnorm(nSites, mean = mean_change_z, sd = 0.1)
      x_mat[, t] <- rho_XZ * z_mat[, t] + sqrt(1 - rho_XZ^2) * rnorm(nSites)
    }
    
    # Simulate occupancy
    logit_psi <- alpha_occ + beta_z * z_mat
    psi <- plogis(logit_psi)
    Y <- matrix(rbinom(nSites * nYears, size = 1, prob = psi), nSites, nYears)
    
    # Simulate sample inclusion
    logit_pi <- alpha_sample + beta_z_sample * z_mat
    pi <- plogis(logit_pi)
    S <- matrix(rbinom(nSites * nYears, 1, pi), nSites, nYears)
    
    # Compute summary stats
    occ_year1 <- mean(Y[, 1])
    occ_diff <- mean(Y[, 2]) - mean(Y[, 1])
    occ_logit_diff <- logit(mean(Y[, 2])) - logit(mean(Y[, 1]))
    bias_year1 <- cor(S[, 1], Y[, 1])
    bias_year2 <- cor(S[, 2], Y[, 2])
    bias_change <- bias_year2 - bias_year1
    sample_size1 <- mean(S[, 1])
    sample_size2 <- mean(S[, 2])
    sample_size_change <- sample_size2 - sample_size1
    
    # Store both parameters and summaries
    results <- rbind(results, data.frame(
      alpha_occ, beta_z, mean_change_z, alpha_sample, beta_z_sample, rho_XZ,
      occ_year1, occ_diff, occ_logit_diff,
      bias_year1, bias_year2, bias_change,
      sample_size1, sample_size2, sample_size_change
    ))
  }
  
  return(results)
}

# Helper function
logit <- function(p) log(p / (1 - p))

# Example usage:
set.seed(123)
scenario_data <- simulate_scenarios_from_priors(nDraws = 10000)
str(scenario_data)
# Required libraries
library(GGally)
library(ggplot2)
library(dplyr)

# Load your data
# Replace with your actual path or data frame
df <- scenario_data

# Optionally sample a subset if plotting all 10,000 is too slow
set.seed(123)
df_sample <- df %>% sample_n(10000)

# Select the columns for visualization
vars_to_plot <- c(
  "occ_year1", "occ_diff", "occ_logit_diff",
  "bias_year1", "bias_year2", "bias_change",
  "sample_size1", "sample_size2", "sample_size_change"
)

# Plot using ggpairs
ggpairs(df_sample[, vars_to_plot]) +
  theme_minimal() +
  ggtitle("Pairwise Relationships Between Summary Statistics")


set.seed(123)  # for reproducibility

library(dplyr)

# Number of bins per stratification variable
n_bins <- 2

# Variables to stratify by
strat_vars <- c("occ_year1", "occ_logit_diff", "bias_year1", 
                "bias_change", "sample_size1", "sample_size_change")

# Create binned versions of each stratification variable
scenario_stratified <- scenario_data %>%
  mutate(across(all_of(strat_vars),
                ~ cut(.x, breaks = quantile(.x, probs = seq(0, 1, length.out = n_bins + 1), 
                                            na.rm = TRUE), include.lowest = TRUE, labels = FALSE),
                .names = "{.col}_bin"))

# Create a unique stratum ID by crossing all bin combinations
scenario_stratified <- scenario_stratified %>%
  mutate(stratum_id = interaction(across(ends_with("_bin")), drop = TRUE))

# Take exactly one sample per stratum
sampled_params <- scenario_stratified %>%
  group_by(stratum_id) %>%
  sample_n(size = 1) %>%
  ungroup()

# Select only the parameter columns
param_set_sample <- sampled_params %>%
  select(alpha_occ, beta_z, mean_change_z, alpha_sample, beta_z_sample, rho_XZ)

# Check result
nrow(param_set_sample)  # should be exactly 64
head(param_set_sample)

library(ggplot2)
library(tidyr)

# Add "source" labels
full_data_long <- scenario_data %>%
  mutate(source = "All") %>%
  select(all_of(strat_vars), source) %>%
  pivot_longer(-source)

sampled_data_long <- sampled_params %>%
  mutate(source = "Sampled") %>%
  select(all_of(strat_vars), source) %>%
  pivot_longer(-source)

combined_long <- bind_rows(full_data_long, sampled_data_long)

# Density plots for each stratification variable
ggplot(combined_long, aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ name, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = "Coverage of Stratification Variables",
       x = "Value", y = "Density")
