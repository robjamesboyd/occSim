simulate_scenarios_from_priors <- function(nDraws = 100, nSites = 10000, nYears = 2) {
  results <- data.frame()
  
  for (i in 1:nDraws) {
    # Sample from priors
    alpha_occ <- runif(1, -2.5, 0)
    beta_z <- -1  # fixed
    mean_change_z <- runif(1, -0.3, 0.3)
    alpha_sample <- runif(1, -7, -4)
    beta_z_sample <- runif(1, -2, 2) 
    beta_year_sample <- 0.2  # fixed
    beta_year_z_sample <- runif(1, -2, 2)  # interaction
    rho_XZ <- runif(1, 0.01, 0.5)

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
    year_mat <- matrix(rep(0:(nYears - 1), each = nSites), nrow = nSites)
    linpred <- alpha_sample + 
      beta_z_sample * z_mat + 
      beta_year_sample * year_mat + 
      beta_year_z_sample * z_mat * year_mat
    pi <- plogis(linpred)
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
    
    # Store everything
    results <- rbind(results, data.frame(
      alpha_occ, beta_z, mean_change_z, 
      alpha_sample, beta_z_sample, beta_year_sample, beta_year_z_sample, rho_XZ,
      occ_year1, occ_diff, occ_logit_diff,
      bias_year1, bias_year2, bias_change,
      sample_size1, sample_size2, sample_size_change
    ))
  }
  
  return(results)
}

# Helper
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
  "sample_size1", "sample_size2", "sample_size_change", "rho_XZ"
)

png("pairs.png", width=10, height= 10, res=500, units = "in")
# Plot using ggpairs
ggpairs(df_sample[, vars_to_plot]) +
  theme_minimal() +
  ggtitle("Pairwise Relationships Between Summary Statistics")
dev.off()

set.seed(123)  # for reproducibility

library(dplyr)

# Number of bins per stratification variable
n_bins <- 2

# Variables to stratify by
strat_vars <- c("occ_year1", "occ_diff", "bias_year1", 
                "bias_change", "sample_size1", "rho_XZ")
names(scenario_data)
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
  dplyr::select(alpha_occ, beta_z, mean_change_z, alpha_sample, beta_z_sample, rho_XZ)

# Check result
nrow(param_set_sample)  # should be exactly 64
head(param_set_sample)
str(param_set_sample)
write.csv(param_set_sample, file = "W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/MAMBO/parameter_sets.csv")
library(ggplot2)
library(tidyr)

# Add "source" labels
full_data_long <- scenario_data %>%
  mutate(source = "All") %>%
  dplyr::select(all_of(strat_vars), source) %>%
  pivot_longer(-source)

sampled_data_long <- sampled_params %>%
  mutate(source = "Sampled") %>%
  dplyr::select(all_of(strat_vars), source) %>%
  pivot_longer(-source)

combined_long <- bind_rows(full_data_long, sampled_data_long)

png("marginal.png", width=10, height= 5, res=500, units = "in")
# Density plots for each stratification variable
ggplot(combined_long, aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ name, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = "Coverage of Stratification Variables",
       x = "Value", y = "Density")
dev.off()
library(ggfortify)
library(ggplot2)

# Run PCA on the full dataset
pca_full <- prcomp(scenario_data[strat_vars], center = TRUE, scale. = TRUE)

# Project sampled points into the same PCA space
sampled_proj <- predict(pca_full, newdata = sampled_params[strat_vars]) %>% as.data.frame()
sampled_proj$group <- "Sampled"

# Get PCA scores for full data (for plotting background)
full_proj <- pca_full$x %>% as.data.frame()
full_proj$group <- "Full"

# Combine
pca_df <- rbind(full_proj, sampled_proj)

png("joint.png", width=5, height= 5, res=500, units = "in")
# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Full" = "grey", "Sampled" = "red")) +
  theme_minimal() +
  labs(title = "PCA of Stratification Space: Sample Coverage")
dev.off()
