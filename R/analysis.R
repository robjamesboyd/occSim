files <- list.files("C:/Users/robboy/OneDrive - UKCEH/Documents/occSim", 
                    pattern = "out.csv")

dat <- lapply(files, read.csv)

dat <- do.call("rbind", dat)
str(dat)
dat$bias_change <- abs(dat$bias_change)
dat$bias_year1 <- abs(dat$bias_year1)

library(dplyr)
library(tidyr)
library(ggplot2)

scenario_summary <- dat %>%
  group_by(scenario) %>%
  summarise(
    # Performance metrics
    mse_I = mean(actual_error_int^2),
    mse_STRS = mean(actual_error_SRS^2),
    mse_NS = mean(actual_error_biased^2),
    
    bias_I = mean(actual_error_int),
    bias_STRS = mean(actual_error_SRS),
    bias_NS = mean(actual_error_biased),
    
    var_I = var(est_int),
    var_STRS = var(est_SRS),
    var_NS = var(est_biased),
    
    coverage_I = mean(covered_int),
    coverage_STRS = mean(covered_SRS),
    coverage_NS = mean(covered_biased),
    
    power_I = mean(trend_detected_int),
    power_STRS = mean(trend_detected_SRS),
    power_NS = mean(trend_detected_biased),
    
    # Scenario conditions (mean summary stats)
    bias_year1 = mean(bias_year1),
    bias_change = mean(bias_change),
    sample_size1 = mean(sample_size1),
    sample_size_change = mean(sample_size_change),
    occ_year1 = mean(occ_year1),
    occ_diff = mean(occ_diff),
    occ_logit_diff = mean(occ_logit_diff),
    
    .groups = "drop"
  )

str(scenario_summary)

scenario_summary <- scenario_summary %>%
  mutate(
    best_estimator_MSE = case_when(
      mse_I < mse_STRS & mse_I < mse_NS ~ "I",
      mse_STRS < mse_I & mse_STRS < mse_NS ~ "STRS",
      TRUE ~ "NS"
    ),
    best_estimator_bias = case_when(
      abs(bias_I) < abs(bias_STRS) & abs(bias_I) < abs(bias_NS) ~ "I",
      abs(bias_STRS) < abs(bias_I) & abs(bias_STRS) < abs(bias_NS) ~ "STRS",
      TRUE ~ "NS"
    ),
    best_estimator_variance = case_when(
      var_I < var_STRS & var_I < var_NS ~ "I",
      var_STRS < var_I & var_STRS < var_NS ~ "STRS",
      TRUE ~ "NS"
    ),
    best_estimator_coverage = case_when(
      coverage_I > coverage_STRS & coverage_I > coverage_NS ~ "I",
      coverage_STRS > coverage_I & coverage_STRS > coverage_NS ~ "STRS",
      TRUE ~ "NS"
    )
  )

# Pivot longer for plotting
scenario_long <- scenario_summary %>%
  pivot_longer(
    cols = starts_with("best_estimator"),
    names_to = "metric",
    values_to = "best_estimator"
  ) %>%
  mutate(
    metric = gsub("best_estimator_", "", metric),
    best_estimator = factor(best_estimator, levels = c("I", "STRS", "NS"))
  )

png("best_by_metric.png", width = 8, height = 4, units = "in", res = 500)
# Plot
ggplot(scenario_long, aes(x = best_estimator, fill = best_estimator)) +
  geom_bar(width = 0.7) +
  facet_wrap(~ metric, nrow = 1) +
  labs(
    title = "Best Estimator by Performance Metric",
    x = "Estimator",
    y = "Number of Scenarios"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

dev.off()
## model to work out what determines which estimator is best 

library(ggplot2)
library(tidyr)
library(dplyr)

# First, pivot the relevant columns to long format
library(dplyr)
library(tidyr)
library(ggplot2)

# Create long-format data with bias, variance, and MSE
scenario_long_metrics <- scenario_summary %>%
  dplyr::select(scenario, 
         bias_I = bias_I, bias_STRS = bias_STRS, bias_NS = bias_NS,
         var_I = var_I, var_STRS = var_STRS, var_NS = var_NS,
         mse_I = mse_I, mse_STRS = mse_STRS, mse_NS = mse_NS) %>%
  pivot_longer(
    cols = -scenario,
    names_to = c("metric", "estimator"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    estimator = factor(estimator, levels = c("I", "STRS", "NS")),
    metric = factor(metric, levels = c("bias", "var", "mse"))
  ) %>%
  filter(!is.na(estimator))  # optional: drop any NA rows just in case

# Plot density curves
ggplot(scenario_long_metrics, aes(x = value, fill = estimator, color = estimator)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ metric, scales = "free", nrow = 1) +
  labs(
    title = "Distribution of Bias, Variance, and MSE Across Estimators",
    x = "Value", y = "Density"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 14)

## decision tree
library(rpart)
library(rpart.plot)
library(dplyr)

# Convert best_estimator_MSE to a factor
scenario_summary <- scenario_summary %>%
  mutate(best_estimator_MSE = factor(best_estimator_MSE, 
                                     levels = c("I", "STRS", "NS")))

# Fit decision tree
tree_mse <- rpart(
  best_estimator_MSE ~ bias_year1 + bias_change + sample_size1 + sample_size_change +
    occ_year1 + occ_diff,
  data = scenario_summary,
  method = "class",      # because it's classification
  control = rpart.control(cp = 0.0001)  # tweak cp to control tree complexity
)
summary(tree_mse)

png("tree.png", width = 8, height = 7, units = "in", res =500)
rpart.plot(tree_mse, type = 2, extra = 104, fallen.leaves = TRUE)
dev.off()
rpart.rules(tree_mse)
printcp(tree_mse)
plotcp(tree_mse)
