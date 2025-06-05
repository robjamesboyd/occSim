nSites <- 10000
nYears <- 2

w <- rnorm(nSites)
z <- rnorm(nSites)

manySamples <- lapply(1:100, FUN=oneSample,
                      nSites = nSites,
                      nYears = nYears,
                      z = z,
                      w = w,
                      alpha_occ = -1.387053,
                      beta_z_occ = 2,
                      alpha_biased = -4,
                      beta_z_biased = 1,
                      beta_year_biased = 0.5,
                      beta_year_z_biased = 1.5,
                      beta_detection = 1, 
                      nVisits_biased = 3,
                      nVisits_SRS = 1,
                      sample_inclusion_prob_SRS = 0.01,
                      target_error_biased_sample = 0.1,
                      rho = 0.4)

## combine summary stats across samples
manySamples <- do.call("rbind", manySamples)
str(manySamples)
## calculate error metrics (turn into a data.frame)
manySamples
metrics <- data.frame(
  model = c("integrated", "PS", "NS"),
  
  mse = c(
    mean(manySamples$actual_error_int^2),
    mean(manySamples$actual_error_SRS^2),
    mean(manySamples$actual_error_biased^2)
  ),
  
  bias_squared = c(
    mean(manySamples$actual_error_int)^2,
    mean(manySamples$actual_error_SRS)^2,
    mean(manySamples$actual_error_biased)^2
  ),
  
  variance = c(
    var(manySamples$est_int),
    var(manySamples$est_SRS),
    var(manySamples$est_biased)
  ),
  
  coverage = c(
    mean(manySamples$covered_int),
    mean(manySamples$covered_SRS),
    mean(manySamples$covered_biased)
  ),
  
  power = c(
    mean(manySamples$trend_detected_int),
    mean(manySamples$trend_detected_SRS),
    mean(manySamples$trend_detected_biased)
  )
)

print(metrics)
?occAssess::assessRarityBias

library(ggplot2)
library(dplyr)
library(tidyr)

# Start from metrics table
# First reshape MSE into components
mse_parts <- metrics %>%
  select(model, bias_squared, variance) %>%
  pivot_longer(cols = c(bias_squared, variance),
               names_to = "component",
               values_to = "value") %>%
  mutate(metric = "MSE")

# Reshape power and coverage into long format
other_metrics <- metrics %>%
  select(model, power, coverage) %>%
  pivot_longer(cols = c(power, coverage),
               names_to = "metric",
               values_to = "value") %>%
  mutate(component = "total")

# Combine all together
plot_data <- bind_rows(mse_parts, other_metrics)

# Set factor levels for better ordering
plot_data$metric <- factor(plot_data$metric, levels = c("MSE", "power", "coverage"))
plot_data$component <- factor(plot_data$component, levels = c("bias_squared", "variance", "total"))
plot_data$model <- factor(plot_data$model, levels = c("integrated", "PS", "NS"))

png("scenario.png", width=10, height = 4, units = "in", res = 500)
# Plot
ggplot(plot_data, aes(x = model, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("bias_squared" = "tomato", "variance" = "skyblue")) +
  labs(
    x = "Model",
    y = "",
    fill = "MSE
component",
    title = ""
  ) +
  theme_linedraw(base_size = 14)
dev.off()
?occAssess::assessEnvBias
