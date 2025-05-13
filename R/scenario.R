nSites <- 10000

w <- rnorm(nSites)
z <- rnorm(nSites)

manySamples <- lapply(1:2, FUN=oneSample,
                      nSites = 10000,
                      nYears = 2,
                      z = z,
                      w = w,
                      alpha_initial = 0.03834128,
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
                      beta_detection = 2,
                      rho = 0.5)

## combine summary stats across samples
manySamples <- do.call("rbind", manySamples)

## calculate error metrics
mse <- mean(manySamples$actual_error^2)
bias_squared <- mean(manySamples$actual_error)^2
variance <- var(manySamples$est)
coverage <- nrow(manySamples[manySamples$covered==1,]) /
  nrow(manySamples)
power <- nrow(manySamples[manySamples$trend_detected==1,]) /
  nrow(manySamples)