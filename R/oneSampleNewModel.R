oneSample <- function(sample,
                      nSites,
                      nYears,
                      z,
                      w,
                      alpha_occ,
                      beta_z_occ,
                      alpha_biased,
                      beta_z_biased,
                      beta_year_biased,
                      beta_year_z_biased,
                      beta_detection, 
                      nVisits_biased,
                      nVisits_SRS,
                      sample_inclusion_prob_SRS,
                      target_error_biased_sample,
                      rho) {
  
  cat("sample ", sample, "\n")
  
  # Use simulateOccupancy() to get the true states
  occupancy_matrix <- simulateOccupancyStatic(nYears = nyears,
                                              nSites = nSites, 
                                              z_mat = z,
                                              beta0 = alpha_occ ,
                                              beta_z = beta_z_occ,
                                              plot = TRUE)

  # 2. -------------------------------
  # SIMULATE SAMPLING VISITS
  # --------------------------------
  # simple random sample 
  SRS <- simulateSample(Y = occupancy_matrix$simulated_data,
                        nVisits = nVisits_SRS,
                        z_mat = z,
                        alpha_psi = log(sample_inclusion_prob_SRS / (1 - sample_inclusion_prob_SRS)),
                        beta_psi = 0,
                        beta_year = 0,
                        beta_year_z = 0,
                        fixed_panel = TRUE)
  
  # large biased sample
  S <- simulateSample(Y = occ$simulated_data,
                      nVisits = nVisits_biased,
                      z_mat = z,
                      alpha_psi = alpha_biased,
                      beta_psi = beta_z_biased,
                      beta_year = beta_year_biased,
                      beta_year_z = beta_year_z_biased)

  # 3. -------------------------------
  # SIMULATE DETECTIONS CONDITIONAL ON OCCUPANCY + SAMPLING
  # --------------------------------

  # detection for srs
  D_SRS <- simulateDetection(
    Y = occupancy_matrix$simulated_data,
    S = SRS$S,
    w = w,
    alpha_p = log(0.99/(1-0.99)), # essentially perfect detection
    beta_p = 0,
    plot = FALSE
  )
  D_SRS
  
  # detection for biased sample
  D <- simulateDetection(
    Y = occupancy_matrix$simulated_data,
    S = S$S,
    w = w,
    alpha_p = find_detection_intercept(target_error=target_error_biased_sample, w=w, V=nVisits_biased,beta=beta_detection),
    beta_p = beta_detection,
    plot = FALSE
  )
  
  cat("Data simulated", "\n")
  
  ## create covariate and format data to fit model
  x <- rho * z + sqrt(1 - rho^2) * rnorm(nSites)  # x ~ correlated with z
  x_mat <- matrix(rep(x, nYears), nrow = nSites)

  # Whether species was detected at least once in each year
  obs_biased <- apply(D$D, c(1,2), max)
  obs_SRS  <- apply(D_SRS$D,  c(1,2), max)

  # Year matrix: 0 for year 1, 1 for year 2
  year_mat <- matrix(rep(0:(nYears - 1), each = nSites), nrow = nSites)

  ## fit model
  
  # first specify likelihood
  
  # ===== Likelihood Function =====
  loglik_fn <- function(par) {
    beta0_p   <- par[1]
    beta0_np  <- par[2]
    beta_year <- par[3]  # shared
    beta_x_p  <- par[4]
    beta_x_np <- par[5]
    
    ll <- 0
    
    # Probability sample
    for (i in 1:nSites) {
      for (t in 1:nYears) {
        eta_p <- beta0_p + beta_year * year_mat[i,t] + beta_x_p * x[i]
        psi_p <- plogis(eta_p)
        y <- obs_SRS[i,t]
        ll <- ll + y * log(psi_p) + (1 - y) * log(1 - psi_p)
      }
    }

    # Non-probability sample
    for (i in 1:nSites) {
      for (t in 1:nYears) {
        eta_np <- beta0_np + beta_year * year_mat[i,t] + beta_x_np * x[i]
        psi_np <- plogis(eta_np)
        y <- obs_biased[i,t]
        ll <- ll + y * log(psi_np) + (1 - y) * log(1 - psi_np)
      }
    }
    
    return(-ll)  # negative log-likelihood
  }
  
  ## now optimise
  
  # ===== Optimization =====
  start_vals <- rep(0, 5)
  fit <- optim(start_vals, loglik_fn, method = "BFGS", hessian = TRUE)
  
  # Output
  mle <- fit$par
  names(mle) <- c("beta0_p", "beta0_np", "beta_year", "beta_x_p", "beta_x_np")
  se <- sqrt(diag(solve(fit$hessian)))
  
  ## now fit per-dataset models, which do not involve shared parameters and can
  ## therefore be fitted using GLM
  
  # Reshape to long format for P sample
  df_p <- data.frame(
    y = as.vector(obs_SRS),
    year = as.vector(year_mat),
    x = rep(x, times = nYears)
  )

  # Fit GLM to P sample
  fit_p <- glm(y ~ year + x, data = df_p, family = binomial)

  est_SRS <- coef(fit_p)[2]
  
  Confint <- confint(fit_p)
  
  confint_SRS <- data.frame(CI_Lower = Confint[2,1],
                            CI_Upper = Confint[2,2])

  # Reshape to long format for NP sample
  df_np <- data.frame(
    y = as.vector(obs_biased),
    year = as.vector(year_mat),
    x = rep(x, times = nYears)
  )
  
  # Fit GLM to NP sample
  fit_np <- glm(y ~ year + x, data = df_np, family = binomial)

  est_biased <- coef(fit_np)[2]

  Confint <- confint(fit_np)
  
  confint_biased <- data.frame(CI_Lower = Confint[2,1],
                            CI_Upper = Confint[2,2])
  
  cat("Models fitted", "\n")
  
  ## compute parameter uncertainty 
  
  # Confidence intervals: Wald-type (estimate ± 1.96 * SE)
  ci_lower <- mle - 1.96 * se
  ci_upper <- mle + 1.96 * se
  
  conf_int <- data.frame(
    Parameter = names(mle),
    Estimate = round(mle, 3),
    SE = round(se, 3),
    CI_Lower = round(ci_lower, 3),
    CI_Upper = round(ci_upper, 3)
  ); conf_int

  # Does the confidence interval for the year effect span zero? 1 if yes and 0 otherwise
  trend_detected <- as.integer(conf_int$CI_Lower[3] > 0 | conf_int$CI_Upper[3] < 0) # integrated model
  trend_detected_SRS <- as.integer(confint_SRS$CI_Lower > 0 | confint_SRS$CI_Upper < 0) # SRS
  trend_detected_biased <- as.integer(confint_biased$CI_Lower > 0 | confint_biased$CI_Upper < 0) # nonprob sample

  ## calculate true parameter for comparison 
  occ <- colMeans(occupancy_matrix$simulated_data)
  
  # Log-odds ratio of occupancy (year 1 vs year 0)
  logit <- function(p) log(p / (1 - p))

  true_beta_year <- logit(occ[nYears]) - logit(occ[1])

  # Does the confidence interval for beta_year cover the true value? 
  covered <- as.integer(conf_int$CI_Lower[3] <= true_beta_year & conf_int$CI_Upper[3] >= true_beta_year) # integrated model
  covered_SRS <- as.integer(confint_SRS$CI_Lower <= true_beta_year & confint_SRS$CI_Upper >= true_beta_year) # SRS
  covered_biased <- as.integer(confint_biased$CI_Lower <= true_beta_year & confint_biased$CI_Upper >= true_beta_year) # SRS

  out <- data.frame(true_trend = true_beta_year,
                    trend_detected_int = trend_detected,
                    trend_detected_SRS = trend_detected_SRS,
                    trend_detected_biased = trend_detected_biased,
                    covered_int = covered,
                    covered_SRS = covered_SRS,
                    covered_biased = covered_biased,
                    est_int = conf_int$Estimate[3],
                    est_SRS = est_SRS,
                    est_biased = est_biased,
                    actual_error_int = conf_int$Estimate[3] - true_beta_year,
                    actual_error_SRS = est_SRS - true_beta_year,
                    actual_error_biased = est_biased - true_beta_year)
  
  return(out)
  
}
              
