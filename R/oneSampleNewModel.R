oneSample <- function(sample,
                      nSites,
                      nYears,
                      nStrata,
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
                      rho,
                      mean_change_z,
                      sd_change_z) {
  
  cat("sample ", sample, "\n")
  
  # Initialise x and z matrices
  x_mat <- matrix(NA, nrow = nSites, ncol = nYears)
  z_mat <- matrix(NA, nrow = nSites, ncol = nYears)
  
  # Year 1: z is given, and x is generated with correlation rho to z
  z_mat[, 1] <- z
  x_mat[, 1] <- rho * z + sqrt(1 - rho^2) * rnorm(nSites)
  
  # Generate stratum IDs based on year 1 x values (used for fixed stratified design)
  strat_IDs <- as.integer(cut(x_mat[, 1],
                              breaks = quantile(x_mat[, 1], probs = seq(0, 1, length.out = nStrata + 1)),
                              include.lowest = TRUE, labels = FALSE))
  
  # Years 2+: evolve z and generate new x with same rho correlation to updated z
  for (t in 2:nYears) {
    z_mat[, t] <- z + rnorm(nSites, mean = mean_change_z, sd = sd_change_z)
    x_mat[, t] <- rho * z_mat[, t] + sqrt(1 - rho^2) * rnorm(nSites)
  }
  
  # Use simulateOccupancy() to get the true states
  occupancy_matrix <- simulateOccupancyStatic(nYears = nYears,
                                              nSites = nSites, 
                                              z_mat = z_mat,
                                              beta0 = alpha_occ ,
                                              beta_z = beta_z_occ,
                                              plot = TRUE)

  # 2. -------------------------------
  # SIMULATE SAMPLING VISITS
  # --------------------------------
  # simple stratified random sample with proportional allocation
  SRS <- simulateStratifiedSample(Y = occupancy_matrix$simulated_data,
                                  nVisits = nVisits_SRS,
                                  strata = strat_IDs,
                                  prop_sample = sample_inclusion_prob_SRS)

  
  # large biased sample
  S <- simulateSample(Y = occupancy_matrix$simulated_data,
                      nVisits = 1,
                      z_mat = z_mat,
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
  D$D
  
  cat("Data simulated", "\n")

  # Whether species was detected at least once in each year
  obs_biased <- apply(D$D, c(1,2), max)
  obs_SRS  <- apply(D_SRS$D,  c(1,2), max)

  # Year matrix: 0 for year 1, 1 for year 2
  year_mat <- matrix(rep(0:(nYears - 1), each = nSites), nrow = nSites)

  ## fit integrated model
  
  # first estimate sample inclusion probabilities for the nonprobability sample.
  # the inverses of these probabilities are used weights for the nonprobability sample
  # part of the likelihood
  
  # Reshape NS inclusion info to long format
  detected_once_long <- data.frame(
    R = as.vector(S$sampled_once),              # response: sampled at least once
    year = rep(0:(nYears - 1), each = nSites),  # 0-based coding
    X = as.vector(x_mat)
  )

  # Fit logistic regression for sample inclusion
  inclusion_model <- glm(R ~ year * X, data = detected_once_long, family = binomial)

  # Predict inclusion probabilities for each site-year (avoid extreme probabilities)
  p_hat <- predict(inclusion_model, type = "response",newdata = detected_once_long)
  p_hat <- pmin(pmax(p_hat, 0.01), 0.99)  # Truncate for numerical stability

  # Reshape to [nSites x nYears] matrix
  weights_mat <- matrix(1 / p_hat, nrow = nSites, ncol = nYears)
  
  # first specify likelihood
  
  loglik_fn <- function(par) {
    beta0_p   <- par[1]
    beta0_np  <- par[2]
    beta_year <- par[3]  # shared year effect
    
    ll <- 0  # total log-likelihood
    
    # Probability sample (STRS) — no weights
    for (i in 1:nSites) {
      for (t in 1:nYears) {
        eta_p <- beta0_p + beta_year * year_mat[i, t]
        psi_p <- plogis(eta_p)
        y <- obs_SRS[i, t]
        ll <- ll + y * log(psi_p) + (1 - y) * log(1 - psi_p)
      }
    }
    
    # Non-probability sample (NS) — apply estimated weights
    for (i in 1:nSites) {
      for (t in 1:nYears) {
        eta_np <- beta0_np + beta_year * year_mat[i, t]
        psi_np <- plogis(eta_np)
        y <- obs_biased[i, t]
        w <- weights_mat[i, t]
        ll <- ll + w * (y * log(psi_np) + (1 - y) * log(1 - psi_np))
      }
    }
    
    return(-ll)  # negative log-likelihood for minimization
  }
  
  
  ## now optimise
  
  start_vals <- rep(0, 3)
  
  fit <- optim(
    par = start_vals,
    fn = loglik_fn,
    method = "BFGS",
    hessian = TRUE
  )
  
  mle <- fit$par
  names(mle) <- c("beta0_p", "beta0_np", "beta_year")  
  se <- sqrt(diag(solve(fit$hessian)))
  mle
  
  ## now fit per-dataset models, which do not involve shared parameters and can
  ## therefore be fitted using GLM
  
  # Reshape to long format for P sample
  df_p <- data.frame(
    y = as.vector(obs_SRS),
    year = as.vector(year_mat),
    x = rep(x, times = nYears)
  )

  # Fit GLM to P sample
  fit_p <- glm(y ~ year, data = df_p, family = binomial)

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
  fit_np <- glm(y ~ year, data = df_np, family = binomial)

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
  
  ## calculate true parameter for comparison 
  occ <- colMeans(occupancy_matrix$simulated_data)
  print(paste("Real % occupancy change =", ((occ[2] - occ[1])/occ[1]) * 100))
        
  # Log-odds ratio of occupancy (year 1 vs year 0)
  logit <- function(p) log(p / (1 - p))
  
  true_beta_year <- logit(occ[nYears]) - logit(occ[1])

  # Does the confidence interval for the year effect exclude zero and agree in sign with the true value?
  trend_detected <- as.integer(
    (conf_int$CI_Lower[3] > 0 | conf_int$CI_Upper[3] < 0) &
      (sign(conf_int$Estimate[3]) == sign(true_beta_year))
  ) # integrated model
  
  trend_detected_SRS <- as.integer(
    (confint_SRS$CI_Lower > 0 | confint_SRS$CI_Upper < 0) &
      (sign(est_SRS) == sign(true_beta_year))
  ) # SRS
  
  trend_detected_biased <- as.integer(
    (confint_biased$CI_Lower > 0 | confint_biased$CI_Upper < 0) &
      (sign(est_biased) == sign(true_beta_year))
  ) # nonprobability sample

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


