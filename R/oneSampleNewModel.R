oneSample <- function(sample,
                      nSites,
                      nYears,
                      z,
                      w,
                      alpha_initial,
                      beta_initial,
                      alpha_colonisation,
                      beta_colonisation,
                      alpha_extinction,
                      beta_extinction,
                      target_f_biased_sample,
                      beta_sample_inclusion,
                      beta_detection, 
                      nVisits_biased,
                      nVisits_SRS,
                      sample_inclusion_prob_SRS,
                      target_error_biased_sample,
                      rho) {
  
  
  # Use simulateOccupancy() to get the true states
  occupancy_matrix <- simulateOccupancy(
    nYears = nYears,
    nSites = nSites,
    explanatoryVar = z,
    discreteExplanatoryVar = FALSE,
    beta0_extinction = alpha_extinction,
    beta0_colonization = alpha_colonisation,
    beta0_initial = alpha_initial,
    beta_extinction = beta_extinction,
    beta_colonization = beta_colonisation,
    beta_initial = beta_initial,
    randomInitial = FALSE,
    randomInitialPsi = 0.5,
    nBins = 5,
    plot = FALSE
  )$simulated_data
  
  
  # 2. -------------------------------
  # SIMULATE SAMPLING VISITS
  # --------------------------------
  
  # small simple random sample
  SRS <- simulateSample(
    Y = occupancy_matrix,
    nVisits = nVisits_SRS,
    z = z,
    alpha_psi = log(sample_inclusion_prob_SRS / (1 - sample_inclusion_prob_SRS)),
    beta_psi = 0,
    fixed_panel = TRUE,
    plot = FALSE
  )
  
  # large biased sample
  S <- simulateSample(
    Y = occupancy_matrix,
    nVisits = nVisits_biased,
    z = z,
    alpha_psi = find_sampling_intercept(target_f = target_f_biased_sample, Z=z, V=nVisits_biased,beta=beta_sample_inclusion),
    beta_psi = beta_sample_inclusion,
    plot = FALSE
  )
  
  # 3. -------------------------------
  # SIMULATE DETECTIONS CONDITIONAL ON OCCUPANCY + SAMPLING
  # --------------------------------

  # detection for srs
  D_SRS <- simulateDetection(
    Y = occupancy_matrix,
    S = SRS$S,
    w = w,
    alpha_p = log(0.99/(1-0.99)), # essentially perfect detection
    beta_p = 0,
    plot = FALSE
  )
  D_SRS
  
  # detection for biased sample
  D <- simulateDetection(
    Y = occupancy_matrix,
    S = S$S,
    w = w,
    alpha_p = find_detection_intercept(target_error=target_error_biased_sample, w=w, V=nVisits_biased,beta=beta_detection),
    beta_p = beta_detection,
    plot = FALSE
  )
  
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
  trend_detected <- as.integer(conf_int$CI_Lower[3] > 0 | conf_int$CI_Upper[3] < 0)

  ## calculate true parameter for comparison 
  occ <- colMeans(occupancy_matrix)
  
  # Log-odds ratio of occupancy (year 1 vs year 0)
  logit <- function(p) log(p / (1 - p))

  true_beta_year <- logit(occ[nYears]) - logit(occ[1])
  
  cat("True beta_year (from simulated occupancy change):", round(true_beta_year, 3), "\n")

  # Does the confidence interval for beta_year cover the true value? 
  covered <- as.integer(conf_int$CI_Lower[3] <= true_beta_year & conf_int$CI_Upper[3] >= true_beta_year)

  out <- data.frame(trend_detected = trend_detected,
                    covered = covered,
                    est = conf_int$Estimate[3],
                    actual_error = conf_int$Estimate[3] - true_beta_year)
  
  return(out)
  
}
              
