#' Simulate Sampling Visits Across Sites and Years with Time-Varying Effects
#'
#' This function simulates whether a site is sampled on each of several visits per year,
#' where the probability of sampling is determined by a logistic function of a site-level
#' covariate that may vary across time. The effect of the covariate may also vary with time
#' through an interaction with year. Optionally, a fixed panel design can be used, where the
#' same sites are sampled in all years.
#'
#' @param Y Matrix [nSites x nYears] of true occupancy states.
#' @param nVisits Integer. Number of visit occasions per year.
#' @param z_mat Matrix [nSites x nYears] of site-level covariates affecting sampling probability.
#' @param alpha_psi Numeric. Intercept for the sampling probability model (logit scale).
#' @param beta_psi Numeric. Coefficient for the effect of z on sampling probability.
#' @param beta_year Numeric. Coefficient for the effect of year on sampling probability.
#' @param beta_year_z Numeric. Coefficient for the interaction between z and year.
#' @param fixed_panel Logical. If TRUE, the same sites sampled in year 1 are sampled in all subsequent years.
#' @param plot Logical. If TRUE, display summary plots of occupancy, sampling, and data defect.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{S}: Sampling array [nSites x nYears x nVisits]
#'   \item \code{sampled_once}: Matrix [nSites x nYears], 1 if sampled at least once
#'   \item \code{summary_df}: Data frame with one row per year:
#'     \code{Year, sampling_fraction, mean_occ_all, mean_occ_sampled, data_defect}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_linedraw scale_color_manual
#' @importFrom patchwork wrap_plots
#' @export
simulateSample <- function(Y, nVisits, z_mat, alpha_psi, beta_psi, beta_year, beta_year_z, fixed_panel = FALSE, plot = TRUE) {
  if (!is.matrix(Y)) stop("Y must be a matrix of dimensions [nSites, nYears]")
  
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  
  if (!is.matrix(z_mat) || !all(dim(z_mat) == c(nSites, nYears))) {
    stop("z_mat must be a matrix with dimensions [nSites x nYears]")
  }
  
  S <- array(0, dim = c(nSites, nYears, nVisits))  # sampling visits
  psi_mat <- matrix(NA, nSites, nYears)            # sampling probabilities
  
  # Compute sampling probabilities (logistic model) for each site-year
  for (t in 1:nYears) {
    year_val <- t - 1  # use 0-based coding for year
    z_t <- z_mat[, t]  # covariate values for year t
    linpred <- alpha_psi + beta_psi * z_t + beta_year * year_val + beta_year_z * z_t * year_val
    psi_mat[, t] <- plogis(linpred)  # inverse logit
  }
  
  if (fixed_panel) {
    # Sample site-visit combinations in year 1 and repeat across all years
    sampled_first_year <- matrix(0, nrow = nSites, ncol = nVisits)
    for (v in 1:nVisits) {
      sampled_first_year[, v] <- rbinom(nSites, 1, psi_mat[, 1])
    }
    for (t in 1:nYears) {
      S[, t, ] <- sampled_first_year
    }
  } else {
    # Sample independently each year using year-specific probabilities
    for (t in 1:nYears) {
      for (v in 1:nVisits) {
        S[, t, v] <- rbinom(nSites, 1, psi_mat[, t])
      }
    }
  }
  
  # Determine whether each site was sampled at least once in each year
  sampled_once <- apply(S, c(1, 2), max)
  
  # --- Summary statistics per year ---
  summary_df <- data.frame(
    Year = 1:nYears,
    sampling_fraction = NA_real_,
    mean_occ_all = NA_real_,
    mean_occ_sampled = NA_real_,
    data_defect = NA_real_
  )
  
  for (t in 1:nYears) {
    sampled <- sampled_once[, t]
    occupancy <- Y[, t]
    summary_df$sampling_fraction[t] <- mean(sampled)
    summary_df$mean_occ_all[t] <- mean(occupancy)
    if (any(sampled == 1)) {
      summary_df$mean_occ_sampled[t] <- mean(occupancy[sampled == 1])
      summary_df$data_defect[t] <- cor(occupancy, sampled)
    }
  }
  
  # --- Optional diagnostic plots ---
  if (plot) {
    df_occ <- data.frame(
      Year = rep(1:nYears, 2),
      Occupancy = c(summary_df$mean_occ_all, summary_df$mean_occ_sampled),
      Group = rep(c("All Sites", "Sampled Sites"), each = nYears)
    )
    
    p_occ <- ggplot2::ggplot(df_occ, ggplot2::aes(x = Year, y = Occupancy, color = Group)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::scale_color_manual(values = c("black", "red")) +
      ggplot2::labs(title = "Occupancy Trends",
                    y = "Mean Occupancy", x = "Year", color = "") +
      ggplot2::theme_linedraw()
    
    p_frac <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Year, y = sampling_fraction)) +
      ggplot2::geom_line(size = 1.2, color = "blue") +
      ggplot2::labs(title = "Sampling Fraction",
                    y = "Proportion of Sites Sampled", x = "Year") +
      ggplot2::theme_linedraw()
    
    p_corr <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Year, y = data_defect)) +
      ggplot2::geom_line(size = 1.2, color = "darkgreen") +
      ggplot2::labs(title = "Data Defect Correlation",
                    y = "r", x = "Year") +
      ggplot2::theme_linedraw()
    
    print(patchwork::wrap_plots(p_occ, p_frac, p_corr, ncol = 1))
  }
  
  return(list(
    S = S,
    sampled_once = sampled_once,
    summary_df = summary_df
  ))
}

