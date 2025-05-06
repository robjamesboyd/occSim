#' Simulate Sampling Visits Across Sites and Years
#'
#' This function simulates whether a site is sampled on each of several visits per year,
#' with sampling probability determined by a logistic function of a site-level explanatory variable.
#' Optionally, the same sites may be sampled in all years (a fixed panel design).
#'
#' @param Y Matrix [nSites x nYears] of true occupancy states.
#' @param nVisits Integer. Number of visit occasions per year.
#' @param z Numeric vector. Site-level covariate for sampling probability.
#' @param alpha_psi, beta_psi Numeric. Intercept and slope for sampling probability (on logit scale).
#' @param fixed_panel Logical. If TRUE, the same sites sampled in year 1 are sampled in all subsequent years.
#' @param plot Logical. If TRUE, display 3 time-series: (1) occupancy at sampled vs. all sites,
#' (2) the sampling fraction, and (3) correlation between sample inclusion and occupancy.
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
simulateSample <- function(Y, nVisits, z, alpha_psi, beta_psi, fixed_panel = FALSE, plot = TRUE) {
  if (!is.matrix(Y)) stop("Y must be a matrix of dimensions [nSites, nYears]")
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  
  if (length(z) != nSites) {
    stop(paste0("Length of z (", length(z), ") must match number of sites (", nSites, ")"))
  }
  
  S <- array(0, dim = c(nSites, nYears, nVisits))
  psi <- plogis(alpha_psi + beta_psi * z)
  
  if (fixed_panel) {
    # Sample only in year 1
    sampled_first_year <- matrix(0, nrow = nSites, ncol = nVisits)
    for (v in 1:nVisits) {
      sampled_first_year[, v] <- rbinom(nSites, 1, psi)
    }
    for (t in 1:nYears) {
      S[, t, ] <- sampled_first_year
    }
  } else {
    for (t in 1:nYears) {
      for (v in 1:nVisits) {
        S[, t, v] <- rbinom(nSites, 1, psi)
      }
    }
  }
  
  sampled_once <- apply(S, c(1, 2), max)
  
  # --- Summary stats ---
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
  
  # --- Optional plots ---
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
