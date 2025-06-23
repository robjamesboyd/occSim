#' Simulate Stratified Random Sampling with Proportional Allocation
#'
#' This function implements stratified random sampling with proportional allocation.
#' Sites are selected in year 1 based on stratified sampling, and the same panel
#' is used across all years. All selected sites are visited on every occasion.
#'
#' @param Y Matrix [nSites x nYears] of true occupancy states.
#' @param nVisits Integer. Number of visit occasions per year.
#' @param strata Numeric vector of stratum IDs for each site (length = nSites).
#' @param prop_sample Numeric. Proportion of total sites to sample (between 0 and 1).
#' @param plot Logical. If TRUE, plot summary of sampling design and effects.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{S}: Sampling array [nSites x nYears x nVisits]
#'   \item \code{sampled_once}: Matrix [nSites x nYears], 1 if sampled at least once
#'   \item \code{summary_df}: Data frame summarising year-wise sampling and bias
#' }
#'
#' @export
simulateStratifiedSample <- function(Y, nVisits, strata, prop_sample = 0.1, plot = TRUE) {
  if (!is.matrix(Y)) stop("Y must be a matrix [nSites x nYears]")
  if (!is.numeric(strata) || length(strata) != nrow(Y)) stop("strata must be a numeric vector of length nSites")
  
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  
  S <- array(0, dim = c(nSites, nYears, nVisits))  # sampling array
  
  # ----- STRATIFIED RANDOM SAMPLING -----
  unique_strata <- unique(strata)
  sampled_sites <- c()
  
  for (s in unique_strata) {
    stratum_sites <- which(strata == s)
    n_stratum <- length(stratum_sites)
    
    # Proportional allocation: sample a fraction from each stratum
    n_to_sample <- round(prop_sample * n_stratum)
    sampled_in_stratum <- sample(stratum_sites, n_to_sample)
    
    sampled_sites <- c(sampled_sites, sampled_in_stratum)
  }
  
  # ----- FIXED PANEL: SAME SITES SAMPLED IN ALL YEARS -----
  for (t in 1:nYears) {
    for (v in 1:nVisits) {
      S[sampled_sites, t, v] <- 1  # full visit coverage for selected sites
    }
  }
  
  # ----- SUMMARISE -----
  sampled_once <- apply(S, c(1, 2), max)  # at least one visit per site-year
  
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
  
  # ----- OPTIONAL PLOTS -----
  if (plot) {
    library(ggplot2)
    library(patchwork)
    
    df_occ <- data.frame(
      Year = rep(1:nYears, 2),
      Occupancy = c(summary_df$mean_occ_all, summary_df$mean_occ_sampled),
      Group = rep(c("All Sites", "Sampled Sites"), each = nYears)
    )
    
    p_occ <- ggplot(df_occ, aes(x = Year, y = Occupancy, color = Group)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("black", "red")) +
      labs(title = "Occupancy Trends", y = "Mean Occupancy", x = "Year") +
      theme_linedraw()
    
    p_frac <- ggplot(summary_df, aes(x = Year, y = sampling_fraction)) +
      geom_line(size = 1.2, color = "blue") +
      labs(title = "Sampling Fraction", y = "Proportion of Sites Sampled", x = "Year") +
      theme_linedraw()
    
    p_corr <- ggplot(summary_df, aes(x = Year, y = data_defect)) +
      geom_line(size = 1.2, color = "darkgreen") +
      labs(title = "Data Defect Correlation", y = "r", x = "Year") +
      theme_linedraw()
    
    print(p_occ / p_frac / p_corr)
  }
  
  return(list(
    S = S,
    sampled_once = sampled_once,
    summary_df = summary_df
  ))
}
