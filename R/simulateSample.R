#' Simulate Sampling Visits Across Sites and Years
#'
#' This function simulates whether a site is sampled on each of several visits per year,
#' with sampling probability determined by a logistic function of a site-level explanatory variable.
#' It also returns a matrix indicating whether each site was sampled at least once per year,
#' and plots occupancy trends and sampling fraction over time.
#'
#' @param Y A matrix of true occupancy states for each site and year, dimensions \code{[nSites, nYears]}.
#' @param nVisits Integer. Number of visit occasions per year.
#' @param z Numeric vector. Site-level covariate affecting sampling probability.
#' @param alpha_psi Numeric. Intercept on logit scale for sampling probability.
#' @param beta_psi Numeric. Slope for covariate effect on sampling probability.
#' @param plot Logical. If TRUE, display two-panel ggplot: one comparing occupancy at sampled (at least once in a year) vs. non-sampled sites,
#' and one showing the sampling fraction (number of sampled sites / total number of sites) per year.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{S}: Sampling array [nSites x nYears x nVisits]
#'   \item \code{sampled_once}: Binary matrix [nSites x nYears] for any-year sampling
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_linedraw scale_color_manual
#' @importFrom reshape2 melt
#' @importFrom patchwork wrap_plots
#' @export
simulateSample <- function(Y, nVisits, z, alpha_psi, beta_psi, plot = TRUE) {
  if (!is.matrix(Y)) stop("Y must be a matrix of dimensions [nSites, nYears]")
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  
  if (length(z) != nSites) {
    stop(paste0("Length of z (", length(z), ") must match number of sites (", nSites, ")"))
  }
  
  S <- array(0, dim = c(nSites, nYears, nVisits))
  psi <- plogis(alpha_psi + beta_psi * z)
  
  for (t in 1:nYears) {
    for (v in 1:nVisits) {
      S[, t, v] <- rbinom(nSites, size = 1, prob = psi)
    }
  }
  
  sampled_once <- apply(S, c(1, 2), max)
  
  if (plot) {
    # Panel 1: Occupancy
    mean_occ_all <- colMeans(Y)
    mean_occ_sampled <- sapply(1:nYears, function(t) {
      idx <- which(sampled_once[, t] == 1)
      if (length(idx) > 0) mean(Y[idx, t]) else NA
    })
    
    df_occ <- data.frame(
      Year = rep(1:nYears, 2),
      Occupancy = c(mean_occ_all, mean_occ_sampled),
      Group = rep(c("All Sites", "Sampled Sites"), each = nYears)
    )
    
    p_occ <- ggplot2::ggplot(df_occ, ggplot2::aes(x = Year, y = Occupancy, color = Group)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::scale_color_manual(values = c("black", "red")) +
      ggplot2::labs(title = "Occupancy Trends",
                    y = "Mean Occupancy", x = "Year", color = "") +
      ggplot2::theme_linedraw()
    
    # Panel 2: Sampling Fraction
    sampling_fraction <- colMeans(sampled_once)
    
    df_frac <- data.frame(
      Year = 1:nYears,
      FractionSampled = sampling_fraction
    )
    
    p_frac <- ggplot2::ggplot(df_frac, ggplot2::aes(x = Year, y = FractionSampled)) +
      ggplot2::geom_line(size = 1.2, color = "blue") +
      ggplot2::labs(title = "Sampling Fraction",
                    y = "Proportion of Sites Sampled", x = "Year") +
      ggplot2::theme_linedraw()
    
    # Combine plots using patchwork
    print(patchwork::wrap_plots(p_occ, p_frac, ncol = 1))
  }
  
  return(list(S = S, sampled_once = sampled_once))
}

