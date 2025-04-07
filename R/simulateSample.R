#' Simulate Sampling Visits Across Sites and Years
#'
#' This function simulates whether a site is sampled on each of several visits per year,
#' with sampling probability determined by a logistic function of a site-level explanatory variable.
#' The number of sites and years is inferred from the dimensions of a provided occupancy matrix.
#'
#' @param Y A matrix of true occupancy states for each site and year, with dimensions \code{[nSites, nYears]}.
#' This is used to determine the number of sites and years to simulate.
#' @param nVisits Integer. The number of sampling visits per year (assumed to be the same for all site-year combinations).
#' @param z Numeric vector. A site-level covariate influencing the probability of being sampled on each visit.
#' Must be of length equal to the number of rows (sites) in \code{Y}.
#' @param alpha_psi Numeric. The intercept on the logit scale for the sampling probability.
#' @param beta_psi Numeric. The slope on the logit scale for the effect of the covariate \code{z} on the sampling probability.
#'
#' @return A 3-dimensional array of sampling states with dimensions \code{[nSites, nYears, nVisits]}.
#'
#' @export
simulateSample <- function(Y, nVisits, z, alpha_psi, beta_psi) {
  # --- Check dimensions ---
  if (!is.matrix(Y)) stop("Y must be a matrix of dimensions [nSites, nYears]")
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  
  if (length(z) != nSites) {
    stop(paste0("Length of z (", length(z), ") must match number of sites (", nSites, ")"))
  }
  
  # --- Initialise output array ---
  S <- array(0, dim = c(nSites, nYears, nVisits))
  
  # --- Compute site-level sampling probabilities ---
  psi <- plogis(alpha_psi + beta_psi * z)
  
  # --- Loop to simulate visit occurrence ---
  for (t in 1:nYears) {
    for (v in 1:nVisits) {
      S[, t, v] <- rbinom(nSites, size = 1, prob = psi)
    }
  }
  
  return(S)
}

