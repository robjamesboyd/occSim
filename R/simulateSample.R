#' Simulate Sampling Visits Across Sites and Years
#'
#' This function simulates whether a site is sampled on each of several visits per year,
#' with sampling probability determined by a logistic function of a site-level explanatory variable.
#' The number of sites and years is inferred from the dimensions of a provided occupancy matrix.
#'
#' @param Y A matrix of true occupancy states for each site and year, with dimensions \code{nSites × nYears}.
#' This is used to determine the number of sites and years to simulate.
#' @param nVisits Integer. The number of sampling visits per year (assumed to be the same for all site-year combinations).
#' @param z Numeric vector. A site-level covariate influencing the probability of being sampled on each visit.
#' Must be of length equal to the number of rows (sites) in \code{Y}.
#' @param alpha_psi Numeric. The intercept on the logit scale for the sampling probability.
#' @param beta_psi Numeric. The slope on the logit scale for the effect of the covariate \code{z} on the sampling probability.
#'
#' @return A 3-dimensional array of sampling states with dimensions \code{[nSites, nYears, nVisits]}.
#' Each element \code{S[i, t, v]} is a binary indicator equal to 1 if site \code{i} was sampled on visit \code{v} in year \code{t}, and 0 otherwise.
#'
#' @examples
#' # Simulate an occupancy matrix first
#' set.seed(123)
#' nSites <- 100
#' nYears <- 5
#' Y <- matrix(rbinom(nSites * nYears, 1, 0.5), nrow = nSites, ncol = nYears)
#'
#' # Simulate site-level sampling covariate
#' z <- rnorm(nSites)
#'
#' # Simulate 3 visits per year using a moderate sampling bias
#' S <- simulateSample(Y, nVisits = 3, z = z, alpha_psi = -1, beta_psi = 2)
#'
#' # Check the dimensions and a sample of visit outcomes
#' dim(S)
#' S[1, , ]
#'
#' @export
simulateSample <- function(Y, nVisits, z, alpha_psi, beta_psi) {
  # --- Get dimensions from occupancy matrix ---
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  
  # --- Initialise output array ---
  S <- array(0, dim = c(nSites, nYears, nVisits))
  
  # --- Compute sampling probabilities per site ---
  # psi[i] is the probability that site i is sampled on any one visit
  psi <- plogis(alpha_psi + beta_psi * z)
  
  # --- Loop over years and visits ---
  for (t in 1:nYears) {
    for (v in 1:nVisits) {
      # For each visit in each year, sample visit occurrence
      S[, t, v] <- rbinom(nSites, size = 1, prob = psi)
    }
  }
  
  # Return sampling array
  return(S)
}
