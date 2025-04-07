#' Simulate Detection Histories Conditional on Occupancy and Sampling
#'
#' This function simulates detection/non-detection data over multiple visits per year,
#' conditional on both true occupancy and visit occurrence (i.e., whether a site was sampled).
#' Detection probability is defined by a logistic function of a site-level covariate.
#'
#' @param Y A matrix of true occupancy states for each site and year, with dimensions \code{[nSites, nYears]}.
#' @param S A 3D array of sampling events with dimensions \code{[nSites, nYears, nVisits]}.
#' Each element \code{S[i, t, v]} should be 1 if a visit occurred and 0 otherwise.
#' @param w Numeric vector. A site-level covariate influencing detection probability. Must be of length equal to the number of sites.
#' @param alpha_p Numeric. Intercept on the logit scale for detection probability.
#' @param beta_p Numeric. Slope on the logit scale for the effect of \code{w} on detection probability.
#'
#' @return A 3D array of detection outcomes \code{D[i, t, v]}, where each element is:
#' \itemize{
#'   \item 1 if a species was detected (i.e., site was sampled, species was present, and detected)
#'   \item 0 if site was sampled but species was either absent or not detected
#'   \item 0 if no visit occurred (i.e., \code{S[i, t, v] = 0})
#' }
#'
#' @examples
#' # Create dummy occupancy and sampling data
#' set.seed(123)
#' Y <- matrix(rbinom(100 * 5, 1, 0.4), nrow = 100, ncol = 5)
#' S <- simulateSample(Y, nVisits = 3, z = rnorm(100), alpha_psi = -0.5, beta_psi = 1)
#' w <- rnorm(100)  # site-level detection covariate
#'
#' # Simulate detection
#' D <- simulateDetection(Y, S, w, alpha_p = -1, beta_p = 2)
#'
#' # Look at a site's detection history
#' D[1, , ]
#'
#' @export
simulateDetection <- function(Y, S, w, alpha_p, beta_p) {
  # --- Dimensions ---
  nSites <- dim(S)[1]
  nYears <- dim(S)[2]
  nVisits <- dim(S)[3]
  
  # --- Initialise detection array ---
  D <- array(0, dim = c(nSites, nYears, nVisits))
  
  # --- Compute site-level detection probabilities ---
  p <- plogis(alpha_p + beta_p * w)
  
  # --- Loop over sites, years, visits ---
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      for (v in 1:nVisits) {
        if (S[i, t, v] == 1 && Y[i, t] == 1) {
          # Species present and visit occurred → simulate detection
          D[i, t, v] <- rbinom(1, size = 1, prob = p[i])
        }
        # Else remains 0: no visit or species absent
      }
    }
  }
  
  # Return detection array
  return(D)
}
