#' Simulate Detection Histories Conditional on Occupancy and Sampling
#'
#' @param Y A matrix of true occupancy states for each site and year, dimensions \code{[nSites, nYears]}.
#' @param S A 3D array of visit indicators, dimensions \code{[nSites, nYears, nVisits]}.
#' @param w Numeric vector of detection covariates per site, length equal to \code{nSites}.
#' @param alpha_p Numeric. Intercept on the logit scale for detection probability.
#' @param beta_p Numeric. Slope for the detection covariate.
#'
#' @return A 3D array of detections with dimensions \code{[nSites, nYears, nVisits]}.
#'
#' @export
simulateDetection <- function(Y, S, w, alpha_p, beta_p) {
  # --- Check dimensions ---
  if (!is.matrix(Y)) stop("Y must be a matrix")
  if (!is.array(S) || length(dim(S)) != 3) stop("S must be a 3D array")
  
  nSites <- nrow(Y)
  nYears <- ncol(Y)
  nVisits <- dim(S)[3]
  
  if (dim(S)[1] != nSites || dim(S)[2] != nYears) {
    stop("Dimensions of S must match dimensions of Y: [nSites, nYears, nVisits]")
  }
  
  if (length(w) != nSites) {
    stop(paste0("Length of w (", length(w), ") must match number of sites (", nSites, ")"))
  }
  
  # --- Initialise detection array ---
  D <- array(0, dim = c(nSites, nYears, nVisits))
  
  # --- Detection probability per site ---
  p <- plogis(alpha_p + beta_p * w)
  
  # --- Simulate detections conditional on occupancy and sampling ---
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      if (Y[i, t] == 1) {
        for (v in 1:nVisits) {
          if (S[i, t, v] == 1) {
            D[i, t, v] <- rbinom(1, 1, p[i])
          }
        }
      }
    }
  }
  
  return(D)
}

