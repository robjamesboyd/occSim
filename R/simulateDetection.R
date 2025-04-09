#' Simulate Detection Histories Conditional on Occupancy and Sampling
#'
#' @param Y A matrix of true occupancy states for each site and year, dimensions \code{[nSites, nYears]}.
#' @param S A 3D array of visit indicators, dimensions \code{[nSites, nYears, nVisits]}.
#' @param w Numeric vector of detection covariates per site, length equal to \code{nSites}.
#' @param alpha_p Numeric. Intercept on the logit scale for detection probability.
#' @param beta_p Numeric. Slope for the detection covariate.
#' @param plot Logical. If TRUE, plots mean sampled-but-missed error and correlation of error with sampling (at least one visit=1; no visit =0).
#'
#' @return A list:
#' \itemize{
#'   \item \code{D}: 3D array of detections [nSites x nYears x nVisits]
#'   \item \code{error_summary}: data frame with \code{Year}, \code{mean_error}, \code{mean_sampled_missed}, \code{cor_sampling_error}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_linedraw geom_point
#' @importFrom patchwork wrap_plots
#' @export
simulateDetection <- function(Y, S, w, alpha_p, beta_p, plot = TRUE) {
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
  
  D <- array(0, dim = c(nSites, nYears, nVisits))
  p <- plogis(alpha_p + beta_p * w)
  
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
  
  error_matrix <- matrix(0, nrow = nSites, ncol = nYears)
  sampled_missed_matrix <- matrix(NA, nrow = nSites, ncol = nYears)
  sampled_matrix <- apply(S, c(1, 2), max)  # binary: sampled at least once
  
  cor_sampling_error <- rep(NA, nYears)
  
  for (t in 1:nYears) {
    for (i in 1:nSites) {
      is_present <- Y[i, t] == 1
      was_sampled <- sampled_matrix[i, t] == 1
      was_detected <- sum(D[i, t, ]) > 0
      
      if (is_present && !was_detected) {
        error_matrix[i, t] <- -1
        sampled_missed_matrix[i, t] <- if (was_sampled) -1 else NA
      } else {
        error_matrix[i, t] <- 0
        sampled_missed_matrix[i, t] <- if (was_sampled) 0 else NA
      }
    }
    
    # Correlation between sampling (binary) and error (0/-1)
    cor_sampling_error[t] <- suppressWarnings(
      cor(sampled_matrix[, t], error_matrix[, t])
    )
  }
  
  mean_error <- colMeans(error_matrix)
  mean_sampled_missed <- colMeans(sampled_missed_matrix, na.rm = TRUE)
  
  error_summary <- data.frame(
    Year = 1:nYears,
    mean_error = mean_error,
    mean_sampled_missed = mean_sampled_missed,
    cor_sampling_error = cor_sampling_error
  )
  
  if (plot) {
    library(ggplot2)
    library(patchwork)
    
    p1 <- ggplot(error_summary, aes(x = Year, y = mean_sampled_missed)) +
      geom_line(color = "darkorange", size = 1.2) +
      geom_point(color = "darkorange") +
      labs(title = "Mean Sampled-But-Missed Error", y = "Mean Error", x = "Year") +
      theme_linedraw()
    
    p2 <- ggplot(error_summary, aes(x = Year, y = cor_sampling_error)) +
      geom_line(color = "steelblue", size = 1.2) +
      geom_point(color = "steelblue") +
      labs(title = "Data defect correlation", y = "r", x = "Year") +
      theme_linedraw()
    
    print(wrap_plots(p1, p2, ncol = 1))
  }
  
  return(list(D = D, error_summary = error_summary))
}




