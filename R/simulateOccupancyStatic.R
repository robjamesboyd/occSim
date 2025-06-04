#' Simulate Static Species Occupancy Over Time
#'
#' This function simulates species occupancy across multiple years and sites using a static occupancy model.
#' Occupancy is modelled as a Bernoulli random variable whose probability is a logit-linear function of
#' a year effect and a site-level covariate \code{z}. This setup allows investigation of occupancy dynamics
#' driven by changing environmental conditions or sampling covariates.
#'
#' @param nYears Integer. Number of years to simulate.
#' @param nSites Integer. Number of sites to simulate.
#' @param z_mat Numeric matrix. A matrix of site-level covariate values, dimension \code{nSites x nYears}.
#' @param beta0 Numeric. Intercept for the occupancy model.
#' @param beta_z Numeric. Coefficient for the covariate \code{z}.
#' @param plot Logical. Whether to plot the mean occupancy over time.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{occupancy_all_sites}: A data frame with the mean occupancy across sites for each year.
#'   \item \code{simulated_data}: A matrix of simulated occupancy states (1 = occupied, 0 = not occupied), dimensions \code{nSites x nYears}.
#' }
#'
#' @examples
#' nSites <- 500
#' nYears <- 10
#' z_mat <- matrix(rnorm(nSites * nYears), nrow = nSites)
#' sim <- simulateOccupancyStatic(nYears = nYears, nSites = nSites, z_mat = z_mat,
#'                                beta0 = -1, beta_year = 0.5, beta_z = 1)
#' head(sim$occupancy_all_sites)
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_linedraw
#' @export

simulateOccupancyStatic <- function(nYears,
                                    nSites,
                                    z_mat,
                                    beta0,
                                    beta_z,
                                    plot = TRUE) {
  
  # Checks
  if (!is.matrix(z_mat) || nrow(z_mat) != nSites || ncol(z_mat) != nYears) {
    stop("z_mat must be a matrix of dimensions nSites x nYears.")
  }
  
  # Construct year variable
  year_vec <- 0:(nYears - 1)
  year_mat <- matrix(rep(year_vec, each = nSites), nrow = nSites)
  
  # Linear predictor
  eta <- beta0 + beta_z * z_mat
  
  # Apply inverse logit to get probabilities
  psi <- 1 / (1 + exp(-eta))
  
  # Simulate occupancy
  occupancy <- matrix(rbinom(nSites * nYears, 1, as.vector(psi)), nrow = nSites)
  
  # Compute mean occupancy over time
  true_occupancy <- data.frame(
    year = 1:nYears,
    occupancy = colMeans(occupancy)
  )
  
  # Plot if requested
  if (plot) {
    print(
      ggplot2::ggplot(true_occupancy, ggplot2::aes(x = year, y = occupancy)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "", x = "Year", y = "Occupancy") +
        ggplot2::theme_linedraw()
    )
  }
  
  return(list(
    occupancy_all_sites = true_occupancy,
    simulated_data = occupancy
  ))
}
