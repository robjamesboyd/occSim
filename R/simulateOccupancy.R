#' Simulate Species Occupancy Over Time
#'
#' This function simulates species occupancy over time for a given number of sites and years,
#' incorporating explanatory variables and probabilities for extinction and colonization.
#'
#' @param nYears Integer. The number of years over which to simulate occupancy.
#' @param nSites Integer. The number of sites to include in the simulation.
#' @param explanatoryVar Numeric vector. A vector of explanatory variables influencing occupancy, which must have length nSite.
#' @param beta0_extinction Numeric. Intercept for the extinction probability.
#' @param beta0_colonization Numeric. Intercept for the colonization probability.
#' @param beta0_initial Numeric. Intercept for the initial occupancy probability.
#' @param beta_extinction Numeric. Slope for the extinction probability.
#' @param beta_colonization Numeric. Slope for the colonization probability.
#' @param beta_initial Numeric. Slope for the initial occupancy probability.
#' @param randomInitial Logical. Whether to generate initial occupancy randomly. If false, then beta parameters must be set.
#' @param randomInitialPsi Numeric. If `randomInitial` is `TRUE`, the probability of initial occupancy.
#' @param nBins Integer. The number of bins to divide the explanatory variable into.
#' @param plot Character. Specifies the plot type to generate: "occupancy" for overall occupancy or "bin" for occupancy by bin.
#'
#' @return A list containing:
#' \itemize{
#'   \item `occupancy_all_sites`: A data frame of overall occupancy probabilities over time.
#'   \item `occupancy_per_bin`: A data frame of bin-wise occupancy probabilities over time.
#'   \item `simulated_data`: A matrix of simulated occupancy states for each site and year.
#' }
#'
#' @examples
#'sim_data <- simulateOccupancy(
#'  nYears = 10,
#'  nSites = 500,
#'  explanatoryVar = runif(500, min = 0, max = 1),
#'  beta0_extinction = -1,
#'  beta0_colonization = -3,
#'  beta0_initial = NULL,
#'  beta_extinction = -2,
#'  beta_colonization = 2,
#'  beta_initial = NULL,
#'  randomInitial = TRUE,
#'  randomInitialPsi = 0.3,
#'  nBins = 5,
#'  plot = "occupancy"
#')
#' 
#' true_occupancy <- sim_data$occupancy_all_sites
#' occupancy_by_bin <- sim_data$occupancy_per_bin
#'
#' # Further analysis or visualization can follow
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line theme_linedraw labs
#' @export

simulateOccupancy <- function(nYears, 
                              nSites, 
                              explanatoryVar,
                              beta0_extinction,
                              beta0_colonization,
                              beta0_initial = NULL,
                              beta_extinction,
                              beta_colonization,
                              beta_initial = NULL,
                              randomInitial = FALSE,
                              randomInitialPsi = NULL,
                              nBins,
                              plot = "occupancy") {
  
  # Calculate extinction and colonization probabilities
  # These include both intercept and slope components
  epsilon <- 1 / (1 + exp(-(beta0_extinction + beta_extinction * explanatoryVar)))
  gamma <- 1 / (1 + exp(-(beta0_colonization + beta_colonization * explanatoryVar)))

  # initial occupancy
  if (randomInitial) {
    occupancy <- sample(c(0, 1), nSites, replace = TRUE, prob = c(1-randomInitialPsi, randomInitialPsi))
  } else {
    initial_occupancy_prob <- 1 / (1 + exp(-(beta0_initial + beta_initial * explanatoryVar)))
    occupancy <- rbinom(nSites, 1, initial_occupancy_prob)
  }

  # Initialize results matrix
  results <- matrix(nrow = nYears, ncol = nSites)
  results[1, ] <- occupancy

  # Simulation over time using updated probabilities
  for (t in 2:nYears) {
    for (i in 1:nSites) {
      if (results[t-1, i] == 1) {
        # Extinction check
        results[t, i] = rbinom(1, 1, 1 - epsilon[i])
      } else {
        # Colonization check
        results[t, i] = rbinom(1, 1, gamma[i])
      }
    }
  }

  true_occupancy <- data.frame(occupancy = rowMeans(results),
                               year = 1:nYears)

  # Discretize explanatoryVar
  bins <- cut(explanatoryVar, breaks = nBins, labels = FALSE)

  # Prepare a data frame to store results by bin for each time step
  occupancy_by_bin <- data.frame(matrix(ncol = nBins, nrow = nrow(results)))
  
  colnames(occupancy_by_bin) <- paste("Bin", 1:nBins, sep = "_")
  
  # Step 2: Calculate mean occupancy for each bin at each time step
  for (i in 1:nrow(results)) {
    for (j in 1:nBins) {
      # Select columns (sites) that fall into the current bin
      sites_in_bin <- which(bins == j)
      # Calculate mean occupancy for these sites at this time step
      occupancy_by_bin[i, j] <- mean(results[i, sites_in_bin])
    }
  }
  
  # Add a column for the year indices
  occupancy_by_bin$Year <- 1:nrow(occupancy_by_bin)
  occupancy_long <- reshape2::melt(occupancy_by_bin, id.vars = "Year", variable.name = "Bin", value.name = "Occupancy")
  
  # Create the plot
  if (plot == "bin") {
    print(
      ggplot2::ggplot(occupancy_long, ggplot2::aes(x = Year, y = Occupancy, color = Bin)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "",
             x = "Year",
             y = "Occupancy") +
        ggplot2::theme_linedraw()
    )

  } else {
    print(
     ggplot2::ggplot(true_occupancy, ggplot2::aes(x = year, y = occupancy)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "",
             x = "Year",
             y = "Occupancy") +
        ggplot2::theme_linedraw()
    )
  }
  
  return(list(occupancy_all_sites = true_occupancy,
              occupancy_per_bin = occupancy_by_bin,
              simulated_data = results))
}
