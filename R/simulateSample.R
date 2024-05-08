#' Simulate Sampling Process for Occupancy Data
#'
#' This function simulates the sampling process for a given occupancy dataset over multiple years and sites,
#' based on logistic regression parameters. It applies logistic regression to predict sampling probabilities
#' using an explanatory variable. It returns a list that includes a dataframe summarizing the real and observed
#' occupancy rates per year, the matrix of observed occupancy data, and a matrix indicating whether each site
#' was sampled each year.
#'
#' @param occ Output of simulateOccupancy. A list containing at least one element named 'simulated_data', which should be a matrix where
#'   each entry denotes occupancy (1 = occupied, 0 = not occupied) for each site across several years.
#' @param beta0 Intercept (scalar) for the logistic regression model used to calculate sampling probabilities.
#' @param beta Coefficient (scalar) for the logistic regression model associated with the explanatory variable.
#' @param explanatoryVar A numeric vector or scalar of the explanatory variable used in logistic regression,
#'   matching the number of rows in `occ$simulated_data`.
#'
#' @return A list with three elements:
#'   - `summary`: A dataframe with columns `year`, `real`, and `obs` representing the year, the mean real occupancy,
#'     and the mean observed occupancy, respectively.
#'   - `data`: A matrix of observed occupancy values (1 or NA) matching the structure of the input `occ$simulated_data`.
#'   - `sample_inclusion`: A matrix indicating whether a site was sampled (1) or not (0) across all years.
#'
#' @examples
#' # Assuming 'occ' is a list with a matrix of simulated data
#' nYears <- 10
#' nSites <- 5
#' simulated_data <- matrix(rbinom(nYears * nSites, 1, 0.5), nrow = nYears)
#' occ <- list(simulated_data = simulated_data)
#' beta0 <- 0.5
#' beta <- -0.3
#' explanatoryVar <- runif(nYears)
#'
#' results <- simulateSample(occ, beta0, beta, explanatoryVar)
#' print(results$summary)
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line theme_linedraw labs
#' @export

simulateSample <- function(occ,
                           beta0,
                           beta,
                           explanatoryVar) {
  
  # extract simulated data from occ object
  occ <- occ$simulated_data
  
  # get dimensions
  nYears <- nrow(occ)
  nSites <- ncol(occ)

  # Calculate sampling probabilities
  sampling_probability <- 1 / (1 + exp(-(beta0 + beta * explanatoryVar)))

  # Initialize sampling results matrix (0 = not sampled, 1 = sampled)
  sampling_results <- matrix(nrow = nYears, ncol = nSites)

  # Simulate the sampling process over time
  for (t in 1:nYears) {
    # For each site, decide if it is sampled based on the sampling probability
    sampling_results[t, ] <- rbinom(nSites, 1, sampling_probability)
  }

  # Calculate the observed occupancy based on both the true occupancy and whether the site was sampled
  observed_occupancy <- matrix(nrow = nYears, ncol = nSites)

  # Fill in the observed occupancy based on both the true occupancy and whether the site was sampled
  for (t in 1:nYears) {
    for (i in 1:nSites) {
      if (sampling_results[t, i] == 1) {
        observed_occupancy[t, i] = occ[t, i]
      } else {
        observed_occupancy[t, i] = NA  # Mark as missing data
      }
    }
  }

  df <- data.frame(year = 1:nYears,
                   real = rowMeans(occ),
                   obs = rowMeans(observed_occupancy, na.rm = T))

  df_long <- reshape2::melt(df, id.vars = "year", measure.vars = c("real", "obs"),
                  variable.name = "type", value.name = "value")

  print(
    ggplot2::ggplot(data=df_long, ggplot2::aes(x=year, y=value, colour=type)) +
      ggplot2::geom_line() +
      ggplot2::theme_linedraw() +
      ggplot2::labs(x="Year", y="Occupancy", colour="")
  )
  
  return(list(summary=df,
              data=observed_occupancy,
              sample_inclusion=sampling_results))
}
