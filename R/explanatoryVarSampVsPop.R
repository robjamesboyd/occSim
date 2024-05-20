#' Compare distribution of explanatoryVar in sample vs. population
#'
#' This function takes a sample inclusion matrix and an explanatory variable vector
#' and returns a ggplot object that shows the density distribution of the explanatory variable
#' for the population and each sample.
#'
#' @param samp Output of simulateSample.
#' @param explanatoryVar The explanatory variable used to simulate occupancy and the sample.
#' @return A ggplot object showing the density plots.
#' @import ggplot2
#' @export
#' @examples
#' samp <- list(sample_inclusion = matrix(sample(0:1, 100, replace = TRUE), nrow = 10))
#' explanatoryVar <- rnorm(100)
#' plot <- explanatoryVarSampVsPop(samp, explanatoryVar)
#' print(plot)

explanatoryVarSampVsPop <- function(samp, explanatoryVar) {
  n <- nrow(samp$sample_inclusion)
  
  # Create a data frame for the Population
  population_data <- data.frame(
    value = explanatoryVar,
    group = factor(rep("Population", length(explanatoryVar)))
  )
  
  # Initialize a list to store data frames for each sample
  list_of_dataframes <- list(population_data)
  
  # Loop over each row in samp$sample_inclusion to create data frames for samples
  for (i in 1:n) {
    sample_data <- explanatoryVar[samp$sample_inclusion[i,] == 1]
    df <- data.frame(
      value = sample_data,
      group = factor(rep(paste("Sample", i), length(sample_data)))
    )
    list_of_dataframes[[i + 1]] <- df
  }
  
  # Combine all data frames into one
  combined_data <- do.call(rbind, list_of_dataframes)
  
  # Plot the histograms using ggplot2
  plot <- ggplot2::ggplot(combined_data, ggplot2::aes(x = value, colour = group)) +
    ggplot2::geom_density() +
    ggplot2::xlim(-2, 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Comparison of Population and Sample Distributions",
                  x = "Value",
                  y = "Density",
                  colour = "")
  
  return(plot)
}

# Example usage:
# samp <- list(sample_inclusion = matrix(sample(0:1, 100, replace = TRUE), nrow = 10))
# explanatoryVar <- rnorm(100)
# plot <- explanatoryVarSampVsPop(samp, explanatoryVar)
# print(plot)
