#' Title Plot Posterior Distribution from Bootstrap Results
#'
#' @param bootdata A data frame or vector containing the bootstrap estimates of comparator and reference effects.
#' @param ... Additional graphical parameters passed to the plot function.
#'
#' @return a graph
#' @export
#'
plot_posterior <- function(bootdata, ...) {
  # Validate input
  if (!is.data.frame(bootdata) || !("effect_comparator" %in% names(bootdata)) || !("effect_ref_int" %in% names(bootdata))) {
    stop("bootdata must be a data frame containing 'effect_comparator' and 'effect_ref_int' columns.")
  }

  # Extract the relevant columns
  effect_comparator <- bootdata$effect_comparator
  effect_ref_int <- bootdata$effect_ref_int

  # Calculate densities
  density_comparator <- density(effect_comparator)
  density_ref_int <- density(effect_ref_int)

  # Plotting
  plot(density_comparator, col = "blue", lwd = 2, main = "Posterior Distributions of Effects", xlim = range(c(density_comparator$x, density_ref_int$x)), ylim = range(c(density_comparator$y, density_ref_int$y)), xlab = "Effect", ylab = "Density", ...)
  lines(density_ref_int, col = "red", lwd = 2)

  # Adding a legend
  legend("topright", legend = c("Comparator", "Reference"), col = c("blue", "red"), lwd = 2)
}



# Test
plot_posterior(model1$bootdata)


