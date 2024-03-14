#' Title Plot Posterior Distribution from Bootstrap Results
#'
#' @param bootdata A data frame or vector containing the bootstrap estimates of comparator and reference effects.
#' @param ... Additional graphical parameters passed to the plot function.
#'
#' @return a graph
#' @export
#'
plot_posterior <- function(input, effect_type, ...) {
  # Validate input
  if ("bootdata" %in% names(input)) {
    bootdata <- input$bootdata
  } else if (is.data.frame(input)) {
    bootdata <- input
  } else {
    stop("Input must be a data frame or a model object containing a 'bootdata' data frame.")
  }

  if (!is.data.frame(bootdata) || !("effect_comparator" %in% names(bootdata)) || !("effect_ref_int" %in% names(bootdata))) {
    stop("bootdata must be a data frame containing 'effect_comparator' and 'effect_ref_int' columns.")
  }

  if (!is.character(effect_type) || length(effect_type) != 1) {
    stop("effect_type must be a single character string specifying the effect to plot.")
  }

  if (!effect_type %in% c("effect_comparator", "effect_ref_int")) {
    stop("effect_type must be either 'effect_comparator' or 'effect_ref_int'.")
  }

  # Extract the relevant column
  effect <- bootdata[, effect_type, drop = FALSE]

  # Calculate density
  density_effect <- density(effect[[1]])

  # Define titles and colors based on effect_type
  titles <- c(effect_comparator = "Comparator Level", effect_ref_int = "Reference Level")
  colors <- c(effect_comparator = "blue", effect_ref_int = "red")

  # Plotting
  plot(density_effect, main = paste("Average Potential Outcome (APO) of", titles[effect_type]), xlab = "Effect", ylab = "Density", col = colors[effect_type], lwd = 2, ...)

  # Legend
  legend("topright", legend = titles[effect_type], col = colors[effect_type], lwd = 2)
}



# Test
plot_posterior(model1$bootdata, effect_type = "effect_comparator")
plot_posterior(model1, effect_type = "effect_ref_int")
# plot_posterior(model1$bootdata)

