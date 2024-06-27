#' Title Plot Posterior Distribution from Bootstrap Results
#'
#' @param bootdata A data frame or vector containing the bootstrap estimates of comparator and reference effects.
#' @param ... Additional graphical parameters passed to the plot function.
#'
#' @return a graph
#' @export
#'
plot_APO <- function(model, effect_type,
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Average Potential Outcome (APO)",
                     xlab = "Effect", ylab = "Density",
                     xlim = NULL, ylim = NULL, ...) {

  # Validate input
  if ("bootdata" %in% names(model)) {
    bootdata <- model$bootdata
  } else if (is.data.frame(model)) {
    bootdata <- model
  } else {
    stop("Input must be a data frame or a model object containing a 'bootdata' data frame.")
  }

  if (!is.data.frame(bootdata) || !("effect_comparator" %in% names(bootdata)) || !("effect_reference" %in% names(bootdata))) {
    stop("bootdata must be a data frame containing 'effect_comparator' and 'effect_reference' columns.")
  }

  if (!is.character(effect_type) || length(effect_type) != 1) {
    stop("effect_type must be a single character string specifying the effect to plot.")
  }

  if (!effect_type %in% c("effect_comparator", "effect_reference")) {
    stop("effect_type must be either 'effect_comparator' or 'effect_reference'.")
  }

  # Extract the relevant column
  effect <- bootdata[[effect_type]]

  # Calculate density
  density_effect <- density(effect)

  # Define titles and colors based on effect_type
  titles <- c(effect_comparator = "Comparator Level", effect_reference = "Reference Level")
  colors <- c(effect_comparator = "blue", effect_reference = "red")

  # Calculate mean and CI
  mean_effect <- mean(effect)
  ci <- quantile(effect, probs = c(0.025, 0.975))
  density_ci <- density(effect, from = ci[1], to = ci[2])

  # Set layout to allocate space for legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

  # Plot the density
  par(mar = c(5, 4, 4, 0)) # Adjust margins for the plot
  plot(density_effect, col = col_density, main = paste(main, titles[effect_type], sep = ": "), xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(density_effect$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))
  abline(v = mean_effect, col = "purple", lwd = 2, lty = 3)
  abline(v = ci[1], col = "darkgreen", lty = 2)
  abline(v = ci[2], col = "darkgreen", lty = 2)

  # Plot the legend in a separate plot
  par(mar = c(5, 0, 4, 0)) # Adjust margins for the legend
  plot.new()
  legend("center", legend = c(titles[effect_type],
                              paste("Mean:", round(mean_effect, 3)),
                              paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]")),
         col = c(col_density, "purple", "darkgreen"),
         lwd = 2, lty = c(1, 3, 2),
         bty = "n")
}

# Usage example:
# Assume model is the output from bayesmsm
# Plot for effect_comparator
plot_APO(model2, effect_type = "effect_comparator")
# Plot for effect_reference
plot_APO(model2, effect_type = "effect_reference")



# plot_APO <- function(input, effect_type, ...) {
#   # Validate input
#   if ("bootdata" %in% names(input)) {
#     bootdata <- input$bootdata
#   } else if (is.data.frame(input)) {
#     bootdata <- input
#   } else {
#     stop("Input must be a data frame or a model object containing a 'bootdata' data frame.")
#   }
#   if (!is.data.frame(bootdata) || !("effect_comparator" %in% names(bootdata)) || !("effect_reference" %in% names(bootdata))) {
#     stop("bootdata must be a data frame containing 'effect_comparator' and 'effect_reference' columns.")
#   }
#   if (!is.character(effect_type) || length(effect_type) != 1) {
#     stop("effect_type must be a single character string specifying the effect to plot.")
#   }
#   if (!effect_type %in% c("effect_comparator", "effect_reference")) {
#     stop("effect_type must be either 'effect_comparator' or 'effect_reference'.")
#   }
#
#   # Extract the relevant column
#   effect <- bootdata[, effect_type, drop = FALSE]
#
#   # Calculate density
#   density_effect <- density(effect[[1]])
#
#   # Define titles and colors based on effect_type
#   titles <- c(effect_comparator = "Comparator Level", effect_reference = "Reference Level")
#   colors <- c(effect_comparator = "blue", effect_reference = "red")
#
#   # Calculate mean and median
#   mean_effect <- mean(effect[[1]])
#   # median_effect <- median(effect[[1]])
#
#   # Calculate CI
#   ci <- quantile(effect[[1]], probs = c(0.025, 0.975))
#   density_ci <- density(effect[[1]], from = ci[1], to = ci[2])
#
#   # Plot the density curve
#   plot(density_effect, main = paste("Average Potential Outcome (APO) of", titles[effect_type]), xlab = "Effect", ylab = "Density", col = colors[effect_type], lwd = 2, ...)
#
#   # Shade the area under the curve within the 95% CI
#   polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(density_effect$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))
#
#   # Add a vertical lines for the mean and 95% CI bounds
#   abline(v = mean_effect, col = "purple", lty = 3)
#   abline(v = ci[1], col = "darkgreen", lty = 2)
#   abline(v = ci[2], col = "darkgreen", lty = 2)
#
#   # Legend with mean and median
#   legend_text <- c(titles[effect_type],
#                    paste("Mean:", round(mean_effect, 3)),
#                    # ,paste("Median:", round(median_effect, 3)
#                    paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]"))
#
#   # Update the legend text to include the 95% CI
#   # legend_text <- c(legend_text, paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]"))
#
#   legend("topright", legend = legend_text,
#          col = c(colors[effect_type], "purple"
#          # , "darkgrey"
#          , "darkgreen"),
#          lwd = 2, lty = c(1, 3, 2))
# }



# Test
plot_APO(model1$bootdata, effect_type = "effect_comparator")
plot_APO(model1, effect_type = "effect_reference")
# plot_APO(model1$bootdata)


