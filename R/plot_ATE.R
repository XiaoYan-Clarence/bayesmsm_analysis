#' Title Plot Average Treatment Effect (ATE) Density from Bootstrap Results
#'
#' @param input A data frame or vector containing the bootstrap estimates of ATE.
#' @param col_density Color for the density plot (default is "blue").
#' @param fill_density Fill color for the density plot (default is "lightblue").
#' @param main Title of the plot (default is "Density of ATE Estimates").
#' @param xlab X-axis label (default is "ATE").
#' @param ylab Y-axis label (default is "Density").
#' @param xlim Limits for the x-axis (default is NULL).
#' @param ylim Limits for the y-axis (default is NULL).
#' @param ... Additional graphical parameters passed to the plot function.
#'
#' @export
#'
plot_ATE <- function(input,
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Posterior Predictive Distribution of Average Treatment Effect (ATE）",
                     xlab = "ATE", ylab = "Posterior Predictive Distribution",
                     xlim = NULL, ylim = NULL, ...) {
  # Check if input is either a data frame or part of a model object
  if (is.list(input) && "bootdata" %in% names(input)) {
    # If input is a list and has bootdata, check for ATE column within bootdata
    if ("ATE" %in% names(input$bootdata)) {
      ate_values <- input$bootdata$ATE
    } else {
      stop("bootdata within the model object must have an 'ATE' column.")
    }
  } else if (is.data.frame(input) && "ATE" %in% names(input)) {
    ate_values <- input$ATE
  } else if (is.vector(input)) {
    ate_values <- input
  } else {
    stop("input must be a vector of ATE estimates, a data frame, or a model object containing a 'bootdata' data frame with an 'ATE' column.")
  }

  # # Calculate the density of ATE estimates
  # ate_density <- density(ate_values)
  #
  # # Plot the density
  # plot(ate_density, col = col_density, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  # polygon(ate_density, col = fill_density, border = col_density)
  #
  # # Add a vertical line at the mean ATE value
  # abline(v = mean(ate_values), col = "red", lwd = 2, lty = 2)
  #
  # # Add a legend
  # legend("topright", legend = c("ATE Density", "Mean ATE"), fill = c(fill_density, NA), border = c(col_density, "red"), lty = c(NA, 2), lwd = c(NA, 2))

  ate_density <- density(ate_values)
  ci <- quantile(ate_values, probs = c(0.025, 0.975))
  density_ci <- density(ate_values, from = ci[1], to = ci[2])

  plot(ate_density, col = col_density, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  polygon(c(density_ci$x, rev(density_ci$x)), c(rep(min(ate_density$y), length(density_ci$x)), rev(density_ci$y)), col = rgb(0, 0, 1, alpha = 0.3))
  abline(v = mean(ate_values), col = "purple", lwd = 2, lty = 3)
  abline(v = ci[1], col = "darkgreen", lty = 2)
  abline(v = ci[2], col = "darkgreen", lty = 2)

  legend_text <- c("ATE Density",
                   paste("Mean:", round(mean(ate_values), 3)),
                   paste("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]"))

  legend("topright", legend = legend_text,
         col = c(col_density, "purple", "darkgreen"),
         lwd = 2, lty = c(1, 3, 2))
}



# Test
plot_ATE(model2)
# or
plot_ATE(model2$bootdata)
# or
plot_ATE(model2$bootdata$ATE)



# Compared to:
plot(density(model1$bootdata$ATE))

