#' Title Plot Average Treatment Effect (ATE) Density from Bootstrap Results
#'
#' @param bootdata A data frame or vector containing the bootstrap estimates of ATE.
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
plot_ATE <- function(bootdata,
                     col_density = "blue",
                     fill_density = "lightblue",
                     main = "Density of ATE Estimates",
                     xlab = "ATE", ylab = "Density",
                     xlim = NULL, ylim = NULL, ...) {
  # Check if bootdata is a data frame with a specific column or a vector
  if (is.data.frame(bootdata) && "ATE" %in% names(bootdata)) {
    ate_values <- bootdata$ATE
  } else if (is.vector(bootdata)) {
    ate_values <- bootdata
  } else {
    stop("bootdata must be a vector of ATE estimates or a data frame with an 'ATE' column.")
  }

  # Calculate the density of ATE estimates
  ate_density <- density(ate_values)

  # Plot the density
  plot(ate_density, col = col_density, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  polygon(ate_density, col = fill_density, border = col_density)

  # Add a vertical line at the mean ATE value
  abline(v = mean(ate_values), col = "red", lwd = 2, lty = 2)

  # Add a legend
  legend("topright", legend = c("ATE Density", "Mean ATE"), fill = c(fill_density, NA), border = c(col_density, "red"), lty = c(NA, 2), lwd = c(NA, 2))
}



# Test
plot_ATE(model1$bootdata)
# or
plot_ATE(model1$bootdata$ATE)

# Compared to:
plot(density(model1$bootdata$ATE))

