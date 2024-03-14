plot_est_box <- function(input, ...) {
  # Extract bootdata from the model or use the data frame directly
  bootdata <- if (is.data.frame(input)) {
    input
  } else if ("bootdata" %in% names(input)) {
    input$bootdata
  } else {
    stop("Input must be a data frame or a model object containing 'bootdata'.")
  }

  # Validate bootdata
  required_columns <- c("effect_comparator", "effect_ref_int", "ATE")
  if (!all(required_columns %in% names(bootdata))) {
    stop("bootdata must contain 'effect_comparator', 'effect_ref_int', and 'ATE' columns.")
  }

  # Calculate means and standard deviations
  means <- sapply(bootdata[required_columns], mean)
  sds <- sapply(bootdata[required_columns], sd)
  ses <- sds / sqrt(nrow(bootdata))

  # Define the position for each point
  position <- 1:length(means)


  # Plotting
  plot(position, means, ylim = range(means - ses * 1.96, means + ses * 1.96), pch = 19, xaxt = "n",
       xlab = "Treatment Level", ylab = "Effect", main = "Treatment Effect Estimates", ...)
  axis(1, at = position, labels = c("Comparator Level", "Reference Level", "ATE"))

  # Error bars
  arrows(position, means - ses * 1.96, position, means + ses * 1.96, angle = 90, code = 3, length = 0.1, ...)

}



# Test

plot_est_box(model1$bootdata)
plot_est_box(model1)
