plot_est_box <- function(model, ...) {
  # Extract bootdata from the model or use the data frame directly
  bootdata <- if (is.data.frame(model)) {
    model
  } else if ("bootdata" %in% names(model)) {
    model$bootdata
  } else {
    stop("Input must be a data frame or a model object containing 'bootdata'.")
  }

  # Identify if the family is binomial to include RR and OR
  is_binomial <- "RR" %in% names(bootdata) && "OR" %in% names(bootdata)

  # Validate bootdata
  required_columns <- c("effect_comparator", "effect_reference", "RD")
  if (is_binomial) {
    required_columns <- c(required_columns, "RR", "OR")
  }
  if (!all(required_columns %in% names(bootdata))) {
    stop("bootdata must contain required columns.")
  }

  # Calculate means and confidence intervals
  means <- sapply(bootdata[required_columns], mean)
  lowerbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.025))
  upperbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.975))

  # Define the position for each point
  position <- 1:length(means)

  # Define some offsets for text placement
  text_offset <- (max(upperbd) - min(lowerbd)) * 0.05

  # Plotting
  old_par <- par(no.readonly = TRUE) # Save current graphical parameters
  on.exit(par(old_par)) # Restore graphical parameters on exit
  par(mar = c(7, 5, 5, 2) + 0.1) # Adjust margins to avoid "figure margins too large" error

  plot(position, means, ylim = range(lowerbd - text_offset, upperbd + text_offset), pch = 19, xaxt = "n",
       xlab = "Treatment Level", ylab = "Effect", main = "Treatment Effect Estimates", ...)
  axis(1, at = position, labels = if (is_binomial) c("Comparator Level", "Reference Level", "RD", "RR", "OR") else c("Comparator Level", "Reference Level", "RD"))

  # Error bars
  arrows(position, lowerbd, position, upperbd, angle = 90, code = 3, length = 0.1, ...)

  # Adding text for means and CIs
  for (i in seq_along(means)) {
    text(position[i], upperbd[i] + text_offset, labels = paste("Mean:", round(means[i], 2), "\n95% CI: [", round(lowerbd[i], 2), ", ", round(upperbd[i], 2), "]"), cex = 0.8, pos = 3)
  }

  # Check if the input is a model and extract treatment sequences if they exist
  has_treatment_info <- "reference" %in% names(model) && "comparator" %in% names(model)

  # Conditional treatment sequence information below x-axis labels
  if (has_treatment_info) {
    mtext(paste("(", paste(model$reference, collapse = ", "), ")", sep = ""), side = 1, at = position[2], line = 2)
    mtext(paste("(", paste(model$comparator, collapse = ", "), ")", sep = ""), side = 1, at = position[1], line = 2)
  }
}


# Example usage:
# Assume model is the output from bayesmsm
plot_est_box(model2)
plot_est_box(model3)




# plot_est_box <- function(input, ...) {
#   # Extract bootdata from the model or use the data frame directly
#   bootdata <- if (is.data.frame(input)) {
#     input
#   } else if ("bootdata" %in% names(input)) {
#     input$bootdata
#   } else {
#     stop("Input must be a data frame or a model object containing 'bootdata'.")
#   }
#
#   # Validate bootdata
#   required_columns <- c("effect_comparator", "effect_reference", "ATE")
#   if (!all(required_columns %in% names(bootdata))) {
#     stop("bootdata must contain 'effect_comparator', 'effect_reference', and 'ATE' columns.")
#   }
#
#   # Adjust margins if necessary
#   par(mar = c(5, 4, 4, 3) + 0.1) # Adjust the last value if text is plotted outside; bottom, left, top, and right margins
#
#   # Calculate means and standard deviations
#   means <- sapply(bootdata[required_columns], mean)
#   # sds <- sapply(bootdata[required_columns], sd)
#   # ses <- sds / sqrt(nrow(bootdata))
#   lowerbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.025))
#   upperbd <- sapply(bootdata[required_columns], function(x) quantile(x, probs = 0.975))
#
#   # Define the position for each point
#   position <- 1:length(means)
#
#   # Define some offsets for text placement
#   text_offset <- (max(upperbd) - min(lowerbd)) * 0.05
#
#   # Plotting
#   plot(position, means, ylim = range(lowerbd - text_offset, upperbd + text_offset), pch = 19, xaxt = "n", # round down vs round up;
#        xlab = "Treatment Level", ylab = "Effect", main = "Treatment Effect Estimates", ...)
#   axis(1, at = position, labels = c("Comparator Level", "Reference Level", "ATE"))
#
#   # Error bars
#   arrows(position, lowerbd, position, upperbd, angle = 90, code = 3, length = 0.1, ...)
#
#   # Adding text for means and CIs
#   # text(position, lowerbd - 0.1, labels = paste("Mean:", round(means, 2), "\n95% CI:", round(lowerbd, 2), "-", round(upperbd, 2)), adj = c(0,1))
#
#   # for (i in seq_along(means)) {
#   #   text(position[i], upperbd[i] + text_offset, labels = paste("Mean:", round(means[i], 2), "\n95% CI:", round(lowerbd[i], 2), "-", round(upperbd[i], 2)), cex = 0.8, pos = 3)
#   # }
#
#   # Adjust the y-position offset for clarity
#   offset <- 0.1
#   # Adding text for the first and second items with left adjustment
#   text(position[1:2], lowerbd[1:2] - offset, labels = paste("Mean:", round(means[1:2], 2), "\n95% CI: [", round(lowerbd[1:2], 2), ", ", round(upperbd[1:2], 2), "]"), adj = c(0,1), ...)
#   # Adding text for the third item with right adjustment
#   text(position[3], upperbd[3] + offset, labels = paste("Mean:", round(means[3], 2), "\n95% CI: [", round(lowerbd[3], 2), ", ", round(upperbd[3], 2), "]"), adj = c(1,0), ...)
#
#   # Check if the input is a model and extract treatment sequences if they exist
#   has_treatment_info <- "reference" %in% names(input) && "comparator" %in% names(input)
#
#   # Conditional treatment sequence information below x-axis labels
#   if (has_treatment_info) {
#     # mtext(paste("(", paste(input$reference, collapse = ", "), ")", sep = ""), side = 1, at = position[2], line = 2.5, cex = 0.7)
#     # mtext(paste("(", paste(input$comparator, collapse = ", "), ")", sep = ""), side = 1, at = position[1], line = 2.5, cex = 0.7)
#     mtext(paste("(", paste(input$reference, collapse = ", "), ")", sep = ""), side = 1, at = position[2], line = 2)
#     mtext(paste("(", paste(input$comparator, collapse = ", "), ")", sep = ""), side = 1, at = position[1], line = 2)
#   }
# }



# Test

plot_est_box(model2$bootdata) # without reference & comparator information below labels
plot_est_box(model2) # with reference & comparator information below labels

