summary.bayesmsm <- function(model) {
  # Extract bootstrapped data
  bootdata <- model$bootdata

  # Calculate summary statistics for each metric
  summary_stats <- function(metric) {
    mean_val <- mean(bootdata[[metric]])
    sd_val <- sd(bootdata[[metric]])
    quantiles <- quantile(bootdata[[metric]], probs = c(0.025, 0.975))

    return(c(mean = mean_val, sd = sd_val, `2.5%` = quantiles[1], `97.5%` = quantiles[2]))
  }

  # Initialize a list to store summary statistics
  results_list <- list()

  # Add summary statistics for reference and comparator
  results_list$Reference <- summary_stats("effect_reference")
  results_list$Comparator <- summary_stats("effect_comparator")

  # Add summary statistics for RD
  results_list$RD <- summary_stats("RD")

  # Check if RR and OR exist in the bootdata and add them if they do
  if ("RR" %in% colnames(bootdata)) {
    results_list$RR <- summary_stats("RR")
  }

  if ("OR" %in% colnames(bootdata)) {
    results_list$OR <- summary_stats("OR")
  }

  # Convert the list to a data frame for better presentation
  summary_table <- do.call(rbind, results_list)

  # Set the column names explicitly
  colnames(summary_table) <- c("mean", "sd", "2.5%", "97.5%")

  # Return the summary table
  return(summary_table)
}

# Example usage
summary.bayesmsm(model2)
summary.bayesmsm(model3)
