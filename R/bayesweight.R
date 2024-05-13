#' Title
#'
#' @param data
#' @param n.iter
#' @param n.burnin
#' @param n.thin
#' @param formula.list a list of formulas corresponding to each time point with the time-specific treatment variable on the left hand side and pre-treatment covariates to be balanced on the right hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point (i.e., treatments and covariates featured in early formulas should appear in later ones). Interactions and functions of covariates are allowed.
#'
#' @return
#' @export
#'
#' @examples
bayesweight <- function(formula.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                            a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = causaldata,
                        n.iter = 2500,
                        n.burnin = 1500,
                        n.thin = 5){

  # Testing
  formula.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                      a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1)
  causaldata <- read.csv("continuous_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
  data = causaldata
  n.iter = 250
  n.burnin = 50
  n.thin = 5

  extract_variables_list <- function(input) {
    # Function to extract variables from a single formula
    extract_from_formula <- function(formula) {
      # Get the terms of the formula
      formula_terms <- terms(formula)

      # Extract the response variable name (if there is one)
      response_variable <- attr(formula_terms, "response")
      response_name <- if (response_variable > 0) {
        all_vars <- all.vars(formula)
        all_vars[response_variable]
      } else {NA}

      # Extract predictor variable names
      predictor_names <- attr(formula_terms, "term.labels")

      # Return a list of response and predictor variables
      list(response = response_name, predictors = predictor_names)
    }

    # Check if input is a list of formulas
    if (is.list(input)) {
      # Apply the function to each formula in the list
      lapply(input, extract_from_formula)
    } else {
      # Input is a single formula
      extract_from_formula(input)
    }
  }

  # Test `extract_variables_list`
  # extract_variables_list(a_1 ~ w1 + w2 + L1_1 + L2_1)
  # extract_variables_list(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
  #                        a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1))
  # str(extract_variables_list(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
  #                            a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1)))

  # Function to generate and write the JAGS model
  write_jags_model <- function(formula.list) {
    # Extract variable information for each formula
    var_info <- lapply(formula.list, extract_variables_list)

    # Start writing the model string
    model_string <- "model{\n#N = nobs\nfor(i in 1:N){\n"

    # Process each visit
    for (v in seq_along(var_info)) {
      visit <- var_info[[v]]
      response <- visit$response
      predictors <- visit$predictors

      # Marginal treatment assignment model
      model_string <- paste0(model_string,
                             "\n# visit ", v, ";\n",
                             "# marginal treatment assignment model, visit ", v, ";\n",
                             response, "s[i] ~ dbern(p", response, "s[i])\n")

      if (v == 1) {
        model_string <- paste0(model_string, "p", response, "s[i] <- ilogit(bs", v, "0)\n")
      } else {
        previous_response_s = var_info[[v - 1]]$response # Get previous visit's response
        model_string <- paste0(model_string, "p", response, "s[i] <- ilogit(bs", v, "0 + bs", v, "1*", previous_response_s, "s[i])\n")
      }

      # Conditional treatment assignment model
      model_string <- paste0(model_string,
                             "\n# conditional treatment assignment model, visit ", v, ";\n",
                             response, "[i] ~ dbern(p", response, "[i])\n",
                             "p", response, "[i] <- ilogit(b", v, "0")

      for (p in seq_along(predictors)) {
        model_string <- paste0(model_string, " + b", v, p, "*", predictors[p], "[i]")
      }
      model_string <- paste0(model_string, ")\n")
    }

    # Close loop and add priors
    model_string <- paste0(model_string, "\n# export quantity in full posterior specification;\n",
                           "w[i] <- (", paste(sapply(seq_along(var_info), function(x) paste0("p", var_info[[x]]$response, "s[i]")), collapse = "*"),
                           ")/(", paste(sapply(seq_along(var_info), function(x) paste0("p", var_info[[x]]$response, "[i]")), collapse = "*"), ")\n}\n\n#prior;\n")

    # Add priors for all parameters
    for (v in seq_along(var_info)) {
      visit <- var_info[[v]]
      predictors <- visit$predictors
      num_preds <- length(predictors) + 1  # +1 for intercept

      # bs parameters
      if (v == 1) {
        model_string <- paste0(model_string, "bs", v, "0~dnorm(0,.01)\n")
      } else {
        model_string <- paste0(model_string, "bs", v, "0~dnorm(0,.01)\nbs", v, "1~dnorm(0,.01)\n")
      }

      # b parameters
      for (p in 0:(num_preds - 1)) {
        model_string <- paste0(model_string, "b", v, p, "~dnorm(0,.01)\n")
      }
    }

    # Add the closing brace for the model block
    model_string <- paste(model_string, "}\n", sep="")

    # Example initializations for the example; this part needs to be your actual model code generation
    all_parameters <- c()

    # Process each visit
    for (v in seq_along(formula.list)) {
      var_info <- extract_variables(formula.list[[v]])
      response <- var_info$response
      predictors <- var_info$predictors

      # Add bs parameters for marginal models
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        all_parameters <- c(all_parameters, sprintf("bs%d1", v))
      }

      # Add b parameters for conditional models
      all_parameters <- c(all_parameters, sprintf("b%d0", v))  # intercept
      for (p in seq_along(predictors)) {
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Assume model_string is being built here...
      # This is where you'd actually generate the JAGS model code
    }

    # Example writing the model to a file
    cat(model_string, file = "treatment_model.txt")

    return(unique(all_parameters))
  }

  # Test `write_jags_model`
  formula.list <- list(
    a_1 ~ w1 + w2 + L1_1 + L2_1,
    a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1
  )
  #
  # formula.list3 <- list(
  #   a_1 ~ w1 + w2 + L1_1 + L2_1,
  #   a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  #   a_3 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2
  # )
  #
  # formula.list4 <- list(
  #   a_1 ~ w1 + w2 + L1_1 + L2_1,
  #   a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  #   a_3 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2,
  #   a_4 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2 + L1_4 + L2_4 + a_3
  # )

  # Generate and write the JAGS model
  write_jags_model(formula.list)
  # write_jags_model(formula.list3)
  # write_jags_model(formula.list4)

  extract_variables <- function(formula) {
    # Get the terms of the formula
    formula_terms <- terms(formula)

    # Extract the response variable name (if there is one)
    response_variable <- attr(formula_terms, "response")
    response_name <- if (response_variable > 0) {
      all_vars <- all.vars(formula)
      all_vars[response_variable]
    } else {NA}

    # Extract predictor variable names
    predictor_names <- attr(formula_terms, "term.labels")

    # Return a list of response and predictor variables
    list(response = response_name, predictors = predictor_names)
  }

  # Prepare data and parameters for JAGS

  prepare_jags_data <- function(data, formula.list) {
    # Extract variable information from formulas
    variable_info <- lapply(formula.list, extract_variables)

    # Initialize the list to pass to JAGS
    jags.data <- list(N = nrow(data))

    # Loop over each formula to add necessary variables
    for (info in variable_info) {
      response <- info$response
      predictors <- info$predictors

      # Add response and its 's' version
      jags.data[[response]] <- data[[response]]
      jags.data[[paste0(response, "s")]] <- data[[response]]  # Assume treatment assignment variable is same as response

      # Add predictors
      for (predictor in predictors) {
        jags.data[[predictor]] <- data[[predictor]]
      }
    }

    return(jags.data)
  }

  # Use this function to prepare your JAGS data
  jags.data <- prepare_jags_data(data, formula.list)
  # jags.data <- as.list(data)
  # jags.data$N <- nrow(data)
  jags.params <- write_jags_model(formula.list)

  # Run JAGS model
  set.seed(890123) # Is this needed?

  jagsfit <- jags(data = jags.data,
                  parameters.to.save = jags.params,
                  n.iter = n.iter,
                  model.file = "treatment_model.txt",
                  n.chains = 1,
                  n.burnin = n.burnin,
                  n.thin = n.thin)

  # Extract MCMC output
  out.mcmc <- as.mcmc(jagsfit)








  ################################# THE CODE BELOW IS NOT YET DONE #################################

  # Step 1: Parameter matrix (use appropriate parameter range)
  # parameter_matrix <- as.matrix(out.mcmc[[1]][, grep("^b[0-9]+[0-9]+$", colnames(out.mcmc[[1]]))])

  parameter_names <- colnames(out.mcmc[[1]])  # Directly pull parameter names
  parameter_matrix <- as.matrix(out.mcmc[[1]][, parameter_names])  # Select all parameters directly

  # Step 2: Build the Design Matrix from data based on model formulas
  # Let's assume `formula.list` is a list of formulas used for each model component.
  prepare_design_matrix <- function(data, formula.list) {
    # Extract unique predictor names from all formulas
    variable_names <- unique(unlist(lapply(formula.list, function(f) {
      all_vars <- all.vars(f)
      # Remove the response variable if needed, assuming it might be on the left of ~ in any formula
      return(all_vars[!all_vars %in% all.vars(formula(f, rhs = 1))])
    })))

    # Check if all variables are present in the dataframe
    if (!all(variable_names %in% names(data))) {
      missing_vars <- variable_names[!variable_names %in% names(data)]
      stop("The following required variables are missing from the data: ", paste(missing_vars, collapse = ", "))
    }

    # Build the design matrix with an intercept added
    design_formula <- as.formula(paste("~", paste(variable_names, collapse = " + "), "+ 0"))
    design_matrix <- model.matrix(design_formula, data)
    intercept <- matrix(1, nrow = nrow(design_matrix), ncol = 1)  # Create intercept column
    design_matrix <- cbind(intercept, design_matrix)  # Add the intercept column manually

    return(design_matrix)
  }

  design_matrix <- prepare_design_matrix(data, formula.list)

  # Step 3: Calculate the Weight Matrix
  # Ensure transposition is correctly done to match dimensions
  weight_matrix <- parameter_matrix %*% t(design_matrix)  # Correctly use matrix multiplication

  # Step 4: Calculate Mean Weights across all posteriors
  individual_weights <- rowMeans(weight_matrix)  # Get mean across rows for each posterior

  # Return the calculated individual weights
  return(individual_weights)

}


# Example
causaldata <- read.csv("continuous_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
bayesweight(formula.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
            data = causaldata,
            n.iter = 250,
            n.burnin = 50,
            n.thin = 5)
