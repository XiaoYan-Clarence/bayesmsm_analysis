#' Title
#'
#' @param data
#' @param n.iter
#' @param n.burnin
#' @param n.thin
#' @param trtmodel.list a list of formulas corresponding to each time point with the time-specific treatment variable on the left hand side and pre-treatment covariates to be balanced on the right hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point (i.e., treatments and covariates featured in early formulas should appear in later ones). Interactions and functions of covariates are allowed.
#'
#' @return
#' @export
#'
#' @examples
bayesweight <- function(trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = causaldata,
                        n.chains = 1,
                        n.iter = 25000,
                        n.burnin = 15000,
                        n.thin = 5,
                        seed = 890123,
                        parallel = FALSE){

  # Testing
  # trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
  #                     a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1)
  # data <- read.csv("continuous_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
  # n.chains = 2
  # n.iter = 250
  # n.burnin = 50
  # n.thin = 5
  # parallel = TRUE
  # seed = 890123

  # Load all the required R packages;
  if (!require(R2jags)){
    install.packages("R2jags",repos="http://cran.r-project.org")
    library(R2jags)
  }
  if (!require(coda)){
    install.packages("coda",repos="http://cran.r-project.org")
    library(coda)
  }
  if (!require(parallel)){
    install.packages("parallel",repos="http://cran.r-project.org")
    library(parallel)
  }

  # Test `write_jags_model` with different trtmodel.list
  #
  # trtmodel.list <- list(
  #   a_1 ~ w1 + w2 + L1_1 + L2_1,
  #   a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  #   a_3 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2
  # )
  # trtmodel.list_s <- list(
  # a_1 ~ 1,
  # a_2 ~ a_1,
  # a_3 ~ a_2 + a_1
  # )
  #
  # trtmodel.list <- list(
  #   a_1 ~ w1 + w2 + L1_1 + L2_1,
  #   a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
  #   a_3 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2,
  #   a_4 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2 + L1_4 + L2_4 + a_3
  # )

  # Test with different data
  # causaldata$a_3 <- rbinom(n=length(causaldata$y),1,p=0.4)
  # causaldata$a_4 <- rbinom(n=length(causaldata$y),1,p=0.6)
  # causaldata$L1_3 <- rbinom(n=length(causaldata$y),1,p=0.35)
  # causaldata$L2_3 <- rbinom(n=length(causaldata$y),1,p=0.65)
  # causaldata$L1_4 <- rbinom(n=length(causaldata$y),1,p=0.45)
  # causaldata$L2_4 <- rbinom(n=length(causaldata$y),1,p=0.55)
  # data = causaldata

  create_marginal_treatment_models <- function(trtmodel.list) {
    # Initialize the list for the marginal treatment models
    trtmodel.list_s <- list()

    # Loop through each model in the original list
    for (i in seq_along(trtmodel.list)) {
      # Extract the response variable (treatment variable) from each model
      response_var <- all.vars(trtmodel.list[[i]])[1]  # assuming the response is the first variable on the LHS

      # Create the marginal model formula
      if (i == 1) {
        # The first treatment model does not depend on any previous treatments
        formula_s <- as.formula(paste(response_var, "~ 1"))
      } else {
        # Subsequent treatment models depend on all previous treatments
        previous_treatments <- sapply(seq_len(i-1), function(j) {
          all.vars(trtmodel.list[[j]])[1]
        })
        formula_s <- as.formula(paste(response_var, "~", paste(previous_treatments, collapse = " + ")))
      }

      # Append the new formula to the list
      trtmodel.list_s[[i]] <- formula_s
    }

    return(trtmodel.list_s)
  }

  # Generate trtmodel.list_s
  trtmodel.list_s <- create_marginal_treatment_models(trtmodel.list)

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
  # str(extract_variables_list(list(a_1 ~ w1 + w2 + L1_1 + L2_1,
  #                            a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1)))

  trtmodel <- extract_variables_list(trtmodel.list)
  trtmodel_s <- extract_variables_list(trtmodel.list_s)

  # Function to generate and write the JAGS model
  write_jags_model <- function(trtmodel.list) {
    # Extract variable information for each formula
    var_info <- lapply(trtmodel.list, extract_variables_list)

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

      # Build the logistic model including all previous treatments
      if (v == 1) {
        model_string <- paste0(model_string, "p", response, "s[i] <- ilogit(bs", v, "0)\n")
      } else {
        bs_terms <- paste0("bs", v, "0")
        for (j in 1:(v-1)) {
          prev_response_s <- var_info[[j]]$response
          bs_terms <- paste(bs_terms, " + bs", v, j, "*", prev_response_s, "s[i]", sep = "")
        }
        model_string <- paste0(model_string, "p", response, "s[i] <- ilogit(", bs_terms, ")\n")
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
      model_string <- paste0(model_string, "bs", v, "0~dnorm(0,.01)\n")
      if (v > 1) {
        for (j in 1:(v-1)) {
          model_string <- paste0(model_string, "bs", v, j, "~dnorm(0,.01)\n")
        }
      }

      # b parameters
      for (p in 0:(num_preds - 1)) {
        model_string <- paste0(model_string, "b", v, p, "~dnorm(0,.01)\n")
      }
    }

    # Add the closing brace for the model block
    model_string <- paste(model_string, "}\n", sep="")

    # Write the finalized model string to a file
    cat(model_string, file = "treatment_model.txt")

    # Assume model_string is being built here...
    # Example initializations for the example; this part needs to be your actual model code generation
    all_parameters <- c()

    # Process each visit
    for (v in seq_along(trtmodel.list)) {
      var_info <- extract_variables_list(trtmodel.list[[v]])
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

  # Generate and write the JAGS model
  write_jags_model(trtmodel.list)

  # Prepare data and parameters for JAGS
  prepare_jags_data <- function(data, trtmodel.list) {
    # Extract variable information from formulas
    variable_info <- lapply(trtmodel.list, extract_variables_list)

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



  # Use this function to prepare JAGS data
  jags.data <- prepare_jags_data(data, trtmodel.list)
  jags.params <- write_jags_model(trtmodel.list)

  # Run JAGS model
  # set.seed(seed) # Is this needed?
  # jagsfit <- jags(data = jags.data,
  #                 parameters.to.save = jags.params,
  #                 n.iter = n.iter,
  #                 model.file = "treatment_model.txt",
  #                 n.chains = 1,
  #                 n.burnin = n.burnin,
  #                 n.thin = n.thin,
  #                 jags.seed = seed)
  #
  # # Extract MCMC output
  # out.mcmc <- as.mcmc(jagsfit)
  #
  # posterior <- as.matrix(ldply(out.mcmc, data.frame))
  #
  # parallel = TRUE or FALSE; n.chains = 1 by default
  # if FALSE, you run your old jagsfit
  # if parallel = TRUE, #user specification;
  # #check
  # if n.chains == 1 then return error message "Parallel MCMC require at least 2 chains, computing is running on 1 core per chain"
  # if n.chains > 1 then check total number of available cores on user computer,
  #
  # library(parallel)
  # available_core = detectCores(logical = FALSE)
  # if n.chains >= available_core return error message
  #
  # print(paste("Parallel MCMC require 1 core per chain", "you have", available_core, "We recommand using",available_core-2))
  #
  #
  # only run jags.parallel if nchains > 1 & nchains < available_core
  #
  # jagsfit <-jags.parallel(data=,
  #               parameters.to.save=,
  #               model.file = "xxx",
  #               n.chains = n.chains,
  #               n.iter = n.iter,
  #               n.burnin = n.burin,
  #               n.thin = n.thin,
  #               n.cluster= n.chains,
  #               jags.seed = 123)
  #
  # if you do parallel, as.mcmc(jagsfit) will be a list with posterior data per chain,
  # if you have 2 chains
  # the combined posterior data set is full_post <- rbind(out.mcmc[[1]], out.mcmc[[2]]) the dimension is row = n.posterior * 2 (number of chains), column = number of parameters;

  # set.seed(seed) # Not needed

  # Check if parallel computing is requested
  if (parallel == TRUE) {
    if (n.chains == 1) {
      stop("Parallel MCMC requires at least 2 chains. Computing is running on 1 core per chain.")
    }
    available_cores <- detectCores(logical = FALSE)
    if (n.chains >= available_cores) {
      stop(paste("Parallel MCMC requires 1 core per chain. You have", available_cores, "cores. We recommend using", available_cores - 2, "cores."))
    }
    # Run JAGS model in parallel
    library(doParallel)
    cl <- makeCluster(n.chains)
    registerDoParallel(cl)
    jags.model.wd <- paste(getwd(), '/treatment_model.txt',sep='')

    posterior <- foreach(i=1:n.chains, .packages=c('R2jags'),
                                   .combine='rbind') %dopar%{

                                     jagsfit <- jags(data = jags.data,
                                                     parameters.to.save = jags.params,
                                                     model.file = jags.model.wd,
                                                     n.chains = 1,
                                                     n.iter = n.iter,
                                                     n.burnin = n.burnin,
                                                     n.thin = n.thin,
                                                     jags.seed = seed+i)
                                     # Combine MCMC output from multiple chains
                                     out.mcmc <- as.mcmc(jagsfit)
                                     return(do.call(rbind, lapply(out.mcmc, as.matrix)))
                                   }

    stopCluster(cl)

  } else if (parallel == FALSE) {

    if (n.chains != 1) {
      stop("Non-parallel MCMC requires 1 chain.")
    }

    # Run JAGS model without parallel computing
    jagsfit <- jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = "treatment_model.txt",
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = seed)

    # Extract MCMC output
    out.mcmc <- as.mcmc(jagsfit)
    posterior <- as.matrix(out.mcmc[[1]])
  }

  diagnostics <- geweke.diag(out.mcmc)

  # Check diagnostics for convergence issues
  significant_indices <- which(abs(diagnostics[[1]]$z) > 1.96)
  if (length(significant_indices) > 0) {
    warning("Some parameters have not converged with Geweke index > 1.96. More iterations may be needed.")
  }

  # n_posterior = (n.iter - n.burnin)/n.thin
  # posterior <- as.matrix(out.mcmc[[1]])

  # pa_1[i] <- ilogit(b10 + b11*w1[i] + b12*w2[i] + b13*L1_1[i] + b14*L2_1[i]), so #covariates = 4

  # number of parameters for this model is 5 and design matrix is 1 variables
  # looping throught each treatment model 1 visit at a time;
  expit <- function(x){exp(x) / (1+exp(x))}

  # psc three dimention 1 dimention for specifying visit, 1 dimention for nposterior, 1 dimention for nobs;
  # psc means conditional model with confounders and previous treatment;
  # psm means marginal model with previous treatment only;

  # Idea: (using n.iter=250, n.burnin=50, n.thin=5, so #posteriors=(250-50)/5=40)
  # Step 1: Parameter matrix with dim(#posteriors, 1 + #covariates); in this case, dim(40, 5)
  # Step 2: Build design matrix from data based on model formulas with dim(1 + #covariates, #obs); in this case, dim(5, 1000)
  # Step 3: Calculate weight matrix as parameter matrix %*% design matrix with dim(#posteriors, #obs); in this case, dim(40, 1000)
  # Step 4: Calculate mean weights across all posteriors using colMean()
  #
  # for ( nvisit in 1:length(trtmodel)){
  #   #treatment probability by visit;
  #   # dim( expit(posterior[,1:(1+length(trtmodel[[nvisit]]$predictors))] %*% t(cbind(1,data[,c(trtmodel[[1]]$predictors)]))))
  #   # dim(n_posterior, dim(data)[1])
  #   psc[nvisit, , ] <- expit(posterior[,1:(1+length(trtmodel[[nvisit]]$predictors))] %*% t(cbind(1,data[,c(trtmodel[[1]]$predictors)])))
  #   psm[nvisit, , ] <- expit(posterior[,1:(1+length(trtmodel_s[[nvisit]]$predictors))] %*% t(cbind(1,data[,c(trtmodel_s[[1]]$predictors)])))
  # }
  # Xiao's note: the way to extract columns from "posterior", posterior[,1:(1+length(trtmodel[[nvisit]]$predictors))], is wrong.
  # This is because the b's and bs's in "posterior" are not in correct order, and correct columns need to match the predictors in the treatment models.

  # [a,b,c] [d,e,f]
  # a*d + b*e + c*f
  #
  # [a*d, b*e, c*f]

  # weight calculation;
  # numerator: psm[1,,]*psm[2,,]*...*psm[nvisit,,]
  # denominator: psc[1,,]*psc[2,,]*...*psc[nvisit,,]
  # wmean <- colMean(numertor/denominator)

  n_visits <- length(trtmodel)
  n_posterior <- dim(posterior)[1]
  n_obs <- nrow(data)

  # Initialize arrays for storing probabilities
  psc <- array(dim = c(n_visits, n_posterior, n_obs))
  psm <- array(dim = c(n_visits, n_posterior, n_obs))

  # # Calculate probabilities for each visit
  # for (nvisit in 1:n_visits) {
  #   # Extract predictors for the conditional model
  #   predictors_c <- trtmodel[[nvisit]]$predictors
  #   predictors_s <- trtmodel_s[[nvisit]]$predictors
  #
  #   # Design matrices for conditional and marginal models
  #   design_matrix_c <- cbind(1, data[, predictors_c, drop = FALSE])
  #   design_matrix_s <- cbind(1, data[, predictors_s, drop = FALSE])
  #
  #   for (j in 1:n_posterior) {
  #     psc[nvisit, j, ] <- expit(posterior[j, 1:(1+length(predictors_c))] %*% t(design_matrix_c))
  #     psm[nvisit, j, ] <- expit(posterior[j, 1:(1+length(predictors_s))] %*% t(design_matrix_s))
  #   }
  # }
  # Xiao's note: the way to extract columns from "posterior", posterior[j, 1:(1+length(predictors_c))], is wrong.
  # This is because the b's and bs's in "posterior" are not in correct order, and correct columns need to match the predictors in the treatment models.

  # Calculate probabilities for each visit
  parameter_map <- colnames(posterior)

  # Complete data: length N3
  # Subset data: non-missing

  for (nvisit in 1:n_visits) {
    predictors_c <- trtmodel[[nvisit]]$predictors
    predictors_s <- trtmodel_s[[nvisit]]$predictors

    design_matrix_c <- cbind(1, data[, predictors_c, drop = FALSE]) # Code with complete data "data"
    design_matrix_s <- cbind(1, data[, predictors_s, drop = FALSE])

    beta_indices_c <- match(c(sprintf("b%d0", nvisit),
                              sapply(1:length(predictors_c), function(p) sprintf("b%d%d", nvisit, p))),
                            parameter_map) # This produces the correct indices that corresponds to the predictors in the treatment models

    beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
                              if (nvisit > 1) sprintf("bs%d1", nvisit) else NULL),
                            parameter_map)

    for (j in 1:n_posterior) {
      psc[nvisit, j, ] <- expit(posterior[j, beta_indices_c] %*% t(design_matrix_c))
      psm[nvisit, j, ] <- expit(posterior[j, beta_indices_s] %*% t(design_matrix_s))
      # For censoring: need to have two more lines: csc and csm (censoring score marginal) for stablizers
    }
  }

  # Calculate the product of probabilities across visits for each posterior and observation
  numerator <- apply(psm, c(2, 3), prod)  # Apply 'prod' across the first dimension (visits)
  denominator <- apply(psc, c(2, 3), prod)  # Same for psc
  # For censoring: two more numerator and denominator

  # Calculate weights by element-wise division
  weights <- numerator / denominator  # Resulting in a 40 (posterior samples) x 1000 (observations) matrix
  # For censoring: weights <- (numerator_trt * numerator_cen) / (denominator_trt * denominator_cen)

  # Mean weight across all observations
  wmean <- colMeans(weights)

  return(wmean)
}



# Testing
# Continuous outcome with a_3 and a_4
testdata <- read.csv("continuous_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
testdata$a_3 <- rbinom(n=length(testdata$y),1,p=0.4)
testdata$a_4 <- rbinom(n=length(testdata$y),1,p=0.6)
testdata$L1_3 <- rbinom(n=length(testdata$y),1,p=0.35)
testdata$L2_3 <- rbinom(n=length(testdata$y),1,p=0.65)
testdata$L1_4 <- rbinom(n=length(testdata$y),1,p=0.45)
testdata$L2_4 <- rbinom(n=length(testdata$y),1,p=0.55)
start<-Sys.time()
weights1 <- bayesweight(trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1,
                                             a_3 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2,
                                             a_4 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1 + L1_3 + L2_3 + a_2 + L1_4 + L2_4 + a_3),
                       data = testdata,
                       n.iter = 250,
                       n.burnin = 50,
                       n.thin = 5)
Sys.time()-start
weights1



# Continuous outcome; original dataset (same as that on Quarto)
testdata2 <- read.csv("continuous_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
start<-Sys.time()
weights2 <- bayesweight(trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = testdata2,
                        n.iter = 25000,
                        n.burnin = 15000,
                        n.thin = 5,
                        n.chains = 2,
                        seed = 890123,
                        parallel = TRUE)
Sys.time()-start
weights2

mean(weights2)
var(weights2)

# As a comparison:
causaldata_with_bayeswt <- read.csv("R/continuous_outcome_data2.csv", header = TRUE, fileEncoding="UTF-8-BOM")
bayeswt <- causaldata_with_bayeswt$bayeswt
mean(bayeswt)
var(bayeswt)

# Binary outcome; original dataset
testdata3 <- read.csv("binary_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")
start<-Sys.time()
weights3 <- bayesweight(trtmodel.list = list(a_1 ~ w1 + w2 + L1_1 + L2_1,
                                             a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
                        data = testdata3,
                        n.iter = 25000,
                        n.burnin = 15000,
                        n.thin = 5,
                        n.chains = 2,
                        seed = 890123,
                        parallel = TRUE)
Sys.time()-start
weights3
