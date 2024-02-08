#' Bayesian Marginal Structural Model Bootstrap Estimation
#'
#' This function performs Bayesian non-parametric bootstrap to calculate causal effects
#' for Bayesian marginal structural models.
#'
#' @param wloglik_function The weighted log-likelihood function to optimize.
#' @param testdata A data frame containing the outcome variable 'y' and treatment variables 'a_1', 'a_2', etc.
#' @param wmean A vector of weights to be used in the bootstrap procedure.
#' @param nboot The number of bootstrap iterations (default is 1000).
#' @param inits1 Initial values for the optimization algorithm (default is c(rep(0.1, length(treatment_vars)+1),4).
#' @param optim_method The optimization method to be used in the optim() function (default is 'BFGS').
#'
#' @return A list of summary measures (mean, sd, quantile)
#' @export
#'
bayesm_bootstrap <- function(outcomevar = Y, # Do we really need this? outcomevar is variables$response extracted from data
                             ymodel = "myoutcome ~ thing1 + a1^2 + thing3",
                             intervention = list(ref_int = c(rep(0,n_visits)), comparator = c(rep(1,n_visits))), # just an example for 2 vectors of treatments
                             nvisit = 10,
                             ref_int = c(rep(0,n_visits)), # An example of never treated
                             family = "gaussian", # "gaussian" or "binomial"
                             data = testdata,
                             wmean = rep(1, 1000),
                             nboot = 1000,
                             optim_method = 'BFGS') {

  # check
  # return error message if the input weight vector has different length comparing to the outcome Y;
  if (length(wmean) != nrow(data)) {
    stop("The length of the weight vector does not match the length of Y.")
  }

  extract_variables <- function(formula) {
    # Get the terms of the formula
    formula_terms <- terms(formula)

    # Extract the response variable name (if there is one)
    response_variable <- attr(formula_terms, "response")
    response_name <- if (response_variable > 0) {
      all_vars <- all.vars(formula)
      all_vars[response_variable]
    } else {
      NA
    }

    # Extract predictor variable names
    predictor_names <- attr(formula_terms, "term.labels")

    # Return a list of response and predictor variables
    list(response = response_name, predictors = predictor_names)
  }

  variables <- extract_variables(ymodel) # Extract variable names from the formula
  Y_name <- variables$response
  A_name <- variables$predictors

  Y <- data[[Y_name]]
  A <- cbind(1, data[, A_name]) # How to deal with interaction terms like thing1:thing2?



  wloglik_normal<-function(param,
                           Y,
                           A,
                           weight){
    #number of observations;
    n <- length(Y)
    theta <- param[1:dim(A)[2]] #causal parameters on the mean
    #number of parameter is determined by number of treatment variables, plus intercept;
    sigma <- param[(dim(A)[2]+1)] # the remaining the parameter represent the standard deviation;
    mmat <- as.matrix(A) #design matrix of the causal outcome model, e.g., A = cbind(1, a_1, a_2);
    logl<- -0.5*log(sigma^2) - 0.5*((Y - mmat%*%theta)^2)/(sigma^2)
    wlogl<-sum(weight*logl)

    return(wlogl)
  }



  if (family == "gaussian"){
    wfn = wloglik_normal
  } else if (family == "binomial"){
    wfn = wloglik_binomial
  } else if (!family %in% c("gaussian","binomial")){
    stop("Current version only handles continuous (gaussian) and binary (binomial) outcomes.")
  }

  bootest <- numeric(nboot)

  for (j in 1:nboot) {
    alpha <- as.numeric(rdirichlet(1, rep(1.0, length(Y))))

    inits1 <- c(rep(0.1, length(A_name)), 4)  # Default initial values

    maxim <- optim(inits1,
                   fn = wfn,
                   Y = Y,
                   A = A,
                   weight = alpha * wmean,
                   control = list(fnscale = -1),
                   method = optim_method,
                   hessian = FALSE)

    effect_ref_int <- sum(maxim$par[intervention$ref_int]) # Example: ref_int = c(0,0), never treated for 2 visits
    effect_comparator <- sum(maxim$par[intervention$comparator]) # Example: comparator = c(1,1), always treated for 2 visits
    bootest[j] <- effect_comparator - effect_ref_int # Not sure if this is correct; for example always treated vs never treated
    # Sum the effects for all treatment variables

    if (j %% 100 == 0) {
      print(j)
    }
  }

  return(list(
    mean = mean(bootest),
    sd = sqrt(var(bootest)),
    quantile = quantile(bootest, probs = c(0.025, 0.975))
  ))
}



