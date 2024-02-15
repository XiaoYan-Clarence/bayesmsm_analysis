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
bayesm_bootstrap <- function(ymodel = "myoutcome ~ thing1 + a1^2 + thing3",
                             nvisit = 10,
                             intervention = list(ref_int = c(rep(0,n_visits)),
                                                 comparator = c(rep(1,n_visits))), # just an example for 2 vectors of treatments
                             ref_int = c(rep(0,n_visits)), # An example of never treated
                             family = "gaussian", # "gaussian" or "binomial"
                             data = testdata,
                             wmean = rep(1, 1000),
                             nboot = 1000,
                             optim_method = 'BFGS'){

  #testing;
  ymodel = y ~ a_1+a_2+a_1*a_2*a_3 + a_2*a_3;
  nvisit = 3;
  intervention = list(ref_int = c(rep(0,nvisit)),
                      comparator = c(rep(1,nvisit))); # just an example for 2 vectors of treatments
  ref_int = c(rep(0,nvisit)); # An example of never treated
  family = "gaussian";
  data = testdata;
  wmean = rep(1, 1000);
  nboot = 1000;
  optim_method = 'BFGS';

  # first thing first load all the required R packages;
  library(MCMCpack)


  # check
  # return error message if the input weight vector has different length comparing to the outcome Y;
  if (length(wmean) != nrow(data)) {
    stop("The length of the weight vector does not match the length of Y.")
  }

  # # return error message if ref_int or intervention n_visit differenct from nvisit;
  # if (xxx) {
  #   stop("The potential outcome treatment sequence does not match with nvisit in length")
  # }

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

  Y <- data[[Y_name]]
  A_base <- data.frame(matrix(data = NA,
                   nrow = nrow(data),
                   ncol = length(variables$predictors)))
  for (i in 1:length(variables$predictors)){
    initial_vector <- variables$predictors[i]
    split_vector <- strsplit(initial_vector, ":")
    new_vector <- unlist(split_vector)
    if (length(new_vector)==1){
      A_base[,i] <-  data[, new_vector]
    } else if (length(new_vector)>1){
      A_base[,i] <-  apply(data[, new_vector],1,prod)
    }
  }

  A <- cbind(1, A_base)
  colnames(A)[2:ncol(A)]<- variables$predictors

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

    j = 1;
    alpha <- as.numeric(rdirichlet(1, rep(1.0, length(Y))))

    inits1 <- c(rep(0.1, length(A)), 4)  # Default initial values

    maxim <- optim(inits1,
                   fn = wfn,
                   Y = Y,
                   A = A,
                   weight = alpha * wmean,
                   control = list(fnscale = -1),
                   method = optim_method,
                   hessian = FALSE)
    #Feb 8th, 2024;
    maxim$par
    [1]  2.40706022 -1.33106703 -2.31545410 -0.09949307  0.58522148  0.32738170  0.14387898 -0.22532452
    [9] -1.14763173
    y= intercept+ "a_1"+"a_2"+"a_3" + "a_1:a_2"+"a_1:a_3"+"a_2:a_3"+"a_1:a_2:a_3"
    ref_int <- c(0,1,1)
    # y given a_1 = 0, a_2=1, a_3=0 what is that value given the parameter estimates =


    # Function to calculate the effect of an intervention given the parameter estimates and intervention levels
    calculate_effect <- function(intervention_levels, variables, param_estimates) {
      # Start with the intercept term
      effect <- param_estimates[1]

      # Go through each predictor and add its contribution
      for (i in 1:length(variables$predictors)) {
        term <- variables$predictors[i]
        term_variables <- unlist(strsplit(term, ":"))
        term_index <- which(names(param_estimates) == term)

        # Calculate the product of intervention levels for the interaction term
        term_contribution <- param_estimates[term_index]
        for (term_variable in term_variables) {
          var_index <- which(variables$predictors == term_variable)
          term_contribution <- term_contribution * intervention_levels[var_index]
        }

        # Add the term contribution to the effect
        effect <- effect + term_contribution
      }

      return(effect)
    }

    # Testing with our previous results
    maxim$par <- c(2.40706022, -1.33106703, -2.31545410, -0.09949307, 0.58522148, 0.32738170, 0.14387898, -0.22532452, -1.14763173)
    # names(maxim$par) <- c("(Intercept)", "a_1", "a_2", "a_3", "a_1:a_2", "a_1:a_3", "a_2:a_3", "a_1:a_2:a_3")
    names(maxim$par) <- c("(Intercept)", variables$predictors)

    # Treatment history
    ref_int <- c(0, 1, 1)  # Example: a_1 = 0, a_2 = 1, a_3 = 1
    comparator <- c(1, 1, 1)  # Example: a_1 = 1, a_2 = 1, a_3 = 1

    # Calculate the effects
    effect_ref_int <- calculate_effect(ref_int, variables, param_estimates=maxim$par)
    effect_comparator <- calculate_effect(comparator, variables, param_estimates=maxim$par)

    # effect_ref_int
    # (Intercept)
    # 0.135992
    # which is exactly 2.40706022-2.31545410-0.09949307+0.14387898!!!

    # Calculate the ATE
    bootest <- effect_comparator - effect_ref_int




      # adding empty initation similar to bootest outside the loop;
    # effect_ref_int[j] <- sum(maxim$par[intervention$ref_int]) # Example: ref_int = c(0,0), never treated for 2 visits
    # effect_comparator[j] <- sum(maxim$par[intervention$comparator]) # Example: comparator = c(1,1), always treated for 2 visits
    # bootest[j] <- effect_comparator - effect_ref_int # Not sure if this is correct; for example always treated vs never treated
    # Sum the effects for all treatment variables

    if (j %% 100 == 0) {
      print(j)
    }
  }

  return(list(
    mean = mean(bootest),
    sd = sqrt(var(bootest)),
    quantile = quantile(bootest, probs = c(0.025, 0.975)),
    data.frame(effect_ref_int, effect_comparator, bootest)
  ))

  #also save output dataframe this will have 3 columns each corrsponding to effect_ref_int,effect_comparator, difference;
}


#test;
testdata <- readr::read_csv("R/continuous_outcome_data.csv")
testdata$a_3 <- rbinom(n=length(testdata$y),1,p=0.4)
bayesm_bootstrap <- function(ymodel = "y ~ a_1+a_2+a_1*a_2*a_3 + a_2*a_3",
                             nvisit = 10,
                             intervention = list(ref_int = c(rep(0,n_visits)),
                                                 comparator = c(rep(1,n_visits))), # just an example for 2 vectors of treatments
                             ref_int = c(rep(0,n_visits)), # An example of never treated
                             family = "gaussian", # "gaussian" or "binomial"
                             data = testdata,
                             wmean = rep(1, 1000),
                             nboot = 1000,
                             optim_method = 'BFGS')


