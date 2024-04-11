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

bayesm_bootstrap <- function(ymodel = y ~ a_1*a_2*a_3*a_4,
                             nvisit = 10,
                             reference = c(rep(0,nvisits)), # An example of never treated
                             comparator = c(rep(1,nvisits)),
                             family = "gaussian", # "gaussian" or "binomial"
                             data = testdata,
                             wmean = rep(1, 1000),
                             nboot = 1000,
                             optim_method = 'BFGS',
                             estimand = "RD",
                             parallel = TRUE,
                             ncore = 4){

  # #testing;
  # ymodel = y ~ a_1*a_2;
  # nvisit = 2;
  # reference = c(rep(0,2));
  # comparator = c(rep(1,2));
  # family = "binomial";
  # data = testdata;
  # wmean = rep(1, 1000);
  # nboot = 1000;
  # optim_method = "BFGS";

  # first thing first load all the required R packages;
  if (!require(foreach)){
    install.packages("foreach",repos="http://cran.r-project.org")
    library(foreach)
  }
  if (!require(doParallel)){
    install.packages("doParallel",repos="http://cran.r-project.org")
    library(doParallel)
  }
  if (!require(MCMCpack)){
    install.packages("MCMCpack",repos="http://cran.r-project.org")
    library(MCMCpack)
  }

  # check
  # return error message if the input weight vector has different length comparing to the outcome Y;
  if (length(wmean) != nrow(data)) {
    stop("The length of the weight vector does not match the length of Y.")
  }

  # # return error message if reference or intervention n_visit differenct from nvisit;
  # if (xxx) {
  #   stop("The potential outcome treatment sequence does not match with nvisit in length")
  # }

  # loading utility functions;
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

  source("R/rdirichlet.R")
  source("R/calculate_effect.R")

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

  wloglik_binomial <- function(param,
                               Y,
                               A,
                               weight){
    # number of observations;
    n <- length(Y)
    beta <- param[1:dim(A)[2]] # causal parameters on the log-odds scale (no sigma for binomial?)
    mmat <- as.matrix(A)
    eta<-mmat %*% beta # linear predictor
    # p <- 1 / (1 + exp(-eta))
    logl <- Y*eta - log(1+exp(eta))
      # Y * log(p + 0.0001) + (1 - Y) * log(1 - p + 0.0001)
    wlogl<-sum(weight*logl)
    return(wlogl)
  }

  expit <- function(x){exp(x) / (1+exp(x))}

  if (family == "gaussian"){
    wfn = wloglik_normal
    inits1 <- c(rep(0.1, length(A)), 4)  # Default initial values, 4 is for the SD;
  } else if (family == "binomial"){
    wfn = wloglik_binomial
    inits1 <- c(rep(1, length(A)))
  } else if (!family %in% c("gaussian","binomial")){
    stop("Current version only handles continuous (gaussian) and binary (binomial) outcomes.")
  }


  # parallel computing only for this bootstrap step;
  if (parallel == TRUE){
  numCores <- ncore
  registerDoParallel(cores = numCores)

  results <- foreach(i=1:nboot, .combine = 'rbind') %dopar% {

    results.it <- matrix(NA, 1, 3) #result matrix, three columns for bootest, effect_ref, and effect_comp;

    alpha <- as.numeric(rdirichlet(1, rep(1.0, length(Y))))

    maxim <- optim(inits1,
                   fn = wfn,
                   Y = Y,
                   A = A,
                   weight = alpha * wmean,
                   control = list(fnscale = -1),
                   method = optim_method,
                   hessian = FALSE)

    names(maxim$par) <- c("(Intercept)", variables$predictors)

    # Calculate the effects
    results.it[1,1] <- calculate_effect(reference, variables, param_estimates=maxim$par)
    results.it[1,2] <- calculate_effect(comparator, variables, param_estimates=maxim$par)

    # Calculate the ATE
    # results.it[1,3] <- results.it[1,1] - results.it[1,2]
    if (family == "binomial") { # binary outcomes
      if (estimand == "RD") { # Risk Difference
        bootest[j] <- expit(results.it[1,1]) - expit(results.it[1,2])
      } else if (estimand == "RR") {
        bootest[j] <- expit(results.it[1,1]) / expit(results.it[1,2])
      } else if (estimand == "OR") {
        bootest[j] <- (expit(results.it[1,1]) / (1 - expit(results.it[1,1]))) /
                      (expit(results.it[1,2]) / (1 - expit(results.it[1,2])))
      }

    } else if (family == "gaussian"){ # continuous outcomes
      if (estimand == "RD") { # Risk Difference
        results.it[1,3] <- results.it[1,1] - results.it[1,2]
      } else if (estimand == "RR") {
        results.it[1,3] <- results.it[1,1] / results.it[1,2]
      } else if (estimand == "OR") {
        bootest[j] <- (results.it[1,1] / (1 - results.it[1,1])) /
                      (results.it[1,2] / (1 - results.it[1,2]))
      }

    # combining parallel results;
    cbind(i,results.it) #end of parallel;
  }

  #saving output for the parallel setting;
  return(list(
    mean = mean(results[,4]),
    sd = sqrt(var(results[,4])),
    quantile = quantile(results[,4], probs = c(0.025, 0.975)),
    bootdata <- data.frame(results[,-1]),
    reference = reference,
    comparator = comparator
  ))

  }

  else if (parallel == FALSE) {

    bootest <- numeric(nboot)
    effect_reference <- numeric(nboot)
    effect_comparator <- numeric(nboot)

    for (j in 1:nboot) {
      alpha <- as.numeric(rdirichlet(1, rep(1.0, length(Y))))

      maxim <- optim(inits1,
                     fn = wfn,
                     Y = Y,
                     A = A,
                     weight = alpha * wmean,
                     control = list(fnscale = -1),
                     method = optim_method,
                     hessian = FALSE)

      # maxim <- optimx(inits1,
      #                fn = wfn,
      #                Y = Y,
      #                A = A,
      #                weight = alpha * wmean,
      #                control = list(fnscale = -1),
      #                method = optim_method)

      # maxim <- nlminb(start = inits1,
      #                 objective = wfn,
      #                 Y = Y,
      #                 A = A,
      #                 weight = alpha * wmean)

      # maxim <- nloptr(x0 = inits1,
      #                 eval_f = function(param) -wfn(param, Y = Y, A = A, weight = alpha * wmean),  # Note the negation for minimization
      #                 opts = list(algorithm = "NLOPT_LN_NELDERMEAD",
      #                             xtol_rel = 1.0e-8,
      #                             maxeval = 1500))

      # optim()
      names(maxim$par) <- c("(Intercept)", variables$predictors)

      # Calculate the effects
      effect_reference[j] <- calculate_effect(reference, variables, param_estimates=maxim$par)
      effect_comparator[j] <- calculate_effect(comparator, variables, param_estimates=maxim$par)

      # #optimx()
      # param_vector = unlist(lapply(1:length(inits1), function(i) maxim[[paste0("p", i)]]))
      # names(param_vector) <- c("(Intercept)", variables$predictors)
      #
      # # Calculate the effects
      # effect_reference[j] <- calculate_effect(reference, variables, param_estimates=param_vector)
      # effect_comparator[j] <- calculate_effect(comparator, variables, param_estimates=param_vector)

      #nloptr()
      # names(maxim$solution) <- c("(Intercept)", variables$predictors)
      #
      # # Calculate the effects
      # effect_reference[j] <- calculate_effect(reference, variables, param_estimates=maxim$solution)
      # effect_comparator[j] <- calculate_effect(comparator, variables, param_estimates=maxim$solution)

      # Calculate the ATE
      if (family == "binomial") { # binary outcomes
        if (estimand == "RD") { # Risk Difference
          bootest[j] <- expit(effect_comparator[j]) - expit(effect_reference[j])
        } else if (estimand == "RR") {
          bootest[j] <- expit(effect_comparator[j]) / expit(effect_reference[j])
        } else if (estimand == "OR") {
          bootest[j] <- (expit(effect_comparator[j]) / (1 - expit(effect_comparator[j]))) /
                        (expit(effect_reference[j]) / (1 - expit(effect_reference[j])))
        }

      } else if (family == "gaussian"){ # continuous outcomes
        if (estimand == "RD") { # Risk Difference
          bootest[j] <- effect_comparator[j] - effect_reference[j]
        } else if (estimand == "RR") {
          bootest[j] <- effect_comparator[j] / effect_reference[j]
        } else if (estimand == "OR") {
          bootest[j] <- (effect_comparator[j] / (1 - effect_comparator[j])) /
                        (effect_reference[j] / (1 - effect_reference[j]))
        }

      }


    }
    #saving output for the non-parallel setting;
    return(list(
      mean = mean(bootest),
      sd = sqrt(var(bootest)),
      quantile = quantile(bootest, probs = c(0.025, 0.975)),
      bootdata = data.frame(effect_reference, effect_comparator, ATE=bootest),
      reference = reference,
      comparator = comparator
    ))

  }
}


#test;
# Continuous outcome
# setwd("C:/Users/YanXi/Downloads")
testdata <- readr::read_csv("R/continuous_outcome_data.csv")
testdata$a_3 <- rbinom(n=length(testdata$y),1,p=0.4)
testdata$a_4 <- rbinom(n=length(testdata$y),1,p=0.6)
start<-Sys.time()
model1 <- bayesm_bootstrap(ymodel = y ~ a_1+a_2+a_3+a_4,
                 nvisit = 4,
                 reference = c(rep(0,4)),
                 comparator = c(rep(1,4)),
                 family = "gaussian",
                 data = testdata,
                 wmean = rep(1, 1000),
                 nboot = 1000,
                 optim_method = "BFGS",
                 estimand = "RD",
                 parallel = FALSE,
                 ncore = 10)
Sys.time()-start
bootoutput = model1$bootdata

mean(bootoutput$ATE)
var(bootoutput$ATE)
sqrt(var(bootoutput$ATE))
quantile(bootoutput$ATE, probs=c(0.025,0.975))

# Continuous outcome; original dataset (same as that on Quarto)
# setwd("C:/Users/YanXi/Downloads")
testdata2 <- readr::read_csv("R/continuous_outcome_data.csv")
start<-Sys.time()
model2 <- bayesm_bootstrap(ymodel = y ~ a_1+a_2,
                           nvisit = 2,
                           reference = c(rep(0,2)),
                           comparator = c(rep(1,2)),
                           family = "gaussian",
                           data = testdata2,
                           wmean = rep(1, 1000),
                           nboot = 1000,
                           optim_method = "BFGS",
                           estimand = "RD",
                           parallel = FALSE,
                           ncore = 6)
Sys.time()-start
bootoutput2 = model2$bootdata

mean(bootoutput2$ATE)
var(bootoutput2$ATE)
sqrt(var(bootoutput2$ATE))
quantile(bootoutput2$ATE, probs=c(0.025,0.975))
# Result
# > mean(bootoutput2$ATE)
# [1] -3.164354
# > var(bootoutput2$ATE)
# [1] 0.009445813
# > sqrt(var(bootoutput2$ATE))
# [1] 0.09718957
# > quantile(bootoutput2$ATE, probs=c(0.025,0.975))
# 2.5%     97.5%
#   -3.356886 -2.992084

# Binary outcome; original dataset
# setwd("C:/Users/YanXi/Downloads")
testdata3 <- readr::read_csv("R/binary_outcome_data.csv")
start<-Sys.time()
model3 <- bayesm_bootstrap(ymodel = y ~ a_1+a_2,
                           nvisit = 2,
                           reference = c(rep(0,2)),
                           comparator = c(rep(1,2)),
                           family = "binomial",
                           data = testdata3,
                           wmean = rep(1, 1000),
                           nboot = 1000,
                           optim_method = "BFGS",
                           estimand = "RD",
                           parallel = FALSE,
                           ncore = 6)
Sys.time()-start
bootoutput3 = model3$bootdata

mean(bootoutput3$ATE)
var(bootoutput3$ATE)
sqrt(var(bootoutput3$ATE))
quantile(bootoutput3$ATE, probs=c(0.025,0.975))
# Result
# > mean(bootoutput3$ATE)
# [1] -0.1044813
# > var(bootoutput3$ATE)
# [1] 0.001627635
# > sqrt(var(bootoutput3$ATE))
# [1] 0.04034396
# > quantile(bootoutput3$ATE, probs=c(0.025,0.975))
# 2.5%       97.5%
#   -0.18221276 -0.02746802

