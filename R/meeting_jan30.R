#' @param wmean A vector of weights to be used in the bootstrap procedure.
#' @param nboot The number of bootstrap iterations (default is 1000).
#' @param inits1 Initial values for the optimization algorithm (default is c(rep(0.1, length(treatment_vars)+1),4).
#' @param optim_control List of control parameters to pass to the optim() function (default is list(fnscale = -1)).
#' @param optim_method The optimization method to be used in the optim() function (default is 'BFGS').
#' @param optim_hessian Logical; argument of optim() function (default is FALSE).
#'
#' @return A list of summary measures (mean, sd, quantile)'
#' @param param Numeric vector containing parameters (coefficients for treatment and standard deviation).
#' @param Y Numeric vector of dependent variable (outcome).
#' @param A Matrix of independent variables (treatment variables).
#' @param weight Numeric vector of treatment weights for each observation.
#' @return Numeric value of the weighted log-likelihood.
#' @export
#'
bayesm_bootstrap <- function(
  outcomevar = Y,
  ymodel = "myoutcome ~ thing1 + a1^2 + thing3",
  intervention = (2 vector of treatments),
  nvisit = 10,
  ref_int = ,
  family = , data = mydate,
  wmean = rep(1, 1000),
  nboot = 1000,
  optim_method = 'BFGS') {
  # check
  # return error message if the input weight vector has different length comparing to the outcome Y;

  #define user specificed outcome model variables;
  Y_name = extracting outcome variable name from the formula, trick is to use string before "~"
  A_name = c(1, List of covariates), this can be a vector = extracting strings after "~" and seperate by "+"


  # if people specify initial value, return error message if the initial value length is different from the total number of parameters;
  # also return error message for normal case, if initial value of standiviation is negative.
  # in your help file, you can tell user that the first few parameters for the normal case is mean and the last one is variance;

  next step is using these names to map it to the data

  Y <- from data get outcome variable with name Y_name
  A <- from data get design matrix

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
  } else if (family !%in% c("gaussian","binomial")){
    return error message saying "current version only handle continous and binary outcome"

  }

  bootest <- numeric(nboot)
  treatment_vars <- names(testdata)[grepl("^a_", names(testdata))]



  for (j in 1:nboot) {
    alpha <- as.numeric(rdirichlet(1, rep(1.0, length(testdata$y))))

    # we fix initial value for user and add this info to helpful incase people are curious;
    # if people didn't specificy initial value, they you use defult,
    # normal: which is 0.1 for all the mean and 1 for variance this for normal,
    # binomial: if we have binary all parameters are log odds which can be negative, we can still use 0.1 as logdds


    maxim <- optim(inits1,
                   fn = wfn,
                   Y = Y,
                   A = A,
                   weight = alpha * wmean,
                   control = list(fnscale = -1),
                   method = optim_method,
                   hessian = FALSE)
    # this line is returning causal estimand we want it to reflect user specification;
    # we let people to specify 2 intervantion sequence for comparision;
    # intervention = list (first thing, second thing)
    # we ask people to indicate the reference intervention;
    # for example ref_int is = c(rep(0,nvisit))
    # a_j = 0 or 1, for j=1, ..., nvisit based on user specification
    # Y = vector of theta*vector of a ordered by visits
    # effectunderreferenceintervention = sum of subset of mean parameters
    # effectundercomparater = sum of another subset of mean parameters
    # bootest[j] should be the difference between effectundercomparater-effectunderreferenceintervention

    bootest[j] <- sum(maxim$par[2:length(treatment_vars)+1]) # Sum the effects for all treatment variables

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



