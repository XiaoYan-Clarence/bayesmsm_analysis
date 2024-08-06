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
bayesweight_cen <- function(trtmodel.list = list(A1 ~ L11 + L21,
                                                 A2 ~ L11 + L21 + L12 + L22 + A1,
                                                 A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
                            cenmodel.list = list(C1 ~ L11 + L21,
                                                 C2 ~ L11 + L21 + A1,
                                                 C3 ~ L11 + L21 + A1 + L12 + L22 + A2),
                            data,
                            n.iter = 25000,
                            n.burnin = 15000,
                            n.thin = 5,
                            parallel = FALSE,
                            n.chains = 1,
                            seed = 890123) {


  # Testing
  # simulating causal data; (see sim_causal.R)
  # with censoring;
  sim.rc <- function(samplesize = 500)
  {
    set.seed(123)
    expit <- function(x){exp(x)/(1+exp(x))}
    # visit 1;
    L11 <- rbinom(n=samplesize, size=1, prob=0.5)
    L21 <- rnorm(n=samplesize, mean=0, sd=1)

    C1prob <- expit(-2-0.1*L11-0.1*L21) #right-censoring also known as lost to followup require non missing baseline. If there is missing baseline, ask user to impute missing baseline variables or remove missing observations;
    C1 <- rbinom(n=samplesize, size=1, prob=C1prob)

    A1prob <- expit(0.5*L11-0.2*L21)
    A1 <- rbinom(n=samplesize, size=1, prob=A1prob)

    # visit 2;
    C2prob <- expit(-2-0.1*L11-0.1*L21 - 0.1*A1)
    C2 <- rbinom(n=samplesize, size=1, prob=C2prob)

    L12prob <- expit(0.5*A1+0.5*L11-0.2*L21)
    L12 <- rbinom(n=samplesize, size=1, prob=L12prob)
    meanL22 <- 0.5*L21-0.5*A1-0.2*L11
    L22 <- rnorm(n=samplesize, mean=meanL22, sd=1)

    A2prob <- expit(0.5*L12-0.2*L22+0.2*A1)
    A2 <- rbinom(n = samplesize, size = 1, prob = A2prob)

    # visit 3;
    C3prob <- expit(-2-0.1*L12-0.1*L22 - 0.1*A2)
    C3 <- rbinom(n=samplesize, size=1, prob=C3prob)

    L13prob <- expit(0.5*A2+0.5*L12-0.2*L22)
    L13 <- rbinom(n = samplesize, size = 1, prob = L13prob)
    meanL23 <- 0.5*L22-0.5*A2-0.2*L12
    L23 <- rnorm(n=samplesize, mean=meanL23, sd=1)

    A3prob <- expit(0.5*L13-0.2*L23+0.2*A2)
    A3 <- rbinom(n = samplesize, size = 1, prob = A3prob)

    # end-of-study outcome;
    Yprob <- expit(0.3*A3+0.1*A2-0.1*A1+0.1*L13-0.2*L23)
    Y <- rbinom(n = samplesize, size = 1, prob = Yprob)
    dat <- cbind(L11, L21, A1, L12, L22, A2, L13, L23, A3, C1, C2, C3, Y)
    dat <- data.frame(dat)
    return(dat)
  }
  simdat <- sim.rc(samplesize = 500)

  library(tidyverse)
  simdat_cen <- simdat %>%
    mutate(A1 = ifelse(C1==1, NA, A1),
           C2 = ifelse(C1==1, NA, C2),
           L12 = ifelse(C1==1, NA, L12),
           L22 = ifelse(C1==1, NA, L22),
           A2 = ifelse(C1==1, NA, A2),
           L13 = ifelse(C1==1, NA, L13),
           L23 = ifelse(C1==1, NA, L23),
           A3 = ifelse(C1==1, NA, A3),
           C3 = ifelse(C1==1, NA, C3),
           Y = ifelse(C1==1, NA, Y)) %>%
    mutate(L12 = ifelse(C2==1, NA, L12),
           L22 = ifelse(C2==1, NA, L22),
           A2 = ifelse(C2==1, NA, A2),
           L13 = ifelse(C2==1, NA, L13),
           L23 = ifelse(C2==1, NA, L23),
           A3 = ifelse(C2==1, NA, A3),
           C3 = ifelse(C2==1, NA, C3),
           Y = ifelse(C2==1, NA, Y)) %>%
    mutate(L13 = ifelse(C3==1, NA, L13),
           L23 = ifelse(C3==1, NA, L23),
           A3 = ifelse(C3==1, NA, A3),
           Y = ifelse(C3==1, NA, Y))

  # Testing
  trtmodel.list = list(A1 ~ L11 + L21,
                       A2 ~ L11 + L21 + L12 + L22 + A1,
                       A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2)
  cenmodel.list = list(C1 ~ L11 + L21,
                       C2 ~ L11 + L21 + A1,
                       C3 ~ L11 + L21 + A1 + L12 + L22 + A2)
  data = simdat_cen
  n.iter = 25000
  n.burnin = 15000
  n.thin = 5
  parallel = TRUE
  n.chains = 2
  seed = 890123



  # Load all the required R packages;
  if (!require(R2jags)){
    install.packages("R2jags",repos="http://cran.r-project.org")
    library(R2jags)
  }
  if (!require(parallel)){
    install.packages("parallel",repos="http://cran.r-project.org")
    library(parallel)
  }

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
  cenmodel.list_s <- create_marginal_treatment_models(cenmodel.list)

  # Extract variables from treatment and censoring models
  extract_variables_list <- function(input) {
    extract_from_formula <- function(formula) {
      formula_terms <- terms(formula)
      response_variable <- attr(formula_terms, "response")
      response_name <- if (response_variable > 0) {
        all_vars <- all.vars(formula)
        all_vars[response_variable]
      } else {NA}
      predictor_names <- attr(formula_terms, "term.labels")
      list(response = response_name, predictors = predictor_names)
    }
    if (is.list(input)) {
      lapply(input, extract_from_formula)
    } else {
      extract_from_formula(input)
    }
  }

  trtmodel <- extract_variables_list(trtmodel.list)
  trtmodel_s <- extract_variables_list(trtmodel.list_s)
  cenmodel <- extract_variables_list(cenmodel.list)
  cenmodel_s <- extract_variables_list(cenmodel.list_s)

  # Define JAGS model for treatment and censoring
  write_jags_model <- function(trtmodel.list, cenmodel.list) {
    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    model_string <- "model{\n"

    all_parameters <- c()

    for (v in seq_along(var_info_trt)) {
      visit_trt <- var_info_trt[[v]]
      response_trt <- visit_trt$response
      predictors_trt <- visit_trt$predictors
      visit_cen <- var_info_cen[[v]]
      response_cen <- visit_cen$response
      predictors_cen <- visit_cen$predictors

      model_string <- paste0(model_string, "\nfor (i in 1:N", v, ") {\n")

      # Conditional treatment assignment model
      model_string <- paste0(model_string,
                             "\n# conditional model;\n",
                             response_trt, "[i] ~ dbern(p", v, "[i])\n",
                             "logit(p", v, "[i]) <- b", v, "0")
      for (p in seq_along(predictors_trt)) {
        model_string <- paste0(model_string, " + b", v, p, "*", predictors_trt[p], "[i]")
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }
      model_string <- paste0(model_string, "\n")

      # Censoring model
      model_string <- paste0(model_string,
                             response_cen, "[i] ~ dbern(cp", v, "[i])\n",
                             "logit(cp", v, "[i]) <- s", v, "0")
      for (p in seq_along(predictors_cen)) {
        model_string <- paste0(model_string, " + s", v, p, "*", predictors_cen[p], "[i]")
        all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
      }
      model_string <- paste0(model_string, "\n")

      # Marginal treatment assignment model
      model_string <- paste0(model_string,
                             "\n# marginal model;\n",
                             response_trt, "s[i] ~ dbern(p", v, "s[i])\n",
                             "logit(p", v, "s[i]) <- bs", v, "0")
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + bs", v, j, "*", prev_response_trt, "s[i]")
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }
      model_string <- paste0(model_string, "\n")

      # Marginal censoring model
      model_string <- paste0(model_string,
                             response_cen, "s[i] ~ dbern(cp", v, "s[i])\n",
                             "logit(cp", v, "s[i]) <- ts", v, "0")
      all_parameters <- c(all_parameters, sprintf("ts%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          model_string <- paste0(model_string, " + ts", v, j, "*", prev_response_trt, "s[i]")
          all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
        }
      }
      model_string <- paste0(model_string, "\n}\n")
    }

    # Priors section
    model_string <- paste0(model_string, "\n# Priors\n")
    for (v in seq_along(var_info_trt)) {
      num_preds_trt <- length(var_info_trt[[v]]$predictors)
      num_preds_cen <- length(var_info_cen[[v]]$predictors)

      # Treatment priors
      for (p in 0:num_preds_trt) {
        model_string <- paste0(model_string, "b", v, p, " ~ dunif(-10, 10)\n")
        all_parameters <- c(all_parameters, sprintf("b%d%d", v, p))
      }

      # Censoring priors
      for (p in 0:num_preds_cen) {
        model_string <- paste0(model_string, "s", v, p, " ~ dunif(-10, 10)\n")
        all_parameters <- c(all_parameters, sprintf("s%d%d", v, p))
      }

      # Marginal treatment priors
      model_string <- paste0(model_string, "bs", v, "0 ~ dunif(-10, 10)\n")
      all_parameters <- c(all_parameters, sprintf("bs%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "bs", v, j, " ~ dunif(-10, 10)\n")
          all_parameters <- c(all_parameters, sprintf("bs%d%d", v, j))
        }
      }

      # Marginal censoring priors
      model_string <- paste0(model_string, "ts", v, "0 ~ dunif(-10, 10)\n")
      all_parameters <- c(all_parameters, sprintf("ts%d0", v))
      if (v > 1) {
        for (j in 1:(v - 1)) {
          model_string <- paste0(model_string, "ts", v, j, " ~ dunif(-10, 10)\n")
          all_parameters <- c(all_parameters, sprintf("ts%d%d", v, j))
        }
      }
    }

    # Add the closing brace for the model block
    model_string <- paste0(model_string, "}\n")

    # Write the finalized model string to a file
    cat(model_string, file = "censoring_model.txt")

    return(unique(all_parameters))
  }

  write_jags_model(trtmodel.list, cenmodel.list)

  # Prepare data for JAGS
  prepare_jags_data <- function(data, trtmodel.list, cenmodel.list) {
    variable_info_trt <- lapply(trtmodel.list, extract_variables_list)
    variable_info_cen <- lapply(cenmodel.list, extract_variables_list)

    # Collect all variables from trtmodel and cenmodel
    all_vars <- unique(unlist(lapply(c(variable_info_trt, variable_info_cen), function(info) {
      c(info$response, info$predictors)
    })))

    # Initialize the list to pass to JAGS
    jags.data <- list()

    # Add each column of data to the jags.data list, with NA values removed
    for (col in all_vars) {
      if (!is.na(col)) {
        jags.data[[col]] <- data[[col]][!is.na(data[[col]])]
      }
    }

    # Determine N for each visit based on the length of non-NA values of the response variable
    for (v in seq_along(variable_info_trt)) {
      response_trt <- variable_info_trt[[v]]$response
      jags.data[[paste0("N", v)]] <- sum(!is.na(data[[response_trt]]))
      jags.data[[paste0("A", v, "s")]] <- data[[paste0("A", v)]][!is.na(data[[paste0("A", v)]])] # Remove the NAs
      jags.data[[paste0("C", v, "s")]] <- data[[paste0("C", v)]][!is.na(data[[paste0("C", v)]])]
    }

    return(jags.data)
  }

  jags.data <- prepare_jags_data(data, trtmodel.list, cenmodel.list)
  jags.params <- write_jags_model(trtmodel.list, cenmodel.list)

  # Run JAGS model
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
    jags.model.wd <- paste(getwd(), '/censoring_model.txt',sep='')

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
      stop("Non-parallel MCMC requires exactly 1 chain.")
    }
    jagsfit <- jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = "censoring_model.txt",
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = seed)

    out.mcmc <- as.mcmc(jagsfit)
    posterior <- as.matrix(out.mcmc[[1]])
  }

  diagnostics <- geweke.diag(out.mcmc)

  # Check diagnostics for convergence issues
  significant_indices <- which(abs(diagnostics[[1]]$z) > 1.96)
  if (length(significant_indices) > 0) {
    warning("Some parameters have not converged with Geweke index > 1.96. More iterations may be needed.")
  }

  expit <- function(x){exp(x) / (1+exp(x))}


  # Initialize arrays for storing probabilities
  n_visits <- length(trtmodel)
  n_posterior <- dim(posterior)[1]
  # Calculate the number of observations for the last visit (i.e. complete data; in this case N3)
  # last_visit_response <- trtmodel[[n_visits]]$response
  # non_missing_indices <- which(!is.na(data[[last_visit_response]]))
  # n_obs <- length(non_missing_indices)
  n_obs <- dim(data)[1]

  psc <- array(dim = c(n_visits, n_posterior, n_obs))
  psm <- array(dim = c(n_visits, n_posterior, n_obs))
  csc <- array(dim = c(n_visits, n_posterior, n_obs))
  csm <- array(dim = c(n_visits, n_posterior, n_obs))

  # Calculate probabilities for each visit
  parameter_map <- colnames(posterior)

  for (nvisit in 1:n_visits) {
    predictors_c <- trtmodel[[nvisit]]$predictors
    predictors_s <- trtmodel_s[[nvisit]]$predictors
    predictors_cen_c <- cenmodel[[nvisit]]$predictors
    predictors_cen_s <- cenmodel_s[[nvisit]]$predictors

    # design_matrix_c <- cbind(1, data[non_missing_indices, predictors_c, drop = FALSE]) # It should deal with the whole dataset! not just complete data
    # design_matrix_s <- cbind(1, data[non_missing_indices, predictors_s, drop = FALSE])
    # design_matrix_cen_c <- cbind(1, data[non_missing_indices, predictors_cen_c, drop = FALSE])
    # design_matrix_cen_s <- cbind(1, data[non_missing_indices, predictors_cen_s, drop = FALSE])

    design_matrix_c <- cbind(1, data[predictors_c]) # It should deal with the whole dataset! not just complete data
    design_matrix_s <- cbind(1, data[predictors_s])
    design_matrix_cen_c <- cbind(1, data[predictors_cen_c])
    design_matrix_cen_s <- cbind(1, data[predictors_cen_s])

    beta_indices_c <- match(c(sprintf("b%d0", nvisit),
                              sapply(1:length(predictors_c), function(p) sprintf("b%d%d", nvisit, p))),
                            parameter_map)  #> beta_indices_c [1] 1 2 3

    # beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
    #                           if (nvisit > 1) sprintf(paste0("bs%d", "1"), nvisit) ... sprintf(paste0("bs%d", nvisit-1), nvisit) else NULL),
    #                         parameter_map)  #> beta_indices_s [1] 19

    beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
                              if (nvisit > 1) sapply(1:length(predictors_s), function(p) sprintf("bs%d%d", nvisit, p)) else NULL),
                            parameter_map)

    beta_indices_cen_c <- match(c(sprintf("s%d0", nvisit),
                                sapply(1:length(predictors_cen_c), function(p) sprintf("s%d%d", nvisit, p))),
                              parameter_map)  #> beta_indices_cen_c [1] 26 27 28

    # beta_indices_cen_s <- match(c(sprintf("ts%d0", nvisit),
    #                               if (nvisit > 1) sprintf("ts%d1", nvisit) else NULL),
    #                             parameter_map)  #> beta_indices_cen_s [1] 40

    beta_indices_cen_s <- match(c(sprintf("ts%d0", nvisit),
                                  if (nvisit > 1) sapply(1:length(predictors_cen_s), function(p) sprintf("ts%d%d", nvisit, p)) else NULL),
                                parameter_map)

    for (j in 1:n_posterior) {
      psc[nvisit, j, ] <- expit(posterior[j, beta_indices_c] %*% t(design_matrix_c)) # Dimension: 40x3 %*% 3x364 3x500 (nvisit=1, j=1)
      psm[nvisit, j, ] <- expit(posterior[j, beta_indices_s] %*% t(design_matrix_s)) # Dimension: 40x1 %*% 1x364 1x500
      csc[nvisit, j, ] <- expit(posterior[j, beta_indices_cen_c] %*% t(design_matrix_cen_c)) # Dimension: 40x3 %*% 3x364 3x500
      csm[nvisit, j, ] <- expit(posterior[j, beta_indices_cen_s] %*% t(design_matrix_cen_s)) # Dimension: 40x1 %*% 1x364 1x500
    }
  }

  numerator_trt <- apply(psm, c(2, 3), prod)
  denominator_trt <- apply(psc, c(2, 3), prod)
  numerator_cen <- apply(csm, c(2, 3), prod)
  denominator_cen <- apply(csc, c(2, 3), prod)

  weights <- (numerator_trt / denominator_trt) * (numerator_cen / denominator_cen)
  wmean <- colMeans(weights)

  return(wmean)
 #  > summary(wmean)
 #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
 # 0.3037  0.9072  1.5633  2.3929  3.0790 34.6315     136

}




# Testing
# simulating causal data; (see sim_causal.R)
# with censoring;
sim.rc <- function(samplesize = 500)
{
  set.seed(123)
  expit <- function(x){exp(x)/(1+exp(x))}
  # visit 1;
  L11 <- rbinom(n=samplesize, size=1, prob=0.5)
  L21 <- rnorm(n=samplesize, mean=0, sd=1)

  C1prob <- expit(-2-0.1*L11-0.1*L21) #right-censoring also known as lost to followup require non missing baseline. If there is missing baseline, ask user to impute missing baseline variables or remove missing observations;
  C1 <- rbinom(n=samplesize, size=1, prob=C1prob)

  A1prob <- expit(0.5*L11-0.2*L21)
  A1 <- rbinom(n=samplesize, size=1, prob=A1prob)

  # visit 2;
  C2prob <- expit(-2-0.1*L11-0.1*L21 - 0.1*A1)
  C2 <- rbinom(n=samplesize, size=1, prob=C2prob)

  L12prob <- expit(0.5*A1+0.5*L11-0.2*L21)
  L12 <- rbinom(n=samplesize, size=1, prob=L12prob)
  meanL22 <- 0.5*L21-0.5*A1-0.2*L11
  L22 <- rnorm(n=samplesize, mean=meanL22, sd=1)

  A2prob <- expit(0.5*L12-0.2*L22+0.2*A1)
  A2 <- rbinom(n = samplesize, size = 1, prob = A2prob)

  # visit 3;
  C3prob <- expit(-2-0.1*L12-0.1*L22 - 0.1*A2)
  C3 <- rbinom(n=samplesize, size=1, prob=C3prob)

  L13prob <- expit(0.5*A2+0.5*L12-0.2*L22)
  L13 <- rbinom(n = samplesize, size = 1, prob = L13prob)
  meanL23 <- 0.5*L22-0.5*A2-0.2*L12
  L23 <- rnorm(n=samplesize, mean=meanL23, sd=1)

  A3prob <- expit(0.5*L13-0.2*L23+0.2*A2)
  A3 <- rbinom(n = samplesize, size = 1, prob = A3prob)

  # end-of-study outcome;
  Yprob <- expit(0.3*A3+0.1*A2-0.1*A1+0.1*L13-0.2*L23)
  Y <- rbinom(n = samplesize, size = 1, prob = Yprob)
  dat <- cbind(L11, L21, A1, L12, L22, A2, L13, L23, A3, C1, C2, C3, Y)
  dat <- data.frame(dat)
  return(dat)
}
simdat <- sim.rc(samplesize = 500)

library(tidyverse)
simdat_cen <- simdat %>%
  mutate(A1 = ifelse(C1==1, NA, A1),
         C2 = ifelse(C1==1, NA, C2),
         L12 = ifelse(C1==1, NA, L12),
         L22 = ifelse(C1==1, NA, L22),
         A2 = ifelse(C1==1, NA, A2),
         L13 = ifelse(C1==1, NA, L13),
         L23 = ifelse(C1==1, NA, L23),
         A3 = ifelse(C1==1, NA, A3),
         C3 = ifelse(C1==1, NA, C3),
         Y = ifelse(C1==1, NA, Y)) %>%
  mutate(L12 = ifelse(C2==1, NA, L12),
         L22 = ifelse(C2==1, NA, L22),
         A2 = ifelse(C2==1, NA, A2),
         L13 = ifelse(C2==1, NA, L13),
         L23 = ifelse(C2==1, NA, L23),
         A3 = ifelse(C2==1, NA, A3),
         C3 = ifelse(C2==1, NA, C3),
         Y = ifelse(C2==1, NA, Y)) %>%
  mutate(L13 = ifelse(C3==1, NA, L13),
         L23 = ifelse(C3==1, NA, L23),
         A3 = ifelse(C3==1, NA, A3),
         Y = ifelse(C3==1, NA, Y))


# Censored data
start<-Sys.time()
weights <- bayesweight_cen(trtmodel.list = list(A1 ~ L11 + L21,
                                                 A2 ~ L11 + L21 + L12 + L22 + A1,
                                                 A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2),
                            cenmodel.list = list(C1 ~ L11 + L21,
                                                 C2 ~ L11 + L21 + A1,
                                                 C3 ~ L11 + L21 + A1 + L12 + L22 + A2), # kitchen sink
                            data = simdat_cen,
                            n.iter = 250,
                            n.burnin = 50,
                            n.thin = 5)
Sys.time()-start
weights




######################################################################

# Calculating frequentist weights R code;

psc1m <- glm(A1 ~ L11 + L21, family = binomial(link = "logit") , data = simdat_cen)

psc2m <- glm(A2 ~ L11 + L21 + L12 + L22 + A1, family = binomial(link = "logit") , data = simdat_cen)

psc3m <- glm(A3 ~ L11 + L21 + L12 + L22 + A1 + L13 + L23 + A2, family = binomial(link = "logit") , data = simdat_cen)

psc1 <- predict(psc1m, newdata = simdat_cen, type="response")

psc2 <- predict(psc2m, newdata = simdat_cen, type="response")

psc3 <- predict(psc3m, newdata = simdat_cen, type="response")



psm1m <- glm(A1 ~ 1, family = binomial(link = "logit") , data = simdat_cen)

psm2m <- glm(A2 ~ A1, family = binomial(link = "logit") , data = simdat_cen)

psm3m <- glm(A3 ~ A1 + A2, family = binomial(link = "logit") , data = simdat_cen)

psm1 <- predict(psm1m, newdata = simdat_cen, type="response")

psm2 <- predict(psm2m, newdata = simdat_cen, type="response")

psm3 <- predict(psm3m, newdata = simdat_cen, type="response")



weight_trt1 <- (psm1/psc1)

weight_trt2 <- (psm2/psc2)

weight_trt3 <- (psm3/psc3)



weight_trt<- weight_trt1*weight_trt2*weight_trt3



csc1m <- glm(C1 ~ L11 + L21, family = binomial(link = "logit") , data = simdat_cen)

csc2m <- glm(C2 ~ L11 + L21 + A1, family = binomial(link = "logit") , data = simdat_cen)

csc3m <- glm(C3 ~ L11 + L21 + L12 + L22 + A1 + A2, family = binomial(link = "logit") , data = simdat_cen)

csc1 <- predict(csc1m, newdata = simdat_cen, type="response")

csc2 <- predict(csc2m, newdata = simdat_cen, type="response")

csc3 <- predict(csc3m, newdata = simdat_cen, type="response")



csm1m <- glm(C1 ~ 1, family = binomial(link = "logit") , data = simdat_cen)

csm2m <- glm(C2 ~ A1, family = binomial(link = "logit") , data = simdat_cen)

csm3m <- glm(C3 ~ A1 + A2, family = binomial(link = "logit") , data = simdat_cen)

csm1 <- predict(csm1m, newdata = simdat_cen, type="response")

csm2 <- predict(csm2m, newdata = simdat_cen, type="response")

csm3 <- predict(csm3m, newdata = simdat_cen, type="response")



weight_cen1 <- (csm1/csc1)

weight_cen2 <- (csm2/csc2)

weight_cen3 <- (csm3/csc3)



weight_cen <- weight_cen1*weight_cen2*weight_cen3



weight_freq <- weight_trt*weight_cen

summary(weight_freq)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's

#  0.2356  0.8193  1.1823  1.5666  1.9554 13.8196     136
