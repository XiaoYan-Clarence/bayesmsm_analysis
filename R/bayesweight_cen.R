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
                                                 A3 ~ L12 + L22 + L13 + L23 + A2),
                            cenmodel.list = list(C1 ~ L11 + L21,
                                                 C2 ~ L11 + L21 + A1,
                                                 C3 ~ L12 + L22 + A2),
                            data,
                            n.iter = 25000,
                            n.burnin = 1500,
                            n.thin = 5,
                            parallel = FALSE,
                            n.chains = 1,
                            seed = 890123) {

  # Load all the required R packages;
  if (!require(R2jags)){
    install.packages("R2jags",repos="http://cran.r-project.org")
    library(R2jags)
  }
  if (!require(parallel)){
    install.packages("parallel",repos="http://cran.r-project.org")
    library(parallel)
  }

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
  censmodel <- extract_variables_list(cenmodel.list)

  # Prepare data for JAGS
  prepare_jags_data <- function(data, trtmodel.list, cenmodel.list) {
    variable_info_trt <- lapply(trtmodel.list, extract_variables_list)
    variable_info_cen <- lapply(cenmodel.list, extract_variables_list)

    jags.data <- list(N = nrow(data))

    for (info in variable_info_trt) {
      response <- info$response
      predictors <- info$predictors
      jags.data[[response]] <- data[[response]]
      for (predictor in predictors) {
        jags.data[[predictor]] <- data[[predictor]]
      }
    }

    for (info in variable_info_cen) {
      response <- info$response
      predictors <- info$predictors
      jags.data[[response]] <- data[[response]]
      for (predictor in predictors) {
        jags.data[[predictor]] <- data[[predictor]]
      }
    }

    return(jags.data)
  }

  jags.data <- prepare_jags_data(data, trtmodel.list, cenmodel.list)

  # Define JAGS model for treatment and censoring
  write_jags_model <- function(trtmodel.list, cenmodel.list) {
    var_info_trt <- lapply(trtmodel.list, extract_variables_list)
    var_info_cen <- lapply(cenmodel.list, extract_variables_list)

    model_string <- "model{\n#N = nobs\nfor(i in 1:N){\n"

    for (v in seq_along(var_info_trt)) {
      visit_trt <- var_info_trt[[v]]
      response_trt <- visit_trt$response
      predictors_trt <- visit_trt$predictors
      visit_cen <- var_info_cen[[v]]
      response_cen <- visit_cen$response
      predictors_cen <- visit_cen$predictors

      model_string <- paste0(model_string,
                             "\n# visit ", v, ";\n",
                             "# marginal treatment assignment model, visit ", v, ";\n",
                             response_trt, "s[i] ~ dbern(p", response_trt, "s[i])\n")

      if (v == 1) {
        model_string <- paste0(model_string, "p", response_trt, "s[i] <- ilogit(bs", v, "0)\n")
      } else {
        bs_terms <- paste0("bs", v, "0")
        for (j in 1:(v-1)) {
          prev_response_trt <- var_info_trt[[j]]$response
          bs_terms <- paste(bs_terms, " + bs", v, j, "*", prev_response_trt, "s[i]", sep = "")
        }
        model_string <- paste0(model_string, "p", response_trt, "s[i] <- ilogit(", bs_terms, ")\n")
      }

      model_string <- paste0(model_string,
                             "\n# conditional treatment assignment model, visit ", v, ";\n",
                             response_trt, "[i] ~ dbern(p", response_trt, "[i])\n",
                             "p", response_trt, "[i] <- ilogit(b", v, "0")
      for (p in seq_along(predictors_trt)) {
        model_string <- paste0(model_string, " + b", v, p, "*", predictors_trt[p], "[i]")
      }
      model_string <- paste0(model_string, ")\n")

      model_string <- paste0(model_string,
                             "\n# censoring model, visit ", v, ";\n",
                             response_cen, "[i] ~ dbern(p", response_cen, "[i])\n",
                             "p", response_cen, "[i] <- ilogit(c", v, "0")
      for (p in seq_along(predictors_cen)) {
        model_string <- paste0(model_string, " + c", v, p, "*", predictors_cen[p], "[i]")
      }
      model_string <- paste0(model_string, ")\n")
    }

    model_string <- paste0(model_string, "\n# export quantity in full posterior specification;\n",
                           "w[i] <- (", paste(sapply(seq_along(var_info_trt), function(x) paste0("p", var_info_trt[[x]]$response, "s[i]")), collapse = "*"),
                           ")/(", paste(sapply(seq_along(var_info_trt), function(x) paste0("p", var_info_trt[[x]]$response, "[i]")), collapse = "*"), ")\n}\n\n#prior;\n")

    for (v in seq_along(var_info_trt)) {
      visit_trt <- var_info_trt[[v]]
      predictors_trt <- visit_trt$predictors
      num_preds_trt <- length(predictors_trt) + 1

      model_string <- paste0(model_string, "bs", v, "0~dnorm(0,.01)\n")
      if (v > 1) {
        for (j in 1:(v-1)) {
          model_string <- paste0(model_string, "bs", v, j, "~dnorm(0,.01)\n")
        }
      }
      for (p in 0:(num_preds_trt - 1)) {
        model_string <- paste0(model_string, "b", v, p, "~dnorm(0,.01)\n")
      }

      visit_cen <- var_info_cen[[v]]
      predictors_cen <- visit_cen$predictors
      num_preds_cen <- length(predictors_cen) + 1

      model_string <- paste0(model_string, "c", v, "0~dnorm(0,.01)\n")
      for (p in 0:(num_preds_cen - 1)) {
        model_string <- paste0(model_string, "c", v, p, "~dnorm(0,.01)\n")
      }
    }

    model_string <- paste(model_string, "}\n", sep = "")
    cat(model_string, file = "treatment_model.txt")
    return(model_string)
  }

  write_jags_model(trtmodel.list, cenmodel.list)

  # Run JAGS model
  if (parallel == TRUE) {
    if (n.chains == 1) {
      stop("Parallel MCMC requires at least 2 chains. Computing is running on 1 core per chain.")
    }
    available_cores <- detectCores(logical = FALSE)
    if (n.chains >= available_cores) {
      stop(paste("Parallel MCMC requires 1 core per chain. You have", available_cores, "cores. We recommend using", available_cores - 2, "cores."))
    }
    cl <- makeCluster(n.chains)
    clusterExport(cl, varlist = c("jags.data", "jags.params", "n.iter", "n.burnin", "n.thin", "model.file", "seed"))
    jagsfit <- parLapply(cl, 1:n.chains, function(i) {
      library(R2jags)
      jags.parallel(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = "treatment_model.txt",
                    n.chains = n.chains,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    n.cluster = n.chains,
                    jags.seed = seed)
    })
    stopCluster(cl)

    out.mcmc <- as.mcmc(jagsfit)
    posterior <- do.call(rbind, lapply(out.mcmc, as.matrix))
  } else if (parallel == FALSE) {
    if (n.chains != 1) {
      stop("Non-parallel MCMC requires exactly 1 chain.")
    }
    jagsfit <- jags(data = jags.data,
                    parameters.to.save = jags.params,
                    model.file = "treatment_model.txt",
                    n.chains = 1,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.thin = n.thin,
                    jags.seed = seed)

    out.mcmc <- as.mcmc(jagsfit)
    posterior <- as.matrix(out.mcmc[[1]])
  }

  # Calculate probabilities
  n_visits <- length(trtmodel)
  n_posterior <- dim(posterior)[1]
  n_obs <- nrow(data)

  psc <- array(dim = c(n_visits, n_posterior, n_obs))
  psm <- array(dim = c(n_visits, n_posterior, n_obs))

  parameter_map <- colnames(posterior)

  for (nvisit in 1:n_visits) {
    predictors_c <- trtmodel[[nvisit]]$predictors
    predictors_s <- trtmodel_s[[nvisit]]$predictors

    design_matrix_c <- cbind(1, data[, predictors_c, drop = FALSE])
    design_matrix_s <- cbind(1, data[, predictors_s, drop = FALSE])

    beta_indices_c <- match(c(sprintf("b%d0", nvisit),
                              sapply(1:length(predictors_c), function(p) sprintf("b%d%d", nvisit, p))),
                            parameter_map)

    beta_indices_s <- match(c(sprintf("bs%d0", nvisit),
                              if (nvisit > 1) sprintf("bs%d1", nvisit) else NULL),
                            parameter_map)

    for (j in 1:n_posterior) {
      psc[nvisit, j, ] <- expit(posterior[j, beta_indices_c] %*% t(design_matrix_c))
      psm[nvisit, j, ] <- expit(posterior[j, beta_indices_s] %*% t(design_matrix_s))
    }
  }

  numerator <- apply(psm, c(2, 3), prod)
  denominator <- apply(psc, c(2, 3), prod)

  weights <- numerator / denominator
  wmean <- colMeans(weights)

  return(wmean)
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

