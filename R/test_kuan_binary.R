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

testdata3 <- readr::read_csv("R/binary_outcome_data.csv")

results.it <- matrix(NA, 1, 3) #result matrix, three columns for bootest, effect_ref, and effect_comp;

alpha <- as.numeric(rdirichlet(1, rep(1.0, length(Y))))

maxim <- optim(inits1,
               fn = wloglik_binomial,
               Y = Y,
               A = A,
               weight = alpha * wmean,
               control = list(fnscale = -1),
               method = optim_method,
               hessian = FALSE)

names(maxim$par) <- c("(Intercept)", variables$predictors)




