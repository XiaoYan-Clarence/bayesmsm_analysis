# Compare Step 2 to Frequentist
# Source: https://kuan-liu.github.io/causal_Quarto/section3.html

# Binary outcome
# setwd("C:/Users/YanXi/Downloads")
testdata_binary <- readr::read_csv("R/binary_outcome_data.csv")

library(WeightIt)
Wmsm <- weightitMSM(
  list(a_1 ~ w1 + w2 + L1_1 + L2_1,
       a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
  data = testdata_binary,
  method = "ps",
  stabilize = TRUE)

# summary(Wmsm)

library(survey)
msm_design <- svydesign(~1, weights = Wmsm$weights, data = testdata_binary)
fitMSM <- svyglm(y ~ a_1*a_2, design = msm_design)
# summary(fitMSM)

APO_11 <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1), type = "response")
APO_00 <- predict(fitMSM, newdata = data.frame(a_1=0,a_2=0), type = "response")
APO_11 - APO_00
# link     SE
# 1 -0.10715 0.0366

set.seed(123)
boot.est <- rep(NA, 1000)
for (i in 1:1000){
  boot.idx <- sample(1:dim(testdata_binary)[1], size = dim(testdata_binary)[1], replace = T)
  boot.data <- testdata_binary[boot.idx,]
  msm_design <- svydesign(~1, weights = Wmsm$weights, data = boot.data)
  fitMSM <- svyglm(y ~ a_1*a_2, design = msm_design)
  boot.est[i] <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1), type = "response")[1] - predict(fitMSM, newdata = data.frame(a_1=0,a_2=0), type = "response")[1]
}

# SE of ATE;
sd(boot.est)
# > sd(boot.est)
# [1] 0.04250907

#95% CI
quantile(boot.est, probs = c(0.025, 0.975))
# > #95% CI
#   > quantile(boot.est, probs = c(0.025, 0.975))
# 2.5%       97.5%
#   -0.19420932 -0.03258665


#bayesian
wloglik_binomial <- function(param,
                             Y,
                             A,
                             weight){
  # number of observations;
  n <- length(Y)
  beta <- param[1:dim(A)[2]] # causal parameters on the log-odds scale (no sigma for binomial?)
  mmat <- as.matrix(A)
  eta<-mmat %*% beta # linear predictor
  logl <- Y*eta - log(1+exp(eta))
  wlogl<-sum(weight*logl)
  return(wlogl)
}


inits1<-c(0.1,0.1,0.1,0.1) #three mean parameters + 1 variance parameter
nboot <- 1000
bootest<-numeric(nboot)

Wmean<-Wmsm$weights
# Wmean<-rep(1:dim(testdata_continuous)[1])
expit <- function(x){exp(x) / (1+exp(x))}

for (j in 1:nboot) {
  alpha <- as.numeric(rdirichlet(1, rep(1.0, length(testdata$y))))
  maxim <- optim(inits1,
                 fn=wloglik_binomial,
                 Y=testdata_binary$y,
                 A=cbind(1,testdata_binary$a_1, testdata_binary$a_2, testdata_binary$a_1*testdata_binary$a_2), #three mean parameters (intercept + coefficient for a_1 and coefficient for a_2);
                 weight=alpha*Wmean,
                 control=list(fnscale=-1), method='BFGS', hessian=F)
  bootest[j] <- expit(maxim$par[1]+maxim$par[2]+maxim$par[3]+maxim$par[4]) - expit(maxim$par[1]) #difference on the mean of Y between always treated and never treated;

  if (j %% 100 == 0) {
    print(j)
  }
}

mean(bootest)
sd(bootest)










# Continuous outcome (same as in example_BMSMs_for_non_repeated_continuous_Y.R)
testdata_continuous <- readr::read_csv("R/continuous_outcome_data.csv")

library(WeightIt)
Wmsm <- weightitMSM(
  list(a_1 ~ w1 + w2 + L1_1 + L2_1,
       a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
  data = testdata_continuous,
  method = "ps",
  stabilize = TRUE)

# summary(Wmsm)

library(survey)
msm_design <- svydesign(~1, weights = Wmsm$weights, data = testdata_continuous)
fitMSM <- svyglm(y ~ a_1*a_2, design = msm_design)
# summary(fitMSM)

APO_11 <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1))
APO_00 <- predict(fitMSM, newdata = data.frame(a_1=0,a_2=0))
APO_11 - APO_00
# link     SE
# -3.1134 0.0888


fitMSM2 <- glm(y ~ a_1*a_2, data=testdata_continuous)
summary(fitMSM)
APO_11 <- predict(fitMSM2, newdata = data.frame(a_1=1,a_2=1))
APO_00 <- predict(fitMSM2, newdata = data.frame(a_1=0,a_2=0))
APO_11 - APO_00


# using bootstrap to obtain SE and confidence interval of the ATE;
set.seed(123)
boot.est <- rep(NA, 1000)
for (i in 1:1000){
  boot.idx <- sample(1:dim(testdata_continuous)[1], size = dim(testdata_continuous)[1], replace = T)
  boot.data <- testdata_continuous[boot.idx,]
  msm_design <- svydesign(~1, weights = Wmsm$weights, data = boot.data)
  fitMSM <- svyglm(y ~ a_1*a_2, design = msm_design)
  boot.est[i] <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1))[1] - predict(fitMSM, newdata = data.frame(a_1=0,a_2=0))[1]
}

# SE of ATE;
sd(boot.est)
# > sd(boot.est)
# [1] 0.1000114

#95% CI
quantile(boot.est, probs = c(0.025, 0.975))
# > #95% CI
#   > quantile(boot.est, probs = c(0.025, 0.975))
# 2.5%     97.5%
#   -3.314099 -2.910470



#bayesian

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

inits1<-c(0.1,0.1,0.1,0.1,4) #three mean parameters + 1 variance parameter
nboot <- 1000
bootest<-numeric(nboot)

Wmean<-Wmsm$weights
Wmean<-rep(1:dim(testdata_continuous)[1])

for (j in 1:nboot) {
  alpha <- as.numeric(rdirichlet(1, rep(1.0, length(testdata$y))))

  maxim <- optim(inits1,
                 fn=wloglik_normal,
                 Y=testdata_continuous$y,
                 A=cbind(1,testdata_continuous$a_1, testdata_continuous$a_2, testdata_continuous$a_1*testdata_continuous$a_2), #three mean parameters (intercept + coefficient for a_1 and coefficient for a_2);
                 weight=alpha*Wmean,
                 control=list(fnscale=-1), method='BFGS', hessian=F)
  bootest[j] <- maxim$par[2]+maxim$par[3]+maxim$par[4] #difference on the mean of Y between always treated and never treated;

  if (j %% 100 == 0) {
    print(j)
  }
}

mean(bootest)
sd(bootest)

> mean(bootest)
[1] -3.090517
> sd(bootest)
[1] 0.1094727

