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

summary(Wmsm)



library(survey)
msm_design <- svydesign(~1, weights = Wmsm$weights, data = testdata_binary)
fitMSM <- svyglm(y ~ a_1+a_2, design = msm_design)
summary(fitMSM)

APO_11 <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1))
APO_00 <- predict(fitMSM, newdata = data.frame(a_1=0,a_2=0))
APO_11 - APO_00

# link     SE
# 1 -0.095954 0.0316

set.seed(123)
boot.est <- rep(NA, 1000)
for (i in 1:1000){
  boot.idx <- sample(1:dim(testdata_binary)[1], size = dim(testdata_binary)[1], replace = T)
  boot.data <- testdata_binary[boot.idx,]
  msm_design <- svydesign(~1, weights = Wmsm$weights, data = boot.data)
  fitMSM <- svyglm(y ~ a_1*a_2, design = msm_design)
  boot.est[i] <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1))[1] - predict(fitMSM, newdata = data.frame(a_1=0,a_2=0))[1]
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



# Compared to the result of BMSM, step 2:
# > mean(bootoutput3$ATE)
# [1] -0.4507853
# > var(bootoutput3$ATE)
# [1] 0.03255215
# > sqrt(var(bootoutput3$ATE))
# [1] 0.1804221
# > quantile(bootoutput3$ATE, probs=c(0.025,0.975))
# 2.5%      97.5%
#   -0.8202015 -0.1112203















# Continuous outcome (same as in example_BMSMs_for_non_repeated_continuous_Y.R)
testdata_continuous <- readr::read_csv("R/continuous_outcome_data.csv")

library(WeightIt)
Wmsm <- weightitMSM(
  list(a_1 ~ w1 + w2 + L1_1 + L2_1,
       a_2 ~ w1 + w2 + L1_1 + L2_1 + L1_2 + L2_2 + a_1),
  data = testdata_continuous,
  method = "ps",
  stabilize = TRUE)

summary(Wmsm)



library(survey)
msm_design <- svydesign(~1, weights = Wmsm$weights, data = testdata_continuous)
fitMSM <- svyglm(y ~ a_1+a_2, design = msm_design)
summary(fitMSM)

APO_11 <- predict(fitMSM, newdata = data.frame(a_1=1,a_2=1))
APO_00 <- predict(fitMSM, newdata = data.frame(a_1=0,a_2=0))
APO_11 - APO_00
# link     SE
# 1 -3.1611 0.0758

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



# Compared to the result of BMSM, step 2:
# > mean(bootoutput2$ATE)
# [1] -3.164354
# > var(bootoutput2$ATE)
# [1] 0.009445813
# > sqrt(var(bootoutput2$ATE))
# [1] 0.09718957
# > quantile(bootoutput2$ATE, probs=c(0.025,0.975))
# 2.5%     97.5%
#   -3.356886 -2.992084
# super similar!
