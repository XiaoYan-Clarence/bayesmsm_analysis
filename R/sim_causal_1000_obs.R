#simulating causal data;
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

# full length of data is
dim(simdat_cen)[1]
# [1] 500
# number of observed observations at visit 1 that has L11 and L21 measured are
sum(simdat_cen$C1==0)
# [1] 445
# number of observed observations at visit 2
sum(simdat_cen$C2==0, na.rm = T)
# [1] 400
# number of observed observations at visit 3
sum(simdat_cen$C3==0, na.rm = T)
# [1] 364

# no censoring;
sim.r <- function(samplesize = 500)
{
  set.seed(123)
  expit <- function(x){exp(x)/(1+exp(x))}
  # visit 1;
  L11 <- rbinom(n=samplesize, size=1, prob=0.5)
  L21 <- rnorm(n=samplesize, mean=0, sd=1)
  A1prob <- expit(0.5*L11-0.2*L21)
  A1 <- rbinom(n=samplesize, size=1, prob=A1prob)
  # visit 2;
  # c2prob <-
  L12prob <- expit(0.5*A1+0.5*L11-0.2*L21)
  L12 <- rbinom(n=samplesize, size=1, prob=L12prob)
  meanL22 <- 0.5*L21-0.5*A1-0.2*L11
  L22 <- rnorm(n=samplesize, mean=meanL22, sd=1)
  A2prob <- expit(0.5*L12-0.2*L22+0.2*A1)
  A2 <- rbinom(n = samplesize, size = 1, prob = A2prob)
  # visit 3;
  # c3prob <-
  L13prob <- expit(0.5*A2+0.5*L12-0.2*L22)
  L13 <- rbinom(n = samplesize, size = 1, prob = L13prob)
  meanL23 <- 0.5*L22-0.5*A2-0.2*L12
  L23 <- rnorm(n=samplesize, mean=meanL23, sd=1)
  A3prob <- expit(0.5*L13-0.2*L23+0.2*A2)
  A3 <- rbinom(n = samplesize, size = 1, prob = A3prob)
  # end-of-study outcome;
  Yprob <- expit(0.3*A3+0.1*A2-0.1*A1+0.1*L13-0.2*L23)
  Y <- rbinom(n = samplesize, size = 1, prob = Yprob)
  dat <- cbind(L11, L21, A1, L12, L22, A2, L13, L23, A3, Y)
  dat <- data.frame(dat)
  return(dat)
}
simdat <- sim.r(samplesize = 500)


