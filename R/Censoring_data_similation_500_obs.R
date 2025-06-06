simc.r <- function(samplesize = 500)

{

  set.seed(123)

  expit <- function(x){exp(x)/(1+exp(x))}

  # visit 1;

  L11 <- rbinom(n=samplesize, size=1, prob=0.5)

  L21 <- rnorm(n=samplesize, mean=0, sd=1)

  A1prob <- expit(0.5*L11-0.1*L21)

  A1 <- rbinom(n=samplesize, size=1, prob=A1prob)

  # visit 2;

  L12prob <- expit(0.5*A1+0.5*L11-0.2*L21)

  L12 <- rbinom(n=samplesize, size=1, prob=L12prob)

  meanL22 <- 0.5*L21-0.5*A1-0.2*L11

  L22 <- rnorm(n=samplesize, mean=meanL22, sd=1)

  A2prob <- expit(0.5*L12-0.1*L22+0.2*A1)

  A2 <- rbinom(n = samplesize, size = 1, prob = A2prob)

  # visit 3;

  Cprob <- expit(-2-0.1*L12-0.1*L22 - 0.1*A2)

  C <- rbinom(n=samplesize, size=1, prob=Cprob)

  L13prob <- expit(0.5*A2+0.5*L12-0.2*L22)

  L13 <- rbinom(n = samplesize, size = 1, prob = L13prob)

  meanL23 <- 0.5*L22-0.5*A2-0.2*L12

  L23 <- rnorm(n=samplesize, mean=meanL23, sd=1)

  A3prob <- expit(0.5*L13-0.1*L23+0.2*A2)

  A3 <- rbinom(n = samplesize, size = 1, prob = A3prob)

  # end-of-study outcome;

  Yprob <- expit(0.3*A3+0.1*A2-0.1*A1+0.1*L13-0.2*L23)

  Y <- rbinom(n = samplesize, size = 1, prob = Yprob)

  dat <- cbind(L11, L21, A1, L12, L22, A2, L13, L23, A3, C, Y)

  dat <- data.frame(dat)

  return(dat)

}



simdat <- simc.r(samplesize = 500)

library(tidyverse)

simdat_cen <- simdat %>%

  mutate(L13 = ifelse(C==1, NA, L13),

         L23 = ifelse(C==1, NA, L23),

         A3 = ifelse(C==1, NA, A3),

         Y = ifelse(C==1, NA, Y))


# Save simdat_cen to sim.csv and sim_causal.csv
write.csv(simdat, file = "sim_causal.csv", row.names = FALSE)
write.csv(simdat_cen, file = "sim_causal_cen.csv", row.names = FALSE)
