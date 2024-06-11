cat("
model{
for(i in 1:N){

# visit 1;

a_1s[i] ~ dbern(pa1s[i])
pa1s[i] <- ilogit(bs11)

a_1[i] ~ dbern(pa1[i])
pa1[i] <- ilogit(b11 + b12*w1[i] + b13*w2[i] + b14*L1_1[i] +b15*L2_1[i])

# visit 2;
a_2s[i] ~ dbern(pa2s[i])
pa2s[i] <- ilogit(bs21 + bs22*a_1s[i])

a_2[i] ~ dbern(pa2[i])
pa2[i] <- ilogit(b21+b22*w1[i]+b23*w2[i]+b24*L1_1[i]+b25*L2_1[i]+b26*a_1[i]+b27*L1_2[i]+b28*L2_2[i])

# export quantity in full posterior specification;
w[i] <- (pa1s[i]*pa2s[i])/(pa1[i]*pa2[i])
}

#prior;
#all parameters here use dnorm(0,.01);
bs11~dnorm(0,.01)
bs21~dnorm(0,.01)
bs22~dnorm(0,.01)

b11~dnorm(0,.01)
b12~dnorm(0,.01)
b13~dnorm(0,.01)
b14~dnorm(0,.01)
b15~dnorm(0,.01)

b21~dnorm(0,.01)
b22~dnorm(0,.01)
b23~dnorm(0,.01)
b24~dnorm(0,.01)
b25~dnorm(0,.01)
b26~dnorm(0,.01)
b27~dnorm(0,.01)
b28~dnorm(0,.01)

}
", fill = TRUE, file = "treatment_model_2visit.txt")

#Bayesian PS;
library(tidyverse)
library(R2jags)
options(scipen = 999)

causaldata <- read.csv("continuous_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")

jags.data<-list(w1=causaldata$w1, w2=causaldata$w2,
                L1_1 = causaldata$L1_1,L2_1 = causaldata$L2_1,
                L1_2 = causaldata$L1_2,L2_2 = causaldata$L2_2,
                a_1 = causaldata$a_1,  a_1s = causaldata$a_1,
                a_2 = causaldata$a_2,  a_2s = causaldata$a_2,
                N = dim(causaldata)[1])

jags.params<-c("w")

set.seed(890123)
jagsfit<- jags(data = jags.data,
               parameters.to.save=jags.params,
               n.iter = 25000,
               model.file = "treatment_model_2visit.txt",
               n.chains = 1,
               n.burnin = 15000,
               n.thin = 5)

out.mcmc <- as.mcmc(jagsfit)
wmean <- data.frame(id = gsub(".*[[]([^.]+)[]].*", "\\1", colnames(out.mcmc[[1]])[-1]), weight = colMeans(out.mcmc[[1]])[-1])
wmean <- wmean %>% arrange(as.numeric(id))

# write.csv(wmean, file = "bayesweight_ycont.csv",row.names = F)

causaldata$bayeswt<- wmean$weight
write.csv(causaldata, file = "continuous_outcome_data2.csv",row.names = F)

# binary case;
causaldata <- read.csv("binary_outcome_data.csv", header = TRUE, fileEncoding="UTF-8-BOM")

jags.data<-list(w1=causaldata$w1, w2=causaldata$w2,
                L1_1 = causaldata$L1_1,L2_1 = causaldata$L2_1,
                L1_2 = causaldata$L1_2,L2_2 = causaldata$L2_2,
                a_1 = causaldata$a_1,  a_1s = causaldata$a_1,
                a_2 = causaldata$a_2,  a_2s = causaldata$a_2,
                N = dim(causaldata)[1])

jags.params<-c("w")

set.seed(890123)
jagsfit<- jags(data = jags.data,
               parameters.to.save=jags.params,
               n.iter = 25000,
               model.file = "treatment_model_2visit.txt",
               n.chains = 1,
               n.burnin = 15000,
               n.thin = 5)

out.mcmc <- as.mcmc(jagsfit)
wmean <- data.frame(id = gsub(".*[[]([^.]+)[]].*", "\\1", colnames(out.mcmc[[1]])[-1]), weight = colMeans(out.mcmc[[1]])[-1])
wmean <- wmean %>% arrange(as.numeric(id))

causaldata$bayeswt<- wmean$weight
write.csv(causaldata, file = "binary_outcome_data2.csv",row.names = F)
