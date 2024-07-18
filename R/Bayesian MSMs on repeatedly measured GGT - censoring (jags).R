
###########################################################################################

################################# Primary outcome 1 #######################################

# Bayesian MSMs on repeatedly measured GGT

###########################################################################################



# Bayesian MSM JAGS model;

cat( " model {



     for (i in 1:N2) {



     #among patients who had none-missing 1 year primary outcome and who were OVT free between 0 and 1;

     # conditional model among those who were OVT free between 0 to 1 year;

     z1[i] ~ dbern(p2[i])

     p2[i] <- (z0[i]==0)*ilogit(b0 + inprod(b , bx1[i,]))+(z0[i]==1)



     # marginal model among those who were OVT free between 0 to 1 year;

     z1s[i] ~ dbern(p2s[i])

     logit(p2s[i]) <- (z0[i]==0)*ilogit(bs1) +(z0[i]==1)



     }



     for (i in 1:N1) { #among patients who had none-missing baseline (all patients here);



     # conditional model;

     z0[i] ~ dbern(p1[i])

     logit(p1[i]) <- a0 + inprod(a , ax0[i,])



     c1[i] ~ dbern(cp1[i])

     logit(cp1[i])<- c0 + inprod(c , cx1[i,])



     # marginal model;

     z0s[i] ~ dbern(p1s[i])

     logit(p1s[i]) <- as0



     c1s[i] ~ dbern(cp1s[i])

     logit(cp1s[i]) <- cs0 + cs1*z0[i]



    }



     # Priors;

     a0 ~ dunif(-10,10)

     as0 ~ dunif(-10,10)

     b0 ~ dunif(-10,10)

     bs1 ~ dunif(-10,10)

     c0 ~ dunif(-10,10)

     cs0 ~ dunif(-10,10)

     cs1 ~ dunif(-10,10)



     for (j1 in 1:14){

     a[j1] ~ dunif(-10,10)

     b[j1] ~ dunif(-10,10)

     }



     for (j2 in 1:15){

     c[j2] ~ dunif(-10,10)

     }



     }",

     file = "model_unif_bmsm.txt")



#setting up bayesian data;

mypsc<-data.frame("id"=lab_wide3$id, ovt_t0, ovt_t12,

                  lnapri0, ggt0, lnapri1, ggt1, lnapri2, ggt2,

                  sex, age, largeduct,

                  ibd0, aih0, urso0, stero0, immun0, bio0, FASA0, fibro0, hepa0,

                  ibd1, aih1, urso1, stero1, immun1, bio1, FASA1, fibro1, hepa1,

                  ibd2, aih2, urso2, stero2, immun2, bio2, FASA2, fibro2, hepa2)



#creating missing indicator;

mypsc$C1<- ifelse(is.na(mypsc$ggt1) | is.na(mypsc$lnapri1),1,0)

mypsc<-mypsc[order(ovt_t0),]



# complete data N3 = 323 and among those N2= 294 were OVT free between 0 to 1 year;

mypsc2<-subset(mypsc, C1==0)



ax0 <- mypsc[,c("lnapri0", "ggt0", "sex", "age", "largeduct", "ibd0", "aih0", "urso0","stero0", "immun0", "bio0", "FASA0", "fibro0", "hepa0")]

cx1 <- cbind(mypsc$ovt_t0,ax0)

bx1 <- mypsc2[,c("lnapri1","ggt1", "sex", "age", "largeduct", "ibd1", "aih1", "urso1","stero1", "immun1", "bio1", "FASA1", "fibro1", "hepa1")]



jags.data<-list(N1=401, N2=323,

                z0=mypsc$ovt_t0, z0s=mypsc$ovt_t0, z1=mypsc2$ovt_t12, z1s=mypsc2$ovt_t12,

                c1=mypsc$C1, c1s=mypsc$C1,

                ax0=ax0, bx1=bx1, cx1=cx1)





jags.params<-c("a0","as0","b0","bs1","c0","cs0","cs1","a","b","c")





jagsfit<- jags(data = jags.data, parameters.to.save=jags.params, n.iter = 50000,

               model.file = "model_unif_bmsm.txt",jags.seed = 123, n.chains = 1,n.burnin = 25000, n.thin = 5)



# use as.mcmmc to convert rjags object into mcmc.list

jags.mcmc <- as.mcmc(jagsfit)

out.mcmc <- as.matrix(jags.mcmc[[1]])

geweke.diag(out.mcmc)



samplesize <- 5000

ntot<-401



obs_prob1<-matrix(NA, samplesize, ntot)

exp_prob1<-matrix(NA, samplesize, ntot)

obs_prob2<-matrix(NA, samplesize, ntot)

exp_prob2<-matrix(NA, samplesize, ntot)



obs_cprob1<-matrix(NA, samplesize, ntot)

exp_cprob1<-matrix(NA, samplesize, ntot)

obs_cprob2<-matrix(NA, samplesize, ntot)

exp_cprob2<-matrix(NA, samplesize, ntot)



# mcmcweight_ws<-array(NA, dim = c(ntot,3,samplesize))



#calculating the MCMC weights;

for (i2 in 1:(samplesize)){

  for (j2 in 1:(ntot)){



    exp_prob1[i2,j2] <- (exp(mypsc$ovt_t0[j2]*out.mcmc[i2,16]))/(1.0+exp(out.mcmc[i2,16]))

    obs_prob1[i2,j2] <- (exp(mypsc$ovt_t0[j2]*(out.mcmc[i2,15] + out.mcmc[i2,c(1,7:14,2:6)]%*%t(mypsc[j2,c(4:5,10:21)]))))/(1.0+exp(out.mcmc[i2,15] + out.mcmc[i2,c(1,7:14,2:6)]%*%t(mypsc[j2,c(4:5,10:21)])))



    exp_cprob1[i2,j2] <- (exp(mypsc$C1[j2]*(out.mcmc[i2,49]+out.mcmc[i2,50]*ovt_t0[j2])))/(1.0+exp(out.mcmc[i2,49]+out.mcmc[i2,50]*ovt_t0[j2]))

    obs_cprob1[i2,j2] <- (exp(mypsc$C1[j2]*(out.mcmc[i2,48]+out.mcmc[i2,c(33,40:47,34:39)]%*%t(mypsc[j2,c(2, 4:5,10:21)]))))/(1.0+exp(out.mcmc[i2,48] + out.mcmc[i2,c(33,40:47,34:39)]%*%t(mypsc[j2,c(2,4:5,10:21)])))



    exp_prob2[i2,j2] <- ifelse(mypsc$ovt_t0[j2]==1,exp_prob1[i2,j2],exp_prob1[i2,j2]*(exp(mypsc$ovt_t12[j2]*(out.mcmc[i2,32])))/(1.0+exp(out.mcmc[i2,32])))

    obs_prob2[i2,j2] <- ifelse(mypsc$ovt_t0[j2]==1,obs_prob1[i2,j2],obs_prob1[i2,j2]*(exp(mypsc$ovt_t12[j2]*(out.mcmc[i2,31]+out.mcmc[i2,c(17,23:30, 18:22)]%*%t(mypsc[j2,c(4:5,10:12,22:30)]))))/(1.0+exp(out.mcmc[i2,31]+out.mcmc[i2,c(17,23:30,18:22)]%*%t(mypsc[j2,c(4:5,10:12,22:30)]))))



  }



  if (i2 %% 100 == 0) {

    print(i2)

  }

}



#visit-specific Bayesian bootstrap weights;

wmean1_s <- rep(1, ntot)

wmean2_s <- (colSums(exp_prob1)/colSums(obs_prob1))

wmean3_s <- (colSums(exp_prob2)/colSums(obs_prob2))*(colSums(exp_cprob1)/colSums(obs_cprob1))
