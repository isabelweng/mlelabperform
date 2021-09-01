
###User Instruction
library("devtools")
install_github("isabelweng/mlelabperform",auth_token = "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")
library(mlelabperform)
n_test<-c(6535,2816,2128,1929,1699,1399,1288, 41051,9273,8172,7227,6273,5430,4780)

lofu_function(n_test)

###Simulated dataset

###Test the bootstrapping functions
participant<-rmultinom(1:14,4000,c(0.05,0.02*(6:1),0.45,c(0.03,rep(0.01,5))))
weights<-rnorm(4000,2,0.5)

sample_data<-data.frame()
for (type in 1:nrow(participant)){
  temp_dat<-data.frame(type=rep(type,participant[type,1]))
  sample_data<-rbind(sample_data,temp_dat)
}

sample_data$id<-1:4000
sample_data$weight<-weights

###Point estimate
input_dat<-data.frame(table(sample_data$type))$Freq
point_estimates<-em_2parameters(nk=input_dat,
               pi_init=0.03,p_init=0.01,
               se=0.9, sp=0.96,
               maxit=1000,tol=1e-6,verbose=F)

point_estimates

###MLE confidence interval
mle_ci<-
  mle_ci(pi= point_estimates$pi,p= point_estimates$p,nk=input_dat,se=0.9,sp=0.96)
mle_ci

###Bootstrapped CI
library(plyr)
boot_results<-boot_function(data=sample_data,
                            n_boot=500,
                            weight_yn=F,
                            weight_var="weight",
                            test_type_var="type",
                            record_id_var="id",
                            pi_init=0.03,
                            p_init=0.01,
                            se=0.9,
                            sp=0.96,
                            maxit=1000,
                            tol=1e-6,
                            seed=1412414)

summary(boot_results$pi)
quantile(boot_results$pi,c(0.025,0.5,0.975))
summary(boot_results$p)
quantile(boot_results$p,c(0.025,0.5,0.975))


