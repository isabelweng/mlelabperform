
###User Instruction
library("devtools")
install_github("isabelweng/mlelabperform",auth_token = "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")
library(mlelabperform)

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

####Combined the Sensitivity and Specificity of the 4 Lab Assays
###simulate prevalence bootstrap
##v1
##5/6/2021
simulate_data<-function(
  n_pat=4000,
  p_shc=0.5,
  pcr_obs=0.02,
  ab_spike_obs_pcrpos=0.8,
  ab_ncap_obs_pcrpos=0.8,
  ab_spike_obs_pcrneg=0.01,
  ab_ncap_obs_pcrneg=0.01,
  ab_spike_obs_c=0.7,# conditional probability of having a confirmatory positive spike ab given n-cap is positive
  ab_ncap_obs_c=0.5,# conditional probability of having a confirmatory positive ncap given spike is positive
  ab_neu_obs_con=0.9,# conditional probability of having a confirmatory positive neutraizing ab given n-cap/spike is positive
  seed=1314124){
  set.seed(seed)
  dat_pat<-data.frame(record_id=1:n_pat,
                      shc=rbinom(n_pat,1,p_shc),
                      pcr_pos=rbinom(n_pat,1,pcr_obs)
  )
  dat_pat$ab_spike<-NA
  dat_pat$ab_ncap<-NA
  dat_pat$ab_neu<-NA
  ###original ab tests
  n_shc<-length(dat_pat$ab_spike[dat_pat$shc==1])
  n_shc_pcrpos<-length(dat_pat$ab_spike[dat_pat$shc==1 & dat_pat$pcr_pos==1])
  dat_pat$ab_spike[dat_pat$shc==1 & dat_pat$pcr_pos==1]<-rbinom(n_shc_pcrpos,1, ab_spike_obs_pcrpos)
  dat_pat$ab_spike[dat_pat$shc==1 & dat_pat$pcr_pos==0]<-rbinom((n_shc-n_shc_pcrpos),1, ab_spike_obs_pcrneg)

  n_ucsf<-length(dat_pat$ab_spike[dat_pat$shc==0])
  n_ucsf_pcrpos<-length(dat_pat$ab_spike[dat_pat$shc==0 & dat_pat$pcr_pos==1])
  dat_pat$ab_ncap[dat_pat$shc==0 &  dat_pat$pcr_pos==1]<-rbinom(n_ucsf_pcrpos,1, ab_ncap_obs_pcrpos)
  dat_pat$ab_ncap[dat_pat$shc==0 &  dat_pat$pcr_pos==0]<-rbinom((n_ucsf-n_ucsf_pcrpos),1, ab_ncap_obs_pcrneg)

  ##confirmatory ab tests
  dat_pat$ab_ncap[dat_pat$shc==1]<-rbinom(n_shc,1, ab_ncap_obs_c)
  dat_pat$ab_spike[dat_pat$shc==0]<-rbinom(n_ucsf,1, ab_spike_obs_c)
  ##neutralizing ab tests
  n_initial_pos<-length(  dat_pat$ab_neu[ (dat_pat$ab_spike==1&dat_pat$shc==1 ) | (dat_pat$ab_ncap==1 & dat_pat$shc==0)])
  dat_pat$ab_neu[ (dat_pat$ab_spike==1&dat_pat$shc==1 ) | (dat_pat$ab_ncap==1 & dat_pat$shc==0)]<-rbinom(n_initial_pos,1, ab_neu_obs_con)
  dat_pat$ab_ncap[dat_pat$shc==1 & dat_pat$ab_spike==0]<-NA
  dat_pat$ab_spike[dat_pat$shc==0 & dat_pat$ab_ncap==0]<-NA

  dat_pat$ab_pos<-apply( dat_pat[,c("ab_spike","ab_ncap","ab_neu")],1,sum,na.rm=T)
  dat_pat$ab_pos<-ifelse( dat_pat$ab_pos>=2,1,0)
  dat_pat$pos_status<-ifelse( dat_pat$ab_pos==1|dat_pat$pcr_pos==1,1,0)
  dat_pat
}
dat_simulated<-simulate_data()
table(dat_simulated$shc)
table(dat_simulated$pos_status)
table(dat_simulated$pcr_pos, dat_simulated$ab_pos)
table(dat_simulated$pcr_pos)
table(dat_simulated$ab_pos)
head(dat_simulated)

library(epiR)
pos <- nrow(dat_simulated[dat_simulated$pos_status==1,])
pop <-  nrow(dat_simulated)
dat <- as.matrix(cbind(pos, pop))
epi.conf(dat, ctype = "prevalence", method = "exact", N = 1000,
         design = 1, conf.level = 0.95) * 100
sesp_adj_prevalence(dat=dat_simulated,ab_nen_hyp_obs=0.03,r1=0.99,s1=0.9999,r2=0.9,s2=0.999,r3=0.9,s3=0.999,r4=0.9,s4=0.999)

