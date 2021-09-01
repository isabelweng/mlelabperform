#' Bootstrapping for CIs Prevalence and Incidence
#'
#' This function calculates Confidence Intervals of MLE-EM estimates of Prevalence and Incidence
#' @param se sensitivity/PPA
#' @param sp specificity/NPA
#' @param pi_init prevalence to initialization
#' @param p_init incidence to initialization
#' @param data input data
#' @param seed starting seed
#' @param weight_yn whether used weighted data
#' @param weight_var variable for weights
#' @param test_type_var variable for scenario indicator
#' @param record_id_var patient id variable
#' @param n_boot # of bootstrapped samples
#' @param maxit maximum number iteriteration allowed
#' @param tol convergence criteria
#' @param verbose show intermediate outputs
#' @return Bootstrapped Confidence intervals.
#' @export
#' @examples
#' boot_function(data=sample_data,n_boot=10000,weight_yn=F,weight_var="weight",test_type_var="type",record_id_var="id",pi_init=0.03,p_init=0.01,se=0.99, sp=0.999,maxit=1000,tol=1e-6,seed=1412414)

library(plyr)
boot_function<-function(data=sample_data,
                        ##sample data contains participant level data and their categories
                        n_boot=10000,
                        weight_yn=F,
                        weight_var="weight",
                        test_type_var="type",
                        record_id_var="id",
                        pi_init=0.03,
                        p_init=0.01,
                        se=0.99,
                        sp=0.999,
                        maxit=1000,
                        tol=1e-6,
                        seed=1412414

){
  set.seed(seed)
  final_boot_data<-data.frame()
  for (i in 1:n_boot){
    ##sample patients with replacement
    boot_sample_id<-sample(data[[record_id_var]],nrow(data),replace=T)
    boot_sample_id<-data.frame(boot_sample_id)
    boot_sample<-merge(boot_sample_id,data,by.x="boot_sample_id",by.y=record_id_var,all.x=T)
    colnames(boot_sample)[colnames(boot_sample)==weight_var]<-"weight"
    colnames(boot_sample)[colnames(boot_sample)==test_type_var]<-"type"

    if (weight_yn==F){
      input_dat<-data.frame(table(boot_sample$type))$Freq
    }else{
      input_dat<-as.numeric(ddply(boot_sample,~type,summarise,
                                  n_w=sum(weight,na.rm=T))$n_w)
    }

    boot_results<-em_2parameters(nk=input_dat,
                                 pi_init=pi_init,p_init=p_init,
                                 se=se, sp=sp,
                                 maxit=maxit,tol=tol,verbose=F)
    final_boot_data<-rbind(final_boot_data,boot_results)

  }
  return(final_boot_data)
}

