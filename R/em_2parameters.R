#' MLE-EM estimates of Prevalence and Incidence
#'
#' This function calculates MLE-EM estimates of Prevalence and Incidence
#' @param se sensitivity/PPA
#' @param sp specificity/NPA
#' @param pi_init prevalence to initialization
#' @param p_init incidence to initialization
#' @param nk input data
#' @param maxit maximum number iteriteration allowed
#' @param tol convergence criteria
#' @param verbose show intermediate outputs
#' @return MLE-EM estimates of Prevalence and Incidence.
#' @export
#' @examples
#' em_2parameters(nk=c(100, rep(200,7),rep(150,6)),pi_init=0.03,p_init=0.01,se=0.99, sp=0.999,maxit=1000,tol=1e-6,verbose=F)

em_2parameters<-function(nk=c(100, rep(200,7),rep(150,6)),
                         pi_init=0.03,p_init=0.01,
                         se=0.99, sp=0.999,
                         maxit=1000,tol=1e-6,verbose=F
){
  # Estimation of parameter(Initial)
  flag <- 0
  pi_cur <- pi_init
  p_cur<-p_init
  log_likelihood_prev<-0
  for(i in 1:maxit){
    ###calculate the expectation
    s_matrix<-likelihood_function(se=se, sp=sp, pi=pi_cur , p=p_cur,
                                  nk=nk)[[3]]


    ##per truth and in our case nrow(s_matrix)==8

    epiloson_matrix<-as.data.frame(matrix(NA,
                                          ncol=ncol(s_matrix),
                                          nrow=nrow(s_matrix)))
    for (j in 1:ncol(s_matrix)){
      epilson_vector<- nk[j]*(s_matrix[,j]/sum(s_matrix[,j]))

      epiloson_matrix[,j]<-epilson_vector
    }

    m<-apply(epiloson_matrix,1,sum,na.rm=T)

    ###Update pi and p to maximization
    pi_cur<-m[1]/sum(m)
    p_cur<-sum(m[2:7])/(m[2]+2*m[3]+3*m[4]+4*m[5]+5*m[6]+6*m[7]+6*m[8])

    ###Calculate the likelihood
    if (1>pi_cur>0 & 1>p_cur>0  ){
    log_likelihood_cur<-m[1]*log(pi_cur)+sum(m[2:8])*log(1-pi_cur)+sum(m[2:7])*log(p_cur)+(m[3]+2*m[4]+3*m[5]+4*m[6]+5*m[7]+6*m[8])*log(1-p_cur)}

    if ((pi_cur==0|pi_cur==1)) & 1>p_cur>0   ){
      log_likelihood_cur<-sum(m[2:7])*log(p_cur)+(m[3]+2*m[4]+3*m[5]+4*m[6]+5*m[7]+6*m[8])*log(1-p_cur)}

    if ((p_cur==0|p_cur==1)) & 1>pi_cur>0   ){
      log_likelihood_cur<-m[1]*log(pi_cur)+sum(m[2:8])*log(1-pi_cur)}

   if ((p_cur==0|p_cur==1)) & (pi_cur==0|pi_cur==1)) ){
     flag <- 1
     print(paste0("After ",i," literations: Unable to estimate likelihood function with 0/1 prevalence and 0/1 incidence."))
     return(data.frame(literation=i,
                       log_likelihodod=NA,
                       pi=pi_cur,
                       p=p_cur))
     break
   }

    log_likelihood_change<-abs(log_likelihood_cur-log_likelihood_prev)
    ###change to previous
    log_likelihood_prev<-log_likelihood_cur
    if (verbose==T){
      if (i %% 5==0){
        print(paste0(i, " of literations:"))
        print(paste0("Current likelihood:",log_likelihood_cur,sep=""))
        print(paste0("Current pi:",pi_cur,sep=""))
        print(paste0("Current p:",p_cur,sep=""))
      }
    }

    if( log_likelihood_change < tol  ){
      flag <- 1
      print(paste0("After ",i," literations: Model Converged."))
      return(data.frame(literation=i,
                        log_likelihodod=log_likelihood_cur,
                        pi=pi_cur,
                        p=p_cur))
      break

    }

  }

}

