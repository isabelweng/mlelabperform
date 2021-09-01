#' MLE CIs Prevalence and Incidence
#'
#' This function calculates Confidence Intervals of MLE-EM estimates of Prevalence and Incidence based on MLE
#' @param se sensitivity/PPA
#' @param sp specificity/NPA
#' @param nK input data vector
#' @param pi adjusted prevalence
#' @param p adjusted incidence
#' @return MLE confidence interval.
#' @export
#' @examples
#' mle_ci(pi=0.07,p=0.02,nk=,se=0.99,sp=0.999)


mle_ci<-function(pi=0.07,p=0.02,nk,se,sp){
  dpi1<-c(1,-p*((1-p)^c(0:5)),-(1-p)^6)
  dp1<-c(0,-pi, (1-pi)*(1-2*p),(1-pi)*(1-p)*(1-3*p),
         (1-pi)*((1-p)^2)*(1-4*p),
         (1-pi)*((1-p)^3)*(1-5*p),
         (1-pi)*((1-p)^4)*(1-6*p),
         -6*(1-pi)*((1-p)^5))
  dpi2<-rep(0,8)
  dp2<-c(0,0, -2*(1-pi),-2*(1-pi)*(2-3*p),
         -6*(1-pi)*(1-2*p),-4*(1-pi)*((1-p)^2)*(2-5*p),
         -10*(1-pi)*((1-p)^3)*(1-3*p),
         30*(1-pi)*((1-p)^4))
  dpip2<-c(0,-1,1-2*p, (1-p)*(3*p-1),((1-p)^2)*(4*p-1),
           ((1-p)^3)*(5*p-1),((1-p)^4)*(6*p-1),
           6*((1-p)^5))
  x_matrix<-xmatrix_function(se=se, sp=sp,nk=nk)
  q.vec=c(pi, (1-pi)*p, (1-pi)*(1-p)^(1:5)*p, (1-pi)*(1-p)^6)

  ##derivatives of the pk
  dpkpi1<-apply(t(x_matrix*dpi1),1,sum,na.rm=T)
  dpkp1<-apply(t(x_matrix*dp1),1,sum,na.rm=T)
  dpkpi2<-apply(t(x_matrix*dpi2),1,sum,na.rm=T)
  dpkp2<-apply(t(x_matrix*dp2),1,sum,na.rm=T)
  dpkpip2<-apply(t(x_matrix*dpip2),1,sum,na.rm=T)

  ##pk
  pk<-apply(x_matrix*q.vec,2,sum,na.rm=T)

  ##derivatives of the likelihood function
  dLpi2<-sum(nk*(dpkpi2*(1/pk)-(dpkpi1^2)*(1/(pk^2))))
  dLp2<-sum(nk*(dpkp2*(1/pk)-(dpkp1^2)*(1/(pk^2))))
  dLpip2<-sum(nk*(dpkpip2*(1/pk)-(dpkp1)*(dpkpi1)*(1/(pk^2))))

  ##Fisher's information matrix
  I<-matrix(c(-dLpi2,-dLpip2,-dLpip2,-dLp2),nrow=2)
  I_inverse<-solve(I)
  varpi<-I_inverse[1,1]
  varp<-I_inverse[2,2]
  output<-data.frame(parameter=c("prevalence","incidence"),
                     estimate=c(pi,p),
                     lower=c(pi-1.96*sqrt(varpi),
                             p-1.96*sqrt(varp)),
                     upper=c(pi+1.96*sqrt(varpi),
                             p+1.96*sqrt(varp))
  )

  return(output)

}





