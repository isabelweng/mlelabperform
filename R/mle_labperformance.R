

#Lost to follow-up function

lofu_function<-function(n){
  rhat<-sum(n[9:14])/(n[2]+2*n[3]+3*n[4]+4*n[5]+5*n[6]+6*n[7]+6*n[8]+n[10]+2*n[11]+3*n[12]+4*n[13]+5*n[14])
  return(rhat)
}

#likelihood function matrix
likelihood_function<-function(se, sp, pi, p, nk){
  rhat<-lofu_function(nk)
  ##Calculate the probability of the 8 truth scenario
  q1<-pi
  q2<-(1-pi)*p
  q3<-(1-pi)*(1-p)*p
  q4<-(1-pi)*((1-p)^2)*p
  q5<-(1-pi)*((1-p)^3)*p
  q6<-(1-pi)*((1-p)^4)*p
  q7<-(1-pi)*((1-p)^5)*p
  q8<-(1-pi)*((1-p)^6)
  q<-c(q1,q2,q3,q4,q5,q6,q7,q8)
  
  ###Calculate the the matrix of 8*14 
  ##Test scenario 1
  s11<-q1*se
  s12<-q2*(1-sp)
  s13<-q3*(1-sp)
  s14<-q4*(1-sp)
  s15<-q5*(1-sp)
  s16<-q6*(1-sp)
  s17<-q7*(1-sp)
  s18<-q8*(1-sp)
  P1<-sum(c(s11,s12,s13,s14,s15,s16,s17,s18))
  S1<-c(s11,s12,s13,s14,s15,s16,s17,s18)
  
  ##Test scenario 2
  s21<-q1*(1-se)*se
  s22<-q2*sp*se
  s23<-q3*sp*(1-sp)
  s24<-q4*sp*(1-sp)
  s25<-q5*sp*(1-sp)
  s26<-q6*sp*(1-sp)
  s27<-q7*sp*(1-sp)
  s28<-q8*sp*(1-sp)
  P2<-(1-rhat)*sum(c(s21,s22,s23,s24,s25,s26,s27,s28))
  S2<-(1-rhat)* c(s21,s22,s23,s24,s25,s26,s27,s28)

  ##Test scenario 3
  s31<-q1*((1-se)^2)*se
  s32<-q2*sp*(1-se)*se
  s33<-q3*(sp^2)*se
  s34<-q4*(sp^2)*(1-sp)
  s35<-q5*(sp^2)*(1-sp)
  s36<-q6*(sp^2)*(1-sp)
  s37<-q7*(sp^2)*(1-sp)
  s38<-q8*(sp^2)*(1-sp)
  P3<-((1-rhat)^2)*sum(c(s31,s32,s33,s34,s35,s36,s37,s38))
  S3<-((1-rhat)^2)*c(s31,s32,s33,s34,s35,s36,s37,s38)
  
  ##Test scenario 4
  s41<-q1*((1-se)^3)*se
  s42<-q2*sp*((1-se)^2)*se
  s43<-q3*(sp^2)*(1-se)*se
  s44<-q4*(sp^3)*se
  s45<-q5*(sp^3)*(1-sp)
  s46<-q6*(sp^3)*(1-sp)
  s47<-q7*(sp^3)*(1-sp)
  s48<-q8*(sp^3)*(1-sp)
  P4<-((1-rhat)^3)*sum(c(s41,s42,s43,s44,s45,s46,s47,s48))
  S4<-((1-rhat)^3)*c(s41,s42,s43,s44,s45,s46,s47,s48)
  
  ##Test scenario 5
  s51<-q1*((1-se)^4)*se
  s52<-q2*sp*((1-se)^3)*se
  s53<-q3*(sp^2)*((1-se)^2)*se
  s54<-q4*(sp^3)*(1-se)*se
  s55<-q5*(sp^4)*se
  s56<-q6*(sp^4)*(1-sp)
  s57<-q7*(sp^4)*(1-sp)
  s58<-q8*(sp^4)*(1-sp)
  P5<-((1-rhat)^4)*sum(c(s51,s52,s53,s54,s55,s56,s57,s58))
  S5<-((1-rhat)^4)*c(s51,s52,s53,s54,s55,s56,s57,s58)
  
  ##Test scenario 6
  s61<-q1*((1-se)^5)*se
  s62<-q2*sp*((1-se)^4)*se
  s63<-q3*(sp^2)*((1-se)^3)*se
  s64<-q4*(sp^3)*((1-se)^2)*se
  s65<-q5*(sp^4)*(1-se)*se
  s66<-q6*(sp^5)*se
  s67<-q7*(sp^5)*(1-sp)
  s68<-q8*(sp^5)*(1-sp)
  P6<-((1-rhat)^5)*sum(c(s61,s62,s63,s64,s65,s66,s67,s68))
  S6<-((1-rhat)^5)*c(s61,s62,s63,s64,s65,s66,s67,s68)
  
  ##Test scenario 7
  s71<-q1*((1-se)^6)*se
  s72<-q2*sp*((1-se)^5)*se
  s73<-q3*(sp^2)*((1-se)^4)*se
  s74<-q4*(sp^3)*((1-se)^3)*se
  s75<-q5*(sp^4)*((1-se)^2)*se
  s76<-q6*(sp^5)*(1-se)*se
  s77<-q7*(sp^6)*se
  s78<-q8*(sp^6)*(1-sp)
  P7<-((1-rhat)^6)*sum(c(s71,s72,s73,s74,s75,s76,s77,s78))
  S7<-((1-rhat)^6)*c(s71,s72,s73,s74,s75,s76,s77,s78)
  
  ##Test scenario 8
  s81<-q1*((1-se)^7)
  s82<-q2*sp*((1-se)^6)
  s83<-q3*(sp^2)*((1-se)^5)
  s84<-q4*(sp^3)*((1-se)^4)
  s85<-q5*(sp^4)*((1-se)^3)
  s86<-q6*(sp^5)*((1-se)^2)
  s87<-q7*(sp^6)*(1-se)
  s88<-q8*(sp^7)
  P8<-((1-rhat)^6)*sum(c(s81,s82,s83,s84,s85,s86,s87,s88))
  S8<-((1-rhat)^6)*c(s81,s82,s83,s84,s85,s86,s87,s88)
  

  ##Test scenario 9
  s91<-q1*(1-se)
  s92<-q2*sp
  s93<-q3*sp
  s94<-q4*sp
  s95<-q5*sp
  s96<-q6*sp
  s97<-q7*sp
  s98<-q8*sp
  P9<-rhat*sum(c(s91,s92,s93,s94,s95,s96,s97,s98))
  S9<-rhat* c(s91,s92,s93,s94,s95,s96,s97,s98)
  
  
  ##Test scenario 10
  s101<-q1*((1-se)^2)
  s102<-q2*sp*(1-se)
  s103<-q3*(sp^2)
  s104<-q4*(sp^2)
  s105<-q5*(sp^2)
  s106<-q6*(sp^2)
  s107<-q7*(sp^2)
  s108<-q8*(sp^2)
  P10<-(1-rhat)*rhat*sum(c(s101,s102,s103,s104,s105,s106,s107,s108))
  S10<-(1-rhat)*rhat*c(s101,s102,s103,s104,s105,s106,s107,s108)
  
  ##Test scenario 11
  s111<-q1*((1-se)^3)
  s112<-q2*sp*((1-se)^2)
  s113<-q3*(sp^2)*(1-se)
  s114<-q4*(sp^3)
  s115<-q5*(sp^3)
  s116<-q6*(sp^3)
  s117<-q7*(sp^3)
  s118<-q8*(sp^3)
  P11<-((1-rhat)^2)*rhat*sum(c(s111,s112,s113,s114,s115,s116,s117,s118))
  S11<-((1-rhat)^2)*rhat*c(s111,s112,s113,s114,s115,s116,s117,s118)
  
  ##Test scenario 12
  s121<-q1*((1-se)^4)
  s122<-q2*sp*((1-se)^3)
  s123<-q3*(sp^2)*((1-se)^2)
  s124<-q4*(sp^3)*(1-se)
  s125<-q5*(sp^4)
  s126<-q6*(sp^4)
  s127<-q7*(sp^4)
  s128<-q8*(sp^4)
  P12<-((1-rhat)^3)*rhat*sum(c(s121,s122,s123,s124,s125,s126,s127,s128))
  S12<-((1-rhat)^3)*rhat*c(s121,s122,s123,s124,s125,s126,s127,s128)
                          
  ##Test scenario 13
  s131<-q1*((1-se)^5)
  s132<-q2*sp*((1-se)^4)
  s133<-q3*(sp^2)*((1-se)^3)
  s134<-q4*(sp^3)*((1-se)^2)
  s135<-q5*(sp^4)*(1-se)
  s136<-q6*(sp^5)
  s137<-q7*(sp^5)
  s138<-q8*(sp^5)
  P13<-((1-rhat)^4)*rhat*sum(c(s131,s132,s133,s134,s135,s136,s137,s138))
  S13<-((1-rhat)^4)*rhat*c(s131,s132,s133,s134,s135,s136,s137,s138)
  
  
  ##Test scenario 14
  s141<-q1*((1-se)^6)
  s142<-q2*sp*((1-se)^5)
  s143<-q3*(sp^2)*((1-se)^4)
  s144<-q4*(sp^3)*((1-se)^3)
  s145<-q5*(sp^4)*((1-se)^2)
  s146<-q6*(sp^5)*(1-se)
  s147<-q7*(sp^6)
  s148<-q8*(sp^6)
  P14<-((1-rhat)^5)*rhat*sum(c(s141,s142,s143,s144,s145,s146,s147,s148))
  S14<-((1-rhat)^5)*rhat*c(s141,s142,s143,s144,s145,s146,s147,s148)
  
  
  Pcap<-c(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14)
  Scap<-data.frame(s1=S1,
                s2=S2,
                s3=S3,
                s4=S4,
                s5=S5,
                s6=S6,
                s7=S7,
                s8=S8,
                s9=S9,
                s10=S10,
                s11=S11,
                s12=S12,
                s13=S13,
                s14=S14)
  return(list(Pcap,q,Scap))
  
}


###EM function
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
    log_likelihood_cur<-m[1]*log(pi_cur)+sum(m[2:8])*log(1-pi_cur)+sum(m[2:7])*log(p_cur)+(m[3]+2*m[4]+3*m[5]+4*m[6]+5*m[7]+6*m[8])*log(1-p_cur)
    
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
    
    if( log_likelihood_change < tol ){ 
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


###one parameter

em_1parameter<-function(nk=c(100, rep(200,7),rep(150,6)),
                         pi=0.01,p_init=0.01,
                         se=0.99, sp=0.999,
                         maxit=1000,tol=1e-6,verbose=F
){
  # Estimation of parameter(Initial)
  flag <- 0
  p_cur<-p_init
  log_likelihood_prev<-0
  rhat<-lofu_function(nk)
  for(i in 1:maxit){
    ###calculate the expectation
    s_matrix<-likelihood_function(se=se, sp=sp, pi=pi , p=p_cur,
                                  nk=nk)[[3]]
    epiloson_matrix<-as.data.frame(matrix(NA,
                                          ncol=ncol(s_matrix),
                                          nrow=nrow(s_matrix)))
    for (j in 1:ncol(s_matrix)){
      epilson_vector<- nk[j]*(s_matrix[,j]/sum(s_matrix[,j]))
      
      epiloson_matrix[,j]<-epilson_vector
    }
    
    m<-apply(epiloson_matrix,1,sum,na.rm=T)
    
    ###Update pi and p to maximization
    p_cur<-sum(m[2:7])/(m[2]+2*m[3]+3*m[4]+4*m[5]+5*m[6]+6*m[7]+6*m[8])
    
    ###Calculate the likelihood
    log_likelihood_cur<-m[1]*log(pi)+sum(m[2:8])*log(1-pi)+sum(m[2:7])*log(p_cur)+(m[3]+2*m[4]+3*m[5]+4*m[6]+5*m[7]+6*m[8])*log(1-p_cur)
    
    log_likelihood_change<-abs(log_likelihood_cur-log_likelihood_prev)
    ###change to previous
    log_likelihood_prev<-log_likelihood_cur
    if (verbose==T){
    if (i %% 5==0){
      print(paste0(i, " of literations:"))
      print(paste0("Current likelihood:",log_likelihood_cur,sep=""))
      print(paste0("Current p:",p_cur,sep=""))
    }
    }
    
    if( log_likelihood_change < tol ){ 
      flag <- 1
      print(paste0("After ",i," literations: Model Converged."))
      return(data.frame(literation=i,
                        log_likelihodod=log_likelihood_cur,
                        p=p_cur))
      break
      
    }
    
  }
  
}




###Bootstrap the CI

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

###mle of the variance estimate

xmatrix_function<-function(se, sp, nk){
  rhat<-lofu_function(nk)
 
  r.vec=c((1-rhat)^c((0:6), 6), rhat*(1-rhat)^(0:5))
  s.mat=matrix(NA, 14, 8)
  
  
  s.mat[,1]=c(se,    (1-se)*se,   (1-se)^2*se,   (1-se)^3*se,     (1-se)^4*se,       (1-se)^5*se,        (1-se)^6*se,       (1-se)^7,                   sp^0*(1-se)^(1:6))
  
  s.mat[,2]=c(1-sp,  sp*se,       sp*(1-se)*se, sp*(1-se)^2*se, sp*(1-se)^3*se,   sp*(1-se)^4*se,    sp*(1-se)^5*se,    sp*(1-se)^6,              sp^1*(1-se)^(0:5))
  
  s.mat[,3]=c(1-sp,  sp*(1-sp),   sp^2*se,       sp^2*(1-se)*se, sp^2*(1-se)^2*se, sp^2*(1-se)^3*se,  sp^2*(1-se)^4*se,  sp^2*(1-se)^5, sp,       sp^2*(1-se)^(0:4))
  
  s.mat[,4]=c(1-sp,  sp*(1-sp),   sp^2*(1-sp),   sp^3*se,         sp^3*(1-se)*se,   sp^3*(1-se)^2*se,  sp^3*(1-se)^3*se,  sp^3*(1-se)^4, sp^(1:2), sp^3*(1-se)^(0:3))
  
  s.mat[,5]=c(1-sp,  sp*(1-sp)  , sp^2*(1-sp),   sp^3*(1-sp),     sp^4*se,           sp^4*(1-se)*se,    sp^4*(1-se)^2*se,  sp^4*(1-se)^3, sp^(1:3), sp^4*(1-se)^(0:2))
  
  s.mat[,6]=c(1-sp,  sp*(1-sp)  , sp^2*(1-sp),   sp^3*(1-sp),     sp^4*(1-sp),       sp^5*se,            sp^5*(1-se)*se,    sp^5*(1-se)^2, sp^(1:4), sp^5*(1-se)^(0:1))
  
  s.mat[,7]=c(1-sp,  sp*(1-sp)  , sp^2*(1-sp),   sp^3*(1-sp),     sp^4*(1-sp),       sp^5*(1-sp),        sp^6*se,            sp^6*(1-se),   sp^(1:5), sp^6)
  
  s.mat[,8]=c(1-sp,  sp*(1-sp)  , sp^2*(1-sp),   sp^3*(1-sp),     sp^4*(1-sp),       sp^5*(1-sp),        sp^6*(1-sp),        sp^7,           sp^(1:5), sp^6)
  
  s.mat=t(s.mat*r.vec)
  return(s.mat)
  
  
}


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





