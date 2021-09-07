#' Combine The Performance From CZI Assays (1 PCR and 3 Ab Assays)
#'
#' This function calculate the S-matrix of the likelihood.
#' @param r1 sensitivity/PPA of PCR
#' @param s1 specificity/NPA of PCR
#' @param r2 sensitivity/PPA of Ab Spike (Stanford Ab Assay)
#' @param s2 specificity/NPA of Ab Spike (Stanford Ab Assay)
#' @param r3 sensitivity/PPA of Ab N-cap (UCSF Ab Assay)
#' @param s3 specificity/NPA of Ab N-cap (UCSF Ab Assay)
#' @param r4 sensitivity/PPA of Ab Neutralizing
#' @param s4 specificity/NPA of Ab Neutralizing
#' @param ab_nen_hyp_obs Hypothesized Ab Neuralizing Observed Prevalence
#' @param dat simulated baseline data on the patient level
#' @return Calculate the Combined Performance of Case Definition.
#' @export
#' @examples
#' function(dat=dat_simulated,ab_nen_hyp_obs=0.03,r1=0.99,s1=0.9999,r2=0.9,s2=0.999,r3=0.9,s3=0.999,r4=0.9,s4=0.999)


sesp_adj_prevalence<-function(dat=dat_simulated,
                              ab_nen_hyp_obs=0.03,
                              ##PCR, sensitivity specificity
                              r1=0.99,
                              s1=0.9999,
                              ##ab spike,
                              r2=0.9,
                              s2=0.999,
                              ##ab ncap,
                              r3=0.9,
                              s3=0.999,
                              ##ab neutralizing,
                              r4=0.9,
                              s4=0.999
){

  s<-nrow(dat)
  tao1<-mean(dat$pcr_pos,na.rm=T)
  beta<-mean(dat$shc[dat$pcr_pos==0])
  tao2<-mean(dat$ab_spike[dat$shc==1 & dat$pcr_pos==0],na.rm=T)
  tao3<-mean(dat$ab_ncap[dat$shc==1 &dat$pcr_pos==0& dat$ab_spike==1],na.rm=T)
  tao4<-mean(dat$ab_neu[dat$shc==1 & dat$pcr_pos==0& dat$ab_spike==1 & dat$ab_ncap==0],na.rm=T)
  tao5<-mean(dat$ab_ncap[dat$shc==0& dat$pcr_pos==0],na.rm=T)
  tao6<-mean(dat$ab_spike[dat$shc==0& dat$pcr_pos==0 & dat$ab_ncap==1],na.rm=T)
  tao7<-mean(dat$ab_neu[dat$shc==0& dat$pcr_pos==0 & dat$ab_ncap==1 & dat$ab_spike==0],na.rm=T)
  q1<-tao1
  q2<-tao2
  q3<-tao5
  q4<-ab_nen_hyp_obs
  p1<-(q1+s1-1)/(r1+s1-1)
  p2<-(q2+s2-1)/(r2+s2-1)
  p3<-(q3+s3-1)/(r3+s3-1)
  p4<-(q4+s4-1)/(r4+s4-1)
  x1<-r1*p1/(r1*p1+(1-s1)*(1-p1))
  x2<-r2*p2/(r2*p2+(1-s2)*(1-p2))
  x3<-r3*p3/(r3*p3+(1-s3)*(1-p3))
  x4<-r4*p4/(r4*p4+(1-s4)*(1-p4))
  y1<-(s1*(1-p1))/((1-r1)*p1+s1*(1-p1))
  y2<-(s2*(1-p2))/((1-r2)*p2+s2*(1-p2))
  y3<-(s3*(1-p3))/((1-r3)*p3+s3*(1-p3))
  y4<-(s4*(1-p4))/((1-r4)*p4+s4*(1-p4))

  phi1<-s*tao1
  alpha1<-s*tao1*x1
  rau1<-phi1-alpha1

  phi2<-s*(1-tao1)*beta*tao2*tao3
  alpha2<-s*(1-tao1)*beta*tao2*tao3*y1*x2*x3
  rau2<-phi2-alpha2

  phi3<-s*(1-tao1)*beta*(1-tao2)
  alpha3<-s*(1-tao1)*beta*(1-tao2)*y1*y2
  rau3<-phi3-alpha3

  phi4<-s*(1-tao1)*(1-beta)*(1-tao5)
  alpha4<-s*(1-tao1)*(1-beta)*(1-tao5)*y1*y3
  rau4<-phi4-alpha4

  phi5<-s*(1-tao1)*beta*tao2*(1-tao3)*tao4
  alpha5<-s*(1-tao1)*beta*tao2*(1-tao3)*tao4*y1*x2*x3*x4
  rau5<-phi5-alpha5


  phi6<-s*(1-tao1)*beta*tao2*(1-tao3)*(1-tao4)
  alpha6<-s*(1-tao1)*beta*tao2*(1-tao3)*(1-tao4)*y1*x2*y3*y4
  rau6<-phi6-alpha6

  phi7<-s*(1-tao1)*(1-beta)*tao5*tao6
  alpha7<-s*(1-tao1)*(1-beta)*tao5*tao6*y1*x3*x2
  rau7<-phi7-alpha7

  phi8<-s*(1-tao1)*(1-beta)*tao5*(1-tao6)*tao7
  alpha8<-s*(1-tao1)*(1-beta)*tao5*(1-tao6)*tao7*y1*x3*x2*x4
  rau8<-phi8-alpha8


  phi9<-s*(1-tao1)*(1-beta)*tao5*(1-tao6)*(1-tao7)
  alpha9<-s*(1-tao1)*(1-beta)*tao5*(1-tao6)*(1-tao7)*y1*x3*y2*y4
  rau9<-phi9-alpha9

  a<-alpha1+alpha2+alpha5+alpha7+alpha8
  b<-rau1+rau2+rau5+rau7+rau8
  c<-rau3+rau4+rau6+rau9
  d<-alpha3+alpha4+alpha6+alpha9
  aPr<-(a+c)/s
  sen<-a/(a+c)
  sp<-d/(b+d)
  Pr<-mean(dat$pos_status)

  output<-data.frame(phi1=phi1,
                     rau1=rau1,
                     alpha1=alpha1,
                     phi2=phi2,
                     rau2=rau2,
                     alpha2=alpha2,
                     phi3=phi3,
                     rau3=rau3,
                     alpha3=alpha3,
                     phi4=phi4,
                     rau4=rau4,
                     alpha4=alpha4,
                     phi5=phi5,
                     rau5=rau5,
                     alpha5=alpha5,
                     phi6=phi6,
                     rau6=rau6,
                     alpha6=alpha6,
                     phi7=phi7,
                     rau7=rau7,
                     alpha7=alpha7,
                     phi8=phi8,
                     rau8=rau8,
                     alpha8=alpha8,
                     phi9=phi9,
                     rau9=rau9,
                     alpha9=alpha9,
                     s=s,
                     Pr=Pr,
                     sen=sen,
                     sp=sp,
                     aPr=aPr)
  output
}
