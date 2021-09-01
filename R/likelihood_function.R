#' S matrix
#'
#' This function calculate the S-matrix of the likelihood.
#' @param se sensitivity/PPA
#' @param sp sepecificity/NPA
#' @param pi prevalence
#' @param p incidence
#' @param nk input data
#' @return Calculate the S-matrix of the likelihood.
#' @export
#' @examples
#' likelihood_function(se=0.99,sp=0.9,pi=0.05,p=0.03,nk=c(6535,2816,2128,1929,1699,1399,1288, 41051,9273,8172,7227,6273,5430,4780))

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

