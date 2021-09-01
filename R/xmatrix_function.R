#' X-matrix
#'
#' This function calculates x-matrix of the likelihood function
#' @param se sensitivity/PPA
#' @param sp specificity/NPA
#' @param nk input data vector
#' @return x-matrix
#' @export
#' @examples
#' xmatrix_function(0.99,0.9,nk=c(6535,2816,2128,1929,1699,1399,1288, 41051,9273,8172,7227,6273,5430,4780))


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

