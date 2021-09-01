#' Lost to Follow-up Function
#'
#' This function calculates risk of LOFU.
#' @param n input data
#' @return risk of lost to follow-up based on MLE approach.
#' @export
#' @examples
#' lofu_function(c(6535,2816,2128,1929,1699,1399,1288, 41051,9273,8172,7227,6273,5430,4780))

lofu_function<-function(n){
  rhat<-sum(n[9:14])/(n[2]+2*n[3]+3*n[4]+4*n[5]+5*n[6]+6*n[7]+6*n[8]+n[10]+2*n[11]+3*n[12]+4*n[13]+5*n[14])
  return(rhat)
}
