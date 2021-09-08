#' Combine The Performance on Population
#'
#' This function calculate the S-matrix of the likelihood.
#' @param r sensitivity/PPAs of combined Case definition in subpopulation
#' @param s specificity/NPA of combined Case definition in ith subpopulation
#' @param prob probability of each subpopulation in the entire cohort
#' @return Calculate the Combined Performance of Case Definition in the entire populatioj.
#' @export
#' @examples
#' combine_perform_subcohort(r=c(0.95,0.96,0.97,0.95,0.95,0.96),s=c(0.999,0.997,0.994,0.998,0.996,0.997),prob=c(0.1,0.2,0.3,0.2,0.1,0.1))


combine_perform_subcohort<-function(r, s,prob){
  r_combined<-sum(r*prob)
  s_combined<-sum(s*prob)
  data.frame(sen_final=r_combined, sp_final=s_combined)
}
