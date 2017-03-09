#' Ploting results of sensitivity analyisis
#'
#' This function allows you to assess how sensitive your results are to unmeasured variable.
#' @param data A matched sample
#' @param exp A variables defining exposure group
#' @param outcome The outcome variable
#' @param failtime Time to event
#' @param Gamma Bias to be assessed
#' @param Gammainterval interval between two hidden bias to be assessed
#' @param alpha Significance level
#' @keywords Sensitivity
#' @export
#' @references Section 4.4.8. of Rosenbaum PR (2002) Observational Studies, 2nd Edition.

add <- function(x,y){
  x+y
}
