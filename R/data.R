#' Right Heart Catheterization Dataset
#'
#' This dataset was used in Connors et al. (1996): The effectiveness of RHC in the initial care of critically
#'  ill patients. J American Medical Association 276:889-897. The dataset pertains to day 1 of hospitalization, i.e., the 'treatment' variable swang1 is whether or not a patient received a RHC (also called the Swan
#'  -Ganz catheter) on the first day in which the patient qualified for the SUPPORT study (see above)
#'
#' @format A data frame with 5735 rows and 63 variables:
#' \describe{
#'   \item{Age}{Age, in years}
#'   \item{Sex}{Sex at birth}
#'   \item{Income}{Yearly Income, in dollars}
#'   \item{meta}{Metabolic Diagnostic}
#'   .
#'   .
#'   .
#' }
#' @source \url{http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.html}
"rhc"


#' Toy dataset for EPBI500 at Case Western Reserve University
#'
#' The Data Set is 100\% fictional. It has three outcomes, with no missing observations anywhere. We  assume that a logical argument suggests that the square of `covA`, as well as the interactions of `covB` with `covC` and with `covD` should be related to treatment assignment, and thus should be  included in our propensity model.
#'
#' @format A data frame with 200 rows and 6 variables. 70 treated - subjects 131-200 and 130 controls - subjects 1-130:
#' \describe{
#'   \item{CovA}{Age, in years}
#'   \item{CovB}{Sex at birth}
#'   \item{CovC}{Yearly Income, in dollars}
#'   \item{CovD}{Metabolic Diagnostic}
#'   .
#'   .
#'   .
#' }
#' @source \url{https://sites.google.com/a/case.edu/love-500/home/data-and-code}
"toy"


#' Simulated dataset to illustrate how \code{ds_function} in the SensitivityR5 package works
#'
#' The data set has one outcome, one exposure and 9 covariates, with no missing observations anywhere.
#'
#' @format A data frame with 100 rows and 11 variables. 300 treated  and 700 controls :
#' \describe{
#'   \item{V1 - V6, V8}{Continious variables}
#'   \item{V7,V9}{Binary variables}
#'   \item{Y2}{Outcome}
#'   \item{exposure}{Treatment indicator}
#'   .
#'   .
#'   .
#' }
"dsn_data"









