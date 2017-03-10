#' Right Heart Catheterization Dataset
#'
#' This dataset was used in Connors et al. (1996): The effectiveness of RHC in the initial care of critically
#'  ill patients. J American Medical Association 276:889-897. The dataset pertains to day 1 of hospitalizatio#' n, i.e., the 'treatment' variable swang1 is whether or not a patient received a RHC (also called the Swan
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


#' Simulated matched sample for class EPBI 500 class at Case Western Reserve University
#'
#' The Data Set is 100% fictional, and is available as toy.csv on the course website. It #' 70 treated - subjects 131-200 and 130 controls -     #' subjects 1-130  on treatment status. It has three outcomes, with no    #' missing observations anywhere. We assume that a logical argument suggests that the   #' square of covA, as well as the interactions of covB with covC and with covD should be #' related to treatment assignment, and thus should be included in our propensity model.
#'
#' @format A data frame with 200 rows and 6 variables:
#' @source \url{https://sites.google.com/a/case.edu/love-500/home/data-and-code}
"toy"



