#' Sensitivity Analysis for A Simple Comparison for Censored Survival
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

sensSurv = function(data, exp, outcome, failtime,Gamma,alpha,Gammainterval) {
  results = list()

  data1s = subset(data, select = c("matches", exp, outcome, failtime))
  names(data1s)[c(2:4)] = c("exp", "outcome", "failtime")
  data2s = subset(data1s, exp == 0)
  data3s = subset(data1s, exp == 1)
  data4s = subset(data2s, select = c("matches", "exp", "outcome", "failtime"))
  data5s = subset(data3s, select = c("matches", "exp", "outcome", "failtime"))
  data6s = dplyr::full_join(data4s, data5s, by = "matches")
  names(data6s)[c(4, 7)] = c("failtimeNotexp", "failtimeExp")
  data6s$timediff = data6s$failtimeNotexp - data6s$failtimeExp

  wonpairs = sum(data6s$timediff != 0)
  expoutlive = sum(data6s$timediff < 0)

  results$wonpairs = sum(data6s$timediff != 0)
  results$expoutlive = sum(data6s$timediff < 0)

  gamVal = seq(1, Gamma, by = Gammainterval)
  pplus = 1/(1 + gamVal)
  pminus = gamVal/(1 + gamVal)

  table1 = data.frame(cbind(gamVal, pplus, pminus))
  table1$expTplus = wonpairs * table1$pplus
  table1$expTminus = wonpairs * table1$pminus
  table1$sd_expT = sqrt(wonpairs * table1$pplus * (1 - table1$pplus))

  # looping to add p upper bounds on each row

  for (i in 1:length(gamVal)) {
    table1$pupper[i] = round(min(1, 2 * pnorm((expoutlive - table1[i, 5])/table1[i, 6], lower.tail = FALSE)), 4)
  }

  for (i in 1:length(gamVal)) {
    table1$plower[i] = round(min(1, 2 * pnorm((expoutlive - table1[i, 4])/table1[i, 6], lower.tail = FALSE)), 4)
  }




  table1$min = abs(alpha - table1$pupper)

  vrt = table1[table1$min == min(table1$min), ]$gamVal
  hrz = table1[table1$min == min(table1$min), ]$pupper

  plot(table1$pupper ~ table1$gamVal, type = "l", xlab = "Gamma", ylab = "p-val upper bound", main = "Sensitivity plot for survival outcomes")
  segments(x0 = 0, y0 = hrz, x1 = vrt, y1 = hrz, col = "pink", lty = "dashed", lwd = 3)
  segments(x0 = vrt, y0 = 0, x1 = vrt, y1 = hrz, col = "pink", lty = "dashed", lwd = 3)


  results$upperbound_pval = hrz = table1[table1$min == min(table1$min), ]$pupper
  results$Gamma = table1[table1$min == min(table1$min), ]$gamVal

  return(results)



}

