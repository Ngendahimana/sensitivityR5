#' Assessing the three rubin rules
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
#' @examples
#' covlist1=c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")
#' rubinRules(data=toy,Treatment="treated",covlist=covlist1)


rubinRules = function(data,Treatment,matchscore="ps",covlist){

  results=list()

  # Rubin 1
  data1a = subset(data, select = c(Treatment,matchscore))
  names(data1a)[c(1:2)] = c("Treatment", "matchscore")
  data1a$Treatment=as.factor(data1a$Treatment)

  results$RUBIN1<- with(data1a, abs(100*(mean(matchscore[Treatment=="1"])-mean(matchscore[Treatment=="0"]))/sd(matchscore)))

  # Rubin 2

  results$RUBIN2 <- with(data1a, var(matchscore[Treatment=="1"])/var(matchscore[Treatment=="0"]))


  # Rubin 3
  data1d = subset(data, select = c(Treatment,matchscore))
  names(data1d)=c("Treatment","matchscore")
  data1f=as.data.frame(cbind(data,data1d))

  data1f$Treatment=as.factor(data1f$Treatment)
  covlist1 = data1f[covlist]
  covlist2 <- as.matrix(covlist1)
  res <- NA
  for(i in 1:ncol(covlist2)) {
    cov <- as.numeric(covlist2[,i])
    num <- var(resid(lm(cov ~ data1f$matchscore))[data1f$Treatment=="1"])
    den <- var(resid(lm(cov ~ data1f$matchscore))[data1f$Treatment=="0"])
    res[i] <- round(num/den, 3)
  }
  names(res) <- names(covlist1)
  #print(res)

  results$RUBIN3=res

  d <- sort(res)
  low <- min(min(d), 0.45)
  high <- max(max(d), 2.05)

  dotchart(d, pch=15, col="black", main="Rubin's Rules plot", xlab="Residual Variance Ratio", xlim=c(low, high))
  abline(v=1, lty=1)
  abline(v=0.8, lty=2, lwd=2, col="blue")
  abline(v=1.25, lty=2, lwd=2, col="blue")
  abline(v=0.5, lty=2, lwd=2, col="red")
  abline(v=2, lty=2, lwd=2, col="red")

  mtext(paste("Rubin Two :",round(results$RUBIN2,2)),side = 3)
  mtext(paste("Rubin One :",round(results$RUBIN1,2)),side = 3,adj=0)

  return(results)

}


#' Amplifying and visualizing gamma parameter
#'
#' @param gamma  gamma parameter to be amplfied.
#' @param lambda  vector of 3 possible n-fold increase in odds of treatment to be amplified.
#' @return An amplification plot.
#' @examples
#' ampPlot(10, c(12,13,30))

ampPlot <- function (gamma, lambda)
{

  stopifnot(length(gamma) == 1)
  stopifnot(gamma > 1)
  stopifnot(min(lambda) > gamma)
  delta <- (gamma * lambda - 1)/(lambda - gamma)
  #plot(lambda,delta,type = "l")
  sensdata=data.frame(lambda,delta)
  ggplot(sensdata, aes(x=lambda, y=delta)) +geom_line()+theme_bw()+annotate("text",x=sensdata[1,1],y=sensdata[1,2],hjust=-.5,label=paste0("(",sensdata[1,1],",",sensdata[1,2],")"))+annotate("text",x=sensdata[2,1],y=sensdata[2,2],hjust=-.5,label=paste0("(",sensdata[2,1],",",sensdata[2,2],")"))+annotate("text",x=sensdata[3,1],y=sensdata[3,2],hjust=1.2,label=paste0("(",sensdata[3,1],",",sensdata[3,2],")"))+geom_point()
}


#' Visualizing sensivity results
#'
#' @param x  Count of pairs where only control has outcome.
#' @param y  Count of pairs where onlytreated has outcome.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' binSensgraph (220, 450)

binSensgraph = function (x, y = NULL, Gamma = 3, GammaInc = .2,alpha = 0.06)
{
  if (length(x) == 1) { # in case your are passing # of discocordant pairs
    ctrl <- x
    trt <- y
  }
  else { # in case you are passing an object from another function
    y.c <- x$mdata$Y[x$mdata$Tr == 0] # select outcome for controls
    y.t <- x$mdata$Y[x$mdata$Tr == 1] #select outcomes for treated
    table(y.t, y.c)
    y.tmp1 <- table(y.t, y.c)[2] # select 1st set of  discocordant pairs i.e treated had outcome, control did not hav outcome
    y.tmp2 <- table(y.t, y.c)[3] # select 2nd set of  discocordant pairs i.e control had outcome, treated did not have outcome
    (if (y.tmp1 >= y.tmp2) {
      trt <- y.tmp1
      ctrl <- y.tmp2
    }                          # code to ensure # treated is always greater than controlls
      else {
        trt <- y.tmp2
        ctrl <- y.tmp1
      })
  }
  gamma <- seq(1, Gamma, by = GammaInc)
  mx <- ctrl + trt
  up <- c() # creating an empty vector to store upper bound
  lo <- c() # creating an empty vector to store lower bound
  series <- seq(trt, mx, by = 1)
  n.it <- length(gamma)
  for (i in 1:n.it) {
    p.plus <- gamma[i]/(1 + gamma[i])
    p.minus <- 1/(1 + gamma[i])
    up.tmp <- sum(dbinom(series, mx, prob = p.plus))
    lo.tmp <- sum(dbinom(series, mx, prob = p.minus))
    up <- c(up, up.tmp)
    lo <- c(lo, lo.tmp)
  }
  pval <- lo[1]
  bounds <- data.frame(gamma, plower = round(lo, 5), pupper = round(up, 5))

  bounds$min = abs(alpha - bounds$pupper)

  vrt = round(bounds[bounds$min == min(bounds$min), ]$gamma,2)
  hrz = round(bounds[bounds$min == min(bounds$min), ]$pupper,2)


 k= ggplot(data = bounds, aes(x = gamma,y = pupper))+geom_line()+geom_point(aes(x=vrt,y=hrz))+ylab("p upper bound")+xlab("gamma (Bias)")+theme_bw()+annotate("text",x=vrt+0.2,y=hrz,label=paste0("(",vrt,",",hrz,")"))+labs(title = "Binary Outcome Sensitivity Plot")

  #plot(bounds$pupper ~ bounds$gamma, type = "l", xlab = "Gamma", ylab = "p-val upper bound", main = "Sensitivity plot for binary outcomes")
 # text(vrt,hrz,paste0("(",vrt,",",hrz,")"),pos = 2)
 # points(vrt,hrz,pch=15)


  colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
  msg <- "Rosenbaum Sensitivity Test \n"
  note <- "Note: Gamma is Odds of Differential Assignment To\n Treatment Due to Unobserved Factors \n"
  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
              msg = msg, bounds = bounds, note = note,plot=k)
  class(Obj) <- c("rbounds", class(Obj))
  Obj

}


#' Creating table One
#'
#' @param x  Count of pairs where only control has outcome.
#' @param y  Count of pairs where onlytreated has outcome.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' edaTable (220, 380)


edaTable = function(baselinevars,outcomeVar,data){
  results = list()
  formula= reformulate(termlabels=baselinevars, response =outcomeVar)
  res1 <- compareGroups(formula, data = dataA, ref = 1)
  table01=createTable(res1,show.p.mul = TRUE,show.p.overall = TRUE)
  try(export2pdf(table01, file='table1.pdf'),silent = T)

}


