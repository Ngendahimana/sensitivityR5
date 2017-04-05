#' Assessing the three rubin rules
#'
#' This function allows you to assess how sensitive your results are to unmeasured variable.
#' @param data data set to be used.
#' @param Treatment A variables defining exposure group.
#' @param matchscore Variable containing matching distance. Default is propensity score
#' @param covlist list of variables to be balanced. Note: All variable should be of numeric type.
#' @export
#' @examples
#' data(toy)
#' psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
#'
#' toy$ps <- psmodel$fitted
#' toy$linps <- psmodel$linear.predictors
#'
#' covlist1=c('covA', 'covB', 'covC', 'covD', 'covE', 'covF.Middle', 'covF.High', 'Asqr','BC', 'BD')
#'
#' rubinRules(data=toy,Treatment='treated',covlist=covlist1)


rubinRules = function(data, Treatment, matchscore = "ps", covlist) {

    results = list()

    # Rubin 1
    data1a = subset(data, select = c(Treatment, matchscore))
    names(data1a)[c(1:2)] = c("Treatment", "matchscore")
    data1a$Treatment = as.factor(data1a$Treatment)

    results$RUBIN1 <- with(data1a, abs(100 * (mean(matchscore[Treatment == "1"]) - mean(matchscore[Treatment == "0"]))/sd(matchscore)))

    # Rubin 2

    results$RUBIN2 <- with(data1a, var(matchscore[Treatment == "1"])/var(matchscore[Treatment == "0"]))


    # Rubin 3
    data1d = subset(data, select = c(Treatment, matchscore))
    names(data1d) = c("Treatment", "matchscore")
    data1f = as.data.frame(cbind(data, data1d))

    data1f$Treatment = as.factor(data1f$Treatment)
    covlist1 = data1f[covlist]
    covlist2 <- as.matrix(covlist1)
    res <- NA
    for (i in 1:ncol(covlist2)) {
        cov <- as.numeric(covlist2[, i])
        num <- var(resid(lm(cov ~ data1f$matchscore))[data1f$Treatment == "1"])
        den <- var(resid(lm(cov ~ data1f$matchscore))[data1f$Treatment == "0"])
        res[i] <- round(num/den, 3)
    }
    names(res) <- names(covlist1)
    # print(res)

    results$RUBIN3 = res

    d <- sort(res)
    low <- min(min(d), 0.45)
    high <- max(max(d), 2.05)

     dotchart(d, pch = 15, col = "black", main = "Rubin's Rules plot", xlab = "Residual Variance Ratio", xlim = c(low, high))
    abline(v = 1, lty = 1)
    abline(v = 0.8, lty = 2, lwd = 2, col = "blue")
    abline(v = 1.25, lty = 2, lwd = 2, col = "blue")
    abline(v = 0.5, lty = 2, lwd = 2, col = "red")
    abline(v = 2, lty = 2, lwd = 2, col = "red")

    mtext(paste("Rubin Two :", round(results$RUBIN2, 2)), side = 3)
    mtext(paste("Rubin One :", round(results$RUBIN1, 2)), side = 3, adj = 0)

    results$plot = plot
    return(results)

}


#' Amplifying and visualizing gamma parameter
#'
#' @param gamma  gamma parameter to be amplfied.
#' @param lambda  vector of 3 possible n-fold increase in odds of treatment to be amplified.
#' @export
#' @return An amplification plot.
#' @examples
#' ampPlot(10, c(12,13,30))

ampPlot <- function(gamma, lambda) {

    stopifnot(length(gamma) == 1)
    stopifnot(gamma > 1)
    stopifnot(min(lambda) > gamma)
    delta <- (gamma * lambda - 1)/(lambda - gamma)
    # plot(lambda,delta,type = 'l')
    sensdata = data.frame(lambda, delta)
    ggplot(sensdata, aes(x = lambda, y = delta)) + geom_line() + theme_bw() + annotate("text", x = sensdata[1, 1], y = sensdata[1,
        2], hjust = -0.5, label = paste0("(", sensdata[1, 1], ",", sensdata[1, 2], ")")) + annotate("text", x = sensdata[2,
        1], y = sensdata[2, 2], hjust = -0.5, label = paste0("(", sensdata[2, 1], ",", sensdata[2, 2], ")")) + annotate("text",
        x = sensdata[3, 1], y = sensdata[3, 2], hjust = 1.2, label = paste0("(", sensdata[3, 1], ",", sensdata[3, 2], ")")) +
        geom_point()
}



#' Visualizing sensivity results
#'
#' @param data  Dataset before matching.
#' @param x  name of the \code{match},\code{matchit} and \code{bimatch} object.
#' @param exposure name of the exposure variable. This name must be in quotation marks.
#' @param outcome name of the outcome variable. This name must be in quotation marks.
#' @param Gamma Sensitivity parameter.
#' @param GammaInc Increamental value to be used in generating Sensitivity table
#' @param alpha Significant value. Default is 0.05
#' @export
#' @examples
#' library(MatchIt);library(Matching);library(ggplot2)
#' data(rhc)
#' rhc1 = rhc[,c('swang1','age', 'female', 'edu', 'income', 'ninsclas', 'race', 'cat1', 'dnr1', 'wtkilo1','hrt1', 'meanbp1', 'resp1', 'temp1', 'card', 'gastr', 'hema', 'meta', 'neuro', 'ortho','renal', 'resp', 'seps', 'trauma' ,'amihx', 'ca','cardiohx' ,'chfhx', 'chrpulhx','dementhx' ,'gibledhx','immunhx', 'liverhx', 'malighx', 'psychhx','renalhx', 'transhx','aps1', 'das2d3pc','scoma1', 'surv2md1','alb1', 'bili1', 'crea1', 'hema1', 'paco21', 'pafi1', 'ph1','pot1', 'sod1','urin1.NA', 'urin1.i', 'wblc1','surv_30')]
#'
#' # fitting score model
#'
#' psmodel <- glm(swang1 ~ age + female + edu + income + ninsclas + race + cat1 + dnr1 + wtkilo1 + hrt1 + meanbp1 + resp1 + temp1 + card + gastr + hema + meta + neuro + ortho + renal + resp + seps + trauma +amihx + ca + cardiohx + chfhx + chrpulhx + dementhx + gibledhx + immunhx + liverhx + malighx + psychhx + renalhx + transhx + aps1 + das2d3pc + scoma1 + surv2md1 + alb1 + bili1 + crea1 + hema1 + paco21 +pafi1 + ph1 + pot1 + sod1 + urin1.NA + urin1.i + wblc1, family=binomial(), data=rhc1)
# rocplot(mod1)
#' rhc1$ps <- psmodel$fitted
#' rhc1$linps <- psmodel$linear.predictors
#'
#' # Creating Match object
#'
#' X <- rhc1$linps
#' Tr <- as.logical(rhc1$swang1)
#' match1 <- Match(Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
#'
#' # Creating Matchit object.
#'
#' m.out1 <- matchit(swang1 ~ age + female + edu + income + ninsclas + race + cat1 + dnr1 + wtkilo1 + hrt1 + meanbp1 + resp1 + temp1 + card + gastr + hema + meta + neuro + ortho + renal + resp + seps + trauma +amihx + ca + cardiohx + chfhx + chrpulhx + dementhx + gibledhx + immunhx + liverhx + malighx + psychhx + renalhx + transhx + aps1 + das2d3pc + scoma1 + surv2md1 + alb1 + bili1 + crea1 + hema1 + paco21 +pafi1 + ph1 + pot1 + sod1 + urin1.NA + urin1.i + wblc1,data = rhc1)
#'
#' # sensitivity analyisis
#'
#' # Using Matchit object
#'
#' sensbin (x=m.out1,data= rhc1,exposure = 'swang1',outcome = 'surv_30',Gamma = 1.5,GammaInc = 0.1)
#'
#' # Using Matching object
#'
#' sensbin (x=match1,data= rhc1,exposure = 'swang1',outcome = 'surv_30',Gamma = 1.5,GammaInc = 0.1)

sensbin <- function(x,y=NULL,data=NULL, exposure, outcome, Gamma, GammaInc, alpha = 0.05) {

    if (class(x) == "Match") {
        extractor = (c(x$index.treated, x$index.control))
        group_id = c(c(1:length(x$index.treated)), c(1:length(x$index.control)))
        match = group_id
        data1 = data[extractor, ]
        x1 = data1[, c(exposure, outcome)]
        x2 = cbind(match, x1)
        X = data.frame(x2[order(match), ])
        names(X) = c("id", "treat", "Y")  # dataframe with treat,outcome and pairID variables


    } else if (class(x) == "bmatch") {
        extractor = c(x$t_id, x$c_id)
        match = x$group_id
        data1 = data[extractor, ]
        x1 = data1[, c(exposure, outcome)]
        x2 = cbind(match, x1)
        X = x2[order(match), ]
        names(x) = c("id", "treat", "Y")  # dataframe with treat,outcome and pairID variables
        # return(x)

    } else if (class(x) == "matchit") {
      t_id = rownames(x$match.matrix)
      c_id = x$match.matrix
      extractor = as.numeric(c(t_id, c_id))
      data1s = x$model$data[extractor, c(exposure,outcome)]

      k2 = length(extractor)/2
      match = rep(1:k2, 2)
      data1s$match = match
      names(x) = c("id", "treat", "Y")  # dataframe with treat,outcome and pairID variables

    } else {
        print("object unknown")
    }

    # return(x)

    y.c <- x$Y[x$treat == 0]
    y.t <- x$Y[x$treat == 1]
    table(y.t, y.c)
    y.tmp1 <- table(y.t, y.c)[2]
    y.tmp2 <- table(y.t, y.c)[3]

    (if (y.tmp1 >= y.tmp2) {
        trt <- y.tmp1
        ctrl <- y.tmp2
    } else {
        trt <- y.tmp2
        ctrl <- y.tmp1
    })

    gamma <- seq(1, Gamma, by = GammaInc)
    mx <- ctrl + trt
    up <- c()  # creating an empty vector to store upper bound
    lo <- c()  # creating an empty vector to store lower bound
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

    vrt = bounds[bounds$min == min(bounds$min), ]$gamma
    hrz = bounds[bounds$min == min(bounds$min), ]$pupper
    vrt1 = round(bounds[bounds$min == min(bounds$min), ]$gamma, 2)
    hrz1 = round(bounds[bounds$min == min(bounds$min), ]$pupper, 2)


    plot = ggplot(data = bounds, aes(x = gamma, y = pupper)) + geom_line() + geom_point(aes(x = vrt, y = hrz)) + ylab("p upper bound") +
        xlab("gamma (Bias)") + theme_bw() + annotate("text", x = vrt1 + 0.1, y = hrz1, label = paste0("(", vrt1, ",", hrz1,
        ")")) + labs(title = "Binary Outcome Sensitivity Plot")

    # plot(bounds$pupper ~ bounds$gamma, type = 'l', xlab = 'Gamma', ylab = 'p-val upper bound', main = 'Sensitivity plot for
    # binary outcomes') text(vrt,hrz,paste0('(',vrt,',',hrz,')'),pos = 2) points(vrt,hrz,pch=15)


    colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
    msg <- "Rosenbaum Sensitivity Test \n"
    note <- "Note: Gamma is Odds of Differential Assignment To\n Treatment Due to Unobserved Factors \n"
    Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval, msg = msg, bounds = bounds[, c(1:3)], note = note, plot = plot)
    # class(Obj) <- c('rbounds', class(Obj))
    return(Obj)


}

#' Sensitivity Test for Matched Time to Event Outcomes
#'
#' This function performs sensitivity analysis for time to event outcomes in observational studies
#'
#' @param x name of the matching object. Should be either a \code{Match} or \code{matchit} object.
#' @param data dataset.
#' @param exp name of treatment variable.
#' @param outcome name of outcome variable.
#' @param failtime name of event time variable.
#' @param Gamma upper bound for sensitivity parameter.
#' @param Gammainterval value by which to increament sensitivity parameter from zero to \code{Gamma}
#' @param plot_title main title of your plot.
#' @export
#' @examples
#' data(toy)
#' psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
#' toy$ps <- psmodel$fitted
#' toy$linps <- psmodel$linear.predictors
#' X <- toy$linps ## matching on the linear propensity score
#' Tr <- as.logical(toy$treated)
#' Y <- toy$out3.time
#' match1 <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)

#' match.it <- matchit(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, data = toy, method='nearest', ratio=1)

#' Survsens(x= match.it,data =toy,exp='treated',outcome = 'out2',failtime = 'out3.time',Gamma=2.4,Gammainterval = 0.01,alpha = 0.05,plot_title = 'Time To Event Outcome Sensitivity Plot')
#'

#' @references Rosenbum 2011
#'
Survsens = function(x, y=NULL,data =NULL, exp=NULL, outcome=NULL, failtime, Gamma, alpha, Gammainterval, plot_title = NULL) {

    results = list()

    if ((class(x) != "Match") & (class(x) != "matchit"))
        {
            trt <- x  # exposure failure times
            ctrl <- y
        }
 else if (class(x) == "Match") {

        extractor = (c(x$index.treated, x$index.control))
        group_id = c(c(1:length(x$index.treated)), c(1:length(x$index.control)))
        data1s = data[extractor, c(exp, outcome, failtime)]
        data1s$match = group_id


    } else if (class(x) == "matchit") {

        # data$rId = row.names(data) data2 = data[, c('rId', exp)] data2$rId = row.names(data2) names(data2)[2] = 'exp'
        t_id = rownames(x$match.matrix)
        c_id = x$match.matrix
        extractor = as.numeric(c(t_id, c_id))
        data1s = x$model$data[extractor, c(exp, outcome, failtime)]

        k2 = length(extractor)/2
        match = rep(1:k2, 2)
        data1s$match = match
    } else {
        t_id = x$t_id
        c_id = x$c_id
        extractor(t_id, c_id)
        extractor = as.numeric(c(t_id, c_id))
        data1s = data[extractor, c(exp, outcome, failtime)]
        data1s$match = x$group_id
    }


    # data1s = subset(data1, select = c('match', exp, outcome, failtime))
    names(data1s)[c(1:3)] = c("exp", "outcome", "failtime")
    data2s = subset(data1s, exp == 0)
    data3s = subset(data1s, exp == 1)
    data4s = subset(data2s, select = c("match", "exp", "outcome", "failtime"))
    data5s = subset(data3s, select = c("match", "exp", "outcome", "failtime"))
    data6s = full_join(data4s, data5s, by = "match")

    names(data6s)[c(3, 4, 6, 7)] = c("ExpOutcome", "failtimeNotexp", "NoExpOutcome", "failtimeExp")
    data6s$timediff = data6s$failtimeNotexp - data6s$failtimeExp
    wonpairs = sum(data6s$timediff != 0)
    expoutlive = sum(data6s$timediff < 0)

    results$wonpairs = sum(data6s$timediff != 0)
    results$expoutlive = sum(data6s$timediff < 0)

    gamVal = seq(1, Gamma, by = Gammainterval)
    pplus = 1/(1 + gamVal)
    pminus = gamVal/(1 + gamVal)

    bounds = data.frame(cbind(gamVal, pplus, pminus))
    bounds$expTplus = wonpairs * bounds$pplus
    bounds$expTminus = wonpairs * bounds$pminus
    bounds$sd_expT = sqrt(wonpairs * bounds$pplus * (1 - bounds$pplus))



    for (i in 1:length(gamVal)) {
        bounds$pupper[i] = round(min(1, 2 * pnorm((expoutlive - bounds[i, 5])/bounds[i, 6], lower.tail = FALSE)), 4)
    }

    for (i in 1:length(gamVal)) {
        bounds$plower[i] = round(min(1, 2 * pnorm((expoutlive - bounds[i, 4])/bounds[i, 6], lower.tail = FALSE)), 4)
    }

    bounds2 = bounds[,c(1,8,7)]
    names(bounds2)=c("Gamma","Lower bound"," Upper bound")
    bounds1 = bounds

    bounds1$min = abs(alpha - bounds1$pupper)
    vrt = bounds1[bounds1$min == min(bounds1$min), ]$gamVal
    hrz = bounds1[bounds1$min == min(bounds1$min), ]$pupper
    vrt1 = round(bounds1[bounds1$min == min(bounds1$min), ]$gamVal, 2)
    hrz1 = round(bounds1[bounds1$min == min(bounds1$min), ]$pupper, 2)


    plot = ggplot(data = bounds1, aes(x = gamVal, y = pupper)) + geom_line() + geom_point(aes(x = vrt, y = hrz)) + ylab("p upper bound") +
        xlab("gamma (Bias)") + ylim(0, 0.06) + theme_bw() + annotate("text", x = vrt + 0.1 * vrt, y = hrz, label = paste0("(",
        vrt1, ",", hrz1, ")")) + labs(title = plot_title, caption = paste("matching done by", class(x), "function")) + theme(plot.title = element_text(hjust = 0.5))

    results$plot = plot
    results$upperbound_pval = hrz = bounds1[bounds1$min == min(bounds1$min), ]$pupper
    results$Gamma = bounds1[bounds1$min == min(bounds1$min), ]$gamVal
    results$bounds = bounds2

    return(results)


}

#' Sensitivity Test for Continous outcomes
#'
#' This function performs sensitivity analysis for continious outcomes in observational studies
#'
#' @param x Treatment group outcomes in same order as treatment group outcomes or an objects from a Match,matchit or bmatch package.
#' @param y Control group outcomes in same order as treatment group outcomes unnecessary when using a Match,matchit or bmatch packages object.
#' @param exp treatment variable. Must be specified for MatchIt objects not required for match objects.
#' @param outcome treatment variable. Must be specified for MatchIt objects not required for match objects.
#' @param Gamma Upper-bound on gamma parameter.
#' @param GammaInc To set user-specified increments for gamma parameter.
#' @param alpha significance level.
#' @param plot_title plot title.
#' @export
#' @examples
#' library(ggplot2);library(Matching);library(MatchIt)
#' data(toy)
#'
#' ## Creating Matching object
#' psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
#' toy$linps = psmodel$fitted.values
#' X <- toy$linps
#' Tr <- as.logical(toy$treated)
#' Y = toy$out1.cost
#' match1 <- Match(Y=Y,Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
#'
#' contSens(x=match1, Gamma=5,GammaInc = 0.1,alpha = 0.05,plot_title = 'Time To Event Outcome Sensitivity Plot') #using mathit object
#'
#' ## Creating MatchIt object
#' match.it <- matchit(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, data = toy, method='nearest', ratio=1)
#'
#'contSens(x=match.it ,exp = 'treated',outcome = 'out1.cost', CausalEst = 15.5,Gamma=5,GammaInc = 0.1,alpha = 0.05,plot_title = 'Time To Event Outcome Sensitivity Plot') #using mathit object
#'

#'
#' @author David Ngendahimana, Case Western Reserve University.
#' @references Rosenbaum, Paul R. (2002) Observational Studies. Springer-Verlag.


contSens = function(x, y = NULL,data =NULL, exp = NULL, outcome = NULL, CausalEst = NULL, Gamma = NULL, GammaInc = NULL, alpha = NULL,
    plot_title = NULL) {
    if ((class(x) != "Match") & (class(x) != "matchit")) {
        trt <- x
        ctrl <- y
    } else if (class(x) == "Match") {
        if (x$est > 0) {
            ctrl <- x$mdata$Y[x$mdata$Tr == 0]
            trt <- x$mdata$Y[x$mdata$Tr == 1]
        } else {
            ctrl <- x$mdata$Y[x$mdata$Tr == 1]
            trt <- x$mdata$Y[x$mdata$Tr == 0]
        }

    } else if (class(x) == "matchit") {

        # data$rId = row.names(data) data2 = data[, c('rId', exp)] data2$rId = row.names(data2) names(data2)[2] = 'exp'
        t_id = rownames(x$match.matrix)
        c_id = x$match.matrix
        extractor = as.numeric(c(t_id, c_id))
        data1 = x$model$data[extractor, c(exp, outcome)]
        names(data1) = c("exp", "outcome")

        if (CausalEst > 0) {
            ctrl <- data1$outcome[data1$exp == 0]
            trt <- data1$outcome[data1$exp == 1]
        } else {
            ctrl <- data1$outcome[data1$exp == 1]
            trt <- data1$outcome[data1$exp == 0]
        }
    } else {

        print("error")
    }

    gamma <- seq(1, Gamma, by = GammaInc)
    m <- length(gamma)
    pvals <- matrix(NA, m, 2)
    diff <- trt - ctrl
    S <- length(diff)
    diff <- diff[diff != 0]
    ranks <- rank(abs(diff), ties.method = "average")
    psi <- as.numeric(diff > 0)
    T <- sum(psi * ranks)
    for (i in 1:m) {
        p.plus <- gamma[i]/(1 + gamma[i])
        p.minus <- 1/(1 + gamma[i])
        E.T.plus <- sum(ranks * p.plus)
        V.T <- sum(ranks^2 * p.plus * (1 - p.plus))
        E.T.minus <- sum(ranks * p.minus)
        z.plus <- (T - E.T.plus)/sqrt(V.T)
        z.minus <- (T - E.T.minus)/sqrt(V.T)
        p.val.up <- 1 - pnorm(z.plus)
        p.val.low <- 1 - pnorm(z.minus)
        pvals[i, 1] <- round(p.val.low, digits = 4)
        pvals[i, 2] <- round(p.val.up, digits = 4)
    }
    pval <- pvals[1, 1]
    bounds <- data.frame(gamma, pvals)
    names(bounds) <- c("Gamma", "Lower bound", "Upper bound")

    bounds1 = bounds
    bounds1$min = abs(alpha - bounds1$`Upper bound`)
    vrt = bounds1[bounds1$min == min(bounds1$min), ]$Gamma
    hrz = bounds1[bounds1$min == min(bounds1$min), ]$`Upper bound`
    vrt1 = round(bounds1[bounds1$min == min(bounds1$min), ]$Gamma, 2)
    hrz1 = round(bounds1[bounds1$min == min(bounds1$min), ]$`Upper bound`, 2)


    plot = ggplot(data = bounds1, aes(x = Gamma, y = `Upper bound`)) + geom_line() + geom_point(aes(x = vrt, y = hrz)) + ylab("p upper bound") +
        xlab("gamma (Bias)") + ylim(0, 0.08) + theme_bw() + annotate("text", x = vrt + 0.1 * vrt, y = hrz, label = paste0("(",
        vrt1, ",", hrz1, ")")) + labs(title = plot_title, caption = paste("matching done by", class(x), "function")) + theme(plot.title = element_text(hjust = 0.5))


    upperbound_pval = hrz = bounds1[bounds1$min == min(bounds1$min), ]$`Upper bound`
    Gamma = bounds1[bounds1$min == min(bounds1$min), ]$Gamma

    msg <- "Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value \n"
    note <- "Note: Gamma is Odds of Differential Assignment To\n Treatment Due to Unobserved Factors \n"
    Obj <- list(Gamma = Gamma, pval = upperbound_pval, alpha = alpha, msg = msg, note = note, plot = plot, bounds = bounds)
    # class(Obj) <- c('rbounds', class(Obj))
    Obj
}


#' Creating table One
#'
#' @param x  Count of pairs where only control has outcome.
#' @param y  Count of pairs where onlytreated has outcome.
#' @export
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' edaTable (220, 380)


edaTable = function(baselinevars, outcomeVar, data) {
    results = list()
    formula = reformulate(termlabels = baselinevars, response = outcomeVar)
    res1 <- compareGroups(formula, data = dataA, ref = 1)
    table01 = createTable(res1, show.p.mul = TRUE, show.p.overall = TRUE)
    try(export2pdf(table01, file = "table1.pdf"), silent = T)

}


