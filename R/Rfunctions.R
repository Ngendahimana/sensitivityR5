

#' Estimating design sensitivity.
#'
#' This function estimates design sensitivity parameter associated with a list of matching algorithms.
#' @param x A list of matched sampling algorithms.
#' @param sf  A fraction used to split the entire matched sample into planning sample (to estimate design sensitivity) and analysis sample (for outcome analysis).
#' @param Card_Dat Dataframe passed to the function implementing cardinality matching. Passed to \code{ds_function} only if cardinality matching is one of the matching algorithms to be evaluated
#' @export
#' @examples
#' # Implementing different types of matching using matchit package
#' m.out_Nrst <- matchit(exposure~V2+V3+V5+V6+V8+V9,data=dat, method = "nearest", distance = "logit",replace = TRUE ,print = FALSE,ratio = 1)
#' m.out_NrstCAL <- matchit(exposure~V2+V3+V5+V6+V8+V9,data=dat, method = "nearest", distance = "logit",replace = TRUE,print = FALSE,ratio = 1,caliper=0.2)
#' m.out_optimal <- matchit(exposure~V2+V3+V5+V6+V8+V9,data=dat, method = "optimal", distance = "logit",replace =TRUE,print = FALSE,ratio = 1)
#' m.out_genetic <- matchit(exposure~V2+V3+V5+V6+V8+V9,data=dat, method = "genetic", distance = "logit",replace = TRUE,print = FALSE,ratio = 1,unif.seed=1945, int.seed=1906)
#'
#'
#' # Implementing cardinality matching
#' # Solver options
#' t_max = 60*5
#' solver = "glpk"
#' approximate = 1
#' solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)
#' Card_Dat = dat[order(dat$exposure==1, decreasing = TRUE), ]
#' Card_Dat$`.distance` <- (glm(exposure~V2+V3+V5+V6+V8+V9, data = Card_Dat ,family = binomial(link="logit")))$fitted.values
#'
#' attach(Card_Dat)
#' t_ind = exposure
#' # Distance matrix
#' dist_mat = NULL
# Subset matching weight
#' subset_weight = 1
#'
#' ## Setting balance criteria
#' mom_covs = cbind(.distance,V2,V3,V5,V6,V8,V9)
#' mom_tols = round(absstddif(mom_covs, t_ind, .1), 2)
#' mom = list(covs = mom_covs, tols = mom_tols)
#'
#' ## Creating cardinality matching object
#' m.out_dsnMatch = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,mom = mom,solver = solver)
#'
#' # Assessing balance
#' t_id1 = m.out_dsnMatch$t_id
#' c_id1 = m.out_dsnMatch$c_id
#' dsn_bal = meantab(mom_covs, t_ind, t_id1, c_id1)[,6]
#' Card_Dat_A = rbind(Card_Dat[t_id1,],Card_Dat[c_id1,])
#' detach(Card_Dat)
#' arg.list1 = list("a"=m.out_Nrst , "b"=m.out_NrstCAL)
#' arg.list2= list("a"=m.out_Nrst , "b"=m.out_NrstCAL, "c"=m.out_dsnMatch, "d"=m.out_optimal)
#' k=ds_function(arg.list2,Card_Dat = Card_Dat)
#' k$designSensitivity


ds_function <-function (x,sf = 1/3,Card_Dat=NULL)
{

  set.seed(1)

  if(sum(unlist(lapply(x,class)) %in%c("matchit","list")) !=length(unlist(lapply(x,class)))){
    stop("All objects must be of class matchit or list")
  }
  else{
    x1=x[unlist(lapply(x,class))=="matchit"]
    unmatchedData = x1$a$model$data
    ds_df1 = list()
    analysisDataframe1 = list()
    #print("... Design Sensitivity Results ...")

    for(i in 1:length(x1)){
      dat = x1[[i]]$model$data
      ctrlUnits = dat[as.numeric(as.vector(x1[[i]]$match.matrix[,1])),]
      trtUnits = dat[as.numeric(row.names(x1[[i]]$match.matrix)),]
      ctrlUnits$pair_id = c(1:nrow(ctrlUnits))
      trtUnits$pair_id = c(1:nrow(trtUnits))

      matchDat = rbind(trtUnits,ctrlUnits) # matched data with pair_ids
      ctrl_outcomes = ctrlUnits[,c("pair_id","Y2")]
      trt_outcomes = trtUnits[,c("pair_id","Y2")]

      matchDat_outcomes= left_join(trt_outcomes,ctrl_outcomes,by = "pair_id")
      names(matchDat_outcomes)[c(2,3)]<-c("trt_outcomes","ctrl_outcomes")
      matchDat_outcomes = na.omit(matchDat_outcomes)
      matchDat_outcomes$diff = matchDat_outcomes$trt_outcomes-matchDat_outcomes$ctrl_outcomes

      matchDat_plan = sample_n(matchDat_outcomes,floor(sf*nrow(matchDat_outcomes)))
      rownames(matchDat_plan) = c(1:nrow(matchDat_plan))
      matchDat_analysis = subset(matchDat,!(pair_id %in% matchDat_plan$pair_id))

      sumlogicP1 = function(df,diff){
        x1 =sample(1:(nrow(df)-2),1)
        x2 =sample((x1+1):(nrow(df)),1)
        sum(df[x1,diff],df[x2,diff])>0
      }

      p1 = mean(replicate(1000,sumlogicP1(matchDat_plan,"diff")))
      #ds_df[i] = p1/(1-p1)
      ds_df1[[i]] = ifelse((p1/(1-p1))<1,1,(p1/(1-p1)))

      analysisDataframe1[[i]] = matchDat_analysis

    }

    # CARDINALITY MATCHING
    ds_df2 = list()
    analysisDataframe2 =list()
    x2=x[unlist(lapply(x,class))=="list"]

    if(length(x2)==0)
    {
      ds_df =ds_df1
      names(ds_df) <- names(x)
      names(analysisDataframe1)=names(x)
      ds_df2=data.frame(matrix(unlist(ds_df), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
      ds_df2$Algorithm = names(x)
      names(ds_df2) = c("Design.Sensitivity","Algorithm")
      ds_df2 = ds_df2[,c("Algorithm","Design.Sensitivity")]
      ds_df2 = ds_df2[order(ds_df2$Design.Sensitivity),]
      dataframes = analysisDataframe1
      names(dataframes)=names(x)


    }

    else
    {

      for(i in 1:length(x2)){

        ctrlUnits = Card_Dat[x2[[i]]$c_id,]
        trtUnits = Card_Dat[x2[[i]]$t_id,]
        ctrlUnits$pair_id = c(1:nrow(ctrlUnits))
        trtUnits$pair_id = c(1:nrow(trtUnits))

        matchDat = rbind(trtUnits,ctrlUnits) # matched data with pair_ids
        ctrl_outcomes = ctrlUnits[,c("pair_id","Y2")]
        trt_outcomes = trtUnits[,c("pair_id","Y2")]

        matchDat_outcomes= left_join(trt_outcomes,ctrl_outcomes,by = "pair_id")
        names(matchDat_outcomes)[c(2,3)]<-c("trt_outcomes","ctrl_outcomes")
        matchDat_outcomes = na.omit(matchDat_outcomes)
        matchDat_outcomes$diff = matchDat_outcomes$trt_outcomes-matchDat_outcomes$ctrl_outcomes

        matchDat_plan = sample_n(matchDat_outcomes,floor(sf*nrow(matchDat_outcomes)))
        rownames(matchDat_plan) = c(1:nrow(matchDat_plan))
        matchDat_analysis = subset(matchDat,!(pair_id %in% matchDat_plan$pair_id))

        sumlogicP1 = function(df,diff){
          x1 =sample(1:(nrow(df)-2),1)
          x2 =sample((x1+1):(nrow(df)),1)
          sum(df[x1,diff],df[x2,diff])>0
        }

        p1 = mean(replicate(1000,sumlogicP1(matchDat_plan,"diff")))
        #ds_df[i] = p1/(1-p1)
        ds_df2[[i]] = ifelse((p1/(1-p1))<1,1,(p1/(1-p1)))
        analysisDataframe2[[i]] = matchDat_analysis

      }
      ds_df = c(ds_df1,ds_df2)
      names(ds_df) <- names(x)
      ds_df2=data.frame(matrix(unlist(ds_df), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
      ds_df2$Algorithm = names(x)
      names(ds_df2) = c("Design.Sensitivity","Algorithm")
      ds_df2 = ds_df2[,c("Algorithm","Design.Sensitivity")]
      ds_df2 = ds_df2[order(ds_df2$Design.Sensitivity),]
      dataframes = c(analysisDataframe1,analysisDataframe2)
      names(dataframes)=names(x)
    }

  }


  Obj <- list(designSensitivity =ds_df2, AnalysisDataFrames = dataframes,UnmatchedData =unmatchedData )
  class(Obj) <- "DSensitivity"
  Obj

}



#' Assessing the three Rubin Rules.
#'
#' This function allows you to assess how sensitive your results are to unmeasured variable.
#' @param data data set to be used.
#' @param Treatment A variables defining exposure group.
#' @param matchscore Variable containing matching distance.Default is propensity score.
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


rubinRules2 = function(data, Treatment, matchscore = "ps", covlist) {

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
  res2 = data.frame(res)
  res2$VarName = rownames(res2)
  rownames(res2) = 1:dim(res2)[1]
  res2 = res2[,c("VarName","res")]
  names(res2)[2] = "Rubin3"

  results$RUBIN3 = res2


  res3 <- res2 %>%
    dplyr::arrange(Rubin3) %>%
    dplyr::mutate(VarName = factor(VarName, levels = .$VarName))
  p0 = ggplot(res3, aes(Rubin3, VarName))+geom_point()+
    geom_vline(xintercept = 1,colour = "black")+
    geom_vline(xintercept = 0.8,colour = "blue",linetype = "dashed")+
    geom_vline(xintercept = 1.25,colour = "blue",linetype = "dashed")+
    geom_vline(xintercept = 0.5,colour = "red",linetype = "dashed")+
    geom_vline(xintercept = 2,colour = "red",linetype = "dashed")+ labs(x = "Residual Variance Ratio",y= " Variables",title = "Rubin's Rules Plot",subtitle = paste0(" Rubin One:",round(results$RUBIN1, 2),"                                                                Rubin Two: ",round(results$RUBIN2, 2)))



  results$plot = p0

  return(results)

}

#' Amplifying and visualizing gamma parameter.
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
Survsens = function(x, y=NULL,data =NULL, exp=NULL, outcome=NULL, failtime, Gamma, alpha=0.05, Gammainterval, plot_title = NULL) {

  results = list()

  if ((class(x) != "Match") & (class(x) != "matchit") & (class(x) !="list"))
  {
    wonpairs <- x  # won pairs
    expoutlive <- y

    results$wonpairs = x
    results$expoutlive = y
  }

  else{

    if (class(x) == "Match") {

      extractor = (c(x$index.treated, x$index.control))
      group_id = c(c(1:length(x$index.treated)), c(1:length(x$index.control)))
      data1s = data[extractor, c(exp, outcome, failtime)]
      data1s$match33 = group_id


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
    data6s = dplyr::full_join(data4s, data5s, by = "match")

    names(data6s)[c(3, 4, 6, 7)] = c("ExpOutcome", "failtimeNotexp", "NoExpOutcome", "failtimeExp")
    data6s$timediff = data6s$failtimeNotexp - data6s$failtimeExp
    wonpairs = sum(data6s$timediff != 0)
    expoutlive = sum(data6s$timediff < 0)

    results$wonpairs = sum(data6s$timediff != 0)
    results$expoutlive = sum(data6s$timediff < 0)

  }

  ## Universal code

  gamVal = seq(1, Gamma, by = Gammainterval)
  pminus = 1/(1 + gamVal)
  pplus = gamVal/(1 + gamVal)

  bounds = data.frame(cbind(gamVal, pplus, pminus))
  bounds$expTplus = wonpairs * bounds$pplus
  bounds$expTminus = wonpairs * bounds$pminus
  bounds$sd_expT = sqrt(wonpairs * bounds$pplus * (1 - bounds$pplus))


  for (i in 1:length(gamVal)) {
    bounds$pupper[i] = round(min(1, 2 * pnorm((expoutlive - bounds[i, 4])/bounds[i, 6], lower.tail = FALSE)), 4)
  }

  for (i in 1:length(gamVal)) {
    bounds$plower[i] = round(min(1, 2 * pnorm((expoutlive - bounds[i, 5])/bounds[i, 6], lower.tail = FALSE)), 4)
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


#' Sensitivity analysis with Matching, MatchIt and designmatch objects for a continous outcome.
#'
#' @param x Treatment group outcomes or an objects from a Match,MatchIt or designmatch.
#' @param y Control group outcomes in same order as treatment group outcomes such that members of a pair occupy the same row in both x and y. Should not be specified x is a Matching, MatchIt and designmatch objects.
#' @param est Treatment effect. Must be specified if x is not a MatchIt object.
#' @param Gamma Upper bound of sensitivity parameter
#' @param GammaInc interval width for increasing gamma from 1 until the specified upper bound of sensitivity parameter is reached.
#' @param data Dataframe used to during matching. You do not have to specify this parameter if x is a MatchIt object
#' @param treat Treatmetn/Exposure variable name.
#' @export
#' @return a table of Rosenbaum bounds
#' @examples
#' ## Loading lalonde data
#' library(Matching);library(MatchIt)
#' data("lalonde",package = "Matching")
#'
#' ## Sensitivity analysis with a matchit object
#' m.out = matchit(treat ~ age  + educ +  black + hisp +married + nodegr + re74  + re75  +
#' u74 + u75, family=binomial, data = lalonde, method = "nearest")
#'
#' ## Estimating treatment effect. Ideally, balance assessement should be done prior to estimating treatment effect

#' mod = lm(re78~age  + educ +  black + hisp +married + nodegr + re74  + re75  +u74 + u75,data = match.data(m.out))
#' pens2(x = m.out, y="re78",Gamma = 2, GammaInc = 0.1,est = 629.7)
#'
#' ## using match object

#' ## Estimate Propensity Score
#' data("lalonde",package = "Matching")
#' DWglm <- glm(treat ~ age + I(age^2) + educ + I(educ^2) + black + hisp +married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +u74 + u75, family=binomial, data=lalonde)
#' ## Specifying outcome and treatment(Exposure)
#' Y  <- lalonde$re78
#' Tr <- lalonde$treat
#' ## Matching using Match object.
#' mDW  <- Match(Y=Y, Tr=Tr, X=DWglm$fitted, replace=FALSE)

#' ## Sensitivity analysis
#' pens2(mDW, Gamma = 2, GammaInc = 0.1)
#'
#' ## Sensitivity analysis with a matchit object

#' data("lalonde",package = "designmatch")
#' library(designmatch)
#' attach(lalonde)

#' ## Treatment indicator
#' t_ind = treatment

#' ## Distance matrix
#' dist_mat = NULL

#' ## Subset matching weight
#' subset_weight = 1

#' ## Moment balance: constrain differences in means to be at most .05 standard deviations apart
#' mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
#' mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
#' mom = list(covs = mom_covs, tols = mom_tols)

#' ## Fine balance
#' fine_covs = cbind(black, hispanic, married, nodegree)
#' fine = list(covs = fine_covs)

#' ## Exact matching
#' exact_covs = cbind(black)
#' exact = list(covs = exact_covs)

#' ## Solver options
#' t_max = 60*5
#' solver = "glpk"
#' approximate = 1
#' solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)
#' ## Match
#' out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,mom = mom, fine = fine, exact = exact, solver = solver)
#'
#' ## Sensitivity analysis with designmatch Object
#' pens2(x = out, y="re78",Gamma = 2, GammaInc = 0.1,est = 234,treat = "treatment",data = lalonde)

#' detach(lalonde)

pens2 =function (x, y = NULL, est = NULL, Gamma = 2, GammaInc = 0.1,
  data = NULL, treat)
{
  if (length(x) == 1) {
    ctrl <- x
    trt <- y
  }
  else {
    if (class(x) == "Match") {
      if (x$est > 0) {
        ctrl <- x$mdata$Y[x$mdata$Tr == 0]
        trt <- x$mdata$Y[x$mdata$Tr == 1]
      }
      else {
        ctrl <- x$mdata$Y[x$mdata$Tr == 1]
        trt <- x$mdata$Y[x$mdata$Tr == 0]
      }
    }
    else if (class(x) == "matchit") {
      if (missing(est)) {
        stop("Est parameter missing")
      }
      else if (est > 0) {
        trt = (x$model$data[as.numeric(row.names(x$match.matrix)),])[[y]]
        ctrl = (x$model$data[as.numeric(x$match.matrix[, 1]), ])[[y]]
      }
      else {
        ctrl = (x$model$data[row.names(x$match.matrix),
          ])[[y]]
        trt = (x$model$data[x$match.matrix[, 1], ])[[y]]
      }
    }
    else if (class(x) == "list") {
      if (missing(est) | missing(data)) {
        warning("Est or Data not specified")
      }
      else {
        data = dplyr::arrange(data, desc(data[[treat]]))
        rownames(data) = 1:dim(data)[1]
        if (est > 0) {
          trt = (data[as.character(x$t_id), ])[[y]]
          ctrl = (data[as.character(x$c_id), ])[[y]]
        }
        else {
          ctrl = (data[as.character(x$t_id), ])[[y]]
          trt = (data[as.character(x$c_id), ])[[y]]
        }
      }
    }
  }
  gamma <- seq(1, Gamma, by = GammaInc)
  m <- length(gamma)
  pvals <- matrix(NA, m, 2)
  diff <- trt - ctrl
  diff <- na.omit(diff)
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
  msg <- "Rosenbaum Sensitivity Test for Wilcoxon Signed Rank P-Value \n"
  note <- "Note: Gamma is Odds of Differential Assignment To\n Treatment Due to Unobserved Factors \n"
  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval,
    msg = msg, bounds = bounds, note = note)
  class(Obj) <- c("rbounds", class(Obj))
  Obj
}


#' Sensitivity analysis with Matching, MatchIt and designmatch objects for a binary outcome.
#'
#' @param x Treatment group outcomes or an objects from a Match,MatchIt or designmatch.
#' @param y Control group outcomes in same order as treatment group outcomes such that members of a pair occupy the same row in both x and y. Should not be specified x is a Matching, MatchIt and designmatch objects.
#' @param Gamma Upper bound of sensitivity parameter
#' @param GammaInc interval width for increasing gamma from 1 until the specified upper bound of sensitivity parameter is reached.
#' @param data Dataframe used to during matching. You do not have to specify this parameter if x is a MatchIt object
#' @param treat Treatmetn/Exposure variable name.
#' @param alpha p-value to define maximum upper bound allowable
#' @export
#' @return a table of Rosenbaum bounds
#' @examples
#'
#' ## Sensitivity analysis with a matchit object
#' library(Matching);library(MatchIt);library(designmatch)
#' data("GerberGreenImai",package = "Matching")

#' ## Estimate Propensity Score
#' pscore.glm <- glm(PHN.C1 ~ PERSONS + VOTE96.1 + NEW +MAJORPTY + AGE + WARD + PERSONS:VOTE96.1 + PERSONS:NEW + AGE2, family = binomial(logit), data = GerberGreenImai)

#' ## save data objects
#' D <- GerberGreenImai$PHN.C1
#' Y <- GerberGreenImai$VOTED98
#' X <- fitted(pscore.glm)

#' ## Match - without replacement
#' m.obj <- Match(Y = Y, Tr = D, X = X, M = 1, replace=FALSE)

#' ## Sensitivity Test
#' binarysens2(m.obj, Gamma=2, GammaInc=.1)

#' ## Sensitivity analysis with a Match object

#' m.out = matchit(PHN.C1 ~ PERSONS + VOTE96.1 + NEW +MAJORPTY + AGE + WARD , family=binomial, data = GerberGreenImai, method = "nearest")

#' mod = lm(VOTED98 ~ PHN.C1+PERSONS + VOTE96.1 + NEW +MAJORPTY + AGE + WARD,data = match.data(m.out))

#' binarysens2(x=m.out,y ="VOTED98", Gamma=2, GammaInc=.1)


#' ## Sensitivity analysis with a designmatch object


#' ## data("GerberGreenImai",package = "Matching")
#' attach(GerberGreenImai)

#' ## Treatment indicator
#' t_ind = PHN.C1

#' ## Distance matrix
#' dist_mat = NULL

#' ## Subset matching weight
#' subset_weight = 1

# Moment balance: constrain differences in means to be at most .05 standard deviations apart
#' mom_covs = cbind(PERSONS,VOTE96.1 ,NEW , MAJORPTY , AGE , WARD)
#' mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
#' mom = list(covs = mom_covs, tols = mom_tols)


#' ## Solver options
#' t_max = 60*5
#' solver = "glpk"
#' approximate = 1
#' solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)
#' ## Match
#' out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,mom = mom,solver = solver)

#' binarysens2(x=out,y ="VOTED98", Gamma=2, GammaInc=.1,treat = "PHN.C1",data = GerberGreenImai)

#' detach(GerberGreenImai)

binarysens2 = function (x, y = NULL, Gamma = 6, GammaInc = 1,data =NULL,treat =NULL,alpha = 0.05)
{
  if (length(x) == 1) {
    ctrl <- x
    trt <- y
  }
  else {
    if (class(x)=="matchit"){
      if(missing(y)){
        warning("y not defined")
      }
      else{
        y.t = (x$model$data[row.names(x$match.matrix),])[[y]]
        y.c = (x$model$data[x$match.matrix[,1],])[[y]]
      }

    }

    else if(class(x) == "Match"){
      y.c <- x$mdata$Y[x$mdata$Tr == 0]
      y.t <- x$mdata$Y[x$mdata$Tr == 1]

    }

    else if(class(x)== "list"){
      data = dplyr::arrange(data,desc(data[[treat]]))
      rownames(data) = 1:dim(data)[1]
      y.t = (data[as.character(x$t_id),])[[y]]
      y.c = (data[as.character(x$c_id),])[[y]]

    }

    else {
      print("Accepts only matchit, match or bmatch objects.")
    }

    table(y.t, y.c)
    y.tmp1 <- table(y.t, y.c)[2]
    y.tmp2 <- table(y.t, y.c)[3]
    if (y.tmp1 >= y.tmp2) {
      trt <- y.tmp1
      ctrl <- y.tmp2
    }
    else {
      trt <- y.tmp2
      ctrl <- y.tmp1
    }

  }

  gamma <- seq(1, Gamma, by = GammaInc)
  mx <- ctrl + trt
  up <- c()
  lo <- c()
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


  colnames(bounds) <- c("Gamma", "Lower bound", "Upper bound")
  msg <- "Rosenbaum Sensitivity Test \n"
  note <- "Note: Gamma is Odds of Differential Assignment To\n Treatment Due to Unobserved Factors \n"
  Obj <- list(Gamma = Gamma, GammaInc = GammaInc, pval = pval, msg = msg, bounds = bounds[, c(1:3)], note = note, plot = plot)
  #class(Obj) <- c("rbounds", class(Obj))
  return(Obj)
}



#' Love plots with `designmatch`,`matchIt`,`Matching` Object
#'
#' ## Love plot with `designmatch` package object
#'
#' @param x A `designmatch`,`matchIt`,`Matching` Object
#' @param data dataset before its preprocessed.
#' @param covList a vector of covariates to be balanced.
#' @param treat Exposure variable name
#' @return a plot of standardized difference.
#' @export
#'
#' @examples
#' data("lalonde",package = "cobalt")
#' attach(lalonde)
#'
#' ## Treatment indicator
#' t_ind =lalonde$treat
#'
#' ## Distance matrix
#' dist_mat = NULL
#'
#' ## Subset matching weight.
#' subset_weight = 1
#'
#' # Moment balance: constrain differences in means to be at most .05 standard deviations apart
#' mom_covs = cbind(age, educ, black, hispan, married, nodegree, re74, re75)
#' mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
#' mom = list(covs = mom_covs, tols = mom_tols)
#'
#' ## Fine balance
#' fine_covs = cbind(black, hispan, married, nodegree)
#' fine = list(covs = fine_covs)
#'
#' ## Exact matching
#' exact_covs = cbind(black)
#' exact = list(covs = exact_covs)
#'
#' ## Solver options
#' t_max = 60*5
#' solver = "glpk"
#' approximate = 1
#' solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)
#'
#' ## Cardinality matching
#' out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, mom = mom, fine = fine, exact = exact, solver = solver)
#'
#' # Indices of the treated units and matched controls
#' t_id = out$t_id
#' c_id = out$c_id
#' detach(lalonde)
#'
#' # Example of cardinality Matching
#' love_plot(X =out, data = lalonde , covList=c("age", "educ", "black", "hispan", "married", "nodegree", "re74", "re75"))
#'
#' # Example with  MatchIt package
#' # data("lalonde") if not yet loaded
#' covs0 <- subset(lalonde, select = -c(treat, re78, nodegree, married))

#' # Nearest neighbor 1:1 matching with replacement
#' library("MatchIt") #if not yet loaded
#' m.out <- matchit(f.build("treat", covs0), data = lalonde, method = "nearest", replace = TRUE)

#' love_plot(X =m.out, data = lalonde , covList=c("age", "educ", "black", "hispan", "married", "nodegree", "re74", "re75"))
#'
#' # Example with Matching package
#'
#' library("Matching")
#' #data("lalonde") #If not yet loaded
#' covs0 <- subset(lalonde, select = -c(treat, re78, nodegree, married))

#' fit <- glm(f.build("treat", covs0), data = lalonde, family = "binomial")
#' p.score <- fit$fitted.values
#' match.out <- Match(Tr = lalonde$treat, X = p.score, estimand = "ATT")

#' std_dif_data =bal.tab(match.out, formula = f.build("treat", covs0), data = lalonde)

#' love_plot(X, data = lalonde , covList=c("age", "educ", "black", "hispan", "re74", "re75"),treat = "treat")

love_plot = function (X,data,covList, legend_position = "topright",treat=NULL)
{
  if (class(X) =="list"){

    if (missing(data)|missing(covList)){

      stop("Data or list of covariate need to be specified")
    }

    else {

      attach(data)
      X_mat = as.matrix(data[,covList])
      detach(data)

      X_mat_t = X_mat[X$t_id, ] # Extract treated observations
      X_mat_c_before = X_mat[-X$t_id, ] # Extract control observations before matching
      X_mat_c_before_mean = apply(X_mat_c_before, 2, mean) # Extract control variable means before matching
      X_mat_t_mean = apply(X_mat_t, 2, mean) # extract mean of treated observations
      X_mat_t_var = apply(X_mat_t, 2, var) # extract variance of treated observations
      X_mat_c_before_var = apply(X_mat_c_before, 2, var) # Extract variance before matching for control obserations
      std_dif_before = (X_mat_t_mean - X_mat_c_before_mean)/sqrt((X_mat_t_var + X_mat_c_before_var)/2) # std_df before matching


      X_mat_c_after = X_mat[X$c_id, ] # extract matched observations
      X_mat_c_after_mean = apply(X_mat_c_after, 2, mean) # extract variable mean for treated observations after matching
      std_dif_after = (X_mat_t_mean - X_mat_c_after_mean)/sqrt((X_mat_t_var +  X_mat_c_before_var)/2) # Compute std_diff after matching

      abs_std_dif_before = data.frame(std_dif_before)
      abs_std_dif_before$VarName =rownames(abs_std_dif_before)
      #abs_std_dif_before$Match = "Before Matching"
      rownames(abs_std_dif_before) = NULL
      names(abs_std_dif_before)[1] = "abs_std_before"

      abs_std_dif_after = data.frame(std_dif_after)
      abs_std_dif_after$VarName =rownames(abs_std_dif_after)
      #abs_std_dif_after$Match = "After Matching"
      rownames(abs_std_dif_after) = NULL
      names(abs_std_dif_after)[1] = "abs_std_after"

      std_dif_dat2 = dplyr::left_join(abs_std_dif_before,abs_std_dif_after,by = c("VarName"))

    }
  }

  else if (class(X) =="matchit"){

    std_dif_dat = bal.tab(X)$Balance[,c("Diff.Un","Diff.Adj")][-1,]
    std_dif_dat$VarName =rownames(std_dif_dat)
    rownames(std_dif_dat) = NULL
    std_dif_dat2 = std_dif_dat
    names(std_dif_dat2) = c("abs_std_before","abs_std_after","VarName")

  }

  else if(class(X) =="Match"){

    if(missing(treat)|missing(data)){
      stop(" Treatment Variable or Data not specified")
    }
    else{

      #bal.tab(X, formula = f.build("treat", covList), data = data)
      std_dif_dat =bal.tab(X, formula = f.build(treat, covList), data = data)$Balance[,c("Diff.Un","Diff.Adj")][-1,]
      std_dif_dat$VarName =rownames(std_dif_dat)
      rownames(std_dif_dat) = NULL
      std_dif_dat2 = std_dif_dat
      names(std_dif_dat2) = c("abs_std_before","abs_std_after","VarName")

    }
  }

  else{

    stop("Object has to be MatchIt, Matching or designMatch object")

  }


  ggplot(std_dif_dat2, aes(x = value, y = VarName, color = Variables)) +
    geom_point(aes(x = abs_std_before, col = "Before Matching")) +
    geom_point(aes(x = abs_std_after, col = "After Matching"))+ theme(legend.title = element_blank())


}

#' Sensitivity Analysis with multiple controls
#'
#' @param X A `MatchIt` object
#' @param outcomeName Name of the outcome.
#' @param Gamma Upper bound of sensitivity parameter
#' @param GammaInc Interval width of sequence between 1 and Gamma.
#' @param n_contrl Number of control matched to each treated observation
#' @export
#'
#' @examples

#' data("lalonde",package ="designmatch")
#' covs0 <- subset(lalonde, select = -c(treat, re78, nodegree, married))
#' m.out <- matchit(f.build("treatment", covs0), data = lalonde, method = "nearest", replace = TRUE,ratio = 2)
#' multiControlSens(X =m.out,outcomeName = "re78",Gamma = 2,GammaInc = 0.1,n_contrl = 2)


multiControlSens = function(X,outcomeName,Gamma,GammaInc,n_contrl){
  results4 = list()
  data1 = get_matches(m.out,model_frame = m.out$model$data)
  treat = data1[which(row.names(data1) %in% rownames(m.out$match.matrix)),][,outcomeName]
  k = list()
  for(i in 1:dim(m.out$match.matrix)[2]){
    k[[i]] = data1[as.vector(m.out$match.matrix[,i]),outcomeName]
    #assign((paste0("c",i)),data1[as.vector(m.out$match.matrix[,i]),treat])
  }
  data2 = as.data.frame(cbind(treat,do.call(cbind, k,n_contrl)))
  names(data2) = c("Treated","Control 1","Control 2")
  results4$data = data2

  gma = NULL
  for(i in seq(1,Gamma,by=GammaInc )){
    gma=append(gma,senmw(data2, gamma = i, method = "t")$pval)
  }

  sensData = data.frame(cbind(seq(1,Gamma,by=GammaInc ),round(gma,2)))
  names(sensData) = c("Gamma","pval Upper Bound")

  results4$pvalues = sensData
  return(results4)
}




