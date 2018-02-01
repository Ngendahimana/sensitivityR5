## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----packages,include=FALSE,echo=FALSE-----------------------------------

library(tidyr);library(ggplot2);library(designmatch);library(MatchIt);library(Matching);library(sensitivitymv);library(sensitivitymw);library(knitr);library(cobalt);library(dplyr)


## ----packageinstall,message=FALSE----------------------------------------
library(devtools)
devtools::install_github("Ngendahimana/SensitivityR5")


## ----rubinRules,fig.width = 5, fig.asp =1,message=FALSE,warning=FALSE----

library(SensitivityR5)

data("toy",package = "SensitivityR5")
psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
toy$ps <- psmodel$fitted
toy$linps <- psmodel$linear.predictors
covlist1=c('covA', 'covB', 'covC', 'covD', 'covE', 'covF.Middle', 'covF.High', 'Asqr','BC', 'BD')
k =rubinRules2(data=toy,Treatment='treated',covlist=covlist1)
k$plot


## ----modifyrubinRules,fig.width = 5, fig.asp =1,message=FALSE,warning=FALSE----
k$plot+theme_bw()+labs(title = "Assessing Rubin (2010) Balance Measures")

## ----individualRubinComponents-------------------------------------------
k$RUBIN3


## ----cardinalityMatch,warning=FALSE,message=FALSE------------------------

data("lalonde",package = "designmatch")
attach(lalonde)
## Treatment indicator
t_ind =lalonde$treat

## Distance matrix
dist_mat = NULL

## Subset matching weight.
subset_weight = 1 

# Moment balance: constrain differences in means to be at most .05 standard deviations apart
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
mom = list(covs = mom_covs, tols = mom_tols)

## Fine balance
#fine_covs = cbind(black, hispan, married, nodegree)
#fine = list(covs = fine_covs)

## Exact matching
#exact_covs = cbind(black)
#exact = list(covs = exact_covs)

## Solver options
t_max = 60*5
solver = "glpk"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)

## Cardinality matching
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, mom = mom,  solver = solver)

# Indices of the treated units and matched controls
t_id = out$t_id
c_id = out$c_id
#detach(lalonde)


## ----removelaldonde_01, echo=FALSE---------------------------------------
detach(lalonde)

## ----LoveplotDesignMatch,fig.width=6,fig.asp=1,message=FALSE,warning=FALSE----
p = love_plot(X =out, data = lalonde , covList=c("age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75"))

p

## ----modifyLoveplotDesignMatch,fig.width=6,fig.asp=1,warning=FALSE,message=FALSE----
p+theme_bw()+geom_vline(xintercept = -0.1)+ geom_vline(xintercept = 0.1) +theme(legend.title = element_blank())


## ----withMatchIt---------------------------------------------------------

library(Matching);library(MatchIt)
data("lalonde",package = "Matching")

## Sensitivity analysis with a matchit object
m.out = matchit(treat ~ age  + educ +  black + hisp +married + nodegr + re74  + re75  +
u74 + u75, family=binomial, data = lalonde, method = "optimal")

## Estimating treatment effect. Ideally, balance assessement should be done prior to estimating treatment effect
mod = lm(re78~age  + educ +  black + hisp +married + nodegr + re74  + re75  +u74 + u75,data = match.data(m.out))

## Sensitivity analysis
sens.out =pens2(x = m.out, y="re78",Gamma = 2, GammaInc = 0.1,est = 629.7)


## ----gammabounds---------------------------------------------------------
kable(sens.out$bounds)

## ----designCont,echo=FALSE,message=FALSE,warning=FALSE-------------------
data("lalonde",package = "designmatch")
attach(lalonde)
## Treatment indicator
t_ind = treatment
## Distance matrix
dist_mat = NULL
## Subset matching weight
subset_weight = 1
## Moment balance: constrain differences in means to be at most .05 standard deviations apart
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
mom = list(covs = mom_covs, tols = mom_tols)
## Fine balance
fine_covs = cbind(black, hispanic, married, nodegree)
fine = list(covs = fine_covs)
## Exact matching
exact_covs = cbind(black) 
exact = list(covs = exact_covs)
## Solver options
t_max = 60*5
solver = "glpk"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)
## Match
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,mom = mom, fine = fine, exact = exact, solver = solver)

## ----sensitiivtyCard-----------------------------------------------------
sens2 = pens2(x = out, y="re78",Gamma = 2, GammaInc = 0.1,est = 234,treat = "treatment",data = lalonde)
kable(sens2$bounds)

## ----removeLalonde01,echo=FALSE------------------------------------------
detach(lalonde)

## ----survSens,warning=FALSE----------------------------------------------
data("toy",package = "SensitivityR5")
psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
toy$ps <- psmodel$fitted
toy$linps <- psmodel$linear.predictors
X <- toy$linps ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
Y <- toy$out3.time
match1 <- Match(Y=Y, Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)

match.it <- matchit(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, data = toy, method='nearest', ratio=1) 

## ----timetovenSens,warning=FALSE-----------------------------------------
res_Surv =Survsens(x= match.it,data =toy,exp='treated',outcome = 'out2',failtime = 'out3.time',Gamma=1.2,Gammainterval = 0.01,alpha ,plot_title = 'Time To Event Outcome Sensitivity Plot')

kable(res_Surv$bounds)

## ----multSens------------------------------------------------------------

data("lalonde",package ="designmatch")
covs0 = c("age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75")
m.out <- matchit(f.build("treatment", covs0), data = lalonde, method = "nearest", replace = TRUE,ratio = 2)
res1 = multiControlSens(X =m.out,outcomeName = "re78",Gamma = 2,GammaInc = 0.1,n_contrl = 2)


## ----multiContrsENS,warning=FALSE----------------------------------------
kable(head(res1$data,5))


## ----pvalueMultiplecontrl------------------------------------------------
kable(res1$pvalues)

## ---- binsens------------------------------------------------------------
data("GerberGreenImai",package = "Matching")
## Sensitivity analysis with a Match object
m.out = matchit(PHN.C1 ~ PERSONS + VOTE96.1 + NEW +MAJORPTY + AGE + WARD , family=binomial, data = GerberGreenImai, method = "nearest")
mod = lm(VOTED98 ~ PHN.C1+PERSONS + VOTE96.1 + NEW +MAJORPTY + AGE + WARD,data = match.data(m.out))
binarysens2(x=m.out,y ="VOTED98", Gamma=2, GammaInc=.1)


## ----match1_matchit------------------------------------------------------
data("lalonde",package = "Matching")
m.out = matchit(treat ~ age  + educ +  black + hisp +married + nodegr + re74  + re75  +
u74 + u75, family=binomial, data = lalonde, method = "nearest")


## ----sens_cont-----------------------------------------------------------
mod = lm(re78~age  + educ +  black + hisp +married + nodegr + re74  + re75  +u74 + u75,data = match.data(m.out))

sens_object = pens2(x = m.out, y="re78",Gamma = 2, GammaInc = 0.1,est = 629.7)


## ----sens_cont_results---------------------------------------------------

kable(sens_object$bounds)

## ----designMatch---------------------------------------------------------

data("lalonde",package = "designmatch")
library(designmatch)
attach(lalonde)
## Treatment indicator
t_ind = treatment
## Distance matrix
dist_mat = NULL
## Subset matching weight
subset_weight = 1
## Moment balance: constrain differences in means to be at most .05 standard deviations apart
mom_covs = cbind(age, education, black, hispanic, married, nodegree, re74, re75)
mom_tols = round(absstddif(mom_covs, t_ind, .05), 2)
mom = list(covs = mom_covs, tols = mom_tols)
## Fine balance
fine_covs = cbind(black, hispanic, married, nodegree)
fine = list(covs = fine_covs)
## Exact matching
exact_covs = cbind(black)
exact = list(covs = exact_covs)
## Solver options
t_max = 60*5
solver = "glpk"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,round_cplex = 0, trace = 0)
## Match
out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,mom = mom, fine = fine, exact = exact, solver = solver)


## ----pens2_cardMatch-----------------------------------------------------

cardMatch = pens2(x = out, y="re78",Gamma = 2, GammaInc = 0.1,est = 234,treat = "treatment",data = lalonde)

kable(head(cardMatch$bounds,5))
detach(lalonde)


