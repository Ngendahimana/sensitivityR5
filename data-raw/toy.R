## Toy example Script
## Libraries used: Epi, survival, sm, Matching, lme4, reshape2, epibasix, Hmisc, lattice and survey
## Assumes toy.csv has been imported (200 observations, 11 variables)

## Preliminary Data Cleanup

summary(toy)

## Re-expressing Binary Variables
toy$treated.f <- factor(toy$treated, levels=c(1,0),labels=c("Treated","Control"))
toy$covB.f <- factor(toy$covB, levels=c(1,0), labels=c("Has B", "No B"))
toy$out2.f <- factor(toy$out2.event, levels=c("Yes","No"), labels=c("Event Occurred", "No Event"))
toy$out2 <- as.numeric(toy$out2.event)-1 # subtracting 1 at the end changes the default 1/2 code to 0/1

## Sanity Checks
table(toy$treated.f, toy$treated)
table(toy$covB.f, toy$covB)
table(toy$out2.f, toy$out2.event)
table(toy$out2, toy$out2.event)
table(toy$out2, toy$out2.f)

## Re-expressing the Multi-Categorical Variable
toy$covF.Low <- as.numeric(toy$covF=="1-Low")
toy$covF.Middle <- as.numeric(toy$covF=="2-Middle")
toy$covF.High <- as.numeric(toy$covF=="3-High")

## Sanity Checks
table(toy$covF, toy$covF.Low)
table(toy$covF, toy$covF.Middle)
table(toy$covF, toy$covF.High)

## We have three transformations to execute for the covariates
## A squared, plus B-C and B-D interactions
## Must use cov B in a numeric (0,1) form to build product terms with C and D
toy$Asqr <- toy$covA^2
toy$BC <- toy$covB*toy$covC
toy$BD <- toy$covB*toy$covD

## Now, toy should contain 200 observations, 21 variables
names(toy)
dim(toy)

## Question A. Ignoring covariate information, what is the estimated effect
## ... on Outcome 1 [a continuous outcome]
by(toy$out1.cost, toy$treated.f, summary)
unadj.out1 <- lm(out1.cost ~ treated, data=toy)
summary(unadj.out1); confint(unadj.out1) ## provides treated effect and confidence interval estimates

## ... on Outcome 2 [a binary outcome]
table(toy$treated.f, toy$out2.f)
library(Epi)
twoby2(table(toy$treated.f, toy$out2.f)) ## provides risk difference and odds ratio estimates

## To get a logistic regression model version of this, we use
unadj.out2 <- glm(out2 ~ treated, data=toy, family=binomial())
summary(unadj.out2)
exp(coef(unadj.out2)) # produces odds ratio estimate
exp(confint(unadj.out2)) # produced 95% CI for odds ratio

## ... on Outcome 3 [a time-to-event outcome]
## Patients with out2.event=No are right-censored, those with out2.event=Yes have their times to event observed
## Fit a simple unadjusted Cox proportional hazards model
## predicting time to event (with event=Yes indicating non-censored cases) based on treatment group (treated)
library(survival)
unadj.out3 <- coxph(Surv(out3.time, out2.event=="Yes") ~ treated, data=toy)
summary(unadj.out3) ## exp(coef) section indicates relative risk estimate and 95% CI
## Check proportional hazards assumption, at least a little bit - would like this p value to be non-significant
cox.zph(unadj.out3)
plot(cox.zph(unadj.out3), var="treated")

## Question B. Fitting the Propensity Score Model, then plotting it simply
psmodel <- glm(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD, family=binomial(), data=toy)
toy$ps <- psmodel$fitted
toy$linps <- psmodel$linear.predictors

by(toy$ps, toy$treated.f, summary)

plot(toy$ps ~ toy$treated.f, ylab="Propensity for Treatment", xlab="")

## for a fancier plot to compare the distributions of PS across treatment groups, we could use
library(sm)
sm.density.compare(toy$ps, toy$treated.f, xlab="Propensity for Treatment", main="Propensity Score Comparison", xlim=c(0,1), col=c("red", "dark green"), lty=1)
legend("topleft", legend=levels(toy$treated.f), lty=c(1,2), lwd=1, col=c("red","dark green"), text.col=c("red","dark green"), bty="n")

## The approach above automatically places the legend on the top left of the plot. To place the legend
## yourself, comment out the legend line above, in favor of ...
## legend(locator(1), legend=levels(toy$treated.f), lty=c(1,2), lwd=1, col=c("red","dark green"), text.col=c("red","dark green"), bty="n")

## Checking Rubin's Rules Before Propensity Score Adjustment

## Rubin's Rule 1 - calculate the (absolute value of the) standardized difference of the linear propensity score
## We want this value to be close to 0, and certainly less than 50 in order to push forward to outcomes analysis
## Before Propensity Adjustment, Rubin's Rule 1 summary is...
rubin1.unadj <- with(toy, abs(100*(mean(linps[treated==1])-mean(linps[treated==0]))/sd(linps)))
rubin1.unadj

## Rubin's Rule 2 - calculate the ratio of variances of the linear propensity score comparing treated to control
## We want this value to be near 1, between 4/5 and 5/4, ideally, but certainly between 1/2 and 2.

## Prior to Propensity Score Adjustment, Rubin's Rule 2 summary is ...
rubin2.unadj <-with(toy, var(linps[treated==1])/var(linps[treated==0]))
rubin2.unadj

## Rubin's Rule 3 - ratio of variances of regression residuals for each covariate included in the propensity model
## comparing treated to control - again looking for near 1, between 4/5 and 5/4, ideally, definitely between 1/2 and 2.

## General function rubin3 to help calculate Rubin's Rule 3
rubin3 <- function(data, covlist, linps) {
  covlist2 <- as.matrix(covlist)
  res <- NA
  for(i in 1:ncol(covlist2)) {
    cov <- as.numeric(covlist2[,i])
    num <- var(resid(lm(cov ~ data$linps))[data$treated==1])
    den <- var(resid(lm(cov ~ data$linps))[data$treated==0])
    res[i] <- round(num/den, 3)
  }
  names(res) <- names(covlist)
  print(res)
}

## Prior to propensity adjustment, Rubin's Rule 3 summaries are:
rubin3.unadj <- rubin3(data=toy, covlist=toy[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])

## Building a dotplot of the Rubin's Rule 3 Variance Ratios:
d0 <- rubin3(data=toy, covlist=toy[c("covA", "covB", "covC", "covD", "covE", "covF.Middle", "covF.High", "Asqr","BC", "BD")])
d <- sort(d0)
low <- min(min(d), 0.45)
high <- max(max(d), 2.05)

dotchart(d, pch=15, col="black", main="Rubin's Rule 3 Results (Unadjusted)", xlab="Residual Variance Ratio", xlim=c(low, high))
abline(v=1, lty=1)
abline(v=0.8, lty=2, lwd=2, col="blue")
abline(v=1.25, lty=2, lwd=2, col="blue")
abline(v=0.5, lty=2, lwd=2, col="red")
abline(v=2, lty=2, lwd=2, col="red")

rm(d, d0, low, high)

## Question C. Use 1:1 greedy matching to match all 70 treated to 70 unique control patients
## on the linear propensity scores. We'll break ties at random, as well.
## Then check balance (and plot standardized differences and variance ratios) appropriately
## including the raw and linear propensity scores in plots (usually unnecessary to show both)
library(Matching)
X <- toy$linps ## matching on the linear propensity score
Tr <- as.logical(toy$treated)
match1 <- Match(Tr=Tr, X=X, M = 1, replace=FALSE, ties=FALSE)
summary(match1)

## Next, we need to assess balance imposed by match1 on covariates (and PS and the linear PS).
mb1 <- MatchBalance(treated ~ covA + covB + covC + covD + covE + covF + Asqr + BC + BD + ps + linps, data=toy, match.out = match1, nboots=500)
covnames <- c("covA", "covB", "covC", "covD", "covE", "covF - Middle", "covF - High", "A^2","B*C", "B*D", "raw PS", "linear PS")

## Extract standardized differences
pre.szd <- NULL; post.szd <- NULL
for(i in 1:length(covnames)) {
  pre.szd[i] <- mb1$BeforeMatching[[i]]$sdiff.pooled
  post.szd[i] <- mb1$AfterMatching[[i]]$sdiff.pooled
}

## Basic Standardized Difference Table
temp <- data.frame(pre.szd, post.szd, row.names=covnames)
print(temp, digits=3)

## Absolute Standardized Difference Plot
temp <- data.frame(pre.szd, post.szd, row.names=covnames)
tempsort <- temp[with(temp, order(abs(pre.szd))),]
high <- max(max(abs(pre.szd)), max(abs(post.szd)), 0.1)

dotchart(abs(tempsort$pre.szd), pch="", xlim=c(0, 1.05*high), labels=row.names(tempsort), main="Absolute Standardized Difference Plot", xlab="Absolute Standardized Difference (%)")
points(abs(tempsort$pre.szd), seq(1:length(tempsort$pre.szd)), pch=15, col="blue", cex=1.2)
points(abs(tempsort$post.szd), seq(1:length(tempsort$post.szd)), pch=19, col="red", cex=1.2)
abline(v=0, lty=1)
abline(v=10, lty=2, col="purple")
legend("bottomright", legend = c("Before Matching", "After Matching"), col=c("blue", "red"), text.col=c("blue", "red"), bty="o", pch = c(15, 19))

## Or, if you prefer, plotting the standardized differences (not absolute values)
temp <- data.frame(pre.szd, post.szd, row.names=covnames)
tempsort <- temp[with(temp, order(pre.szd)), ]
low <- min(min(pre.szd), min(post.szd), -0.1)
high <- max(max(pre.szd), max(post.szd), 0.1)

dotchart(tempsort$pre.szd, xlim=c(1.05*low, 1.05*high), pch="", labels=row.names(tempsort), main="Standardized Difference Plot", xlab="Standardized Difference (%)")
points(tempsort$pre.szd, seq(1:length(tempsort$pre.szd)), pch=15, col="blue", cex=1.2)
points(tempsort$post.szd, seq(1:length(tempsort$post.szd)), pch=19, col="red", cex=1.2)
abline(v=0, lty=1)
abline(v=10, lty=2, col="purple")
abline(v=-10, lty=2, col="purple")
legend("bottomright", legend = c("Before Matching", "After Matching"), col=c("blue", "red"), text.col=c("blue", "red"), bty="o", pch = c(15, 19))

## Extract variance ratios
pre.vratio <- NULL; post.vratio <- NULL
for(i in 1:length(covnames)) {
  pre.vratio[i] <- mb1$BeforeMatching[[i]]$var.ratio
  post.vratio[i] <- mb1$AfterMatching[[i]]$var.ratio
}

## Table of Variance Ratios
temp <- data.frame(pre.vratio, post.vratio, row.names=covnames)
print(temp, digits=2)

## Variance Ratio Plot
temp <- data.frame(pre.vratio, post.vratio, row.names=covnames)
tempsort <- temp[with(temp, order(pre.vratio)), ]
low <- min(min(pre.vratio), min(post.vratio))
high <- max(max(pre.vratio), max(post.vratio))

dotchart(tempsort$pre.vratio, xlim=c(0.95*low, 1.05*high), pch="", labels=row.names(tempsort), main="Plot of Variance Ratios", xlab="Treatment Variance / Control Variance")
points(tempsort$pre.vratio, seq(1:length(tempsort$pre.vratio)), pch=8, col="black", cex=1.2)
points(tempsort$post.vratio, seq(1:length(tempsort$post.vratio)), pch=7, col="magenta", cex=1.2)
abline(v=1, lty=1)
abline(v=3/4, lty=2, col="brown")
abline(v=4/3, lty=2, col="brown")
legend("topleft", legend = c("Before Matching", "After Matching"), col=c("black", "magenta"), text.col=c("black", "magenta"), bty="o", pch = c(8, 7))

## Finally, we'll create a new data frame, containing only the matched sample
matches <- factor(rep(match1$index.treated, 2))
toy.matchedsample <- cbind(matches, toy[c(match1$index.control, match1$index.treated),])

