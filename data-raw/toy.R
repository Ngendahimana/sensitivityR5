## Toy example Script
 library(Epi); library(survival); library(sm); library(Matching); library(lme4); library(reshape2);library(survey)
## Assumes toy.csv has been imported (200 observations, 11 variables)

## Preliminary Data Cleanup
toy = read.csv("data-raw/toy.csv")


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

save(toy, file = 'data/toy.rdata', compress = 'xz')

