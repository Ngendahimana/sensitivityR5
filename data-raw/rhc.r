## REQUIRED PACKAGES##
library(Epi); library(lme4); library(survival)
library(twang); library(tableone); library(survey)
library(Matching); library(Hmisc); library(rms)
library(MatchIt); library(CBPS); library(ebal)
library(cobalt); library(forcats); library(pander)
library(tidyverse);library(designmatch); library(Deducer);library(rbounds)
library(data.table);library(lattice);library(gdata);library(compareGroups);library(psych)
library(knitr)


## DATA CLEAN-UP   ##

column_types_rhc <-
  cols(urin1 = "d", meanbp1 = "d", resp1 = "d",
       swang1 = col_factor(c("RHC", "No RHC")),
       death = col_factor(c("No", "Yes")),
       sex = col_factor(c("Male", "Female")),
       cat1 = col_factor(c("ARF", "CHF", "Cirrhosis", "Colon Cancer", "Coma", "COPD",
                           "Lung Cancer", "MOSF w/Malignancy", "MOSF w/Sepsis")),
       dnr1 = col_factor(c("No", "Yes")),
       card = col_factor(c("No", "Yes")),
       gastr = col_factor(c("No", "Yes")),
       hema = col_factor(c("No", "Yes")),
       meta = col_factor(c("No", "Yes")),
       neuro = col_factor(c("No", "Yes")),
       ortho = col_factor(c("No", "Yes")),
       renal = col_factor(c("No", "Yes")),
       resp = col_factor(c("No", "Yes")),
       seps = col_factor(c("No", "Yes")),
       trauma = col_factor(c("No", "Yes")),
       income = col_factor(c("Under $11k", "$11-$25k", "$25-$50k", "> $50k")),
       ninsclas = col_factor(c("Private", "Private & Medicare", "Medicare",
                               "Medicare & Medicaid", "Medicaid", "No insurance")),
       race = col_factor(c("white", "black", "other")),
       ca = col_factor(c("No", "Yes", "Metastatic"))
  )

rhc.raw <- read_csv("data-raw/rhc.csv", col_types = column_types_rhc)

rhc <- rhc.raw %>% dplyr::rename(subject = X1) %>%
  select(subject, swang1,
         death, sadmdte, dthdte, lstctdte, dschdte,
         age, sex, edu, income, ninsclas, race,
         cat1, dnr1, wtkilo1, hrt1, meanbp1, resp1, temp1,
         card, gastr, hema, meta, neuro, ortho, renal, resp, seps, trauma,
         amihx, ca, cardiohx, chfhx, chrpulhx, dementhx, gibledhx,
         immunhx, liverhx, malighx, psychhx, renalhx, transhx,
         aps1, das2d3pc, scoma1, surv2md1,
         alb1, bili1, crea1, hema1, paco21, pafi1, ph1, pot1, sod1, urin1, wblc1,
         dth30, t3d30, cat2, adld3p, ptid)

rhc$swang1_f <- factor(rhc$swang1, levels=c("RHC", "No RHC"))
rhc$swang1 <- as.numeric(rhc$swang1_f == "RHC")

rhc$death_f <- rhc$death
rhc$death <- as.numeric(rhc$death_f == "Yes")


rhc$dschdte.orig= rhc$dschdte
rhc$dschdte <- impute(rhc$dschdte, 12411)

rhc$hospdays <- as.numeric(rhc$dschdte - rhc$sadmdte)

rhc$survdays <- ifelse(rhc$death == 1, rhc$dthdte - rhc$sadmdte, rhc$lstctdte - rhc$sadmdte)

rhc$urin1.NA <- as.numeric(is.na(rhc$urin1))
rhc$urin1.NA_f <- factor(rhc$urin1.NA == 1, levels = c(F, T), labels = c("No", "Yes"))

rhc <- rhc %>%
  mutate(dnr1_f = dnr1, card_f = card, gastr_f = gastr,
         hema_f = hema, meta_f = meta, neuro_f = neuro,
         ortho_f = ortho, renal_f = renal, resp_f = resp,
         seps_f = seps, trauma_f = trauma) %>%
  mutate(dnr1 = as.numeric(dnr1_f == "Yes"),
         card = as.numeric(card_f == "Yes"),
         gastr = as.numeric(gastr_f == "Yes"),
         hema = as.numeric(hema_f == "Yes"),
         meta = as.numeric(meta_f == "Yes"),
         neuro = as.numeric(neuro_f == "Yes"),
         ortho = as.numeric(ortho_f == "Yes"),
         renal = as.numeric(renal_f == "Yes"),
         resp = as.numeric(resp_f == "Yes"),
         seps = as.numeric(seps_f == "Yes"),
         trauma = as.numeric(trauma_f == "Yes"),
         female = as.numeric(sex == "Female")
  )


rhc <- rhc %>%
  mutate(amihx_f = factor(amihx == "1", levels = c(F, T), labels = c("No", "Yes")),
         cardiohx_f = factor(cardiohx == 1, levels = c(F, T), labels = c("No", "Yes")),
         chfhx_f = factor(chfhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         chrpulhx_f = factor(chrpulhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         dementhx_f = factor(dementhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         gibledhx_f = factor(gibledhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         immunhx_f = factor(immunhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         liverhx_f = factor(liverhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         malighx_f = factor(malighx == 1, levels = c(F, T), labels = c("No", "Yes")),
         psychhx_f = factor(psychhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         renalhx_f = factor(renalhx == 1, levels = c(F, T), labels = c("No", "Yes")),
         transhx_f = factor(transhx == 1, levels = c(F, T), labels = c("No", "Yes"))
  )

save(rhc, file = 'data/rhc.rdata', compress = 'xz')
