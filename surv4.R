library(survival)
library(asaur)
library(lubridate)
library(survminer)

survtime <- c(6, 7, 10, 15, 19, 25)
censor <- c(1, 0, 1, 1, 0, 1)
group <- c( "C", "C", "T", "C", "T", "T")
dtfr <- data.frame(survtime, censor, group)
# skip the 7 and 19 as these are censored
# 4 two by two tables
# FALSE is the failures
# mean is number in control times total failures divided by total
# var is a bit more involved
# statistic is the totals from the survival tables
# as a chi-square with 1 df
with(dtfr, table(survtime >6, group))
with(dtfr[dtfr$survtime >= 10, ]
     , table(survtime >10, group))
with(dtfr[dtfr$survtime >= 15, ]
     , table(survtime >15, group))
with(dtfr[dtfr$survtime >= 25, ]
     , table(survtime >25, group))
# survival function works it out:
with(dtfr
     , survdiff(Surv(survtime
                     , censor) ~ group))

data("pancreatic")
attach(pancreatic)
head(pancreatic)
# convert the text dates into R dates
Progression <- mdy(progression)
OnStudy <- mdy(onstudy)
Death <- mdy(death)
progressionOnly <- Progression - OnStudy
overallSurvival <- Death - OnStudy
pfs <- pmin(progressionOnly, overallSurvival)
pfs[is.na(pfs)] <- overallSurvival[is.na(pfs)]
pfs
# convert to months
pfsmonth <- pfs/30.5
plot(survfit(Surv(pfsmonth)~stage)
     , ylab = "Survival Prob"
     , xlab = "Time in Months"
     , col=c("blue", "red")
     , lwd = 2)
legend("topright", legend=c("Locally advanced"
                            , "Metastatic")
       , col=c("blue", "red")
       , lwd = 2)
# log rank test - are the curves signif diff
survdiff(Surv(pfs) ~ stage, rho=0) # not signif
# log rank test - with Prentice modification
# gives more weight to the early events
survdiff(Surv(pfs) ~ stage, rho=1) # signif
# wilcox test
wilcox.test(as.numeric(pfs)[pancreatic$stage == "M"], as.numeric(pfs)[pancreatic$stage != "M"])
# signif resutls as is - closer to prentice modified.
survdiff(Surv(overallSurvival) ~ stage, rho=0) # not signif
# log rank test - with Prentice modification
# gives more weight to the early events
survdiff(Surv(overallSurvival) ~ stage, rho=1) # signif

attach(pharmacoSmoking)
plot(survfit(Surv(ttr, relapse)~grp)
     , ylab = "Relapse Prob"
     , xlab = "Time in Months"
     , col=c("blue", "red")
     , lwd = 2)
legend("topright", legend=c("Combo"
                            , "Patch Only")
       , col=c("blue", "red")
       , lwd = 2)

survdiff(Surv(ttr, relapse) ~ grp)
table(ageGroup2)
survdiff(Surv(ttr, relapse) ~ grp +
           strata(ageGroup2)) # stratify on age
# as the values are very similar
# we can say that it was not necessary to strat on age

lambda.mutant.0 <- 0.03
lambda.mutant.1 <- 0.03*0.55
lambda.wt.0 <- 0.03*0.2
lambda.wt.1 <- 0.03*0.2*0.55

set.seed(4321)
tt.control.mutant <- rexp(25
                          , rate=lambda.mutant.0)
tt.treat.mutant <- rexp(125
                        , rate=lambda.mutant.1)
tt.control.wt <- rexp(125
                      , rate=lambda.wt.0)
tt.treat.wt <- rexp(25
                    , rate=lambda.wt.1)
ttAll <- c(tt.control.mutant, tt.treat.mutant, tt.control.wt, tt.treat.wt)
status <- rep(1, length(ttAll))
genotype <- c(rep("mutant", 150), rep("wt", 150))
trt <- c(rep(0, 25), rep(1, 125), rep(0, 125), rep(1, 25))
survdiff(Surv(ttAll, status) ~ trt)
# because of presences of mutation, it appears as though survival is reduced by treatment
plot(survfit(Surv(ttAll, status)~trt)
     , ylab = "Survival Prob"
     , xlab = "Time in Months"
     , col=c("blue", "red")
     , lwd = 2)
legend("topright", legend=c("Control"
                            , "Treatment")
       , col=c("blue", "red")
       , lwd = 2)

plot(survfit(Surv(ttAll, status)~strata(genotype)+trt)
     , ylab = "Survival Prob"
     , xlab = "Time in Months"
     , col=c("blue", "red")
     , lwd = 2
     , lty = c(1, 1, 2, 2))
legend("topright", legend=c("Control Mute"
                            , "Treatment Mute"
                            , "Control WT"
                            , "Treatment WT")
       , col=c("blue", "red")
       , lwd = 2
       , lty = c(1, 1, 2, 2))


attach(pharmacoSmoking)
plot(survfit(Surv(ttr, relapse)~grp)
     , ylab = "Relapse Prob"
     , xlab = "Time in Months"
     , col=c("blue", "red")
     , lwd = 2)
legend("topright", legend=c("Combo"
                            , "Patch Only")
       , col=c("blue", "red")
       , lwd = 2)
survdiff(Surv(ttr, relapse) ~ grp, rho=0)
survdiff(Surv(ttr, relapse) ~ grp, rho=1)
# prentice modification seems to make little difference.
survdiff(Surv(ttr, relapse) ~ grp+strata(employment))
survfit(Surv(ttr, relapse) ~ grp+strata(employment))
summary(survfit(Surv(ttr, relapse) ~ grp+strata(employment)))
plot(survfit(Surv(ttr, relapse) ~ grp+strata(employment))
     , ylab = "Relapse Prob"
     , xlab = "Time in Months"
     , col=c("blue", "red")
     , lwd = 2
     , lty = 1:3)
legend("topright", legend=c("Combo FT"
                            , "Patch Only FT"
                            , "Combo other"
                            , "Patch Only other"
                            , "Combo PT"
                            , "Patch Only PT")
       , col=c("blue", "red")
       , lwd = 2
       , lty = 1:3)
ggsurvplot(fit = survfit(Surv(ttr, relapse) ~ strata(employment)+grp)
           , data = pharmacoSmoking
           , break.time.by = 12
           , color = "strata"
           , surv.scale = "percent"
           , censor = T
           , legend = "right"
           , xlab = "Time to relapse - months"
           , ylab = "% relapsed"
           , risk.table = TRUE
           , ggtheme = theme_bw())

View(pancreatic)
