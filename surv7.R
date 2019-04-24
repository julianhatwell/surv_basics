library(survival)
library(asaur)
library(survminer)
library(lubridate)

data("pharmacoSmoking")
pharm <- within(pharmacoSmoking, {
  priorAttemptsT <- priorAttempts;
  priorAttemptsT[priorAttempts > 20] <- 20
  })

result.0.coxph <- coxph(Surv(ttr, relapse) ~ 1, data = pharm)
rr.0 <- residuals(result.0.coxph, type="martingale") 

smoothSEcurve <- function(yy, xx) {
  # use after a call to "plot" 
  # fit a lowess curve and 95% confidence interval curve
  # make list of x values
  xx.list <- seq(min(xx), max(xx), length.out = 101)
  # Then fit loess function through the points (xx, yy) 
  # at the listed values
  yy.xx <- predict(loess(yy ~ xx)
                   , se=TRUE
                   , newdata=data.frame(xx=xx.list))
  
  lines(yy.xx$fit ~ xx.list, lwd=2)
  lines(yy.xx$fit -
          qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list
        , lty=2)
  lines(yy.xx$fit +
          qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list
        , lty=2) }

with(pharm, {
  plot(rr.0 ~ age, data = pharm);
  smoothSEcurve(rr.0, age)
  })
title("Martingale residuals\nversus age")

with(pharm, {
  logAge <- log(age);
  plot(rr.0 ~ logAge);
  smoothSEcurve(rr.0, logAge)
  })
title("Martingale residuals\nversus log age")

with(pharm, {
  plot(rr.0 ~ priorAttemptsT, data = pharm);
  smoothSEcurve(rr.0, priorAttemptsT)
})
title("Martingale residuals\nversus prior attempts")

with(pharm, {
  logpriorAttemptsT <- log(priorAttemptsT + 1);
  plot(rr.0 ~ logpriorAttemptsT);
  smoothSEcurve(rr.0, logpriorAttemptsT)
})
title("Martingale residuals\nversus log prior attempts")

with(pharm, {
  plot(rr.0 ~ longestNoSmoke, data = pharm);
  smoothSEcurve(rr.0, longestNoSmoke)
})
title("Martingale residuals\nversus longest no smoke")

with(pharm, {
  loglongestNoSmoke <- log(longestNoSmoke + 1);
  plot(rr.0 ~ loglongestNoSmoke);
  smoothSEcurve(rr.0, loglongestNoSmoke)
})
title("Martingale residuals\nversus log longest no smoke")

# log longest no smoke will give a linear predictor. The others not
pharm <- within(pharm, {
  logLongestNoSmoke <- log(longestNoSmoke + 1)
})

result.grp.coxph <- coxph(Surv(ttr, relapse) ~ grp, data = pharm)
result.step <- step(result.grp.coxph
                    , scope=list(upper=~ grp + gender +
                                   race + employment +
                                   yearsSmoking +
                                   levelSmoking + age +
                                   priorAttemptsT +
                                   logLongestNoSmoke
                                 , lower=~grp) )

rr.final <- residuals(result.step, type="martingale")
with(pharm, {
  plot(rr.final ~ age);
  smoothSEcurve(rr.final, age)
  })
title("Martingale residuals\nversus age")
with(pharm, {
  plot(rr.final ~ grp);
  })
title("Martingale residuals\nversus grp")
with(pharm, {
  plot(rr.final ~ employment);
})
title("Martingale residuals\nversus employment")

# case deletion residuals (df betas)
result.coxph <- coxph(Surv(ttr, relapse) ~ grp +
                        employment + age
                      , data = pharm)
coef.all <- result.coxph$coef[4]
coef.all 

n.obs <- length(pharm$ttr)
jkbeta.vec <- rep(NA, n.obs)
for (i in 1:n.obs) {
  tt.i <- pharm$ttr[-i]
  delta.i <- pharm$relapse[-i]
  grp.i <- pharm$grp[-i]
  employment.i <- pharm$employment[-i]
  age.i <- pharm$age[-i]
  result.coxph.i <- coxph(Surv(tt.i, delta.i) ~ grp.i +
                            employment.i + age.i
                          , data = pharm)
  coef.i <- result.coxph.i$coef[4]
  jkbeta.vec[i] <- (coef.all - coef.i)
  }
jkbeta.vec

index.obs <- 1:n.obs
plot(jkbeta.vec ~ index.obs
     , type="h"
     , xlab="Observation"
     , ylab="Change in coefficient for age"
     , cex.axis=1.3
     , cex.lab=1.3)
abline(h=0) 
identify(jkbeta.vec ~ index.obs)

# automated approximation
resid.dfbeta <- residuals(result.coxph, type="dfbeta")
plot(resid.dfbeta[,4] ~ index.obs
     , type="h"
     , xlab="Observation"
     , ylab="Change in coefficient")
abline(h=0)
identify(resid.dfbeta[,4] ~ index.obs) 

plot(resid.dfbeta[,3] ~ index.obs
     , type="h"
     , xlab="Observation"
     , ylab="Change in coefficient")
abline(h=0)

# standardised - dfbetas, not dfbeta
resid.dfbeta <- residuals(result.coxph, type="dfbeta")
plot(resid.dfbeta[,4] ~ index.obs
     , type="h"
     , xlab="Observation"
     , ylab="Change in coefficient")
abline(h=0)

# assess proportional odds assumption (should give parallel curves)
data("pancreatic")
attach(pancreatic)
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

# now to check prop hazards assump
result.surv.LA <- survfit(Surv(pfsmonth) ~ stage
                          , subset={stage == "LA"}
                          , data = pancreatic)
time.LA <- result.surv.LA$time
surv.LA <- result.surv.LA$surv
cloglog.LA <- log(-log(surv.LA))
logtime.LA <- log(time.LA)
result.surv.M <- survfit(Surv(pfsmonth) ~ stage
                         , subset={stage == "M"}
                         , data = pancreatic)
time.M <- result.surv.M$time
surv.M <- result.surv.M$surv
cloglog.M <- log(-log(surv.M))
logtime.M <- log(time.M)

plot(cloglog.LA ~ logtime.LA
     , type="s"
     , col="blue"
     , lwd=2)
lines(cloglog.M ~ logtime.M
      , col="red"
      , lwd=2
      , type="s")
legend("bottomright"
       , legend=c("Locally advanced"
                  , "Metastatic")
       , col=c("blue","red")
       , lwd=2)

# schonefeld resids
tt <- c(6, 7, 10, 15, 19, 25)
delta <- c(1, 0, 1, 1, 0, 1)
trt <- c(0, 0, 1, 0, 1, 1)
result.coxph <- coxph(Surv(tt, delta) ~ trt)
result.coxph$coef
residuals(result.coxph, type = "schoenfeld")

# standardising these resids
resid.unscaled <- residuals(result.coxph, type="schoenfeld")
resid.scaled <- resid.unscaled*result.coxph$var*sum(delta)
resid.unscaled
resid.scaled 
