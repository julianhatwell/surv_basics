library(survival)
library(asaur)
library(survminer)
library(forestplot)

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

ttAll.df <- data.frame(ttAll, status, genotype, trt)
ggsurvplot(survfit(Surv(ttAll, status)~strata(genotype)+trt)
           , data = ttAll.df
           , color = "strata")

coxph(Surv(ttAll, status)~trt)
coxph(Surv(ttAll, status)~trt+strata(genotype))
coxph(Surv(ttAll, status)~trt+genotype)

race <- factor(c("black", "black"
                 , "white", "white"
                 , "other", "other"))
age <- c(48, 52, 87, 82, 67, 53)
model.matrix(~ race + age)[,-1]
model.matrix(~ race + age + race:age)[,-1]

age <- runif(n=60, min=40, max=80)
race <- factor(c(rep("white", 20)
                 , rep("black", 20)
                 , rep("other", 20)))
race <- relevel(race, ref="white") 
log.rate.vec <- -4.5 + 
  c(rep(0,20), rep(1,20), rep(2,20)) +
  age * 0.05 
tt <- rexp(n=60, rate=exp(log.rate.vec))
status <- rep(1, 60)
model.matrix(~ race + age + race:age)[,-1]

result.cox <- coxph(Surv(tt, status) ~ race + age)
summary(result.cox)

data("pharmacoSmoking")
levels(pharmacoSmoking$ageGroup4)
levels(pharmacoSmoking$employment)
modelA.coxph <- coxph(Surv(ttr, relapse) ~ ageGroup4
                      , data = pharmacoSmoking)
modelA.coxph
summary(modelA.coxph)

modelB.coxph <- coxph(Surv(ttr, relapse) ~ employment
                      , data = pharmacoSmoking)
modelB.coxph
summary(modelB.coxph)
modelC.coxph <- coxph(Surv(ttr, relapse) ~ ageGroup4 +
                        employment
                      , data = pharmacoSmoking)
modelC.coxph

logLik(modelA.coxph)
logLik(modelB.coxph)
logLik(modelC.coxph)

lik.ratio.test <- function(a, b) as.numeric(2 * (logLik(a) - logLik(b)))
lik.ratio.test(modelC.coxph, modelA.coxph)
pchisq(lik.ratio.test(modelC.coxph, modelA.coxph)
       , df=attr(logLik(modelC.coxph), "df") - attr(logLik(modelA.coxph), "df")
       , lower.tail=FALSE)
lik.ratio.test(modelC.coxph, modelB.coxph)
pchisq(lik.ratio.test(modelC.coxph, modelB.coxph)
       , df=attr(logLik(modelC.coxph), "df") - attr(logLik(modelB.coxph), "df")
       , lower.tail=FALSE)

# need a test of whether age group belongs in the model
model.null.coxph <- coxph(Surv(ttr, relapse) ~ 1
                          , data=pharmacoSmoking)
logLik(model.null.coxph)
lik.ratio.test(modelA.coxph, model.null.coxph)
pchisq(lik.ratio.test(modelA.coxph, model.null.coxph)
       , df=attr(logLik(modelA.coxph), "df") - attr(logLik(model.null.coxph), "df")
       , lower.tail=FALSE)
# result identical to summary from model a alone.
# another way to do a lik.ratio test is to subtract the values given in the summaries.
# best of all is good old anova

anova(modelC.coxph, modelA.coxph)
anova(modelC.coxph, modelB.coxph)

# also the AIC
AIC(modelA.coxph)
AIC(modelB.coxph)
AIC(modelC.coxph)

# also the BIC - favours fewer params
BIC(modelA.coxph)
BIC(modelB.coxph)
BIC(modelC.coxph)


# a backwards stepwise procedure
modelAll.coxph <- coxph(Surv(ttr, relapse) ~ grp + gender +
                          race + employment + yearsSmoking +
                          levelSmoking + ageGroup4 +
                          priorAttempts + longestNoSmoke
                        , data = pharmacoSmoking)
result.step <- step(modelAll.coxph
                    , scope=list(upper =~ grp + gender +
                                   race + employment +
                                   yearsSmoking +
                                   levelSmoking +
                                   ageGroup4 +
                                   priorAttempts +
                                   longestNoSmoke
                                 , lower=~grp))
result.step


trt.f <- factor(veteran$trt, labels=c("standard", "test"))
result <- coxph(Surv(time, status) ~ trt.f + celltype
                , data=veteran)
result

coef.est <- c(NA, NA, 0, 0.198
              , NA, NA, NA, 0
              , 1.096, 1.169, 0.297)
se.est <- c(NA, NA, 0, 0.197
            , NA, NA, NA, 0
            , 0.272, 0.295, 0.286)
lower <- coef.est - 1.96*se.est
upper <- coef.est + 1.96*se.est
label.factors <- matrix(c("Treatment Group"
                          , ""
                          , " standard"
                          , " test"
                          , ""
                          , "Cell Type"
                          , ""
                          , " sqamous"
                          , " smallcell"
                          , " adeno"
                          , " large")
                        , ncol=1) 

forestplot(label.factors
           , coef.est
           , lower=lower
           , upper=upper
           , boxsize=0.4
           , xticks=c(-0.5,0,0.5, 1, 1.5, 2)
           , txt_gp=fpTxtGp(label=gpar(cex=1.5))) 

label.factors <- matrix(c("Treatment"
                          , " combination"
                          , " patch only"
                          , ""
                          , "Employment"
                          , " ft"
                          , " pt"
                          , " other"
                          , ""
                          , "Age Group"
                          , " 21-34"
                          , " 35-49"
                          , " 50-64"
                          , " 65+")
                        , ncol=1) 

coef.est <- c(NA, 0, 0.6564, NA
              , NA, 0, 0.5214, 0.6231, NA
              , NA , 0
              , -0.1119, -1.0233, -0.7071)
se.est <- c(NA, 0, 0.2198, NA
            , NA, 0, 0.3320, 0.2764, NA
            , NA , 0
            , 0.3216, 0.3597, 0.5017)
table.text <- cbind(label.factors, coef.est)
lower <- coef.est - 1.96*se.est
upper <- coef.est + 1.96*se.est
forestplot(table.text
           , coef.est
           , lower=lower
           , upper=upper
           , boxsize=0.4
           , xticks=c(-1,-0.5,0,0.5, 1)
           , clip = c(-1, 1)
           , txt_gp=fpTxtGp(label=gpar(cex=1.5))) 

modelS4.coxph <- coxph(Surv(ttr, relapse) ~ grp +
                         employment +
                         pspline(age, df=4)
                       , data = pharmacoSmoking)
modelS4.coxph # non linear component is not signif
termplot(modelS4.coxph, se=T, terms=3, ylabs="Log hazard") 
# also
termplot(modelS4.coxph, se=T, terms=1, ylabs="Log hazard") 
termplot(modelS4.coxph, se=T, terms=2, ylabs="Log hazard") 

data("hepatoCellular")
smodT <- coxph(Surv(OS, Death) ~ CXCL17T
               , data = hepatoCellular)
smodT
smodP <- coxph(Surv(OS, Death) ~ CXCL17P
               , data = hepatoCellular)
smodP
smodN <- coxph(Surv(OS, Death) ~ CXCL17N
               , data = hepatoCellular)
smodN

AIC(smodT)
AIC(smodP)
AIC(smodN)

smodPsp <- coxph(Surv(OS, Death) ~ pspline(CXCL17P)
               , data = hepatoCellular)
smodPsp
termplot(smodPsp, se=T, terms=1, ylabs="Log hazard") 

smodT <- coxph(Surv(RFS, Recurrence) ~ CXCL17T
               , data = hepatoCellular)
smodT
smodP <- coxph(Surv(RFS, Recurrence) ~ CXCL17P
               , data = hepatoCellular)
smodP
smodN <- coxph(Surv(RFS, Recurrence) ~ CXCL17N
               , data = hepatoCellular)
smodN

AIC(smodT)
AIC(smodP)
AIC(smodN)

smodPsp <- coxph(Surv(RFS, Recurrence) ~ pspline(CXCL17P)
                 , data = hepatoCellular)
smodPsp
termplot(smodPsp, se=T, terms=1, ylabs="Log hazard") 
# worth it to discretize P into 3 parts?

# a backwards stepwise procedure
modelAll.coxph <- coxph(Surv(OS, Death) ~ Age + Gender +
                          HBsAg + Cirrhosis + ALT +
                          AST + AFP + Tumorsize +
                          Tumordifferentiation +
                          Vascularinvasion + 
                          Tumormultiplicity +
                          Capsulation + TNM + BCLC +
                          CXCL17T + CXCL17P + CXCL17N
                        , data = hepatoCellular)
result.step.death <- step(modelAll.coxph
                    , scope=list(upper =~ CXCL17P + CXCL17N +
                                   CXCL17T + Age + Gender +
                                   HBsAg + Cirrhosis + ALT +
                                   AST + AFP + Tumorsize +
                                 Tumordifferentiation +
                                   Vascularinvasion + 
                                   Tumormultiplicity +
                                   Capsulation + TNM + BCLC
                                 , lower=~CXCL17P))

modelAll.coxph <- coxph(Surv(RFS, Recurrence) ~ Age + Gender +
                          HBsAg + Cirrhosis + ALT +
                          AST + AFP + Tumorsize +
                          Tumordifferentiation +
                          Vascularinvasion + 
                          Tumormultiplicity +
                          Capsulation + TNM + BCLC +
                          CXCL17T + CXCL17P + CXCL17N
                        , data = hepatoCellular)
result.step.recur <- step(modelAll.coxph
                    , scope=list(upper =~ CXCL17P + CXCL17N +
                                   CXCL17T + Age + Gender +
                                   HBsAg + Cirrhosis + ALT +
                                   AST + AFP + Tumorsize +
                                   Tumordifferentiation +
                                   Vascularinvasion + 
                                   Tumormultiplicity +
                                   Capsulation + TNM + BCLC
                                 , lower=~CXCL17P))
result.step.death
AIC(result.step.death)
BIC(result.step.death)

result.step.recur
AIC(result.step.recur)
BIC(result.step.recur)
