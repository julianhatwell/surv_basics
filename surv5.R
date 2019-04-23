library(survival)
library(asaur)
library(survminer)
library(numDeriv)

plsimple <- function(beta) {
  psi <- exp(beta)
  result <- log(psi) - log(3*psi + 3) -
    log(3*psi + 1) - log(2*psi + 1)
  result
}
result <- optim(par = 0, fn = plsimple
                , method = "L-BFGS-B"
                , control = list(fnscale = -1)
                , lower = -3, upper = 1)
result$par

fn <- expression(beta - log((3 * exp(beta) + 3)) - 
                   log((3 * exp(beta) + 1)) - 
                   log((2 * exp(beta) + 1)))
grad <- D(fn, "beta")
hess <- D(grad, "beta") # Hessian
beta <- 0
grad <- eval(grad)
hess <- eval(hess)

curve(plsimple, from = -4, to = 1)
points(result$par, result$value
       , pch = 19
       , col = "blue")
segments(x0 = result$par, x1 = result$par
         , y0 = -5.5, y1 = result$value
         , lty = 2, col = "blue")
points(0, plsimple(0)
       , pch = 19
       , col = "red")
segments(x0 = 0, x1 = 0
         , y0 = -5.5, y1 = plsimple(0)
         , lty = 2, col = "red")
segments(x0 = 0, x1 = 1
         , y0 = plsimple(0)
         , y1 = plsimple(0) + grad
         , col = "red")
segments(x0 = 0, x1 = -1
         , y0 = plsimple(0)
         , y1 = plsimple(0) - grad
         , col = "red")

# data from ch 4
tt <- c(6, 7, 10, 15, 19, 25)
status <- c(1, 0, 1, 1, 0, 1) # censor
grp <- c( "C", "C", "T", "C", "T", "T")
dtfr <- data.frame(survtime, censor, group)

# Partial Likelihood Hypothesis Tests
result.cox <- coxph(Surv(tt, status) ~ grp)
summary(result.cox)

grad(func = plsimple, x=0)
hessian(func = plsimple, x=0)

# score test
Z <- grad(func = plsimple, x=0)^2 / 
  (-hessian(func = plsimple, x=0))
pchisq(Z, 1, lower.tail = FALSE)

# wald test
se <- sqrt(1/(-hessian(func = plsimple, x=0)))
Z <- result.cox$coefficients/se
2 * pnorm(abs(Z), lower.tail = FALSE)

# loglik
Z <- 2 * (result$value - plsimple(0))
pchisq(Z, 1, lower.tail = FALSE)

# ties
tt <- c(7, 6, 6, 5, 2, 4, 4, 1, 3, 1)
status <- c(0, 1, 0, 0, 1, 1, 1, 1, 0, 1)
grp <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1) 
result.coxph <- coxph(Surv(tt, status) ~ grp)
result.coxph$coef 
result.coxph <- coxph(Surv(tt, status) ~ grp
                      , ties="exact") # discrete ties, e.g. industrial failure, menstrual cycles
result.coxph$coef 
result.coxph <- coxph(Surv(tt, status) ~ grp
                      , ties="breslow")
result.coxph$coef 
result.coxph <- coxph(Surv(tt, status) ~ grp
                      , ties="efron") # default if absent
result.coxph$coef 

# left truncation - e.g if we have enrolled on study later than diagnosis
tt <- c(6, 7, 10, 15, 19, 25)
status <- c(1, 0, 1, 1, 0, 1)
grp <- c(0, 0, 1, 0, 1, 1)
backtime <-  c(-3, -11, -3, -7, -10, -5) 
coxph(Surv(tt, status) ~ grp) # ingnoring trunc
# problem if diagnosis was not controlled for between groups
tm.enter <- -backtime
tm.exit <- tt - backtime
coxph(Surv(tm.enter, tm.exit, status, type="counting") ~ grp)

data(ChanningHouse, package = "asaur")
Channing68 <- within(ChanningHouse, {
  entryYears <- entry/12
  exitYears <- exit/12
  timeYears <- time/12})

Channing68 <- Channing68[Channing68$exitYears >= 68,]

coxph(Surv(entryYears, exitYears, cens, type="counting") ~ sex
      , data=Channing68)

# exercises
aml_rel <- aml
aml_rel$x <- relevel(aml$x, "Nonmaintained")
result <- coxph(Surv(time, status) ~ x, data=aml_rel)
summary(result)
time.months <- cut(aml$time, breaks=seq(0,161,4), labels=F) 
result <- coxph(Surv(time.months, status) ~ x
                , data=aml_rel
                , ties = "exact")
summary(result)
result <- coxph(Surv(time.months, status) ~ x
                , data=aml_rel
                , ties = "efron")
summary(result)
result <- coxph(Surv(time.months, status) ~ x
                , data=aml_rel
                , ties = "breslow")
summary(result)

# data from ch 4
tt <- c(6, 7, 10, 15, 19, 25)
status <- c(1, 0, 1, 1, 0, 1) # censor
grp <- c( "C", "C", "T", "C", "T", "T")
dtfr <- data.frame(survtime, censor, group)

# Partial Likelihood Hypothesis Tests
result.cox <- coxph(Surv(tt, status) ~ grp)
cxph <- summary(result.cox)

bh <- basehaz(result.cox, centered = FALSE)
tt0 <- c(0, tt[-length(tt)])
tt.diff <- tt - tt0
tt.diff
surv.cont <- exp(-cumsum(bh$hazard*tt.diff))
surv.treat <- exp(-cumsum(bh$hazard*tt.diff*cxph$coefficients[2]))
matplot(c(0, tt), cbind(c(1, surv.cont), c(1, surv.treat))
        , xlab = "time", ylab = "Survival"
        , type = "l")
