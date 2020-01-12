library(survival)
library(muhaz)
library(bshazard)
library(kmconfband)
library(DTDA)

ex1.1 <- data.frame(YoE = c(1990, 1990, 1991, 1991, 1992)
                    , YoX = c(1995, 1995, 1995, 1994, 1993)
                    , delta = c(0, 0, 1, 1, 1))

survex1 <- with(ex1.1, Surv(YoX-YoE, delta))
ex1.km <- survfit(survex1 ~ 1, conf.type = "log-log")
ex1.km
summary(ex1.km)
plot(ex1.km)

ex1.fh <- survfit(survex1 ~ 1, conf.type ="log-log"
               , type = "fh")
summary(ex1.fh)
plot(ex1.fh)


tab1.1 <- data.frame(survtime = c(7,6,6,5,2,4)
                   , status = c(0,1,0,0,1,1))
survtab1 <- with(tab1.1, Surv(survtime, status))
tab1.km <- survfit(survtab1 ~ 1, conf.type = "log-log")
tab1.km
summary(tab1.km)
plot(tab1.km)

tab1.fh <- survfit(survtab1 ~ 1, conf.type ="log-log"
                  , type = "fh")
summary(tab1.fh)
plot(tab1.fh)

timeMonths <- gastricXelox$timeWeeks*7/30.25
delta <- gastricXelox$delta
result.km <- survfit(Surv(timeMonths, delta) ~ 1
                     , conf.type="log-log")
# finding median
result.km$n
med_cl <- quantile(result.km, 0.5)
medn <- med_cl$`quantile`
lcl <- med_cl$lower
ucl <- med_cl$upper

plot(result.km, conf.int=T, mark="|", xlab="Time in months"
     , ylab="Survival probability")
title("Progression-free Survival in Gastric Cancer Patients")
abline(h = 0.5, col = "red", lty = 3)
lines(rep(x = ucl, 2), y = c(0, 0.5), col = "blue", lty = 2)
lines(rep(x = lcl, 2), y = c(0, 0.5), col = "blue", lty = 2)
lines(rep(x = medn, 2), y = c(0, 0.5), col = "green", lty = 2)

# median follow up-time (how long would each subject have been observed if they had survived)
median(timeMonths) # problem is that a lot of early deaths shorten the median
# better to do this by reversing the deltas, treat delta as a censoring event
# effectively the potential amount of total patient study time
delta.followup <- 1 - delta
survfit(Surv(timeMonths, delta.followup) ~ 1)

# kernel function for smoothing
# Epanechnikov kernel
K <- function(u) {
  return(ifelse(u < -1 | (u > 1), 0
                , (3/4) * (1-u^2) ))
}
x <- seq(-2, 2, 0.1)
plot(K(x)~x, type = "l")

hsmooth <- function(b=1, t.vec, cens.vec, t.max=max(t.vec)) {
  t.seq <- seq(0, t.max, length.out = 101)
  cens.vec <- cens.vec[order(t.vec)]
  t.vec <- sort(t.vec)
  d <- cumsum(cens.vec)
  n <- length(t.vec) - d
  kern <- (1/(b * length(t.vec))) * sapply(t.vec, function(x) {K((t.seq-x)/b)}) %*% (d/n)
  return(list(est.grid = t.seq
              , haz.est = kern))
}

t.vec <- c(7,6,6,5,2,4)
cens.vec <- c(0,1,0,0,1,1)
result.man <- hsmooth(b = 2.25, t.vec = t.vec, cens.vec = cens.vec, t.max = 8)
plot(result.man$haz.est~result.man$est.grid, type = "l")

# muhaz library does the same thing
result.simple <- muhaz(t.vec, cens.vec, max.time=8
                       , bw.grid=2.25, bw.method="global", b.cor="none")
plot(result.simple)

# estimating the hazard function using d_i/n_i
# using the gastic data, blocks of 5 months and 1 month superimposed
result.pe5 <- pehaz(timeMonths, delta, width=5, max.time=20)
plot(result.pe5, ylim=c(0,0.15), col="black", lty=3)
result.pe1 <- pehaz(timeMonths, delta, width=1, max.time=20)
lines(result.pe1)
# show the smooth hazard fuction using the kernel method
result.smooth <- muhaz(timeMonths, delta, bw.smooth=10# , bw.grid=10
                       , b.cor="left", max.time=20)
lines(result.smooth, col = "red")
# automatic smoothing
result.smooth2 <- muhaz(timeMonths, delta
                       , b.cor="left", max.time=20
                       , bw.method="local")
lines(result.smooth2, col = "green")

lines(bshazard(Surv(timeMonths, delta) ~ 1, lambda = 100)) # object with fitted values

# getting smooth survival function from the smoothed hazard
haz <- result.smooth$haz.est
times <- result.smooth$est.grid
# numerically evaluate the integral. diff creates the widths
# cumsum adds up the areas of the resulting rectangles
surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))

result.km <- survfit(Surv(timeMonths, delta) ~ 1,
                     conf.type="none")
plot(result.km, conf.int=T, mark="|", xlab="Time in months",
     xlim=c(0,30), ylab="Survival probability")
lines(surv ~ times[1:(length(times) - 1)])

cf <- confband(result.km, 0.1)
lines(result.km$time[1:(nrow(cf))], cf[,1])
lines(result.km$time[1:(nrow(cf))], cf[,2])

# left truncation
ltrunc <- data.frame(diag = -c(2, 5, 3, 3, 2, 5, 4) # when was the diagnosis relative to trial start
                     , survtime = c(7, 6, 6, 5, 2, 4, -2)
                     , status = c(0, 1, 0, 0, 1, 1, 1))

# remove 'patient X' who never entered the trial (died before beginning)
ltrunc <- ltrunc[-7, ]
tm.enter <- -ltrunc$diag
tm.exit <- ltrunc$survtime - ltrunc$diag
ltrunc.surv <- Surv(tm.enter, tm.exit
                    , ltrunc$status
                    , type="counting")
result.ltrunc.km <- survfit(ltrunc.surv ~ 1, conf.type="none")
summary(result.ltrunc.km)

result.ltrunc.fh <- survfit(ltrunc.surv ~ 1, conf.type="none", type="fh")
summary(result.ltrunc.fh)

# data example of left trunc issues - if risk set falls to zero
# km cum.product will drop to zero and remain at zero
# NAA (fh) will drop dramatically and can't recover (survival is monotonic)
# only solution is to start after the problematic risk set
data(ChanningHouse, package = "asaur")
ChanningMales <- within(ChanningHouse[ChanningHouse$sex == "Male",], {
  entryYears <- entry/12
  exitYears <- exit/12})
chan <- with(ChanningMales, Surv(entryYears, exitYears, cens,
             type="counting"))
result.km <- survfit(chan ~ 1)
plot(result.km, xlim=c(64, 101), xlab="Age",
     ylab="Survival probability", conf.int=F)
result.naa <- survfit(chan ~ 1
                      , type="fleming-harrington")
lines(result.naa, col="blue", conf.int=F)
result.km.68 <- survfit(chan ~ 1
                        , start.time=68)
lines(result.km.68, col="green", conf.int=F)
legend("topright", legend=c("KM", "NAA", "KM 68 and older")
       , lty=1, col=c("black", "blue", "green"))

# ex
# find the median
tab1.1 <- data.frame(survtime = c(7,6,6,5,2,4)
                     , status = c(0,1,0,0,1,1))
survtab1 <- with(tab1.1, Surv(survtime, status))
tab1.km <- survfit(survtab1 ~ 1, conf.type = "log-log")
tab1.km
summary(tab1.km)
plot(tab1.km)

med_cl <- quantile(tab1.km, 0.5)
medn <- med_cl$`quantile`
lcl <- med_cl$lower
ucl <- med_cl$upper # undefined because there is no data after

timeMonths <- gastricXelox$timeWeeks*7/30.25
delta <- gastricXelox$delta
result.km <- survfit(Surv(timeMonths, delta) ~ 1
                     , conf.type="log-log")
plot(result.km)
# finding quartiles
med_cl <- quantile(result.km, 0.25)
medn <- med_cl$`quantile`
lcl <- med_cl$lower
ucl <- med_cl$upper
medn; lcl; ucl

med_cl <- quantile(result.km, 0.75)
medn <- med_cl$`quantile`
lcl <- med_cl$lower # undefined because the variance is zero
ucl <- med_cl$upper # undefined because the variance is zero
medn; lcl; ucl

# show the smooth hazard fuction using the kernel method
result.smooth <- muhaz(timeMonths, delta, bw.grid=20
                       , b.cor="left")
plot(result.smooth, col = "red")

# ignoring ltrunc
data(ChanningHouse, package = "asaur")
ChanningMales <- within(ChanningHouse[ChanningHouse$sex == "Male",], {
  entryYears <- entry/12
  exitYears <- exit/12
  timeYears <- time/12})

chan <- with(ChanningMales, Surv(timeYears, cens))
result.km <- survfit(chan ~ 1)
plot(result.km, xlab="Years in Channing House",
     ylab="Survival probability", conf.int=F)
result.naa <- survfit(chan ~ 1
                      , type="fleming-harrington"
                      , data=ChanningMales)
lines(result.naa, col="blue", conf.int=F)
chan68 <- with(ChanningMales[ChanningMales$entryYears >= 68, ], Surv(timeYears, cens))
result.km.68 <- survfit(chan68 ~ 1)
lines(result.km.68, col="green", conf.int=F)
legend("topright", legend=c("KM", "NAA", "KM 68 and older")
       , lty=1, col=c("black", "blue", "green"))
# the data is length biased because...
# the sample are treated as if entering the trial with the same condition
# i.e. there is no accounting for their different ages when they enter
# the trial does not include anyone who died young
