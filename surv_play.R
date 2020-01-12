library(muhaz)
library(bshazard)
library(kmconfband)
library(DTDA)
source("surv1and2_alterntive.R")
# create a sample data set

set.seed(13132)
alpha <- 0.9
lambda <- 1.1
h <- get_weibull_hazard(alpha = alpha, lambda = lambda)

# CDF(t) is 1-S(t)
# theoretical survfunc
S <- function(t) 1 - pweibull(t, shape = alpha, scale = 1/lambda)
S((0:100) * 0.01)

initpop <- 100
pop <- initpop
t <- 0
d <- list()
while(pop > 0) {
  num_d <- sum(runif(pop) > S(t/30.425))
  d[[t + 1]] <- if(num_d > 0) {
    rep(t, num_d)
  } else {0}
  pop <- pop - num_d
  t <- t + 1
}
unlist(d)
wholeyears <- t%/%12
furthermonths <- t%%12
if(furthermonths > 0) wholeyears <- wholeyears + 1
studymonths <- wholeyears * 12
d <- unlist(d)
d <- d[d != 0]
cens <- rbinom(d, 1, 0.5)
fit <- survfit(Surv(unlist(d), cens) ~ 1, conf.type = "log-log")
summary(fit)
plot(fit, mark.time = TRUE)

fit$n
med_cl <- quantile(fit, 0.5)
medn <- med_cl$`quantile`
lcl <- med_cl$lower
ucl <- med_cl$upper

plot(fit, conf.int=T, mark="|", xlab="Time in months"
     , ylab="Survival probability")
title("Survival in Synthetic Cancer Patients")
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

# extracting the hazard function
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
result.man <- hsmooth(b = 3, t.vec = d, cens.vec = cens, t.max = 10)
plot(result.man$haz.est~result.man$est.grid, type = "l")

# muhaz library does the same thing
result.simple <- muhaz(d, cens, bw.grid=3, bw.method="global", b.cor="none", max.time = 10)
plot(result.simple)

# estimating the hazard function using d_i/n_i
# using the gastic data, blocks of 5 months and 1 month superimposed
result.pe5 <- pehaz(d, cens, width=3, max.time=10)
plot(result.pe5, ylim=c(0,0.5), col="black", lty=3)
# show the smooth hazard fuction using the kernel method
result.smooth <- muhaz(d, cens, bw.smooth=3# , bw.grid=10
                       , b.cor="left", max.time=10)
lines(result.smooth, col = "red")

result.smooth2 <- muhaz(d, cens
                        , b.cor="left", max.time=10
                        , bw.method="local")
lines(result.smooth2, col = "green")

lines(bshazard(Surv(d, cens) ~ 1, lambda = 100)) # object with fitted values

# getting smooth survival function from the smoothed hazard
haz <- result.smooth$haz.est
times <- result.smooth$est.grid
# numerically evaluate the integral. diff creates the widths
# cumsum adds up the areas of the resulting rectangles
surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
result.km <- survfit(Surv(d, cens) ~ 1,
                     conf.type="log-log")
plot(result.km, conf.int=T, mark="|", xlab="Time in months",
     xlim=c(0,20), ylab="Survival probability")
lines(surv ~ times[1:(length(times) - 1)])

cf <- confband(result.km, 0.1)
lines(result.km$time[1:(nrow(cf))], cf[,1])
lines(result.km$time[1:(nrow(cf))], cf[,2])
