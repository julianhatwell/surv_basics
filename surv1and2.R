library(survival)
library(asaur)

# parametric - exponential dist
# instantaneous hazard
set_lambda <- function(lambda) {
  h <- function(x) lambda
  return(h)
}
h <- set_lambda(0.5)
h(1) # returns vector of same length
# cumulative hazard
set_hazard <- function(h) {
  H <- function(t) {
    integrand <- function(x) {sapply(x, h)}
    result <- integrate(integrand, 0, t)$value
    return(result)
  }
  return(H)
}
# for exponential dist, H(t) = lambda * t: area under rectangle
H <- set_hazard(h)
H(10)
# Survival is related to cum. hazard by exp(-H(t))
# for exponential, this is -lambda * t
S <- function(t) {exp(-H(t))}
# for exponention dist is it exp(-lambda * t)
S(0.4)
# probability density function is S(t) * h(t)
f <- function(t) {S(t) * h(t)}
f(0.4)

plot_h <- function(x) {sapply(x, h)}
plot_H <- function(x) {sapply(x, H)}
plot_S <- function(x) {sapply(x, S)}
plot_f <- function(x) {sapply(x, f)}

curve(plot_H(x)
      , from=0, to=5, ylab="Hazard", xlab="Time"
      , col="black") 
curve(plot_h(x)
      , from=0, to=5
      , col="red", add = T)
curve(plot_S(x)
      , from=0, to=5
      , col="blue", add=T)
curve(plot_f(x)
      , from=0, to=5
      , col="green", add=T)

# mean: integrate S() or 1/lambda
t_mu <- function(survfunc) {integrate(survfunc, 0, Inf)$value}
t_mu(plot_S) # one day I'll make this better
# t_mu <- 1/h(1) # t is arbitrary here

# median: value at which survival is 50%
# for exponential, Surv is exp(-lambda * t)
# exp(-lambda * t) = 0.5
# exp(lambda * t) = 2  (change sign of exponent = reciprocal)
# lambda * t = log(2)
# t_med = log(2)/lambda
t_med <- log(2)/h(1) # t is arbitrary here
t_med

# weibull
set_alpha_lambda <- function(alpha, lambda) {
  h <- function(t) {
    alpha * lambda^alpha * t^(alpha-1)
  }
  return(h)
}
h <- set_alpha_lambda(0.9, 0.5)
h(0.5)
# by integrating the instantaneous hazard
H <- set_hazard(h)
H(1)
# by simple formula
set_hazard_wf <- function(alpha, lambda) {
  H <- function(t) {
    (lambda * t)^alpha
  }  
  return(H)
}
H <- set_hazard_wf(1.1, 0.5)
H(1)
# deriving the formula
# int alpha * lambda * (lambda * t)^(alpha-1) dt
# int alpha * lambda * lambda^(alpha-1) * t^(alpha-1) dt
# int alpha * lambda^alpha * t^(alpha-1) dt
# alpha * lambda^alpha * int t^(alpha-1) dt
# alpha * lambda^alpha * t^alpha / alpha
# lambda^alpha * t^alpha
# (t * lambda)^alpha ( + C)

S <- function(t) {exp(-H(t))}

f <- function(t) {S(t) * h(t)}
f(0.5)

plot_h <- function(x) {sapply(x, h)}
plot_H <- function(x) {sapply(x, H)}
plot_S <- function(x) {sapply(x, S)}
plot_f <- function(x) {sapply(x, f)}

curve(plot_H(x)
      , from=0, to=3, ylab="Hazard", xlab="Time"
      , col="black") 
curve(plot_h(x)
      , from=0, to=3
      , col="red", add = T)
curve(plot_S(x)
      , from=0, to=3
      , col="blue", add=T)
curve(plot_f(x)
      , from=0, to=3
      , col="green", add=T)

# by integration
t_mu(plot_S)
# by formula for weibull
t_mu <- function(alpha, gamma) {
  gamma(1 + 1/alpha) / gamma # don't be confused by gamma function
}
t_mu(1.1, 0.5)
t_med <- function (alpha, gamma) {
  log(2)^(1/alpha) / gamma
}
t_med(1.1, 0.5)

# using the built in weibull dist
weibSurv <- function(t, shape, scale) {
  pweibull(t, shape=shape, scale=scale, lower.tail=F) 
}


curve(weibSurv(x, shape=1.5, scale=1/0.03)
      , from=0, to=80, ylim=c(0,1)
      , ylab="Survival probability", xlab="Time") 


weibHaz <- function(x, shape, scale) {
  dweibull(x, shape=shape, scale=scale) / 
    pweibull(x, shape=shape, scale=scale, lower.tail=F)
}

curve(weibHaz(x, shape=1.5, scale=1/0.03)
      , from=0, to=80, ylab="Hazard", xlab="Time"
      , col="red")
curve(weibHaz(x, shape=1, scale=1/0.03)
      , from=0, to=80, ylab="Hazard", xlab="Time"
      , col="black", add=T) 
curve(weibHaz(x, shape=0.75, scale=1/0.03)
      , from=0, to=80, ylab="Hazard", xlab="Time"
      , col="blue", add=T)


data("gastricXelox")
View(gastricXelox)
sum(gastricXelox$delta) # 32 death or progression
sum(gastricXelox$timeWeeks) # 2866 person weeks
sum(gastricXelox$delta) / sum(gastricXelox$timeWeeks)
# 0.01116539 events per person week

ex1.1 <- data.frame(YoE = c(1990, 1990, 1991, 1991, 1992)
                    , YoX = c(1995, 1995, 1995, 1994, 1993)
                    , delta = c(0, 0, 1, 1, 1))

personYears <- with(ex1.1, sum(YoX-YoE))
deaths <- with(ex1.1, sum(delta))
deathRate <- deaths/personYears # this also the max. likelihood estim for exponential dist.
var_lambda <- deaths/personYears^2 # this also the max. likelihood estim of the variance

deaths <- with(gastricXelox, sum(delta))
personWeeks <- with(gastricXelox, sum(timeWeeks))
deathRate <- deaths/personWeeks
var_lambda <- deaths/personWeeks^2 # 

weibSurv <- function(t, shape, scale) {
  pweibull(t, shape=shape, scale=scale, lower.tail=F) 
}

curve(weibSurv(x, shape=1.5, scale=1/0.03)
      , from=0, to=80, ylim=c(0,1)
      , ylab="Survival probability", xlab="Time") 


weibHaz <- function(x, shape, scale) {
  dweibull(x, shape=shape, scale=scale) / 
    pweibull(x, shape=shape, scale=scale, lower.tail=F)
}

shape <- 1.5; scale <- 1/0.03

# empirical mean of weib
tt.weib <- rweibull(1000, shape, scale)
mean(tt.weib)
# theoretical mean of weib
gamma(1 + 1/shape) * scale # 1/gamma

# empirical med of weib
tt.weib <- rweibull(1000, shape, scale)
median(tt.weib)
# theoretical mean of weib
(log(2)^(1/shape)) * scale # 1/gamma

curve(weibHaz(x, shape=shape, scale=scale)
      , from=0, to=80, ylab="Hazard", xlab="Time"
      , col="red")
curve(weibHaz(x, shape=shape, scale=scale)
      , from=0, to=80, ylab="Hazard", xlab="Time"
      , col="black", add=T) 
curve(weibHaz(x, shape=shape, scale=scale)
      , from=0, to=80, ylab="Hazard", xlab="Time"
      , col="blue", add=T) 

# gamma dist
gammaHaz <- function(x, shape, scale) {
  dgamma(x, shape=shape, scale=scale) / 
    pgamma(x, shape=shape, scale=scale, lower.tail=F)
} 
curve(gammaHaz(x, shape=shape, scale=scale), from=0, to=80,
      ylab="Hazard", xlab="Time", col="red")

meanlog <- 0; sdlog <- 1
lnormSurv <- function(t, shape, scale) {
  plnorm(t, meanlog = meanlog, sdlog = sdlog, lower.tail=F) 
}

curve(lnormSurv(x, meanlog, sdlog)
      , from=0, to=80, ylim=c(0,1)
      , ylab="Survival probability", xlab="Time") 

lnormHaz <- function(x, meanlog, sdlog) {
  dlnorm(x, meanlog = meanlog, sdlog = sdlog) / 
    plnorm(x, meanlog = meanlog, sdlog = sdlog, lower.tail=F)
}

curve(lnormHaz(x, meanlog, sdlog)
      , from=0, to=80
      , ylab="Survival probability", xlab="Time") 

tm <- c(0, # birth
        1/365, # first day of life
        7/365, # seventh day of life
        28/365, # fourth week of life
        1:106) # subsequent years
hazMale <- survexp.us[,"male","2004"] # 2004 males hazard function
hazFemale <- survexp.us[,"female","2004"] # 2004 females hazard function

tm.diff <- diff(tm)
survMale <- exp(-cumsum(hazMale[-109]*tm.diff)*365.24)
survFemale <- exp(-cumsum(hazFemale[-109]*tm.diff)*365.24)
sum(survMale*tm.diff) # mean age of male death in 2004
sum(survFemale*tm.diff) # mean age of female death in 2004

hazMale1940 <- survexp.us[,"male","1940"] # 2004 males
hazFemale1940 <- survexp.us[,"female","1940"] # 2004 females

hazMale2000 <- survexp.us[,"male","2000"] # 2004 males
hazFemale2000 <- survexp.us[,"female","2000"] # 2004 females

plot(hazMale1940~tm, type="l")
lines(hazFemale1940~tm, col=2)
lines(hazMale2000~tm, col=3)
lines(hazFemale2000~tm, col=4)
legend("topleft"
       , legend=c("Male 1940", "Female 1940", "Male 2000", "Female 2000")
       , col=1:4, lty=1)

survMale1940 <- exp(-cumsum(hazMale1940[-109]*tm.diff)*365.24)
survFemale1940 <- exp(-cumsum(hazFemale1940[-109]*tm.diff)*365.24)
sum(survMale1940*tm.diff) # mean age of male death in 2004
sum(survFemale1940*tm.diff) # mean age of female death in 2004

survMale2000 <- exp(-cumsum(hazMale2000[-109]*tm.diff)*365.24)
survFemale2000 <- exp(-cumsum(hazFemale2000[-109]*tm.diff)*365.24)
sum(survMale2000*tm.diff) # mean age of male death in 2004
sum(survFemale2000*tm.diff) # mean age of female death in 2004

plot(survMale1940~tm[-109], type="l")
lines(survFemale1940~tm[-109], col=2)
lines(survMale2000~tm[-109], col=3)
lines(survFemale2000~tm[-109], col=4)
legend("bottomleft"
       , legend=c("Male 1940", "Female 1940", "Male 2000", "Female 2000")
       , col=1:4, lty=1)

hazWMale1940 <- survexp.usr[,"male", "white","1940"] 
hazBMale1940 <- survexp.usr[,"male", "black","1940"] 

hazWMale2000 <- survexp.usr[,"male", "white","2000"] 
hazBMale2000 <- survexp.usr[,"male", "black","2000"] 

plot(hazWMale1940~tm, type="l")
lines(hazBMale1940~tm, col=2)
lines(hazWMale2000~tm, col=3)
lines(hazBMale2000~tm, col=4)
legend("topleft"
       , legend=c("WMale 1940", "BMale 1940", "WMale 2000", "BMale 2000")
       , col=1:4, lty=1)

# piecewise hazard function
# instantaneous hazard
h <- function(x) {
  ifelse(x <= 5, 0.07, 0.14)
}

# cumulative hazard
set_hazard <- function(h) {
  H <- function(t) {
    integrand <- function(x) {sapply(x, h)}
    result <- integrate(integrand, 0, t)$value
    return(result)
  }
  return(H)
}
H <- set_hazard(h)

H(10)
# Survival is related to cum. hazard by exp(-H(t))
# for exponential, this is -lambda * t
S <- function(t) {exp(-H(t))}
# for exponention dist is it exp(-lambda * t)
S(10)
# probability density function is S(t) * h(t)
f <- function(t) {S(t) * h(t)}
f(10)

plot_h <- function(x) {sapply(x, h)}
plot_H <- function(x) {sapply(x, H)}
plot_S <- function(x) {sapply(x, S)}
plot_f <- function(x) {sapply(x, f)}

curve(plot_H(x)
      , from=0, to=10, ylab="Hazard", xlab="Time"
      , col="black") 
curve(plot_h(x)
      , from=0, to=10
      , col="red", add = T)
curve(plot_S(x)
      , from=0, to=10
      , col="blue", add=T)
curve(plot_f(x)
      , from=0, to=10
      , col="green", add=T)

# mean: integrate S() or 1/lambda
t_mu <- function(survfunc) {integrate(survfunc, 0, Inf)$value}
t_mu(plot_S) # one day I'll make this better

# median: value of t at which survival is 50%
t_med <- function(survfunc) {integrate(survfunc, 0, 5)$value}
t_med(plot_S) # not the right answer

# median calc - area under two rectangles or piecewise constant
a <- 0.07 * 5 # area under first piece
b <- log(2) # target quantity
c <- b - a # area of second rectangle
d <- c / 0.14 # length of c under second rectangle
5 + d # the right answer
