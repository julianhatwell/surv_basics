library(survival)
library(asaur)

# set initial parameter for exponential distrib
set_lambda <- function(lambda, type = "exp") {
  if (type == "exp") { # constant hazard
    h <- function(x) rep(lambda, length(x))
    return(h)  
  }
}
h <- set_lambda(0.5)
h(1) # returns vector of same length
h(1:5)
# cumulative hazard
get_cum_hazard <- function(h) {
  function(t) integrate(h, 0, t)$value
}
H <- get_cum_hazard(h)
H(5)
H(10)
# for exponential dist, H(t) = 
# lambda * t: area under rectangle

# Survival is related to cum. hazard by exp(-H(t))
# for exponential, this is -lambda * t
get_survfunc <- function(H) function(t) exp(-H(t))
S <- get_survfunc(H)
# for exponention dist is it exp(-lambda * t)
S(0.4)
S(4)
# probability density function is S(t) * h(t)
get_pdf <- function(S) function(t) sapply(t, S) * h(t)
pdf <- get_pdf(S)
pdf(0:4)

get_CDF <- function(pdf) function(t) integrate(pdf, 0, t)$value
CDF <- get_CDF(pdf)
CDF(4)
S(4) # complement of CDF

curve_h <- function(x) sapply(x, h)
curve_H <- function(x) sapply(x, H)
curve_S <- function(x) sapply(x, S)
curve_pdf <- function(x) sapply(x, pdf)
curve_CDF <- function(x) sapply(x, CDF)

my_plot <- function(plot_func, t
                    , n = 101
                    , fill = "lightblue"
                    , border = "black"
                    , ylab = "f(t)"
                    , xlab = "t"
                    ) {
  cum_t <- seq(0, t, length.out = n)
  plot(cum_t, plot_func(cum_t), type = "n"
       , ylim = c(0, min(1, max(plot_func(0), plot_func(t))))
       , ylab = ylab, xlab = xlab)
  polygon(c(0, cum_t, t)
          , c(0, plot_func(cum_t), 0)
          , col = fill
          , border = border
          )
  }

my_plot(curve_h, 5, ylab = "hazard")
my_plot(curve_H, 5, ylab = "cumhaz")
my_plot(curve_S, 5, ylab = "survfunc")
my_plot(curve_pdf, 5, ylab = "pdf")
my_plot(curve_CDF, 5, ylab = "CDF")

# mean: integrate S(). For exponential distrib it's 1/lambda
t_mu <- function(survfunc) {integrate(survfunc, 0, Inf)$value}
t_mu(curve_S) # one day I'll make this better
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
H <- set_hazard_wf(0.9, 0.5)
H(1)
# deriving the formula
# int alpha * lambda * (lambda * t)^(alpha-1) dt
# int alpha * lambda * lambda^(alpha-1) * t^(alpha-1) dt
# int alpha * lambda^alpha * t^(alpha-1) dt
# alpha * lambda^alpha * int t^(alpha-1) dt
# alpha * lambda^alpha * t^alpha / alpha
# lambda^alpha * t^alpha
# (t * lambda)^alpha ( + C)

S <- get_survfunc(H)
S(1:4)
pdf <- function(t) S(t) * h(t)
f(0.5)

