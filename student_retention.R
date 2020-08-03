library(survival)
library(asaur)
library(survminer)
library(MASS)
library(car)
library(ggplot2)
set.seed(12345)

expit <- function(o) exp(o)/(1 + exp(o))

n <- 80 # number of students
hs_dip <- sample(1:5, n, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
ethnicity <- factor(sample(1:5, n, replace = TRUE
                      , prob = c(0.6, 0.15, 0.15, 0.05, 0.05))) # guess the background. not important
age <- 16 + rnegbin(n, 10, theta = 4)
gender <- factor(sample(1:4, n, replace = TRUE
                        , prob = c(0.48, 0.47, 0.02, 0.03))) # male, female, non-bin, rather not
clubmemb <- sample(0:3, n, replace = TRUE
                   , prob = c(0.5, 0.1, 0.3, 0.1))
residency <- factor(sample(1:3, n, replace = TRUE
                           , prob = c(0.65, 0.05, 0.3))) # NON-RESIDENT, RESIDENT, INTERNATIONAL
oc <- factor(character(n), levels = c("no", "yes")) # orientation course TRUE FALSE
nresid <- sum(residency == 1)
# probs of non resid doing a orientation course
lodds <- log(0.25) + rnorm(nresid, sd = 0.25) - 0.01 * age[residency == 1]
probs <- expit(lodds)
oc[residency == 1] <- ifelse(runif(nresid) < probs, "yes", "no")
# probs of resid/international doing a orientation course
lodds <- log(0.5) + rnorm(n - nresid, sd = 0.25) - 0.05 * age[residency != 1]
probs <- expit(lodds)
oc[residency != 1] <- ifelse(runif(n - nresid) < probs, "yes", "no")
bridge <- factor(character(n), levels = c("no", "yes"))
bridge[hs_dip %in% 1:3] <- "no"
bridge[hs_dip == 4] <- sample(c("no", "yes"), length(bridge[hs_dip == 4]), replace = TRUE
                              , prob = c(0.8, 0.2))
bridge[hs_dip == 5] <- sample(c("no", "yes"), length(bridge[hs_dip == 5]), replace = TRUE
                              , prob = c(0.6, 0.4))
visitliblogi <- sample(c("no", "yes"), n, replace = TRUE
                       , prob = c(0.3, 0.7)) # for the zero inflation
visitlib <- rep(0, n)
visitlib[visitliblogi == "yes"] <- rpois(sum(visitliblogi == "yes"), 4)

medlogi <- sample(c("no", "yes"), n, replace = TRUE
       , prob = c(0.9, 0.1)) # for the zero inflation
med <- rep(0, n)
med[medlogi == "yes"] <- rpois(sum(medlogi == "yes"), 1) # interactions with student services on medical grounds

perlogi <- sample(c("no", "yes"), n, replace = TRUE
                  , prob = c(0.9, 0.1))# for the zero inflation
per <- rep(0, n)
per[perlogi == "yes"] <- rpois(sum(perlogi == "yes"), 1) # interactions with student services on personal grounds
per[per > 0 & hs_dip > 3] <- ifelse(runif(per[per > 0 & hs_dip > 3]) < 0.5, per[per > 0 & hs_dip > 3] + 1, per[per > 0 & hs_dip > 3])
per[per > 0 & age > 28] <- ifelse(runif(per[per > 0 & age > 28]) < 0.5, per[per > 0 & age > 28] + 1, per[per > 0 & hs_dip > 3])

finlogi <- sample(c("no", "yes"), n, replace = TRUE
                  , prob = c(0.9, 0.1))# for the zero inflation
fin <- rep(0, n)
fin[finlogi == "yes"] <- rpois(sum(finlogi == "yes"), 1) # interactions with student services on financial  grounds
fin[fin > 0 & residency == 3] <- ifelse(runif(fin[fin > 0 & residency == 3]) < 0.5, fin[fin > 0 & hs_dip > 3] + 1, fin[fin > 0 & hs_dip > 3])
fin[fin > 0 & age > 35] <- ifelse(runif(fin[fin > 0 & age > 35]) < 0.5, fin[fin > 0 & age > 35] + 1, fin[fin > 0 & age > 35])

shape_a <-abs(((-((age-20)/15)^2) - ifelse(oc == "no", 3, 1) - per * 2 - fin + rnorm(n, 11))/3)
ft_pt <- factor(ifelse(residency != 1 | age < 20, "ft"
                       , ifelse(rbeta(n, shape_a, 1) < 0.75, "pt", "ft"))
                , levels = c("ft", "pt"))

shape_a <- ifelse(residency == 1, 5, 10) # residential slightly higher attend
shape_a <- ifelse(runif(n) < 0.5, abs(shape_a + (((age-25)/5)^2)*2.5 + rnorm(n)), shape_a)
shape_a <- abs(ifelse(ft_pt == "pt", shape_a - 0.75, shape_a + 0.25) - per*2 - med*3)
attendance <- rbeta(n, shape_a, 1)

week4_test <- expit((logit(
  ifelse(attendance > 0.9 & runif(n) < 0.6, attendance - 0.2, attendance)) +
  runif(n, -0.2, 0.2) + rnorm(n, -1) + visitlib -
  (ifelse(runif(n) < 0.5, age/15, 0)) +
  (ifelse(bridge == "yes" & runif(n) < 0.5, 1, 0)))/5)

first_sem <- expit(logit(rowMeans(cbind(week4_test, attendance))) -
   ifelse(runif(n) < 0.25, fin/2, 0) -
   ifelse(runif(n) < 0.5, per/2, 0) -
   ifelse(runif(n) < 0.75, med/2, 0) +
   runif(n, -0.2, 0.2) + rnorm(n, -0.5))
first_sem <- ifelse(first_sem > 0.5, "pass", "fail")

prog_length <- 52 * 3
weeks <- 1:prog_length

create_events <- function(ev) {
  sapply(ev, function(x) {
    if(x == 0) {
      NA 
    } else {
      sort(sample(weeks, size = x, replace = FALSE, prob = (weeks - mean(weeks))^2))
    }
  })
}

medviz <- create_events(med)
perviz <- create_events(per)
finviz <- create_events(fin)

students <- data.frame(hs_dip, ethnicity, age, gender
                       , clubmemb, residency, oc, bridge, visitlib
                       , med, per, fin, ft_pt
                       , attendance, week4_test, first_sem)

# underlying hazard
# using the built in weibull dist
weibSurv <- function(t, shape, scale) {
  pweibull(t, shape=shape, scale=scale, lower.tail=F) 
}
weibHaz <- function(x, shape, scale) {
  dweibull(x, shape=shape, scale=scale) / 
    pweibull(x, shape=shape, scale=scale, lower.tail=F)
}

# trial and error to get these shape and scale values
shape <- 0.5
scale <- 750

curve(weibSurv(x, shape=shape, scale=scale)
      , from=1, to=prog_length
      , ylim=c(0,1)
      , ylab="Survival probability", xlab="Time"
      , col="blue")

# this is the hazard function I would like as a base:
# drop outs are more likely to happen early on.
curve(weibHaz(x, shape=shape, scale=scale)
      , from=1, to=prog_length
      , ylab="Hazard", xlab="Time"
      , col="blue")

check_events <- function(t, s, event_list) {
  event_t <- event_list[s]
  find_event <- function(t) {
    length(which(unlist(event_t) <= t))
  }
  sapply(t, find_event)  
}

multiplier_s <- function(student, time) {
  (1 - (3 - students[student, "hs_dip"])/100) *
    exp(-(students[student, "attendance"] - 0.9)/5) *
    ifelse(students[student, "gender"] == 1, 1.01, 1) *
    (1 - students[student, "clubmemb"]/10) *
    ifelse(students[student, "residency"] == 1, 1.01, 1) *
    ifelse(students[student, "oc"] == "yes", 0.95, 1) *
    ifelse(students[student, "hs_dip"] > 3 & students[student, "bridge"] == "no", 1.05, 1) *
    (1 - students[student, "visitlib"]/100) *
    (1 + students[student, "med"]/100 + check_events(time, student, medviz)/2) *
    (1 + students[student, "fin"]/100 + check_events(time, student, finviz)/2) *
    (1 + students[student, "per"]/100 + check_events(time, student, perviz)/2) *
    ifelse(students[student, "ft_pt"] == "ft", 1, 1.01) *
    exp(-(students[student, "week4_test"] - 0.7)/20) *
    ifelse(students[student, "first_sem"] == "pass", 1, 1.05) 
}

# undefined for t = 0
# adjusts heavier after med viz and per
indiv_h <- function(t, student) {
  weibHaz(t, shape=shape, scale=scale) *
  multiplier_s(student, t)
}

indiv_S <- function(time, student = student) {
  h <- function(t) indiv_h(t, student)
  return(exp(-integrate(h, 0, time)$value))
}

student_surv_func <- function(s) {
  function(t) {
    sapply(t, indiv_S, student = s)
  }
}

dropout <- integer(n)
for(s in 1:n) { # students
  chances <- runif(prog_length)^prog_length
  surv_func <- student_surv_func(s)
  surv_curve <- surv_func(weeks)
  when <- which(chances > surv_curve)
  dropout[s] <- ifelse(length(when) == 0, 0, when[1])
}

students <- cbind(students
                  , startweek = 1
                  , endweek = ifelse(dropout == 0, prog_length, dropout)
                  , drops = ifelse(dropout == 0, 0, 1))

surv_students <- with(students, Surv(endweek - startweek, drops))
survfit <- with(students, survfit(surv_students ~ gender, conf.type = "log-log"))
survfit
summary(survfit)
plot(survfit, conf.int=FALSE, col=1:5
     , lwd = 2)

ggsurvplot(fit = survfit
       , data = students
       , break.time.by = 52
       , color = "strata"
       , surv.scale = "percent"
       , censor = T
       , legend = "right"
       , xlab = "Time (weeks)"
       , ylab = "Active % of cohort"
       , risk.table = TRUE
       , ggtheme = theme_bw())

s_79 <- student_surv_curve(79)
s_80 <- student_surv_curve(80)

curve(s_79
, from=weeks[1], to=prog_length
, ylim=c(0,1)
, ylab="Survival probability", xlab="Time"
, col="blue") 
# assigment_ext # assignment extensions
# vlm_inactivity

# everyone starts with a weibull for drop out and extend,
# modify by attendance (time dep), 
# first sem effect (very time dep), fin-per-med (per more than the others),
# no bridge (time dep). hs dip, club membership, residence


# status at end of study (response 1)
# dropout (response 2)
# most drops are at the quarters  
# end of first spring quarter is a precipice
# 12 quarters