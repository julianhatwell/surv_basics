library(survival)
library(asaur)
library(survminer)
library(muhaz)
library(MASS)
library(car)
library(foreach)
library(doParallel)
library(ggplot2)
set.seed(12345)

expit <- function(o) exp(o)/(1 + exp(o))

n <- 8000 # number of students

prog_length <- 52 * 3
weeks <- 1:prog_length

semester_ends <- c(13, 26, 40)
semester_ends <- semester_ends + rep(c(0, 53, 105), each = 3)
semester_ends <- c(0, semester_ends)

# underlying hazard
# using the built in weibull dist
weibSurv <- function(t, shape, scale) {
  pweibull(t, shape=shape, scale=scale, lower.tail=F) 
}
weibHaz <- function(x, shape, scale) {
  dweibull(x, shape=shape, scale=scale) / 
    pweibull(x, shape=shape, scale=scale, lower.tail=F)
}

expSurv <- function(t, rate) {
  pexp(t, rate = rate, lower.tail=F) 
}
expHaz <- function(x, rate) {
  dexp(x, rate = rate) / 
    pexp(x, rate = rate, lower.tail=F)
}

# trial and error to get these shape and scale values
shape <- 0.3
scale <- 500000
# this is the hazard function I would like as a base:
# drop outs are more likely to happen early on.
# curve(weibHaz(x, shape=shape, scale=scale)
#       , from=1, to=prog_length
#       , ylab="Hazard", xlab="Time"
#       , col="blue")

# the hazard translate to this survival curve
# curve(weibSurv(x, shape=shape, scale=scale)
#       , from=1, to=prog_length
#       , ylim=c(0,1)
#       , ylab="Survival probability", xlab="Time"
#       , col="blue")

# rate <- 0.01
# curve(expHaz(x, rate = rate)
#       , from=1, to=prog_length
#       , ylab="Hazard", xlab="Time"
#       , col="blue")
# 
# # the hazard translate to this survival curve
# curve(expSurv(x, rate = rate)
#       , from=1, to=prog_length
#       , ylim=c(0,1)
#       , ylab="Survival probability", xlab="Time"
#       , col="blue")

# student synthetic data
hs_dip <- sample(1:5, n, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
ethnicity <- factor(sample(1:5, n, replace = TRUE
                      , prob = c(0.6, 0.15, 0.15, 0.05, 0.05))) # guess the background. not important
age <- 16 + rnegbin(n, 10, theta = 4)
gender <- factor(sample(1:4, n, replace = TRUE
                        , prob = c(0.48, 0.47, 0.01, 0.02))) # male, female, non-bin, rather not
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
lodds <- log(0.5) + rnorm(n - nresid, m = 3, sd = 0.25) - 0.05 * age[residency != 1]
probs <- expit(lodds)
oc[residency != 1] <- ifelse(runif(n - nresid) < probs, "yes", "no")
founda <- factor(character(n), levels = c("no", "yes"))
founda[hs_dip %in% 1:3] <- "no"
founda[hs_dip == 4] <- sample(c("no", "yes"), length(founda[hs_dip == 4]), replace = TRUE
                              , prob = c(0.8, 0.2))
founda[hs_dip == 5] <- sample(c("no", "yes"), length(founda[hs_dip == 5]), replace = TRUE
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
per[per > 0 & hs_dip > 3] <- ifelse(runif(length(per[per > 0 & hs_dip > 3])) < 0.5, per[per > 0 & hs_dip > 3] + 1, per[per > 0 & hs_dip > 3])
per[per > 0 & age > 28] <- ifelse(runif(length(per[per > 0 & age > 28])) < 0.5, per[per > 0 & age > 28] + 1, per[per > 0 & hs_dip > 3])

finlogi <- sample(c("no", "yes"), n, replace = TRUE
                  , prob = c(0.9, 0.1))# for the zero inflation
fin <- rep(0, n)
fin[finlogi == "yes"] <- rpois(sum(finlogi == "yes"), 1) # interactions with student services on financial  grounds
fin[fin > 0 & residency == 3] <- ifelse(runif(length(fin[fin > 0 & residency == 3])) < 0.5, fin[fin > 0 & residency == 3] + 1, fin[fin > 0 & residency == 3])
fin[fin > 0 & age > 35] <- ifelse(runif(length(fin[fin > 0 & age > 35])) < 0.5, fin[fin > 0 & age > 35] + 1, fin[fin > 0 & age > 35])

shape_a <- abs(((-((age-20)/15)^2) - ifelse(oc == "no", 3, 1) - per * 2 - fin + rnorm(n, 11))/3)
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
  (ifelse(founda == "yes" & runif(n) < 0.5, 1, 0)))/5)

first_sem <- expit(logit(rowMeans(cbind(week4_test, attendance))) -
                    (hs_dip/5 - 0.3) -
                    ifelse(runif(n) < 0.25, fin/2, 0) -
                    ifelse(runif(n) < 0.5, per/2, 0) -
                    ifelse(runif(n) < 0.75, med/2, 0) +
                    runif(n, -0.2, 0.2) + rnorm(n))

create_events <- function(ev) {
  sapply(ev, function(x) {
    if(x == 0) {
      NA 
    } else {
      sort(sample(weeks, size = x, replace = FALSE, prob = (weeks - mean(weeks))^2))
    }
  })
}

students <- data.frame(id = 1:n
                       , hs_dip, ethnicity, age, gender
                       , clubmemb, residency, oc, founda, visitlib
                       , med, per, fin, ft_pt
                       , attendance, week4_test, first_sem
                       , attend_disc = factor(ifelse(attendance < 0.50, "poor", ifelse(attendance < 0.75, "ok", "good")), levels = c("good", "ok", "poor"))
                       , first_sem_disc = factor(ifelse(first_sem < 0.50, "fail", ifelse(first_sem < 0.75, "pass", "dist")), levels = c("dist", "pass", "fail"))
                       , week4_test_disc = factor(ifelse(week4_test < 0.50, "fail", ifelse(week4_test < 0.75, "pass", "dist")), levels = c("dist", "pass", "fail"))
                       , age_disc = factor(ifelse(age < 21, "genz", ifelse(age < 35, "millenn", "geny")), levels = c("genz", "millenn", "geny")))

#individual survival curves
# helpers for the individual survival curves
medviz <- create_events(med)
perviz <- create_events(per)
finviz <- create_events(fin)

# causes hazard to knot upwards on event and gradually diminish thereafter
check_events <- function(t, s, event_list) {
  event_t <- event_list[s]
  find_event <- function(t) {
    max(0, length(which(unlist(event_t) <= t)) - t/prog_length)
  }
  sapply(t, find_event)  
}

# the proportional hazard multiplier depends on the student data
multiplier_s <- function(student, time) {
  (1 - (3 - students[student, "hs_dip"])/10) *
    exp(-(students[student, "attendance"] - 0.8)/1.75) *
    ifelse(students[student, "gender"] == 1, 1.02, 1) *
    (1 - students[student, "clubmemb"]/20) *
    ifelse(students[student, "residency"] == 1, 1.2, 0.9) *
    ifelse(students[student, "oc"] == "mo", 1.2, 0.9) *
    ifelse(students[student, "hs_dip"] > 3 & students[student, "founda"] == "no", 1.5, 1) *
    (1 - students[student, "visitlib"]/20) *
    (1 + students[student, "med"]/5 + check_events(time, student, medviz)) *
    (1 + students[student, "fin"]/5 + check_events(time, student, finviz)) *
    (1 + students[student, "per"]/5 + check_events(time, student, perviz)) *
    ifelse(students[student, "ft_pt"] == "ft", 1, 1.2) *
    ifelse(students[student, "week4_test_disc"] == "fail" & time < 20, 1.5, 1) *
    ifelse(students[student, "first_sem_disc"] == "fail", 1.5, 1)
}

# undefined for t = 0
# adjusts heavier after med viz and per
indiv_h <- function(t, student) {
  weibHaz(t, shape=shape, scale=scale) *
  multiplier_s(student, t)
}

indiv_S <- function(time, student = student) {
  h <- function(t) indiv_h(t, student)
  S <- exp(-integrate(h, 0, time)$value)
  return(expit(logit(S) +
                 (((time/ifelse(time < 20, 1, 1 + (time-20)/30)-15.5)^2)-150)/100
               ))
}

# the individual student survival curve generator
student_dropout_prob <- function(s) {
  function(t) {
    sapply(t, indiv_S, student = s)
  }
}

# check indiv surv funcs  
s_1 <- student_dropout_prob(1)
s_2 <- student_dropout_prob(2)
s_3 <- student_dropout_prob(3)
s_4 <- student_dropout_prob(4)
s_5 <- student_dropout_prob(5)
s_6 <- student_dropout_prob(6)
s_7 <- student_dropout_prob(7)
s_8 <- student_dropout_prob(8)
s_9 <- student_dropout_prob(9)
s_10 <- student_dropout_prob(10)
curve(s_2
      , from=weeks[1], to=prog_length
      , ylim=c(0,1)
      , ylab="Dropout probability", xlab="Time"
      , col="blue") 

# parallel data gen
cores <- detectCores() - 2
cl <- makeCluster(cores)
registerDoParallel(cl)

# first round raw drops
raw_dropout <- unlist(foreach(s=1:n, .packages = "car") %dopar% {
  chances <- runif(prog_length)^prog_length
  surv_func <- student_dropout_prob(s)
  surv_curve <- surv_func(weeks)
  when <- which(chances > surv_curve)
  # raw_dropout[s] <- ifelse(length(when) == 0, 0, when[1])
  ifelse(length(when) == 0, 0, when[1])
})

stopCluster(cl)

# semester effect
for(i in seq_along(semester_ends[-length(semester_ends)])) {
  raw_dropout <- ifelse(raw_dropout != 0 &
                      raw_dropout > semester_ends[i] &
                      raw_dropout < semester_ends[i + 1] &
                      runif(length(raw_dropout)) > 0.25
         , semester_ends[i + 1], raw_dropout)
}

endweek <- ifelse(raw_dropout == 0, prog_length, raw_dropout)
drops <- ifelse(raw_dropout == 0, 0, 1)
students <- cbind(students
                  , startweek = 1
                  , endweek
                  , drops)

# smooth haz and survival       
# plot survivals for different covariates
# survdiff for different covariates
# model with cox.ph
# high correl between attend, week4, sem_test, visitlib, PCA?
# select model and diagnose


surv_students <- with(students, Surv(endweek - startweek, drops))
sfit <- with(students
             , survfit(surv_students ~ 1
                       , conf.type = "log-log"))
sfit
summary(sfit)
plot(sfit)

alph <- 0.05
# finding upper quantile
sfit_qtl <- quantile(sfit, 0.05)
medn <- sfit_qtl$`quantile`
lcl <- sfit_qtl$lower
ucl <- sfit_qtl$upper

plot(sfit, conf.int=T, mark="+", xlab="Time in weeks"
     , ylab="Active % of Cohort"
     , xlim=c(0, 156))
abline(h = 0.975, col = "red", lty = 3)
lines(rep(x = ucl, 2), y = c(0, 1-alph), col = "blue", lty = 2)
lines(rep(x = lcl, 2), y = c(0, 1-alph), col = "blue", lty = 2)
lines(rep(x = medn, 2), y = c(0, 1-alph), col = "green", lty = 2)

# muhaz smoothed instantaneous hazard estimation
sfit.smoothhaz <- muhaz(students$endweek, students$drops, max.time=156
                        , bw.grid=35, bw.method="global", b.cor="none")
plot(sfit.smoothhaz)

# getting smooth survival function from the smoothed hazard
haz <- sfit.smoothhaz$haz.est
times <- sfit.smoothhaz$est.grid
# numerically evaluate the integral. diff creates the widths
# cumsum adds up the areas of the resulting rectangles
surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
lines(surv ~ times[1:(length(times) - 1)])

plot(sfit
     , col = "grey"
     , xlab="Time in weeks"
     , ylab="Active % of Cohort"
     , xlim=c(0, 156))
# smooth survival
lines(surv ~ times[1:(length(times) - 1)])

ggsurvevents(surv_students)

sfit.cox <- coxph(surv_students ~ fin +
                    hs_dip + 
                    founda +
                    week4_test + first_sem +
                    ft_pt +
                    visitlib +
                    attend_disc #+ strata(hs_dip)
                  , data = students
                  #, subset = (med < 4)
)
summary(sfit.cox)
# alterntive way to compare groups to cox prop haz
survdiff(surv_students ~ ft_pt
         , data = students, rho = 1) # rho 1 places more emph on early differences
# alterntive way to compare groups to cox prop haz
survdiff(surv_students ~ week4_test_disc + strata(age_disc)
         , data = students, rho = 1) # rho 1 places more emph on early differences

sfit <- with(students
             , survfit(surv_students ~ 1
                       #, subset = (visitlib < 9)
                       , conf.type = "log-log"))
# plot(sfit, conf.int=FALSE, col=1:5
# , lwd = 2)

ggsurvplot(fit = sfit
           , data = students
           #, conf.int = TRUE
           , break.time.by = 52
           , color = "strata"
           , surv.scale = "percent"
           , censor = T
           , legend = "right"
           , xlab = "Time (weeks)"
           , ylab = "Active % of cohort"
           , risk.table = TRUE
           , xlim = c(1, prog_length)
           , ylim = c(0.75, 1.0)
           , ggtheme = theme_bw())


# assigment_ext # assignment extensions
# vlm_inactivity

# most drops are at the quarters  


# chance of returning the following year
# for(yend in c(53, 105)) {
#   dropout <- ifelse(raw_dropout != 0 &
#                       raw_dropout > yend - 52 &
#                       raw_dropout < yend &
#                       runif(length(raw_dropout)) < 0.333
#                     , 0, raw_dropout) 
# }
# dropout <- ifelse(raw_dropout != 0 &
#                     raw_dropout < 80 &
#                     runif(length(raw_dropout)) < 0.333
#                   , 0, raw_dropout)

# restart <- ifelse((raw_dropout - dropout > 0 & raw_dropout - dropout < 30) |
#                     (raw_dropout - dropout > 52 & raw_dropout - dropout < 80)
#                   ,  ((raw_dropout - dropout) %/% 52 + 1) * 52 + 1, NA)

# returning effect (left truncation)
# endweek <- ifelse(raw_dropout == 0, prog_length, raw_dropout)
# drops <- ifelse(dropout == 0, 0, 1)
# students <- cbind(students
#                   , startweek = 1
#                   , endweek
#                   , drops
#                   , restart)
# for (i in 1:n) {
#   if(!is.na(restart[i])) {
#     new_row <- students[i, ]
#     new_row$startweek <- restart[i]
#     new_row$endweek <- 2 * 156 - restart[i]
#     students <- rbind(students, new_row)
#   }
# }

