library(survival)
library(asaur)
library(survminer)
result.heart <- coxph(Surv(futime, fustat) ~ transplant + age + surgery, data = jasa)
summary(result.heart) # wrong because people who lived long enough to get a transplant also survived long

# landmark method is very arbitrary and have to ignore data that don't fit the landmark
ind30 <- jasa$futime >= 30
transplant30 <- ((jasa$transplant == 1) & (jasa$wait.time < 30))
summary(coxph(Surv(futime, fustat) ~ transplant30 + age + surgery, data = jasa, subset = ind30))

# time dependent covariate method
# look at just a few subjects
id <- 1:nrow(jasa)
jasaT <- data.frame(id, jasa)
id.simple <- c(2, 5, 10, 12, 28, 95)
heart.simple <- jasaT[id.simple,c(1, 10, 9, 6, 11)]
summary(coxph(Surv(futime, fustat) ~ transplant, data=heart.simple))
# use tmerge to create the start stop records
sdata <- tmerge(heart.simple, heart.simple, id=id
                , death=event(futime, fustat), transpl=tdc(wait.time))
heart.simple.counting <- sdata[,-(2:5)] # drop columns 2through 5
heart.simple.counting
summary(coxph(Surv(tstart, tstop, death) ~ transpl
              , data=heart.simple.counting))
# for all the data
id <- 1:nrow(jasa)
jasaT <- data.frame(id, jasa)
jasaT <- jasaT
jasaT <- jasaT[, -c(2:4, 11:14)]
# for all the data
jasaT$futime <- pmax(.5, heart$futime)
indx <- ((jasaT$wait.time == jasaT$futime) & !is.na(jasaT$wait.time))
jasaT$wait.time[indx] <- jasaT$wait.time[indx] - .5
# use tmerge to create the start stop records
sdata <- tmerge(jasaT, jasaT, id=id
                , death=event(futime, fustat), transplant=tdc(wait.time))
counting <- sdata[,c(1, 9:12, 4:5)] # drop columns 2through 5
head(counting)
# see that transplant is not signif, and all the data was kept
summary(coxph(Surv(tstart, tstop, death) ~ transplant + surgery + age
              , data=counting))
