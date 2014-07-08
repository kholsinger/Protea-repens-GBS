source("f-statistics.R")

load("full-results.Rsave")

theta.pop   <- fit$BUGSoutput$sims.list$theta.pop
theta.p     <- fit$BUGSoutput$sims.list$theta.p
pi.p        <- fit$BUGSoutput$sims.list$pi.p
theta.locus <- fit$BUGSoutput$sims.list$theta.locus
theta.l     <- fit$BUGSoutput$sims.list$theta.l
pi.l        <- fit$BUGSoutput$sims.list$pi.l

save(theta.pop, theta.p, pi.p, theta.locus, theta.l, pi.l,
     file="thetas.Rsave")

sink(file="full-results.txt")
old <- options(width=180)
print.summary(fit)
options(old)
sink()

