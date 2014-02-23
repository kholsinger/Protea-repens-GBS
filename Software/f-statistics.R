require(R2jags)

DEBUG <- TRUE

if (DEBUG) {
  source("generate-data.R")

  n.pops  <- 20
  n.loci  <- 100
  n.indiv <- 30
  n <- get.data(0.2, n.pops, n.loci, n.indiv)
  N <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    for (j in 1:n.loci) {
      N[i,j] <- sum(n[i,j,])
    }
  }
}

n.sample <- 250000
n.burnin <- 50000
n.thin   <- 50
n.chains <- 5

n.iter   <- n.sample + n.burnin

jags.data <- c("n",
               "N",
               "n.pops",
               "n.loci")
jags.params <- c("f",
                 "theta.pop",
                 "theta.p",
                 "pi.p",
                 "theta.locus",
                 "theta.l",
                 "pi.l")

fit <- jags2(data=jags.data,
            inits=NULL,
            parameters.to.save=jags.params,
            model.file="f-statistics.txt",
            n.chains=n.chains,
            n.iter=n.iter,
            n.burnin=n.burnin,
            n.thin=n.thin,
            clearWD=DEBUG)
old.width <- options(width=180)
print(fit, digits=3)
options(width=old.width$width)
