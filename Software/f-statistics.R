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

n.sample <- 25000
n.burnin <- 5000
n.thin   <- 5
n.chains <- 5

n.iter   <- n.sample + n.burnin

jags.data <- c("n",
               "N",
               "n.pops",
               "n.loci")
jags.params <- c("f",
                 "theta.pop",
                 "theta.locus")

fit <- jags2(data=jags.data,
            inits=NULL,
            parameters.to.save=jags.params,
            model.file="f-statistics.txt",
            n.chains=n.chains,
            n.iter=n.iter,
            n.burnin=n.burnin,
            n.thin=n.thin,
            clearWD=DEBUG)
print(fit)
