require(R2jags)

generate.data <- function(theta=0.2,
                          n.pops,
                          n.loci,
                          n.indiv) {
  source("generate-data.R")

  n.pops  <- 20
  n.loci  <- 100
  n.indiv <- 30
  n <- get.data(theta, n.pops, n.loci, n.indiv)
  N <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    for (j in 1:n.loci) {
      N[i,j] <- sum(n[i,j,])
    }
  }
  list(n=n, N=N)
}

strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

count.genos <- function(x) {
  genos <- numeric(3)
  genos[1] <- sum(x==0, na.rm=TRUE)
  genos[2] <- sum(x==1, na.rm=TRUE)
  genos[3] <- sum(x==2, na.rm=TRUE)
  genos
}

read.data <- function(filename) {
  markers <- read.csv(filename, header=TRUE, na.strings=".")

  markers$pop <- as.factor(strip(rownames(markers)))
  markers <- subset(markers, pop!="EMPTY")
  n.pops <-length(unique(markers$pop))
  n.loci <- ncol(markers)-1

  locus <- colnames(markers)

  n <- array(dim=c(n.pops, n.loci, 3))
  N <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    pop.n <- unique(markers$pop)[i]
    ## x has one row for each individual in the population
    ## each locus is in a different column
    x <- subset(markers, pop==pop.n)
    for (j in 1:n.loci) {
      n[i, j, ] <- count.genos(x[,j])
      N[i,j] <- sum(n[i,j,])
    }
  }
  list(n.pops=n.pops,
       n.loci=n.loci,
       n=n,
       N=N)
}

## n - an n.pops x n.loci x 3 array of genotype counts
## N - an n.pops x n.loci array of sample sizes
## n.pops - number of populations in the sample
## n.loci - number of loci scored
## nu - first parameter of beta for prior on theta.l and theta.p
## omega - second parameter of beta for prior on theta.l and theta.p
##
analyze.data <- function(n, N, n.pops, n.loci,
                         nu=1,
                         omega=3,
                         n.sample=250000,
                         n.burnin=50000,
                         n.thin=50,
                         n.chains=5)
{
  n <- n
  N <- N
  n.pops <- n.pops
  n.loci <- n.loci
  n.sample <- n.sample
  n.thin <- n.thin
  n.chains <- n.chains
  n.iter   <- n.sample + n.burnin

  jags.data <- c("n",
                 "N",
                 "n.pops",
                 "n.loci",
                 "nu",
                 "omega")
  jags.params <- c("f",
                   "theta.pop",
                   "theta.p",
                   "pi.p",
                   "theta.locus",
                   "theta.l",
                   "pi.l")

  fit <- jags(data=jags.data,
              inits=NULL,
              parameters.to.save=jags.params,
              model.file="f-statistics.txt",
              n.chains=n.chains,
              n.iter=n.iter,
              n.burnin=n.burnin,
              n.thin=n.thin)
  fit
}

print.summary <- function(fit) {
  old.width <- options(width=180)
  print(fit, digits=3)
  options(width=old.width$width)
}
