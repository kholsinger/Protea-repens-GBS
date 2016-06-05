## to use
##
## markers <- read.marker.data("filename.csv", shuffle=TRUE/FALSE)
## fit <- analyze.data(markers)
## print.summary(fit)

rm(list=ls())

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

## read marker data from CSV file
## - filename: the name of the CSV file
## - shuffle:  randomly shuffle individuals and population assignment
##             if TRUE
##
read.marker.data <- function(filename, shuffle=FALSE) {
  markers <- read.csv(filename, header=TRUE, na.strings=".")

  markers$pop <- as.factor(strip(as.character(markers$indiv)))
  markers <- subset(markers, pop!="EMPTY")
  if (shuffle) {
    tmp <- markers$pop
    new.pop <- sample(seq(1:length(markers$pop)), replace=TRUE)
    for (i in 1:length(markers$pop)) {
      tmp[i] <- markers$pop[new.pop[i]]
    }
    markers$pop <- tmp
  }
  n.pops <-length(unique(markers$pop))
  n.loci <- ncol(markers)-2

  locus <- colnames(markers)

  n <- array(dim=c(n.pops, n.loci, 3))
  N <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    pop.n <- unique(markers$pop)[i]
    ## x has one row for each individual in the population
    ## each locus is in a different column
    x <- subset(markers, pop==pop.n)
    for (j in 1:n.loci) {
      ## +1 because indiv is in column 1
      n[i, j, ] <- count.genos(x[,j+1])
      N[i,j] <- sum(n[i,j,])
    }
  }
  list(n.pops=n.pops,
       n.loci=n.loci,
       n=n,
       N=N)
}

get.fst <- function(p) {
  fst <- var(p, na.rm=TRUE)/(mean(p, na.rm=TRUE)*(1-mean(p, na.rm=TRUE)))
  if (is.nan(fst)) {
    fst <- NA
  }
  fst
}

get.beta.pars <- function(n, N, tightness) {
  n.pops <- nrow(N)
  n.loci <- ncol(N)
  p <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    for (j in 1:n.loci) {
      p[i,j] <- (2*n[i,j,1] + n[i,j,2])/(2*N[i,j])
    }
  }
  ## population-specific effects
  ##
  fst.ij <- matrix(nrow=n.pops, ncol=n.pops)
  fst.l <- numeric(n.loci)
  for (i in 1:n.pops) {
    for (j in 1:n.pops) {
      if (i == j) {
        fst.ij[i,j] <- NA
      } else {
        for (l in 1:n.loci) {
          fst.l[l] <- get.fst(c(p[i,l], p[j,l]))
        }
        fst.ij[i,j] <- (n.pops/(n.pops-1))*mean(fst.l, na.rm=TRUE)
      }
    }
  }
  sigma.2 <- var(as.vector(fst.ij), na.rm=TRUE)
  mu <- mean(as.vector(fst.ij), na.rm=TRUE)
  theta.pop <- sigma.2/(mu*(1-mu))
  alpha.pop <- ((1-theta.pop)/theta.pop)*mu
  beta.pop  <- ((1-theta.pop)/theta.pop)*(1-mu)
  mu.pop <- alpha.pop/(alpha.pop + beta.pop)
  nu.p.pi <- ((1-tightness)/tightness)*mu.pop
  omega.p.pi <- ((1-tightness)/tightness)*(1-mu.pop)
  ## locus-specific effects
  ##
  fst <- numeric(n.loci)
  for (i in 1:n.loci) {
    fst[i] <- get.fst(p[,i])
  }
  sigma.2 <- var(fst, na.rm=TRUE)
  mu <- mean(fst, na.rm=TRUE)
  theta.locus <- sigma.2/(mu*(1-mu))
  alpha.locus <- ((1-theta.locus)/theta.locus)*mu
  beta.locus  <- ((1-theta.locus)/theta.locus)*(1-mu)
  mu.locus <- alpha.locus/(alpha.locus + beta.locus)
  nu.l.pi <- ((1-tightness)/tightness)*mu.locus
  omega.l.pi <- ((1-tightness)/tightness)*(1-mu.locus)
  list(nu.p.pi=nu.p.pi,
       omega.p.pi=omega.p.pi,
       nu.l.pi=nu.l.pi,
       omega.l.pi=omega.l.pi,
       max.fij=max(fst.ij, na.rm=TRUE),
       min.fij=min(fst.ij, na.rm=TRUE),
       max.fst=max(fst, na.rm=TRUE),
       min.fst=min(fst, na.rm=TRUE),
       fst=fst)
}

## markers - a list containing n, N, n.pops, and n.loc
##   n - an n.pops x n.loci x 3 array of genotype counts
##   N - an n.pops x n.loci array of sample sizes
##   n.pops - number of populations in the sample
##   n.loci - number of loci scored
## nu - first parameter of beta specifying "tightness" of prior for
##      mean locus- and population-specific effect
## omega - second parameter of beta specifying "tightness" of prior for
##      mean locus- and population-specific effect
## digits - number of digits to display in print outs
##
analyze.data <- function(markers,
                         n.sample=250000,
                         n.burnin=50000,
                         n.thin=250,
                         n.chains=5,
                         nu=1,
                         omega=9,
                         digits=3,
                         shuffle=FALSE)
{
  require(R2jags)

  n <- markers$n
  N <- markers$N
  n.pops <- markers$n.pops
  n.loci <- markers$n.loci
  beta.pars <- get.beta.pars(n, N, nu/(nu+omega))
  nu.l.t <- nu
  omega.l.t <- omega
  nu.l.pi <- beta.pars$nu.l.pi
  omega.l.pi <- beta.pars$omega.l.pi
  nu.p.t <- nu
  omega.p.t <- omega
  nu.p.pi <- beta.pars$nu.p.pi
  omega.p.pi <- beta.pars$omega.p.pi
  cat("pi.l: ", round(nu.l.pi/(nu.l.pi+omega.l.pi), digits),
      " (", round(nu.l.pi, digits), ",", round(omega.l.pi, digits),")", "\n",
      "  locus Fst: (", round(beta.pars$min.fst, digits), ",",
      round(beta.pars$max.fst, digits), ")\n",
      "pi.p: ", round(nu.p.pi/(nu.p.pi+omega.p.pi), digits),
      " (", round(nu.p.pi, digits), ",", round(omega.p.pi, digits),")", "\n",
      "  pop Fst:   (", round(beta.pars$min.fij, digits), ",",
      round(beta.pars$max.fij, digits), ")\n",
      sep="")

  n.sample <- n.sample
  n.thin <- n.thin
  n.chains <- n.chains
  n.iter   <- n.sample + n.burnin

  jags.data <- c("n",
                 "N",
                 "n.pops",
                 "n.loci",
                 "nu.l.t",
                 "omega.l.t",
                 "nu.l.pi",
                 "omega.l.pi",
                 "nu.p.t",
                 "omega.p.t",
                 "nu.p.pi",
                 "omega.p.pi")
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
              model.file="f-statistics-with-fij.txt",
              n.chains=n.chains,
              n.iter=n.iter,
              n.burnin=n.burnin,
              n.thin=n.thin)
  fit
}

stan.inits.one <- function(chain_id, n.pops, n.loci) {
  p <- matrix(nrow=n.pops, ncol=n.loci)
  w <- matrix(nrow=n.pops, ncol=n.loci)
  f <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    p[i,] <- runif(n.loci)
    for (j in 1:n.loci) {
      ## sets inital f[i,j] = 0 for all i & j
      ##
      f.min <- max(-p[i,j]/(1-p[i,j]), -(1-p[i,j])/p[i,j])
      w[i,j] <- -f.min/(1 - f.min)
    }
    f[i,] <- runif(n.loci)
  }
  pi <- runif(n.loci)
  thetaLocus <- runif(n.loci)
  thetaPop <- runif(n.pops)
  piL <- var(thetaLocus)/(mean(thetaLocus)*(1-mean(thetaLocus)))
  thetaL <- 0.1
  piP <- var(thetaPop)/(mean(thetaPop)*(1-mean(thetaPop)))
  thetaP <- 0.1
  list(p=p, w=w, f=f, pi=pi, thetaLocus=thetaLocus, thetaPop=thetaPop,
       thetaP=thetaP, piP=piP, thetaL=thetaL, piL=piL)
}

stan.inits <-function(n.pops, n.loci, n.chains) {
  inits.list <- lapply(1:n.chains, stan.inits.one, n.pops, n.loci)
  inits.list
}

## markers - a list containing n, N, n.pops, and n.loc
##   n - an n.pops x n.loci x 3 array of genotype counts
##   N - an n.pops x n.loci array of sample sizes
##   n.pops - number of populations in the sample
##   n.loci - number of loci scored
## nu - first parameter of beta specifying "tightness" of prior for
##      mean locus- and population-specific effect
## omega - second parameter of beta specifying "tightness" of prior for
##      mean locus- and population-specific effect
## digits - number of digits to display in print outs
##
analyze.data.stan <- function(markers,
                              n.sample=5000,
                              n.burnin=1000,
                              n.thin=1,
                              n.chains=5,
                              nu=1,
                              omega=9,
                              digits=3)
{
  require(rstan)

  n <- markers$n
  N <- markers$N
  nPops <- markers$n.pops
  nLoci <- markers$n.loci
  beta.pars <- get.beta.pars(n, N, nu/(nu+omega))
  nuLT <- nu*10
  omegaLT <- omega*10
  nuLPi <- beta.pars$nu.l.pi
  omegaLPi <- beta.pars$omega.l.pi
  nuPT <- nu
  omegaPT <- omega
  nuPPi <- beta.pars$nu.p.pi
  omegaPPi <- beta.pars$omega.p.pi
  cat("piL: ", round(nuLPi/(nuLPi+omegaLPi), digits),
      " (", round(nuLPi, digits), ",", round(omegaLPi, digits),")", "\n",
      "  locus Fst: (", round(beta.pars$min.fst, digits), ",",
      round(beta.pars$max.fst, digits), ")\n",
      "piP: ", round(nuPPi/(nuPPi+omegaPPi), digits),
      " (", round(nuPPi, digits), ",", round(omegaPPi, digits),")", "\n",
      "  pop Fst:   (", round(beta.pars$min.fij, digits), ",",
      round(beta.pars$max.fij, digits), ")\n",
      sep="")

  n.sample <- n.sample
  n.thin <- n.thin
  n.chains <- n.chains
  n.iter   <- n.sample + n.burnin

  data <- list(n=n,
               N=N,
               nPops=nPops,
               nLoci=nLoci,
               nuLT=nuLT,
               omegaLT=omegaLT,
               nuLPi=nuLPi,
               omegaLPi=omegaLPi,
               nuPT=nuPT,
               omegaPT=omegaPT,
               nuPPi=nuPPi,
               omegaPPi=omegaPPi)
  params <- c("f",
              "thetaPop",
              "thetaP",
              "piP",
              "thetaLocus",
              "thetaL",
              "piL")

  fit <- stan(data=data,
              pars=params,
              file="f-statistics.stan",
              init=stan.inits(nPops, nLoci, 1),
              iter=n.iter,
              warmup=n.burnin,
              thin=n.thin,
              chains=n.chains)

  fit
}


print.summary <- function(fit) {
  old.width <- options(width=180)
  print(fit, digits=3)
  options(width=old.width$width)
}
