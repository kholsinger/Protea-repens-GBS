get.freqs <- function(theta, npops=20, p=0.5) {
  alpha <- ((1-theta)/theta)*p
  beta  <- ((1-theta)/theta)*(1-p)
  freqs <- rbeta(npops, alpha, beta)
  freqs
}

get.genos <- function(p, n) {
  x <- numeric(3)
  x[1] <- p^2
  x[2] <- 2*p*(1-p)
  x[3] <- (1-p)^2
  k <- rmultinom(1, n, x)
  k
}

get.data <- function(theta, npops, nloci, nsample) {
  ## set up allele frequency distribution among populations using a
  ## single value of theta (Fst) across all loci
  ##
  pi <- matrix(nrow=nloci, ncol=npops)
  for (locus in 1:nloci) {
    pi[locus,] <- get.freqs(theta, npops)
  }
  ## get allele counts at each locus in each population
  ##
  counts <- array(dim=c(npops, nloci, 3))
  for (i in 1:npops) {
    for (j in 1:nloci) {
      counts[i, j, ] <- get.genos(pi[j, i], nsample)
    }
  }
  counts
}

generate.data <- function(theta=0.2,
                          n.pops=20,
                          n.loci=100,
                          n.indiv=30)
{
  n <- get.data(theta, n.pops, n.loci, n.indiv)
  N <- matrix(nrow=n.pops, ncol=n.loci)
  for (i in 1:n.pops) {
    for (j in 1:n.loci) {
      N[i,j] <- sum(n[i,j,])
    }
  }
  list(n.pops=n.pops,
       n.loci=n.loci,
       n=n,
       N=N)
}
