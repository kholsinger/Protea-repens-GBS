rm(list=ls())

strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

read.marker.data <- function(filename) {
  markers <- read.csv(filename, header=TRUE, na.strings=".")

  markers$pop <- as.factor(strip(as.character(markers$indiv)))
  markers <- subset(markers, pop!="EMPTY")
  markers$indiv <- NULL
  n.pops <-length(unique(markers$pop))
  n.loci <- ncol(markers)-1
  list(markers=markers, n.pops=n.pops, n.loci=n.loci, n.indivs=nrow(markers))
}

get.freqs <- function(x) {
  n.loci <- ncol(x)
  p <- numeric(n.loci)
  for (i in 1:n.loci) {
    k <- sum(x[,i], na.rm=TRUE)
    n <- nrow(x) - length(x[is.na(x[,2]),2])
    p[i] <- k/(2*n)
  }
  p
}

get.ibd <- function(x, y, p) {
  n.loci <- length(x)
  f <- 0.0
  k <- 0
  for (i in 1:n.loci) {
    if (!is.na(x[i]) & !is.na(y[i])) {
      p.1 <- x[i]
      q.1 <- 1.0 - p.1
      p.2 <- y[i]
      q.2 <- 1.0 - p.2
      f <- f + (1.0 - (p.1*q.2 + q.1*p.2)/(2.0*p[i]*(1.0-p[i])))
      k <- k + 1
    }
  }
  f/k
}

dat <- read.marker.data("loci_20.csv")

p <- get.freqs(dat$markers[,1:dat$n.loci])
p.indiv <- matrix(nrow=dat$n.indivs, ncol=dat$n.loci)
cat("Calculating individual allele frequencies at each locus...\n")
markers <- as.matrix(dat$markers)
for (i in 1:dat$n.indivs) {
  for (j in 1:dat$n.loci) {
    p.indiv[i,j] <- as.numeric(markers[i,j])/2.0
  }
}
f.dist <- matrix(nrow=dat$n.indivs, ncol=dat$n.indivs)
cat("Calculating pairwise ibd distances...")
for (i in 1:dat$n.indivs) {
  f.dist[i,i] <- 0.0
  cat("\n from individual ", i, sep="")
  for (j in (i+1):dat$n.indivs) {
    if ((j %% 50) == 0) {
      cat(".", sep="")
    }
    f.dist[i,j] <- 1.0 - get.ibd(p.indiv[i,],
                                 p.indiv[j,],
                                 p)
    f.dist[j,i] <- f.dist[i,j]
  }
}
cat("\n")
print(min(f.dist[,]))
print(max(f.dist[,]))

