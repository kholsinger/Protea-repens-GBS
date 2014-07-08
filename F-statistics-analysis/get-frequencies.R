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

read.marker.data <- function(filename) {
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
  rownames(N) <- unique(markers$pop)
  ## -1 because last column is pop
  ##
  colnames(N) <- colnames(markers)[-1]
  list(n.pops=n.pops,
       n.loci=n.loci,
       n=n,
       N=N)
}

get.frequencies <- function(filename) {
  markers <- read.marker.data(filename)
  p <- matrix(nrow=markers$n.pops, ncol=markers$n.loci)
  for (i in 1:markers$n.pops) {
    for (j in 1:markers$n.loci) {
      p[i,j] = (2*markers$n[i,j,1] + markers$n[i,j,2])/(2*markers$N[i,j])
    }
  }
  rownames(p) <- rownames(markers$N)
  colnames(p) <- colnames(markers$N)
  p
}
