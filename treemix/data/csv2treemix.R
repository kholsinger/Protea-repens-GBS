rm(list=ls())

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

count.alleles <- function(x) {
  genos <- count.genos(x)
  alleles <- numeric(2)
  alleles[1] <- 2*genos[1] + genos[2]
  alleles[2] <- 2*genos[3] + genos[2]
  alleles
}

read.marker.data <- function(filename) {
  markers <- read.csv(filename, header=TRUE, na.strings=".")

  markers$pop <- as.factor(strip(as.character(markers$indiv)))
  markers <- subset(markers, pop!="EMPTY")
  n.pops <-length(unique(markers$pop))
  n.loci <- ncol(markers)-2

  locus <- colnames(markers)

  n <- array(dim=c(n.pops, n.loci, 2))
  names <- character(n.pops)
  for (i in 1:n.pops) {
    pop.n <- unique(markers$pop)[i]
    names[i] <- as.character(pop.n)
    ## x has one row for each individual in the population
    ## each locus is in a different column
    x <- subset(markers, pop==pop.n)
    for (j in 1:n.loci) {
      ## +1 because indiv is in column 1
      n[i, j, ] <- count.alleles(x[,j+1])
    }
  }
  list(n.pops=n.pops,
       names=names,
       n.loci=n.loci,
       n=n)
}

write.treemix <- function(x, outfile) {
  sink(file=outfile)
  for (i in 1:x$n.pops) {
    cat(x$names[i], " ")
  }
  cat("\n")
  for (j in 1:x$n.loci) {
    for (i in 1:x$n.pops) {
      count <- sprintf("%d,%d ", x$n[i, j, 1], x$n[i, j, 2])
      cat(count)
    }
    cat("\n")
  }
  sink()
}

csv2treemix <- function(infile, outfile) {
  allele.count <- read.marker.data(infile)
  write.treemix(allele.count, outfile)
}
