require(HardyWeinberg)

check <- FALSE

strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

get.f <- function(x) {
  p <- sum(x, na.rm=TRUE)/(2*length(x))
  obs <- sum(x == 1, na.rm=TRUE)
  exp <- 2*p*(1-p)*length(x)
  f <- 1.0 - obs/exp
  f
}

count.genos <- function(x) {
  genos <- numeric(3)
  genos[1] <- sum(x==0, na.rm=TRUE)
  genos[2] <- sum(x==1, na.rm=TRUE)
  genos[3] <- sum(x==2, na.rm=TRUE)
  genos
}

markers <- read.csv("loci_20.csv", header=TRUE, na.strings=".")
n.markers <- ncol(markers)

markers$pop <- as.factor(strip(as.character(markers$indiv)))
markers <- subset(markers, pop!="EMPTY")

locus <- colnames(markers)

f.vec <- numeric(0)
pvals <- numeric(0)
ct <- 0
for (pop.n in unique(markers$pop)) {
  x <- subset(markers, pop==pop.n)
  if (nrow(x) > 1) {
    cat(pop.n, ":\n", sep="")
    ## start at column 2 because column 1 has individual label
    ##
    for (l in 2:n.markers) {
      f <- get.f(x[,l])
      if (!is.na(f)) {
        chisq <- suppressWarnings(HWChisq(count.genos(x[,l]), cc=0))
        if (check) {
          if (abs(length(x[,l])*(f^2) - chisq$chisq) > 1.0e-6) {
            cat("f:           ", f, "\n",
                "nf^2:        ", length(x[,l])*(f^2), "\n",
                "Chi-squared: ", chisq$chisq, "\n", sep="")
            print(x[,l])
            stop()
          }
        }
        f.vec <- c(f.vec, f)
        pvals <- c(pvals, chisq$pval)
        cat("  ", locus[l], " - ", f, "\n")
        ct <- ct + 1
      }
    }
  }
}

hist(pvals)
dev.new()
hist(f.vec)

cat("Total no. of p-values calculated: ", ct, "\n", sep="")
