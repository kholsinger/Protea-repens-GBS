require(HardyWeinberg)

check <- TRUE

strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

get.f <- function(x) {
  p <- sum(x)/(2*length(x))
  obs <- sum(x == 1)
  exp <- 2*p*(1-p)*length(x)
  f <- 1.0 - obs/exp
}

count.genos <- function(x) {
  genos <- numeric(3)
  genos[1] <- sum(x==0)
  genos[2] <- sum(x==1)
  genos[3] <- sum(x==2)
  genos
}

markers <- read.csv("loci_90_t.csv", header=TRUE, na.strings=".")
n.markers <- ncol(markers)

markers$pop <- as.factor(strip(rownames(markers)))

locus <- colnames(markers)

f.vec <- numeric(0)
pvals <- numeric(0)
for (pop.n in unique(markers$pop)) {
  x <- subset(markers, pop==pop.n)
  if (nrow(x) > 1) {
    cat(pop.n, ":\n", sep="")
    for (l in 1:n.markers) {
      f <- get.f(x[,l])
      if (!is.na(f)) {
        chisq <- HWChisq(count.genos(x[,l]), cc=0)
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
      }
    }
  }
}

hist(pvals)
dev.new()
hist(f.vec)
