strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

get.beta.pars <- function(x) {
  mu <- mean(x)
  sigma.2 <- var(x)
  alpha <- mu*((mu*(1-mu)/sigma.2) - 1)
  beta  <- (1-mu)*((mu*(1-mu)/sigma.2) - 1)
  list(alpha=alpha,
       beta=beta)
}

get.kld <- function(alpha.0, beta.0, alpha.1, beta.1) {
  log(beta(alpha.1,beta.1)/beta(alpha.0,beta.0)) +
    (alpha.0 - alpha.1)*(digamma(alpha.0) - digamma(alpha.0+beta.0)) +
    (beta.0 - beta.1)*(digamma(beta.0) - digamma(alpha.0+beta.0))
}

get.p <- function(kld) {
  q <- exp(-2.0*(kld - log(0.5)))
  p <- (1 - sqrt(1 - 4.0*q))/2.0
  return(p)
}

markers <- read.csv("loci_20.csv", header=TRUE)
markers$pop <- as.factor(strip(as.character(markers$indiv)))
markers <- subset(markers, pop!="EMPTY")
pops <- unique(markers$pop)
n.pops <-length(unique(markers$pop))

load("thetas-JAGS.Rsave")

sink("outliers-JAGS-revised.txt")
lo <- quantile(pi.l, probs=0.025)
hi <- quantile(pi.l, probs=0.975)
cat("Locus outliers:  ", mean(pi.l), " (", lo, ",", hi, ")...\n", sep="")
bias <- 0.0001
threshold <- log(0.5) - 0.5*log(bias*(1-bias))

alpha.0 <- ((1-mean(theta.l))/mean(theta.l))*mean(pi.l)
beta.0 <-  ((1-mean(theta.l))/mean(theta.l))*(1-mean(pi.l))
for (i in 1:ncol(theta.locus)) {
  beta.pars <- get.beta.pars(theta.locus[,i])
  kld <- get.kld(alpha.0, beta.0, beta.pars$alpha, beta.pars$beta)
  if (kld > threshold) {
    lo <- quantile(theta.locus[,i], probs=0.025)
    hi <- quantile(theta.locus[,i], probs=0.975)
    cat("  ", colnames(markers)[i+1], "[", i, "]: ", mean(theta.locus[,i]),
        " (", lo, ",", hi, ") ",
        "  P=", get.p(kld), "\n",
        sep="")
  }
}
sink()
