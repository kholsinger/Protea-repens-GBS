functions {
  real Max(real x, real y) {
    real z;
    if (x >= y) {
      z <- x;
    } else {
      z <- y;
    }
    return z;
  }
  real Min(real x, real y) {
    real z;
    if (x <= y) {
      z <- x;
    } else {
      z <- y;
    }
    return z;
  }
}
data {
  int<lower=1> nPops;
  int<lower=1> nLoci;
  int<lower=0> n[nPops,nLoci,3];
  real nuPT;
  real omegaPT;
  real nuPPi;
  real omegaPPi;
  real nuLT;
  real omegaLT;
  real nuLPi;
  real omegaLPi;
}
parameters {
  real<lower=0,upper=1> p[nPops,nLoci];
  real<lower=0,upper=1> pi[nLoci];
  real<lower=0,upper=1> thetaLocus[nLoci];
  real<lower=0,upper=1> thetaPop[nPops];
  real<lower=0,upper=1> thetaP;
  real<lower=0,upper=1> piP;
  real<lower=0,upper=1> thetaL;
  real<lower=0,upper=1> piL;
  // used to induce constrained prior on f[i,j]
  //
  real<lower=0,upper=1> w[nPops,nLoci];
  real<lower=0,upper=1> f[nPops,nLoci];
}
transformed parameters {
  vector<lower=0,upper=1>[3] x[nPops,nLoci];
  real<lower=0,upper=1> theta[nPops,nLoci];
  // real<lower=-1,upper=1> f[nPops,nLoci];
  real<lower=-1,upper=1> fMin[nPops,nLoci];
  real<lower=1> alpha[nPops,nLoci];
  real<lower=1> beta[nPops,nLoci];
  real<lower=1> alphaPop;
  real<lower=1> betaPop;
  real<lower=1> alphaLocus;
  real<lower=1> betaLocus;
  for (i in 1:nPops) {
    for (j in 1:nLoci) {
      // genotype frequencies in population i at locus j
      //
      x[i,j,1] <- p[i,j]*p[i,j] + f[i,j]*p[i,j]*(1-p[i,j]);
      x[i,j,2] <- 2*(1-f[i,j])*p[i,j]*(1 - p[i,j]);
      x[i,j,3] <- (1-p[i,j])*(1-p[i,j]) + f[i,j]*p[i,j]*(1-p[i,j]);
      // induced prior on inbreeding coefficient in population i at locus j
      //
      // f[i,j] <- w[i,j]*(1-fMin[i,j]) + fMin[i,j];
      fMin[i,j] <- Max(-p[i,j]/(1-p[i,j]), -(1-p[i,j])/p[i,j]);
      // parameters for prior on allele frequencies
      //
      theta[i,j] <- 1 - (1 - thetaPop[i])*(1 - thetaLocus[j]);
      alpha[i,j] <- Max(1, Min(((1-theta[i,j])/theta[i,j])*pi[j], 1.0e4));
      beta[i,j]  <- Max(1, Min(((1-theta[i,j])/theta[i,j])*(1-pi[j]), 1.0e4));
     }
   }
   // parameters for prior on population effect on theta
   //
   alphaPop <- Max(1, Min(((1-thetaP)/thetaP)*piP, 1.0e4));
   betaPop  <- Max(1, Min(((1-thetaP)/thetaP)*(1-piP), 1.0e4));
   // parameters for prior on locus effect on theta
   //
   alphaLocus <- Max(1, Min(((1-thetaL)/thetaL)*piL, 1.0e4));
   betaLocus  <- Max(1, Min(((1-thetaL)/thetaL)*(1-piL), 1.0e4));
}
model {
  for (i in 1:nPops) {
    for (j in 1:nLoci) {
      // likelihood of sample from population i at locus j
      //
      n[i,j] ~ multinomial(x[i,j]);
      // prior for allele frequency in population i at locus j
      //
      p[i,j] ~ beta(alpha[i,j], beta[i,j]);
    }
  }
  // prior for population effect on theta
  //
  for (i in 1:nPops) {
    thetaPop[i] ~ beta(alphaPop, betaPop);
  }
  thetaP ~ beta(nuPT, omegaPT);
  piP    ~ beta(nuPPi, omegaPPi);
  // prior for locus effect on theta
  //
  for (j in 1:nLoci) {
    thetaLocus[j] ~ beta(alphaLocus, betaLocus);
  }
  thetaL ~ beta(nuLT, omegaLT);
  piL    ~ beta(nuLPi, omegaLPi);
}