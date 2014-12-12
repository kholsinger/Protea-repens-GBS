substitute <- function(x) {
  y <- character(length(x))
  y[x==0] <- "11"
  y[x==1] <- "12"
  y[x==2] <- "22"
  y[is.na(x)] <- NA
  y
}

strip <- function(x) {
  y <- character(length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

## substitute my own hs() for the Hs() in adegenet to ignore missing data
##
hs <- function(x, truenames = TRUE) {
    if (is.genind(x)) {
        x <- genind2genpop(x, quiet = TRUE)
    }
    if (!is.genpop(x)) 
        stop("x is not a valid genpop object")
    if (x@type == "PA") 
        stop("not implemented for presence/absence markers")
    x.byloc <- seploc(x, truenames = truenames)
    lX <- lapply(x.byloc, function(e) makefreq(e, quiet = TRUE, 
        truenames = truenames)$tab)
    ## this is the only change - adding na.rm=TRUE to sum()
    ##
    lres <- lapply(lX, function(X) 1 - apply(X^2, 1, sum, na.rm=TRUE))
    res <- apply(as.matrix(data.frame(lres)), 1, mean)
    return(res)
}


dat <- read.csv("loci_20.csv", header=TRUE, na.strings=".")

df <- sapply(dat[-1], substitute)
pop <- as.character(strip(as.character(dat$indiv)))

markers <- df2genind(df, ind.names=dat$indiv, pop=pop)

hs(markers)
