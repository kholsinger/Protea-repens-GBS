require(ggplot2)

strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

genos <- read.csv("../loci_20.csv", na.strings=".", header=TRUE)
genos$pop <- as.factor(strip(rownames(genos)))
n.pops <- length(unique(genos$pop))
pops <- seq(1:n.pops)

locs <- read.table("repens.loc")

for.plot <- data.frame(pop=genos$pop,
                       x=locs$V7,
                       y=locs$V8)

p <- ggplot(for.plot, aes(x=x, y=y, color=pop, shape=pop)) +
     geom_point() +
     xlab("Axis 1") + ylab("Axis 2") + scale_shape_manual(values=pops)
print(p)
ggsave("SPA-repens.pdf")
