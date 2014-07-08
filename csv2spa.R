genos <- read.csv("loci_20.csv", na.strings=".", header=TRUE)

for (i in 1:nrow(genos)) {
  genos[i,is.na(genos[i,])] <- -1
}

write.table(genos, file="loci_20.gen", row.names=FALSE, col.names=FALSE)

